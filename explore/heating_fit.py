"""Fit topneho modulu na reziduich po bazi + diagnostika (krok 4)."""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import base
import etl
import heating

OUT = Path(__file__).resolve().parent

# paleta (dataviz reference, light mode)
SURFACE = "#fcfcfb"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, BASELINE = "#e1e0d9", "#c3c2b7"
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "axes.edgecolor": BASELINE, "axes.labelcolor": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "text.color": INK, "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "font.family": "sans-serif",
})

df = etl.load()
bparams = base.load()
df["resid"] = df["baseload"] - base.predict(df, bparams)

# fit bez chladici kontaminace: dny s dennim prumerem T <= 17 C
day_temp = df.groupby("date_local")["temp"].transform("mean")
mask = (day_temp <= 17.0).to_numpy()
print(f"radku pro fit: {mask.sum()} z {len(df)} ({mask.mean()*100:.0f} %)")

params = heating.fit(df, df["resid"].to_numpy(), mask)
path = heating.save(params)
print("ulozeno:", path)
print("\nnelinearni parametry:")
print(f"  a (osvit)      = {params['a']:.4f} C/(W/m2)")
print(f"  b (vitr)       = {params['b']:.4f} C/(m/s.C)")
print(f"  T_b            = {params['t_b']:.2f} C")
print(f"  s (sirka)      = {params['s']:.2f} C")
print(f"  w (rychly kanal)= {params['w']:.2f}")
print(f"  tau            = {params['tau']:.1f} h")
print(f"  alpha (COP)    = {params['alpha']:.4f} 1/C")
print(f"in-sample (fit subset): RMSE={params['rmse']:.1f}  R2 rezidua={params['r2']:.3f}")

df["p_heat"] = heating.predict(df, params)
df["resid2"] = df["resid"] - df["p_heat"]
cold = day_temp.to_numpy() < 10
print(f"\nchladne dny (prumer < 10 C): RMSE bez topeni {df.resid[cold].pow(2).mean()**.5:.1f}"
      f" -> s topenim {df.resid2[cold].pow(2).mean()**.5:.1f}")

theta = np.array([params[k] for k in heating.PARAM_NAMES])
ts = heating.t_star(df, params["a"], params["b"], params["w"], params["tau"])

# --- 10. topna krivka: rezidua vs. T* + model ---
fig, ax = plt.subplots(figsize=(8.5, 6))
sub = df[mask].sample(6000, random_state=1)
ax.scatter(ts[mask][df[mask].index.get_indexer(sub.index)], sub.resid,
           s=6, color=BLUE, alpha=0.25, lw=0, label="rezidua po bázi (15min)")
grid_t = np.arange(-12, 20, 0.25)
fake = pd.DataFrame({"temp": grid_t, "sun": 0.0, "wind": 0.0,
                     "tod": 12.0, "daytype": 0})
k_ref = float(np.mean(heating.k_basis(df[mask]) @ params["k_coef"]))
cop = np.maximum(1 + params["alpha"] * grid_t, heating.COP_FLOOR)
curve = k_ref * heating.softplus(params["t_b"] - grid_t, params["s"]) / cop
ax.plot(grid_t, curve, color=ORANGE, lw=2.5, label="model (průměrné k)")
ax.axhline(0, color=BASELINE, lw=0.8)
ax.set_xlabel("efektivní teplota T* [°C]"); ax.set_ylabel("spotřeba − báze")
ax.set_title(f"Topná křivka: T_b={params['t_b']:.1f} °C, s={params['s']:.1f}", loc="left")
ax.legend(frameon=False)
fig.savefig(OUT / "10_topeni_krivka.png", dpi=130, bbox_inches="tight")

# --- 11. k(cas dne) po typech dne ---
tod = np.arange(0, 24, 0.25)
fig, ax = plt.subplots(figsize=(8.5, 4.4))
labels = ["Po–Čt", "pátek/most", "sobota", "neděle/svátek"]
for k, color, lab in [(0, BLUE, labels[0]), (1, ORANGE, labels[1]),
                      (2, AQUA, labels[2]), (3, YELLOW, labels[3])]:
    ax.plot(tod, heating.k_curve(k, params, tod), color=color, lw=2, label=lab)
ax.set_xlabel("čas dne [h]"); ax.set_ylabel("k (výkonový koeficient)")
ax.set_xticks([0, 6, 12, 18, 24])
ax.set_title("Časový koeficient topení k(čas dne, typ dne)", loc="left")
ax.legend(frameon=False)
fig.savefig(OUT / "11_topeni_k.png", dpi=130, bbox_inches="tight")

# --- 12. diagnostika: rezidua vs. teplota + volna impulzni odezva ---
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.5, 4.4))
rd = df[mask].groupby("date_local").agg(r=("resid2", "mean"), t=("temp", "mean"))
a1.scatter(rd.t, rd.r, s=12, color=BLUE, alpha=0.5, lw=0)
a1.axhline(0, color=BASELINE, lw=0.8)
a1.set_xlabel("denní průměrná teplota [°C]"); a1.set_ylabel("denní reziduum po topení")
a1.set_title("Rezidua po topném modulu", loc="left")

# volna impulzni odezva: hodinova rezidua ~ lagovana T_ef anomalie (0-96 h, D2 ridge)
hh = pd.DataFrame({"r": df["resid2"].to_numpy(),
                   "te": heating.t_effective(df, params["a"], params["b"])},
                  index=df.index).resample("1h").mean()
hh = hh - hh.rolling("14D", center=True).mean()
lags = np.arange(0, 97)
Xl = np.column_stack([hh["te"].shift(l).to_numpy() for l in lags])
ok = ~np.isnan(Xl).any(axis=1) & ~hh["r"].isna().to_numpy()
Xl, yl = Xl[ok], hh["r"].to_numpy()[ok]
D2 = np.diff(np.eye(len(lags)), n=2, axis=0)
beta = np.linalg.solve(Xl.T @ Xl + 2e4 * D2.T @ D2, Xl.T @ yl)
a2.plot(lags, beta, color=BLUE, lw=2)
a2.axhline(0, color=BASELINE, lw=0.8)
a2.set_xlabel("zpoždění [h]"); a2.set_ylabel("odezva rezidua na T_ef")
a2.set_title("Volná impulzní odezva zbytku (má být ~0)", loc="left")
fig.savefig(OUT / "12_topeni_diagnostika.png", dpi=130, bbox_inches="tight")

# --- 13. studena vlna: pozorovani vs. baze vs. baze+topeni ---
rd7 = df.groupby("date_local")["temp"].mean().rolling(7).mean()
end = rd7.idxmin()
days = pd.date_range(end=pd.Timestamp(end), periods=10)
sel = df[pd.to_datetime(df["date_local"]).isin(days)]
t_loc = sel.index.tz_convert(etl.TZ)
fig, ax = plt.subplots(figsize=(12, 4.4))
ax.plot(t_loc, sel.baseload, color=INK2, lw=1.0, label="pozorování")
b_pred = base.predict(sel, bparams)
ax.plot(t_loc, b_pred, color=AQUA, lw=1.5, label="báze")
ax.plot(t_loc, b_pred + sel.p_heat, color=ORANGE, lw=1.5, label="báze + topení")
ax.set_ylabel("spotřeba (15min)")
ax.set_title(f"Nejchladnější týden ({days[0].date()} – {days[-1].date()})", loc="left")
ax.legend(frameon=False, ncol=3)
fig.savefig(OUT / "13_studena_vlna.png", dpi=130, bbox_inches="tight")
print("\ngrafy 10-13 ulozeny")
