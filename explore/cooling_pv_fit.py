"""Fit chlazeni + FVE (spolecne, sdileny osvit) a diagnostika (krok 5)."""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import base
import cooling
import etl
import fit
import heating
import pv

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
bp, hp = base.load(), heating.load()
df["p_base"] = base.predict(df, bp)
df["p_heat"] = heating.predict(df, hp)
resid = (df["baseload"] - df["p_base"] - df["p_heat"]).to_numpy()

cp, pp = fit.fit_cooling_pv(df, resid)
fit.save(cp, "cooling_params.npz")
fit.save(pp, "pv_params.npz")

df["p_cool"] = cooling.predict(df, cp)
df["g_pv"] = pv.predict(df, pp)
df["resid"] = resid - df["p_cool"] + df["g_pv"]

print("chlazeni:")
print(f"  a_c (osvit)     = {cp['a_c']:.4f} C/(W/m2)")
print(f"  T_bc            = {cp['t_bc']:.2f} C")
print(f"  s_c (sirka)     = {cp['s_c']:.2f} C")
print(f"  w_c (rychly)    = {cp['w_c']:.2f}")
print(f"  tau_c           = {cp['tau_c']:.1f} h")
print(f"  alpha_c (EER)   = {cp['alpha_c']:.4f} 1/C")
print("FVE:")
print(f"  gamma (derating)= {pp['gamma']:.4f} 1/C")
cap = pv.capacity_curve(pd.DatetimeIndex([df.index[0], df.index[-1]]), pp)
print(f"  kapacita: {cap[0]:.4f} -> {cap[1]:.4f} (jednotky spotreby na W/m2)")
print(f"  spicka pri 845 W/m2: {cap[0]*845:.1f} -> {cap[1]*845:.1f}")
print(f"  rocni vyroba: {df.g_pv.groupby(df.index.year).mean().round(1).to_dict()}")
print(f"spolecny fit: RMSE={cp['rmse']:.1f}  R2 zbytku po bazi+topeni={cp['r2']:.3f}")

day_t = df.groupby("date_local")["temp"].transform("mean").to_numpy()
hot = day_t > 20
print(f"\nteple dny (prumer > 20 C): RMSE {np.sqrt(np.mean(resid[hot]**2)):.1f}"
      f" -> {np.sqrt(np.mean(df.resid.to_numpy()[hot]**2)):.1f}")
print(f"cely rok: RMSE {np.sqrt(np.mean(resid**2)):.1f}"
      f" -> {np.sqrt(np.mean(df.resid**2)):.1f}")

# --- 14. chladici krivka ---
ts_c = cooling.t_star(df, cp["a_c"], cp["w_c"], cp["tau_c"])
fig, ax = plt.subplots(figsize=(8.5, 6))
r_pv = resid + df["g_pv"].to_numpy()  # zbytek po odectu FVE = jen chlazeni
idx = np.random.default_rng(1).choice(len(df), 6000, replace=False)
ax.scatter(ts_c[idx], r_pv[idx], s=6, color=BLUE, alpha=0.25, lw=0,
           label="zbytek po bázi, topení a FVE")
grid_t = np.arange(5, 36, 0.25)
k_ref = float(np.mean(cooling.k_basis(df) @ cp["k_coef"]))
eer = np.maximum(1 - cp["alpha_c"] * (grid_t - cooling.T_REF), cooling.EER_FLOOR)
ax.plot(grid_t, k_ref * heating.softplus(grid_t - cp["t_bc"], cp["s_c"]) / eer,
        color=ORANGE, lw=2.5, label="model (průměrné k)")
ax.axhline(0, color=BASELINE, lw=0.8)
ax.set_xlabel("efektivní teplota chlazení T*_c [°C]")
ax.set_ylabel("spotřeba − báze − topení + FVE")
ax.set_title(f"Chladicí křivka: T_bc={cp['t_bc']:.1f} °C, s={cp['s_c']:.1f}", loc="left")
ax.legend(frameon=False)
fig.savefig(OUT / "14_chlazeni_krivka.png", dpi=130, bbox_inches="tight")

# --- 15. kapacita FVE v case + k_c(cas dne) ---
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.5, 4.4))
days = pd.date_range(df.index[0], df.index[-1], freq="D")
a1.plot(days, pv.capacity_curve(days, pp) * 845, color=BLUE, lw=2)
a1.set_ylabel("výkon FVE při 845 W/m² (špička)")
a1.set_title("Odhadnutá kapacita FVE v čase", loc="left")
tod = np.arange(0, 24, 0.25)
labels = ["Po–Čt", "pátek/most", "sobota", "neděle/svátek"]
for k, color, lab in [(0, BLUE, labels[0]), (1, ORANGE, labels[1]),
                      (2, AQUA, labels[2]), (3, YELLOW, labels[3])]:
    a2.plot(tod, cooling.k_curve(k, cp, tod), color=color, lw=2, label=lab)
a2.set_xlabel("čas dne [h]"); a2.set_ylabel("k_c")
a2.set_xticks([0, 6, 12, 18, 24])
a2.set_title("Časový koeficient chlazení", loc="left")
a2.legend(frameon=False)
fig.savefig(OUT / "15_fve_kapacita.png", dpi=130, bbox_inches="tight")

# --- 16. rezidua po vsech clenech: vs. teplota a vs. osvit ---
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.5, 4.4))
rd = df.groupby("date_local").agg(r=("resid", "mean"), t=("temp", "mean"),
                                  s=("sun", "mean"))
a1.scatter(rd.t, rd.r, s=12, color=BLUE, alpha=0.5, lw=0)
a1.axhline(0, color=BASELINE, lw=0.8)
a1.set_xlabel("denní průměrná teplota [°C]"); a1.set_ylabel("denní reziduum")
a1.set_title("Rezidua po všech členech vs. teplota", loc="left")
a2.scatter(rd.s, rd.r, s=12, color=BLUE, alpha=0.5, lw=0)
a2.axhline(0, color=BASELINE, lw=0.8)
b1, b0 = np.polyfit(rd.s, rd.r, 1)
xx = np.array([rd.s.min(), rd.s.max()])
a2.plot(xx, b0 + b1 * xx, color=INK2, lw=1.2, ls="--")
a2.set_xlabel("denní průměrný osvit [W/m²]"); a2.set_ylabel("denní reziduum")
a2.set_title(f"Rezidua vs. osvit (sklon {b1:+.3f})", loc="left")
fig.savefig(OUT / "16_rezidua_vse.png", dpi=130, bbox_inches="tight")
print(f"sklon rezidui na osvit: {b1:+.3f} na W/m2 (ma byt ~0)")

# --- 17. nejteplejsi tyden: rozklad ---
d7 = df.groupby("date_local")["temp"].mean().rolling(7).mean()
days_hot = pd.date_range(end=pd.Timestamp(d7.idxmax()), periods=8)
sel = df[pd.to_datetime(df["date_local"]).isin(days_hot)]
t_loc = sel.index.tz_convert(etl.TZ)
fig, ax = plt.subplots(figsize=(12, 4.4))
ax.plot(t_loc, sel.baseload, color=INK2, lw=1.0, label="pozorování")
ax.plot(t_loc, sel.p_base, color=AQUA, lw=1.4, label="báze")
ax.plot(t_loc, sel.p_base + sel.p_cool, color=YELLOW, lw=1.4, label="+ chlazení")
ax.plot(t_loc, sel.p_base + sel.p_heat + sel.p_cool - sel.g_pv, color=ORANGE,
        lw=1.6, label="+ chlazení − FVE (model)")
ax.set_ylabel("spotřeba (15min)")
ax.set_title(f"Nejteplejší týden ({days_hot[0].date()} – {days_hot[-1].date()})", loc="left")
ax.legend(frameon=False, ncol=4)
fig.savefig(OUT / "17_horky_tyden.png", dpi=130, bbox_inches="tight")
print("grafy 14-17 ulozeny")
