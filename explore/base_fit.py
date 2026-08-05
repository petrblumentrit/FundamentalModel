"""Fit bazalni slozky z mirnych dnu + diagnostika (krok 3)."""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import base
import etl

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
params = base.fit(df)
path = base.save(params)
daily = base.daily_table(df)
mild = base.mild_dates(df)

print("mirnych dnu:", len(mild), "z", len(daily))
print("po typech dne:", daily.loc[mild].groupby("daytype").size().to_dict())
print(f"in-sample: RMSE={params['rmse']:.1f}  MAE={params['mae']:.1f}  R2={params['r2']:.3f}")
print("ulozeno:", path)

# --- 6. urovnova slozka v case + denni prumery mirnych dnu Po-Ct ---
fig, ax = plt.subplots(figsize=(11, 4.2))
m0 = daily.loc[mild].query("daytype == 0 and not prazdniny")
ax.scatter(m0.index, m0.baseload, s=14, color=MUTED, alpha=0.7, lw=0,
           label="mírné dny Po–Čt (pozorování)")
grid_days = pd.date_range(daily.index.min(), daily.index.max(), freq="D")
ax.plot(grid_days, base.level_curve(grid_days, params), color=BLUE, lw=2,
        label="úroveň báze (P-spline)")
ax.set_ylabel("denní průměrná spotřeba")
ax.set_title("Bazální úroveň v čase — fit z mírných dnů", loc="left")
ax.legend(frameon=False)
fig.savefig(OUT / "06_baze_uroven.png", dpi=130, bbox_inches="tight")

# --- 7. denni profily podle typu dne + prazdninova korekce ---
tod = np.arange(0, 24, 0.25)
ref_level = float(np.median(base.level_curve(pd.DatetimeIndex(mild), params)))
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.4),
                             gridspec_kw={"width_ratios": [3, 2]})
labels = ["Po–Čt", "pátek/most", "sobota", "neděle/svátek"]
for k, color, lab in [(0, BLUE, labels[0]), (1, ORANGE, labels[1]),
                      (2, AQUA, labels[2]), (3, YELLOW, labels[3])]:
    a1.plot(tod, ref_level + base.profile_curve(k, params), color=color, lw=2, label=lab)
a1.set_xlabel("čas dne [h]"); a1.set_ylabel("bazální spotřeba")
a1.set_xticks([0, 6, 12, 18, 24])
a1.set_title("Denní profily báze podle typu dne", loc="left")
a1.legend(frameon=False)
praz = base.profile_curve(0, params, prazdniny=True) - base.profile_curve(0, params)
a2.plot(tod, praz, color=BLUE, lw=2)
a2.axhline(0, color=BASELINE, lw=0.8)
a2.set_xlabel("čas dne [h]"); a2.set_ylabel("Δ spotřeba")
a2.set_xticks([0, 6, 12, 18, 24])
a2.set_title("Prázdninová korekce (7–8)", loc="left")
fig.savefig(OUT / "07_baze_profily.png", dpi=130, bbox_inches="tight")

# --- 8. rezidua na mirnych dnech: vs. cas dne a vs. teplota ---
sel = df[pd.to_datetime(df["date_local"]).isin(mild)].copy()
sel["resid"] = sel["baseload"] - base.predict(sel, params)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.2))
r_tod = sel.groupby("tod")["resid"].agg(["mean", "std"])
a1.fill_between(r_tod.index, r_tod["mean"] - r_tod["std"], r_tod["mean"] + r_tod["std"],
                color=BLUE, alpha=0.15, lw=0)
a1.plot(r_tod.index, r_tod["mean"], color=BLUE, lw=2)
a1.axhline(0, color=BASELINE, lw=0.8)
a1.set_xlabel("čas dne [h]"); a1.set_ylabel("reziduum")
a1.set_xticks([0, 6, 12, 18, 24])
a1.set_title("Rezidua vs. čas dne (průměr ± σ)", loc="left")
rd = sel.groupby("date_local").agg(resid=("resid", "mean"), temp=("temp", "mean"))
a2.scatter(rd.temp, rd.resid, s=14, color=BLUE, alpha=0.6, lw=0)
b1, b0 = np.polyfit(rd.temp, rd.resid, 1)
tt = np.array([rd.temp.min(), rd.temp.max()])
a2.plot(tt, b0 + b1 * tt, color=INK2, lw=1.2, ls="--")
a2.axhline(0, color=BASELINE, lw=0.8)
a2.set_xlabel("denní průměrná teplota [°C]"); a2.set_ylabel("denní reziduum")
a2.set_title(f"Rezidua vs. teplota (sklon {b1:+.1f}/°C)", loc="left")
fig.savefig(OUT / "08_baze_rezidua.png", dpi=130, bbox_inches="tight")
print(f"sklon rezidui v mrtvem pasmu: {b1:+.2f} na °C (ma byt ~0)")

# --- 9. nahled pro krok 4: ocistene denni prumery vs. teplota (vsechna data) ---
all_resid = df["baseload"].to_numpy() - base.predict(df, params)
rd_all = pd.DataFrame({"resid": all_resid, "temp": df["temp"], "date": df["date_local"]}) \
    .groupby("date").agg(resid=("resid", "mean"), temp=("temp", "mean"))
fig, ax = plt.subplots(figsize=(8.5, 6))
ax.scatter(rd_all.temp, rd_all.resid, s=12, color=BLUE, alpha=0.5, lw=0)
ax.axhline(0, color=BASELINE, lw=0.8)
ax.set_xlabel("denní průměrná teplota [°C]")
ax.set_ylabel("spotřeba − báze (denní průměr)")
ax.set_title("Po odečtení báze: topný a chladicí signál", loc="left")
fig.savefig(OUT / "09_baze_signal.png", dpi=130, bbox_inches="tight")
