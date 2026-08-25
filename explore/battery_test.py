"""Bateriovy detekcni test — krok 5 navrhu (NAVRH_MODELU.md).

Baterie se pozna podle mezicasove vazby: vecerni a ponocni spotreba zavisi na
osvitu z **predchozich** hodin. Test je proto cileny na hodiny, kdy FVE uz
nevyrabi (18-05 lokalne) — cokoliv, co tam osvit vysvetluje, neni okamzita
vyroba, ale presunuta energie.

1) volna impulzni odezva vecernich/nocnich rezidui na anomalii osvitu,
   lagy 0-36 h (pres pulnoc) s penalizaci hladkosti,
2) rocni vyvoj vecerni citlivosti — rostouci (zapornejsi) sklon kopirujici
   instalace je dukaz, ze efekt je merittelny a sili.
"""
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
cp = fit.load("cooling_params.npz", cooling.PARAM_NAMES + ["rmse", "r2"])
pp = fit.load("pv_params.npz", ["gamma", "span_days", "knot_days", "rmse", "r2"])
df["resid"] = (df["baseload"] - base.predict(df, bp) - heating.predict(df, hp)
               - cooling.predict(df, cp) + pv.predict(df, pp))

# --- anomalie osvitu: odchylka od 15denniho klouzaveho normalu ve stejnou hodinu ---
loc = df.index.tz_convert(etl.TZ)
h = pd.DataFrame({"r": df["resid"].to_numpy(), "sun": df["sun"].to_numpy()},
                 index=df.index).resample("1h").mean()
hloc = h.index.tz_convert(etl.TZ)
piv = h["sun"].groupby([hloc.date, hloc.hour]).mean().unstack()
norm = piv.rolling(15, center=True, min_periods=5).mean().stack()
# napojeni pres klic (datum, hodina) — rekonstrukce lokalnich casu by narazila na DST
key = pd.MultiIndex.from_arrays([hloc.date, hloc.hour])
h["sun_anom"] = h["sun"].to_numpy() - norm.reindex(key).to_numpy()

# rezidua taky jako anomalie (odstrani pomalou driftovou slozku)
h["r_anom"] = h["r"] - h["r"].rolling(15 * 24, center=True, min_periods=100).mean()

# --- 1) impulzni odezva vecer/v noci na lagovany osvit ---
LAGS = np.arange(0, 37)
X = np.column_stack([h["sun_anom"].shift(l).to_numpy() for l in LAGS])
hour_loc = h.index.tz_convert(etl.TZ).hour
D2 = np.diff(np.eye(len(LAGS)), n=2, axis=0)


def irf(y: np.ndarray, hours: np.ndarray) -> np.ndarray:
    ok = hours & ~np.isnan(X).any(axis=1) & ~np.isnan(y)
    Xn, yn = X[ok], y[ok]
    lam = 2e3 * np.trace(Xn.T @ Xn) / len(LAGS)
    return np.linalg.solve(Xn.T @ Xn + lam * D2.T @ D2, Xn.T @ yn), int(ok.sum())


night = (hour_loc >= 18) | (hour_loc <= 5)
beta, n_night = irf(h["r_anom"].to_numpy(), night)
print(f"vecerni/nocni hodiny v testu: {n_night}")
print(f"odezva v lagu 0: {beta[0]:+.4f}, minimum {beta.min():+.4f} v lagu {LAGS[beta.argmin()]} h")
print(f"soucet odezvy pres lagy 1-36: {beta[1:].sum():+.3f}")

# pozitivni kontrola: stejny odhad na dennich hodinach a reziduu BEZ odectene FVE
# musi najit znamy efekt velikosti kapacity FVE -> test ma silu
ctrl = df["resid"] - pv.predict(df, pp)
ctrl = ctrl.resample("1h").mean()
ctrl = (ctrl - ctrl.rolling(15 * 24, center=True, min_periods=100).mean()).to_numpy()
beta_c, n_day = irf(ctrl, (hour_loc >= 9) & (hour_loc <= 15))
print(f"kontrola (den, bez odectene FVE): lag 0 = {beta_c[0]:+.4f}"
      f" (ocekavana kapacita FVE {-pv.capacity_curve(pd.DatetimeIndex([df.index[-1]]), pp)[0]:+.4f})")

fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.5, 4.4))
a1.plot(LAGS, beta_c, color=MUTED, lw=1.6, ls="--",
        label="kontrola: den, FVE neodečtena")
a1.plot(LAGS, beta, color=BLUE, lw=2, label="test: večer/noc (18–05)")
a1.axhline(0, color=BASELINE, lw=0.8)
a1.set_xlabel("zpoždění osvitu [h]"); a1.set_ylabel("odezva rezidua [na W/m²]")
a1.set_xticks([0, 6, 12, 18, 24, 30, 36])
a1.set_title("Odezva na dřívější osvit", loc="left")
a1.legend(frameon=False)

# --- 2) rocni vyvoj vecerni citlivosti ---
dloc = pd.Series(hour_loc, index=h.index)
day_key = pd.Series(h.index.tz_convert(etl.TZ).normalize(), index=h.index)
# solarni den: vecer 18-23 se pari s osvitem tehoz dne 9-15
ev = h.loc[(hour_loc >= 18) & (hour_loc <= 23), "r_anom"].groupby(
    day_key[(hour_loc >= 18) & (hour_loc <= 23)]).mean()
sa = h.loc[(hour_loc >= 9) & (hour_loc <= 15), "sun_anom"].groupby(
    day_key[(hour_loc >= 9) & (hour_loc <= 15)]).mean()
tab = pd.concat([ev.rename("ev"), sa.rename("sun")], axis=1).dropna()
years, slopes, errs = [], [], []
for y, d in tab.groupby(tab.index.year):
    if len(d) < 100:
        continue
    b1, b0 = np.polyfit(d["sun"], d["ev"], 1)
    pred = b0 + b1 * d["sun"]
    se = np.sqrt(((d["ev"] - pred) ** 2).sum() / (len(d) - 2)
                 / ((d["sun"] - d["sun"].mean()) ** 2).sum())
    years.append(y); slopes.append(b1); errs.append(se)
    print(f"  {y}: sklon vecerni citlivosti {b1:+.4f} +- {se:.4f} (n={len(d)})")

a2.errorbar(years, slopes, yerr=1.96 * np.array(errs), color=BLUE, lw=2,
            marker="o", ms=7, capsize=4)
a2.axhline(0, color=BASELINE, lw=0.8)
a2.set_xlabel("rok"); a2.set_ylabel("sklon: večer vs. denní osvit")
a2.set_xticks(years)
a2.set_title("Roční vývoj večerní citlivosti (± 95 %)", loc="left")
fig.savefig(OUT / "18_baterie_test.png", dpi=130, bbox_inches="tight")
print("graf 18 ulozen")
