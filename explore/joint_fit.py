"""Zaverecny joint fit + ocistena spotreba (krok 6)."""
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
import normal
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
bp0, hp0 = base.load(), heating.load()
cp0 = fit.load("cooling_params.npz", cooling.PARAM_NAMES + ["rmse", "r2"])
pp0 = fit.load("pv_params.npz", ["gamma", "span_days", "knot_days", "rmse", "r2"])

staged = (df["baseload"].to_numpy() - base.predict(df, bp0) - heating.predict(df, hp0)
          - cooling.predict(df, cp0) + pv.predict(df, pp0))
print(f"pred joint fitem: RMSE={np.sqrt(np.mean(staged**2)):.2f}")

bp, hp, cp, pp = fit.joint(df, bp0, hp0, cp0, pp0)
for prm, name in [(bp, "base_joint.npz"), (hp, "heating_joint.npz"),
                  (cp, "cooling_joint.npz"), (pp, "pv_joint.npz")]:
    fit.save(prm, name)
print(f"po joint fitu:    RMSE={hp['rmse']:.2f}  R2={hp['r2']:.4f}\n")

print(f"{'parametr':<16}{'staged':>10}{'joint':>10}")
for k in heating.PARAM_NAMES:
    print(f"  topeni {k:<8}{hp0[k]:>10.4f}{hp[k]:>10.4f}")
for k in cooling.PARAM_NAMES:
    print(f"  chlaz. {k:<8}{cp0[k]:>10.4f}{cp[k]:>10.4f}")
print(f"  FVE    gamma   {pp0['gamma']:>10.4f}{pp['gamma']:>10.4f}")
cap0 = pv.capacity_curve(pd.DatetimeIndex([df.index[-1]]), pp0)[0]
cap1 = pv.capacity_curve(pd.DatetimeIndex([df.index[-1]]), pp)[0]
print(f"  FVE    kapacita{cap0:>10.4f}{cap1:>10.4f}  (spicka {cap0*845:.0f} -> {cap1*845:.0f})")

df["p_base"] = base.predict(df, bp)
df["p_heat"] = heating.predict(df, hp)
df["p_cool"] = cooling.predict(df, cp)
df["g_pv"] = pv.predict(df, pp)
df["model"] = df.p_base + df.p_heat + df.p_cool - df.g_pv
df["resid"] = df.baseload - df.model

rd = df.groupby("date_local").agg(r=("resid", "mean"), t=("temp", "mean"), s=("sun", "mean"))
b_sun = np.polyfit(rd.s, rd.r, 1)[0]
b_t = np.polyfit(rd.t, rd.r, 1)[0]
print(f"\nsklon rezidui: osvit {b_sun:+.4f}/(W/m2) (bylo -0.017), teplota {b_t:+.3f}/C")
mon = df.groupby(df.index.tz_convert(etl.TZ).month)["resid"].apply(lambda x: np.sqrt((x**2).mean()))
print("RMSE po mesicich:", {int(k): round(v, 1) for k, v in mon.items()})

# --- ocistena spotreba: dve definice z navrhu ---
dfn = normal.weather(df)
p_heat_n = heating.predict(dfn, hp)
p_cool_n = cooling.predict(dfn, cp)
g_pv_n = pv.predict(dfn, pp)
df["y_resid"] = df.baseload - df.p_heat - df.p_cool + df.g_pv
df["y_norm"] = df.baseload - (df.p_heat + df.p_cool - df.g_pv) + (p_heat_n + p_cool_n - g_pv_n)
yr = df.groupby(df.index.tz_convert(etl.TZ).year)[["baseload", "y_resid", "y_norm"]].mean()
print("\nrocni prumery (mereno / rezidualne ocisteno / na normal):")
print(yr.round(1).to_string())

# --- 19. uroven baze a kapacita FVE: staged vs. joint ---
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.5, 4.4))
days = pd.date_range(df.index[0], df.index[-1], freq="D")
a1.plot(days, base.level_curve(days, bp0), color=MUTED, lw=1.6, ls="--", label="staged")
a1.plot(days, base.level_curve(days, bp), color=BLUE, lw=2, label="joint")
a1.set_ylabel("bazální úroveň")
a1.set_title("Úroveň báze: staged vs. joint", loc="left")
a1.legend(frameon=False)
a2.plot(days, pv.capacity_curve(days, pp0) * 845, color=MUTED, lw=1.6, ls="--", label="staged")
a2.plot(days, pv.capacity_curve(days, pp) * 845, color=BLUE, lw=2, label="joint")
a2.set_ylabel("výkon FVE při 845 W/m²")
a2.set_title("Kapacita FVE: staged vs. joint", loc="left")
a2.legend(frameon=False)
fig.savefig(OUT / "19_joint_srovnani.png", dpi=130, bbox_inches="tight")

# --- 20. mesicni rozklad prispevku ---
m = df.groupby(df.index.tz_convert(etl.TZ).month)[["p_heat", "p_cool", "g_pv"]].mean()
fig, ax = plt.subplots(figsize=(9, 4.4))
x = m.index.to_numpy()
ax.bar(x - 0.26, m.p_heat, 0.25, color=ORANGE, label="topení")
ax.bar(x, m.p_cool, 0.25, color=BLUE, label="chlazení")
ax.bar(x + 0.26, -m.g_pv, 0.25, color=YELLOW, label="FVE (odečítá se)")
ax.axhline(0, color=BASELINE, lw=0.8)
ax.set_xlabel("měsíc"); ax.set_ylabel("průměrný příspěvek")
ax.set_xticks(x)
ax.set_title("Měsíční rozklad povětrnostních příspěvků", loc="left")
ax.legend(frameon=False)
fig.savefig(OUT / "20_rozklad.png", dpi=130, bbox_inches="tight")

# --- 21. ocistena spotreba: obe definice z navrhu ---
fig, ax = plt.subplots(figsize=(11.5, 4.6))
d = df.groupby("date_local")[["baseload", "y_norm", "y_resid"]].mean()
d.index = pd.to_datetime(d.index)
ax.plot(d.index, d.baseload, color=GRID, lw=1.0, label="měřená spotřeba")
ax.plot(d.index, d.y_norm.rolling(7, center=True).mean(), color=ORANGE, lw=1.8,
        label="normálové podmínky (sezónnost zůstává)")
ax.plot(d.index, d.y_resid.rolling(7, center=True).mean(), color=BLUE, lw=2,
        label="reziduální (bez povětrnostních členů)")
ax.set_ylabel("denní průměrná spotřeba (7denní průměr)")
ax.set_title("Očištěná spotřeba — dvě definice", loc="left")
ax.legend(frameon=False, ncol=3)
fig.savefig(OUT / "21_ocistena_spotreba.png", dpi=130, bbox_inches="tight")
print("\ngrafy 19-21 ulozeny")
