"""Vanoce 2022-2025: jak se lisi od modelu (bez vanocniho bloku) a jak konzistentne.

Rezidua plneho joint fitu (models/*_joint.npz) — model zatim nema vanocni
blok, takze reziduum je primo "vanocni efekt" (plus sum). Obdobi 18. 12. -
8. 1. pro sezony 2022/23 az 2025/26 (+ zacatek ledna 2022).

1. denni odchylka po kalendarnich dnech a sezonach, konzistence mezi roky
2. hodinovy tvar odchylky pro skupiny dnu (Stedry den, svatky, mezi svatky,
   Silvestr, Novy rok, do Tri kralu), pracovni dny vs vikendy
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
BLUE, ORANGE, AQUA, YELLOW, MAGENTA = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "axes.edgecolor": BASELINE, "axes.labelcolor": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "text.color": INK, "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "font.family": "sans-serif",
})
DOW = "Po Út St Čt Pá So Ne".split()

df = etl.load()
bp = base.load(base.MODEL_DIR / "base_joint.npz")
hp = heating.load(heating.MODEL_DIR / "heating_joint.npz")
cp = fit.load("cooling_joint.npz", cooling.PARAM_NAMES + ["rmse", "r2"])
pp = fit.load("pv_joint.npz", ["gamma", "span_days", "knot_days", "rmse", "r2"])
model = base.predict(df, bp) + heating.predict(df, hp) + cooling.predict(df, cp) - pv.predict(df, pp)
loc = df.index.tz_convert(etl.TZ)
d = pd.DataFrame({"y": df["baseload"].to_numpy(), "m": model, "den": pd.to_datetime(df["date_local"]).to_numpy(),
                  "tod": df["tod"].to_numpy(), "dt": df["daytype"].to_numpy()}, index=df.index)
d["e"] = d.y - d.m
d["mes"], d["dd"], d["dow"] = loc.month, loc.day, loc.dayofweek
win = ((d.mes == 12) & (d.dd >= 18)) | ((d.mes == 1) & (d.dd <= 8))
v = d[win].copy()
v["sezona"] = np.where(v.mes == 12, v.den.dt.year, v.den.dt.year - 1)
v["klic"] = v.den.dt.strftime("%d.%m.")

# --- 1. denni odchylka po kalendarnich dnech ---------------------------------
day = v.groupby(["sezona", "den"]).agg(e=("e", "mean"), m=("m", "mean"), dow=("dow", "first"),
                                       klic=("klic", "first")).reset_index()
day["e_pct"] = 100 * day.e / day.m
order = [f"{x:02d}.12." for x in range(18, 32)] + [f"{x:02d}.01." for x in range(1, 9)]
tab = day.pivot_table(index="klic", columns="sezona", values="e_pct").reindex(order)
dows = day.pivot_table(index="klic", columns="sezona", values="dow", aggfunc="first").reindex(order)
pd.set_option("display.width", 200)
print("1) denni odchylka od modelu [% modelu] po kalendarnich dnech (den v tydnu):")
out = tab.round(1).astype(str) + " " + dows.map(lambda x: f"({DOW[int(x)]})" if pd.notna(x) else "")
out["prumer"] = tab.drop(columns=[2021], errors="ignore").mean(axis=1).round(1)
out["rozptyl mezi roky"] = tab.drop(columns=[2021], errors="ignore").std(axis=1).round(1)
print(out.to_string())

fig, ax = plt.subplots(figsize=(12, 4.4))
for s, col in zip([2022, 2023, 2024, 2025], [MUTED, BLUE, AQUA, ORANGE]):
    if s in tab.columns:
        ax.plot(range(len(order)), tab[s], color=col, lw=1.8, marker="o", ms=3, label=f"{s}/{s + 1 - 2000:02d}")
ax.axhline(0, color=BASELINE, lw=1)
ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=90, fontsize=8)
ax.set_ylabel("denní odchylka od modelu [% modelu]")
ax.set_title("Vánoce: odchylka skutečnosti od modelu bez vánočního bloku", loc="left")
ax.legend(frameon=False, ncol=4)
fig.savefig(OUT / "31_vanoce_dny.png", dpi=130, bbox_inches="tight")

# --- 2. hodinovy tvar po skupinach dnu ----------------------------------------
def group(r) -> str:
    md = (r.mes, r.dd)
    if md == (12, 24): return "24.12. Štědrý den"
    if md in ((12, 25), (12, 26)): return "25.-26.12. svátky"
    if md == (12, 31): return "31.12. Silvestr"
    if md == (1, 1): return "1.1. Nový rok"
    if r.mes == 12 and 27 <= r.dd <= 30:
        return "27.-30.12. " + ("pracovní" if r.dow < 5 else "víkend")
    if r.mes == 1 and 2 <= r.dd <= 6:
        return "2.-6.1. " + ("pracovní" if r.dow < 5 else "víkend")
    if r.mes == 12 and 21 <= r.dd <= 23:
        return "21.-23.12. " + ("pracovní" if r.dow < 5 else "víkend")
    return "jiné"


keys = v.groupby("den").agg(mes=("mes", "first"), dd=("dd", "first"), dow=("dow", "first"))
keys["skupina"] = keys.apply(group, axis=1)
v["skupina"] = v.den.map(keys["skupina"])
v["h"] = v.tod.astype(int)
prof = v.groupby(["skupina", "h"]).e.mean().unstack()
n = keys.groupby("skupina").size()
print("\n2) hodinova odchylka od modelu po skupinach (prumer pres vsechny roky), 3h bloky:")
blk = v.groupby(["skupina", v.h // 3 * 3]).e.mean().unstack().round(0)
blk["dni"] = n
blk["denni prumer"] = v.groupby("skupina").e.mean().round(1)
print(blk.to_string())

groups = ["21.-23.12. pracovní", "24.12. Štědrý den", "25.-26.12. svátky", "27.-30.12. pracovní",
          "27.-30.12. víkend", "31.12. Silvestr", "1.1. Nový rok", "2.-6.1. pracovní"]
fig, axes = plt.subplots(2, 4, figsize=(15, 6.5), sharey=True)
for ax, g in zip(axes.flat, groups):
    for s, col in zip([2022, 2023, 2024, 2025], [MUTED, BLUE, AQUA, ORANGE]):
        x = v[(v.skupina == g) & (v.sezona == s)].groupby("tod").e.mean()
        if len(x):
            ax.plot(x.index, x.values, color=col, lw=1.2, label=f"{s}/{s + 1 - 2000:02d}")
    ax.axhline(0, color=BASELINE, lw=1)
    ax.set_title(f"{g} ({n.get(g, 0)} dní)", loc="left", fontsize=9)
    ax.set_xticks(range(0, 25, 6))
axes[0, 0].legend(frameon=False, fontsize=7)
for a in axes[:, 0]:
    a.set_ylabel("skutečnost − model")
fig.suptitle("Vánoce: hodinový průběh odchylky od modelu po skupinách dnů a letech", x=0.01, ha="left")
fig.tight_layout()
fig.savefig(OUT / "32_vanoce_profily.png", dpi=130, bbox_inches="tight")
print("\ngrafy 31-32 ulozeny")
