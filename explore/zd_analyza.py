"""Srovnani portfolia s domacnostmi v PRE (zbytkovy diagram OTE), 2022 - 5/2026.

Zbytkovy diagram (ZD, verze 2) = spotreba v distribucni soustave minus merene
odbery (A, B) — domacnosti a verejne osvetleni; sloupec PRE (sit 0031) je
stejny region jako portfolio. KZD = ZD / ocekavani podle TDD prepoctenych na
skutecnou teplotu — tvar KZD je tak tvar domacnosti ocisteny o pocasi a
typove profily. Data nacita src/ote.py (Analyza/, mimo repo); vse hodinove.

1. Dlouhodoby drift tvaru dne (Po-Ct): podil noci (0-6 h) a odpoledne a
   vecera (15-21 h) na dennim prumeru po letech a sezonach — portfolio
   ocistene o pocasi (model) vs PRE (KZD a surovy ZD).
2. Topna citlivost po zimach: sklon denniho prumeru na teplote (dny -5..10 C)
   pro PRE a portfolio, srovnani s trendem topne citlivosti v modelu.
3. KZD v case: uroven a denni tvar — domacnosti se meni rychleji, nez TDD?
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
import cooling
import etl
import fit
import heating
import ote
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

# --- data: portfolio hodinove + PRE ------------------------------------------
df = etl.load()
hp = heating.load(heating.MODEL_DIR / "heating_joint.npz")
cp = fit.load("cooling_joint.npz", cooling.PARAM_NAMES + ["rmse", "r2"])
pp = fit.load("pv_joint.npz", ["gamma", "span_days", "knot_days", "rmse", "r2"])
df["ocist"] = (df["baseload"] - heating.predict(df, hp) - cooling.predict(df, cp)
               + pv.predict(df, pp))
hr = df[["baseload", "ocist", "temp"]].groupby(df.index.floor("h")).mean()
hr = hr.join(ote.zd_pre().rename("zd"), how="inner").join(ote.kzd_pre().rename("kzd"), how="inner")
cal = df[["daytype", "po_volnu", "date_local"]].groupby(df.index.floor("h")).first()
d = hr.join(cal)
loc = d.index.tz_convert(etl.TZ)
d["h"], d["mes"], d["rok"], d["dd"] = loc.hour, loc.month, loc.year, loc.day
d["den"] = pd.to_datetime(d["date_local"])
print(f"spolecne obdobi: {d.den.min().date()} - {d.den.max().date()} ({d.den.nunique()} dni)")
xmas = ((d.mes == 12) & (d.dd >= 20)) | ((d.mes == 1) & (d.dd <= 6))
reg = (d.daytype == 0) & ~d.po_volnu & ~xmas
SEAS = {12: "zima", 1: "zima", 2: "zima", 3: "jaro", 4: "jaro", 5: "jaro",
        6: "leto", 7: "leto", 8: "leto", 9: "podzim", 10: "podzim", 11: "podzim"}
d["sezona"] = d.mes.map(SEAS)
d["srok"] = np.where(d.mes == 12, d.rok + 1, d.rok)    # zima = prosinec + leden-unor dalsiho roku
for c in ("ocist", "zd", "kzd"):
    d[c + "_rel"] = d[c] / d.groupby("den")[c].transform("mean")

# --- 1. dlouhodoby drift tvaru -------------------------------------------------
def share(c: str) -> pd.DataFrame:
    g = d[reg].groupby(["sezona", "srok"])
    night = g.apply(lambda x: 100 * (x.loc[x.h < 6, c + "_rel"].mean() - 1))
    eve = g.apply(lambda x: 100 * (x.loc[(x.h >= 15) & (x.h < 21), c + "_rel"].mean() - 1))
    n = g["den"].nunique()
    return pd.DataFrame({"noc": night, "vecer": eve, "dni": n})

print("\n1) tvar dne Po-Ct po letech: noc 0-6 h / odpoledne a vecer 15-21 h, % nad/pod dennim prumerem")
S = {c: share(c) for c in ("ocist", "kzd", "zd")}
for s in ("zima", "jaro", "leto", "podzim"):
    rows = []
    for y in sorted(S["ocist"].loc[s].index):
        if S["ocist"].loc[(s, y), "dni"] < 20:
            continue
        rows.append([y] + [f"{S[c].loc[(s, y), 'noc']:+5.1f} / {S[c].loc[(s, y), 'vecer']:+5.1f}" for c in S])
    print(f"   {s}:")
    print(pd.DataFrame(rows, columns=["rok", "portfolio ocist.", "PRE KZD", "PRE ZD surovy"]).to_string(index=False))

fig, axes = plt.subplots(1, 4, figsize=(15, 3.8), sharey=True)
for ax, s in zip(axes, ("zima", "jaro", "leto", "podzim")):
    for c, col, lab in (("ocist", BLUE, "portfolio (očištěno o počasí)"), ("kzd", ORANGE, "PRE: KZD (očištěno o TDD)")):
        t = S[c].loc[s]; t = t[t.dni >= 20]
        ax.plot(t.index, t.noc, color=col, lw=2, marker="o", ms=4, label=f"{lab}, noc")
        ax.plot(t.index, t.vecer, color=col, lw=1.2, ls="--", marker="s", ms=3, label=f"{lab}, 15–21 h")
    ax.set_title(s, loc="left"); ax.set_xticks(sorted(S["ocist"].loc[s].index))
axes[0].set_ylabel("% nad/pod denním průměrem")
axes[0].legend(frameon=False, fontsize=7)
fig.suptitle("Tvar dne Po–Čt po letech: noc a odpoledne/večer (portfolio vs domácnosti PRE)", x=0.01, ha="left")
fig.tight_layout()
fig.savefig(OUT / "29_zd_drift_tvaru.png", dpi=130, bbox_inches="tight")

# --- 2. topna citlivost po zimach ----------------------------------------------
day = d[(d.daytype == 0) & ~xmas].groupby("den").agg(T=("temp", "mean"), y=("baseload", "mean"),
                                                        zd=("zd", "mean"), mes=("mes", "first"))
day["zima"] = np.where(day.mes >= 10, day.index.year, day.index.year - 1)
cold = day[(day["T"] >= -5) & (day["T"] <= 10) & day.mes.isin([11, 12, 1, 2, 3])]
rows = []
for w, g in cold.groupby("zima"):
    if len(g) < 25:
        continue
    rows.append({"zima": f"{w}/{w + 1 - 2000:02d}", "dni": len(g),
                 "portfolio": -np.polyfit(g["T"], g.y, 1)[0], "PRE": -np.polyfit(g["T"], g.zd, 1)[0] / 1000})
w = pd.DataFrame(rows).set_index("zima")
for c in ("portfolio", "PRE"):
    w[c + " index"] = 100 * w[c] / w[c].iloc[0]
days_ = pd.date_range(df.index[0], df.index[-1], freq="D")
k_mean = float(np.mean(heating.k_basis(df) @ hp["k_coef"]))
rel = pd.Series(100 + 100 * heating.trend_curve(days_, hp) / k_mean, index=days_.tz_convert(etl.TZ))
mt = []
for z in w.index:
    y0 = int(z[:4])
    sel = rel[((rel.index.year == y0) & (rel.index.month >= 11)) | ((rel.index.year == y0 + 1) & (rel.index.month <= 3))]
    mt.append(sel.mean() if len(sel) else np.nan)
w["model trend index"] = 100 * np.array(mt) / mt[0]
print("\n2) topna citlivost (Po-Ct, dny -5..10 C, XI-III): sklon na 1 C; index = 100 v prvni zime")
print("   portfolio v jednotkach spotreby, PRE v MWh/h")
print(w.round(2).to_string())

# --- 3. KZD v case --------------------------------------------------------------
k = d[reg].groupby(["rok", "mes"])["kzd"].mean().unstack()
print("\n3) KZD PRE — mesicni prumer (Po-Ct):")
print(k.round(3).to_string())

fig, ax = plt.subplots(figsize=(10, 4))
x = w.index
ax.plot(x, w["portfolio index"], color=BLUE, lw=2, marker="o", label="portfolio (sklon z dat)")
ax.plot(x, w["PRE index"], color=ORANGE, lw=2, marker="o", label="domácnosti PRE (sklon z dat)")
ax.plot(x, w["model trend index"], color=MUTED, lw=1.5, ls="--", marker="s", label="trend topné citlivosti v modelu")
ax.set_ylabel("index topné citlivosti (první zima = 100)")
ax.set_title("Topná citlivost po zimách", loc="left")
ax.legend(frameon=False)
fig.savefig(OUT / "30_zd_topna_citlivost.png", dpi=130, bbox_inches="tight")
print("\ngrafy 29-30 ulozeny")
