"""Test: opozdene zatopeni po teplem obdobi (pozorovani z praxe).

Hypoteza (uzivatel): kdyz se po teplem vikendu nebo tydnu ochladi, topeni v
administrativnich budovach jeste nebezi — zatopi se az behem dne nebo dalsi
den. Model (dvoukanalovy filtr, pomaly kanal tau ~82 h) by pak prvni chladny
den spotrebu nadhodnocoval: reziduum (skutecnost - model) zaporne.

Dva pohledy na rezidua plneho joint fitu (models/*_joint.npz), mimo Vanoce:

A. udalosti — pracovni den s prumernou teplotou <= Tc po k dnech s prumerem
   >= Tw; reziduum v pracovni dobe (7-17 h) a v noci (0-5 h), totez den pote;
B. stav ustredniho topeni podle pravidla vyhlasky 194/2007 (zahajeni pri
   prumerne denni teplote < 13 C dva dny po sobe, preruseni pri > 13 C dva dny
   po sobe, otopne obdobi 1. 9. - 31. 5.) — rezidua chladnych dni ve stavu
   "nebezi" proti "bezi".
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import base
import cooling
import etl
import fit
import heating
import pv

RULE_TEMP = 13.0          # [C] prumerna denni teplota pro zahajeni/preruseni vytapeni
WORK = (7, 17)            # pracovni doba [h]
NIGHT = (0, 5)

pd.set_option("display.width", 200)

df = etl.load()
bp = base.load(base.MODEL_DIR / "base_joint.npz")
hp = heating.load(heating.MODEL_DIR / "heating_joint.npz")
cp = fit.load("cooling_joint.npz", cooling.PARAM_NAMES + ["rmse", "r2"])
pp = fit.load("pv_joint.npz", ["gamma", "span_days", "knot_days", "rmse", "r2"])
df["p_heat"] = heating.predict(df, hp)
df["resid"] = (df["baseload"] - base.predict(df, bp) - df["p_heat"]
               - cooling.predict(df, cp) + pv.predict(df, pp))
df["den"] = pd.to_datetime(df["date_local"])


def se(s: pd.Series) -> float:
    return float(s.std() / np.sqrt(len(s))) if len(s) > 1 else np.nan


work = df[(df.tod >= WORK[0]) & (df.tod < WORK[1])]
d = df.groupby("den").agg(T=("temp", "mean"), topeni=("p_heat", "mean"), rez=("resid", "mean"),
                          dt=("daytype", "first"), most=("most", "first"))
d["rez_prace"] = work.groupby("den")["resid"].mean()
d["rez_noc"] = df[(df.tod >= NIGHT[0]) & (df.tod < NIGHT[1])].groupby("den")["resid"].mean()
d["prac"] = (d.dt <= 1) & ~d.most
xmas = ((d.index.month == 12) & (d.index.day >= 18)) | ((d.index.month == 1) & (d.index.day <= 6))

# stav ustredniho topeni: rozhodnuti z predchozich dvou dni
on = np.zeros(len(d), bool)
state = True
T = d["T"].to_numpy()
for i, day in enumerate(d.index):
    if day.month in (6, 7, 8):
        state = False
    elif i >= 2:
        if T[i - 1] < RULE_TEMP and T[i - 2] < RULE_TEMP:
            state = True
        elif T[i - 1] > RULE_TEMP and T[i - 2] > RULE_TEMP:
            state = False
    on[i] = state
d["bezi"] = on
d = d[~xmas]

print("A. prvni chladny pracovni den po teplem obdobi (reziduum = skutecnost - model)")
rows = []
for k, tw, tc in [(2, 15, 12), (3, 15, 12), (3, 14, 11), (3, 16, 13), (5, 15, 12)]:
    warm = d["T"].shift(1).rolling(k).mean()
    ev = d[(warm >= tw) & (d["T"] <= tc) & d.prac]
    nxt = d[d.index.isin(ev.index + pd.Timedelta(days=1)) & d.prac]
    rows.append({"dni tepla": k, "T pred >=": tw, "T dnes <=": tc, "n": len(ev),
                 "prace": ev.rez_prace.mean(), "+-": se(ev.rez_prace), "noc": ev.rez_noc.mean(),
                 "topeni": ev.topeni.mean(), "den pote n": len(nxt), "den pote prace": nxt.rez_prace.mean()})
print(pd.DataFrame(rows).round(1).to_string(index=False))

print("\nB. chladne dny (topna slozka 15-120) podle stavu ustredniho topeni")
c = d[d.topeni.between(15, 120)]
print(c.groupby("bezi").agg(n=("rez", "size"), rez=("rez", "mean"), se=("rez", se),
                            topeni=("topeni", "mean"), T=("T", "mean")).round(2).to_string())
for name, sel in (("pracovni", c.prac), ("volne", ~c.prac)):
    s = c[sel]
    print(f"  {name}: " + ", ".join(
        f"{'bezi' if st else 'nebezi'} {g.rez.mean():+.1f} +- {se(g.rez):.1f} (n={len(g)})"
        for st, g in s.groupby("bezi")))
aut = c[c.index.month.isin([9, 10, 11]) & c.prac]
print("  podzim, pracovni doba: " + ", ".join(
    f"{'bezi' if st else 'nebezi'} {g.rez_prace.mean():+.1f} +- {se(g.rez_prace):.1f} (n={len(g)})"
    for st, g in aut.groupby("bezi")))
