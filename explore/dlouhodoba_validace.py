"""Dopredna validace dlouhodobe predikce (src/longterm.py) — rok dopredu.

Pro kazdy fold se model odhadne jen z dat pred cutoffem (cely retezec, hodinova
agregace) a predikuje nasledujicich 365 dni. Dve otazky:

1. pravidlo pro pomale slozky — `drzet` vs `posledni` trend — se skutecnym
   pocasim (oddeli chybu trendu od nejistoty pocasi);
2. rezim pocasi — `normal` vs `scenare` (jen z let pred cutoffem) — tedy
   realna chyba predikce na rok dopredu vcetne nejistoty pocasi.

Fity foldu se ukladaji do models/dlouhodoba_validace.pkl (mimo git);
`--refit` je prepocita.
"""
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import base
import etl
import longterm
import validate

CUTOFFS = ["2023-08-05", "2024-08-05", "2025-08-05"]   # mistni pulnoc; data od 1/2022 do 4. 8. 2026
DAYS = 365
CACHE = base.MODEL_DIR / "dlouhodoba_validace.pkl"
VARIANTS = [("skutecne", "drzet"), ("skutecne", "posledni"), ("skutecne", "uroven"), ("skutecne", "topeni,fve"),
            ("skutecne", "posledni:0.5"), ("normal", "posledni"), ("scenare", "posledni"),
            ("normal", "topeni,fve"), ("scenare", "topeni,fve")]

pd.set_option("display.width", 220)
df = etl.load()
dh = etl.hourly(df)
y = dh["baseload"]


def utc(day: str) -> pd.Timestamp:
    return pd.Timestamp(day, tz=etl.TZ).tz_convert("UTC")


fits = {}
if CACHE.exists() and "--refit" not in sys.argv:
    fits = pickle.loads(CACHE.read_bytes())
for c in CUTOFFS:
    if c not in fits:
        t = time.time()
        fits[c] = validate.fit_masked(dh, (dh.index < utc(c)))
        print(f"fit do {c} za {time.time() - t:.0f} s", flush=True)
        CACHE.write_bytes(pickle.dumps(fits))


def score(pred: pd.Series) -> dict:
    j = pd.concat([y.rename("y"), pred.rename("p")], axis=1, join="inner").dropna()
    e = j.y - j.p
    mon = e.groupby(j.index.tz_convert(etl.TZ).strftime("%Y-%m")).mean()
    return {"RMSE": float(np.sqrt(np.mean(e**2))), "MAPE %": float(np.mean(np.abs(e) / j.y) * 100),
            "bias": float(e.mean()), "bias %": float(e.mean() / j.y.mean() * 100),
            "mesicni bias RMS": float(np.sqrt(np.mean(mon**2))),
            "mes. objem MAPE %": float(np.mean(np.abs(mon) / j.y.groupby(j.index.tz_convert(etl.TZ).strftime("%Y-%m")).mean()) * 100)}


rows, monthly = [], {}
for c in CUTOFFS:
    cut = utc(c)
    for weather, trend in VARIANTS:
        f = longterm.forecast(df, fits[c], cut, DAYS, weather=weather, trend=trend, hourly=True)
        s = score(f["predikce"])
        extra = {}
        if weather == "scenare":
            sc = longterm.scenario_columns(f)
            extra = {"scenaru": len(sc)}
            # pokryti pasma: jen pocasi vs. pocasi + chyba modelu; a pro mesicni prumery
            for label, sd in (("jen pocasi", 0.0), ("s chybou modelu", longterm.MODEL_SD)):
                lo, hi = longterm.band(f[sc].to_numpy(), sd)
                j = pd.DataFrame({"y": y.reindex(f.index), "lo": lo, "hi": hi}).dropna()
                extra[f"v pasmu % ({label})"] = float(((j.y >= j.lo) & (j.y <= j.hi)).mean() * 100)
                extra[f"sirka ({label})"] = float((j.hi - j.lo).mean())
            m = longterm.monthly(f)
            ym = y.reindex(f.index).groupby(f.index.tz_convert(etl.TZ).strftime("%Y-%m")).mean()
            extra["mesicu v pasmu"] = f"{int(((ym >= m.pasmo_dolni) & (ym <= m.pasmo_horni)).sum())}/{len(m)}"
        rows.append({"cutoff": c, "pocasi": weather, "trend": trend, **s, **extra})
        e = (y.reindex(f.index) - f["predikce"]).dropna()
        monthly[(c, weather, trend)] = e.groupby(e.index.tz_convert(etl.TZ).month).mean()

res = pd.DataFrame(rows)
print(res.round(2).to_string(index=False))
print("\nprumer pres foldy:")
print(res.groupby(["pocasi", "trend"], sort=False)[["RMSE", "MAPE %", "bias", "bias %", "mesicni bias RMS", "mes. objem MAPE %"]]
      .agg(lambda v: np.mean(np.abs(v)) if v.name.startswith("bias") else np.mean(v)).round(2).to_string())
print("\nbias po mesicich (skutecnost - predikce):")
print(pd.DataFrame(monthly).T.round(1).to_string())
