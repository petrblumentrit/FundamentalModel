"""Teplotni koeficient FVE (gamma) na horni mezi 0,01 — co prebira?

    uv run python explore/fve_teplota_test.py

Fyzikalne ma vykon panelu klesat o ~0,4 %/C teploty clanku; model ma gamma na
teplotu vzduchu a fit ze vsech dat ji tlaci na mez 1 %/C. Test:

1. profil: plny fit (tvar + linearni cast, hodinova data) s gamma drzenou na
   rade hodnot a s uvolnenou horni mezi — jak se meni chyba a ostatni parametry,
2. rezidua za dne podle teploty a osvitu pro fyzikalni gamma = 0,004 a pro
   volnou gamma — kde fyzikalni hodnota nesedi (horko = chlazeni od slunce,
   zima = sklon panelu / nizke slunce).

Warm start z provozniho stavu (provoz/parametry.pkl).
"""
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
import backtest
import config
import etl
import explain
import pv
import validate

GRID = (0.0, 0.002, 0.004, 0.006, 0.01, 0.015, 0.02, 0.03)
T_BINS = [-30, 0, 10, 18, 24, 28, 50]
S_BINS = [100, 300, 600, 2000]

dh = etl.hourly(etl.load())
mask = np.ones(len(dh), bool)
warm = pd.read_pickle(ROOT / "provoz" / "parametry.pkl")["params"]
spec = config.model()["fve"]["parametry"]["gamma"]
y = dh["baseload"].to_numpy()
end = pd.DatetimeIndex([dh.index[-1]])


def fit(gamma: float | None):
    spec["pevna"] = gamma
    t = time.time()
    p = backtest.refit(dh, mask, warm)
    res = y - validate.predict(dh, p)
    row = {"gamma": p[3]["gamma"], "rmse": np.sqrt(np.mean(res**2)),
           "rmse_den": np.sqrt(np.mean(res[dh["sun"].to_numpy() > 100]**2)),
           "fve_spicka": pv.capacity_curve(end, p[3])[0] * explain.PV_REF_SUN,
           "a": p[1]["a"], "t_b": p[1]["t_b"], "tau": p[1]["tau"],
           "a_c": p[2]["a_c"], "t_bc": p[2]["t_bc"], "alpha_c": p[2]["alpha_c"], "s": time.time() - t}
    return row, res


def table(res: np.ndarray) -> pd.DataFrame:
    """Prumerne reziduum (skutecnost - model) podle teploty a osvitu, hodiny s osvitem > 100."""
    d = pd.DataFrame({"r": res, "T": pd.cut(dh["temp"], T_BINS), "I": pd.cut(dh["sun"], S_BINS)}).dropna()
    out = d.pivot_table(index="T", columns="I", values="r", aggfunc="mean", observed=False)
    out["hodin"] = d.groupby("T", observed=False).size()
    return out


rows, resid = [], {}
for g in GRID:
    r, res = fit(g)
    rows.append(r)
    resid[g] = res
    print(f"gamma {g:.3f}: RMSE {r['rmse']:.3f} (den {r['rmse_den']:.3f}), spicka FVE {r['fve_spicka']:.1f}, "
          f"a_c {r['a_c']:.4f}, t_bc {r['t_bc']:.2f}, {r['s']:.0f} s", flush=True)
pv.UPPER[0] = 0.05        # uvolnena horni mez
r, res = fit(None)
rows.append(r)
resid["volna"] = res
print(f"volna gamma (mez 0,05): {r['gamma']:.4f}, RMSE {r['rmse']:.3f}", flush=True)

pd.set_option("display.width", 200)
print("\nprofil pres gamma (posledni radek = volna, horni mez 0,05):")
print(pd.DataFrame(rows).round(4).to_string(index=False))
for key, label in ((0.004, "gamma = 0,004 (fyzikalni)"), ("volna", "volna gamma")):
    print(f"\nprumerne reziduum za dne, {label}; radky = teplota [C], sloupce = osvit [W/m2]:")
    print(table(resid[key]).round(1).to_string())
