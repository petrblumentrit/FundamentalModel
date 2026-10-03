"""Simulace intraday predikce za posledni rok (viz src/intraday.py).

Po kazdem 15min mereni se vyda predikce od konce dat do konce dne D+1: model
s cerstvym pocasim (mereni do konce dat, dal predpoved) + korekce zavisla na
horizontu. Faze 0 a 1 (tvar) jako explore/backtest.py, faze 2 pro kazdy den
linearni prefit a 96 vydani. Vysledky:

    simulace/intraday/intraday.npz      predikce a vstupy korekce (vydani x horizont)
    simulace/intraday/horizonty.csv     chyba po horizontech, model a s korekci
    simulace/intraday/koeficienty.csv   koeficienty korekce po horizontech (konec simulace)

`--meteo=predpoved` predpoved pocasi z archivu (jinak namerene pocasi jako
dokonala predpoved), `--out=slozka`, `--workers=N`, `--jen-korekce` prepocita
korekci a tabulky z ulozeneho npz.
"""
import importlib.util
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
import workspace

# sdilene pomucky simulace D+1 (pocet workeru, predpoved pocasi na 15min ose);
# modul se jmenuje stejne jako src/backtest.py, proto pod jinym jmenem
_spec = importlib.util.spec_from_file_location("backtest_d1", Path(__file__).with_name("backtest.py"))
d1 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(d1)

OUT = workspace.ROOT / "simulace" / "intraday"
for _a in sys.argv:
    if _a.startswith("--out="):
        OUT = workspace.ROOT / _a.split("=", 1)[1]
SHOW = (0.25, 0.5, 1, 2, 3, 6, 12, 24, 36)   # [h] horizonty v tabulce

_df = None
_fc = None


def _init(meteo_until=None):
    global _df, _fc
    import etl
    _df = etl.load()
    if meteo_until is not None:
        _fc = d1._forecast_meteo(_df, meteo_until)


def _shape(args):
    import backtest
    day, warm = args
    t = time.time()
    params = backtest.shape_fit(_df, day, warm)
    print(f"tvar k {day.date()} za {time.time() - t:.0f} s", flush=True)
    return day, params


def _days(args):
    import intraday
    days, shapes = args
    parts = []
    for d in days:
        t = time.time()
        sd = max(k for k in shapes if k <= d)
        parts.append(intraday.day_issues(_df, d, shapes[sd], meteo=_fc))
        print(f"{d.date()}  {len(parts[-1]['s'])} vydani za {time.time() - t:.0f} s", flush=True)
    return parts


def simulate():
    import numpy as np
    import pandas as pd

    import backtest
    import etl
    import intraday
    df = etl.load()
    loc = df.index.tz_convert(etl.TZ)
    counts = pd.Series(1, index=loc.date).groupby(level=0).size()
    last_target = pd.Timestamp(counts.index[counts >= 92][-1])
    targets = pd.date_range(last_target - pd.Timedelta(days=364), last_target, freq="D")
    issue_days = list(targets - pd.Timedelta(days=1))
    print(f"data do {loc[-1]}; vydani po 15 min od {issue_days[0].date()} 09:00 "
          f"({len(issue_days)} dni)", flush=True)

    sdays = backtest.shape_days(issue_days)
    t = time.time()
    p0 = backtest.shape_fit(df, sdays[0], None)
    print(f"faze 0 hotova za {time.time() - t:.0f} s", flush=True)
    for v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[v] = str(d1.THREADS)
    nw = d1.n_workers()
    print(f"workeru: {nw} (volna pamet {d1._free_gb():.1f} GB)", flush=True)
    meteo_until = None
    if d1.METEO == "predpoved":
        meteo_until = backtest.cutoff(sdays[0]).tz_convert("UTC")
        print(f"meteo: predpoved z archivu, kalibrace osvitu do {meteo_until}", flush=True)
    with ProcessPoolExecutor(nw, initializer=_init, initargs=(meteo_until,)) as ex:
        shapes = {sdays[0]: p0}
        shapes.update(dict(ex.map(_shape, [(d, p0) for d in reversed(sdays[1:])])))
        print(f"faze 1 ({len(sdays)} tvaru) hotova za {time.time() - t:.0f} s", flush=True)
        chunks = [list(c) for c in np.array_split(np.array(issue_days, dtype=object), 3 * nw)]
        res = list(ex.map(_days, [(c, shapes) for c in chunks]))
    print(f"simulace hotova za {(time.time() - t) / 60:.1f} min", flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT / "intraday.npz", **intraday.stack([p for r in res for p in r]))


def evaluate():
    import numpy as np
    import pandas as pd

    import intraday
    res = dict(np.load(OUT / "intraday.npz"))
    corr, coefs = intraday.online(res)
    tab = intraday.by_horizon(res, corr)
    tab.to_csv(OUT / "horizonty.csv", sep=";", decimal=",", float_format="%.3f")
    # koeficienty z posledniho vydani, ktere dany horizont melo
    last = np.array([coefs[np.flatnonzero(~np.isnan(coefs[:, h, 0]))[-1], h] for h in range(coefs.shape[1])])
    last = pd.DataFrame(last, columns=intraday.available(res), index=tab.index)
    (last.to_csv(OUT / "koeficienty.csv", sep=";", decimal=",", float_format="%.4f"))

    print(f"\nvydani: {len(res['s'])}  ({res['s'][0]} .. {res['s'][-1]} UTC)")
    print("horizont = konec ciloveho intervalu za koncem dat\n")
    print(tab.loc[[h for h in SHOW if h in tab.index],
                  ["n", "rmse_model", "rmse_kor", "mape_model", "mape_kor", "bias_kor"]].round(2).to_string())
    irr = intraday.irregular_targets(res)
    reg = {k: (np.where(irr, np.nan, v) if v.ndim == 2 else v) for k, v in res.items()}
    print("\nbez svatku, mostu a Vanoc v cili:")
    print(intraday.by_horizon(reg, corr).loc[[h for h in SHOW if h in tab.index],
                                              ["rmse_model", "rmse_kor"]].round(2).to_string())
    print("\nkoeficienty korekce na konci simulace:")
    print(last.loc[[h for h in SHOW if h in tab.index]].round(2).to_string())

    # kontrola proti simulaci D+1: konec dat D 09:00, cile = cely den D+1
    s = intraday.local_time(res)
    rows = np.flatnonzero((s.hour == 9) & (s.minute == 0))
    e, ek = [], []
    for i in rows:
        tgt = (s[i] + pd.to_timedelta(np.arange(intraday.HORIZON) * 15, unit="min")).date
        sel = np.asarray(tgt == s[i].date() + pd.Timedelta(days=1))
        d = res["skut"][i, sel].astype(float) - res["pred"][i, sel]
        e.append(d)
        ek.append(d - corr[i, sel])
    e, ek = np.concatenate(e), np.concatenate(ek)
    print(f"\nD+1 z dat do D 09:00: model RMSE {np.sqrt(np.nanmean(e**2)):.2f}, "
          f"s intraday korekci {np.sqrt(np.nanmean(ek**2)):.2f}")


if __name__ == "__main__":
    if "--jen-korekce" not in sys.argv:
        simulate()
    evaluate()
