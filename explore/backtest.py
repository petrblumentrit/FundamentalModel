"""Simulace provozni predikce D+1 za posledni rok platnych dat (viz src/backtest.py).

V den D v 10:00 se z dat do D 09:00 prefituje model a predikuje cely den D+1.
Dny vydani se deli do souvislych bloku, ktere bezi paralelne (v bloku warm
start ze dne na den). Vysledky:

    simulace/predikce_D1.csv   15min predikce D+1 se slozkami a skutecnosti,
                               korekce z chyb drivejsich predikci (src/correction.py)
    simulace/parametry.csv     parametry modelu po dnech vydani
    simulace/predikce_D1.html  interaktivni offline graf (zoom), otevre se v prohlizeci

`--plot-only` jen prepocita korekci a prekresli HTML z ulozeneho CSV.
"""
import json
import os
import sys
import time
import webbrowser
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
OUT = ROOT / "simulace"
N_WORKERS = 10
THREADS = 3

_df = None


def _init():
    global _df
    import etl
    _df = etl.load()


def _work(days):
    import backtest
    t = time.time()
    pred, pars = backtest.run_days(_df, days, log=lambda s: print(s, flush=True))
    print(f"blok {days[0].date()}..{days[-1].date()} hotov za {time.time() - t:.0f} s",
          flush=True)
    return pred, pars


def simulate():
    import numpy as np
    import pandas as pd

    import etl
    df = etl.load()
    loc = df.index.tz_convert(etl.TZ)
    counts = pd.Series(1, index=loc.date).groupby(level=0).size()
    full_days = counts.index[counts >= 92]           # DST dny maji 92/100
    last_target = pd.Timestamp(full_days[-1])
    targets = pd.date_range(last_target - pd.Timedelta(days=364), last_target, freq="D")
    issue_days = list(targets - pd.Timedelta(days=1))
    print(f"data do {loc[-1]}; predikce D+1 pro {targets[0].date()} .. {targets[-1].date()} "
          f"({len(targets)} dni)", flush=True)

    chunks = [list(c) for c in np.array_split(np.array(issue_days, dtype=object), N_WORKERS)]
    for v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[v] = str(THREADS)
    t = time.time()
    with ProcessPoolExecutor(N_WORKERS, initializer=_init) as ex:
        res = list(ex.map(_work, chunks))
    print(f"simulace hotova za {(time.time() - t) / 60:.1f} min", flush=True)

    pred = pd.concat([r[0] for r in res]).sort_index()
    pars = pd.concat([r[1] for r in res]).sort_index()
    OUT.mkdir(exist_ok=True)
    pred.index.name = "timestamp_utc"
    pred.insert(0, "timestamp_mistni", pred.index.tz_convert(etl.TZ))
    pred.to_csv(OUT / "predikce_D1.csv", sep=";", decimal=",", float_format="%.3f")
    pars.to_csv(OUT / "parametry.csv", sep=";", decimal=",", float_format="%.5f")


def correct():
    """Online korekce z chyb drivejsich predikci; sloupce se pridaji do CSV."""
    import pandas as pd

    import correction
    path = OUT / "predikce_D1.csv"
    pred = pd.read_csv(path, sep=";", decimal=",", index_col="timestamp_utc")
    pred["korekce"], coefs = correction.online(pred)
    pred["predikce_kor"] = pred["predikce"] + pred["korekce"]
    pred.to_csv(path, sep=";", decimal=",", float_format="%.3f")
    coefs.to_csv(OUT / "korekce_koeficienty.csv", sep=";", decimal=",", float_format="%.4f")


def load_results():
    import pandas as pd
    pred = pd.read_csv(OUT / "predikce_D1.csv", sep=";", decimal=",")
    pred["timestamp_mistni"] = pd.to_datetime(pred["timestamp_mistni"], utc=True).dt.tz_convert("Europe/Prague")
    return pred


def scores(e, y):
    import numpy as np
    return {"rmse": float(np.sqrt(np.mean(e**2))), "mae": float(np.mean(np.abs(e))),
            "bias": float(np.mean(e)), "mape": float(np.mean(np.abs(e) / y) * 100),
            "r2": float(1 - np.var(e) / np.var(y))}


def plot():
    import pandas as pd

    pred = load_results()
    y, f, k = (pred[c].to_numpy() for c in ("skutecnost", "predikce", "predikce_kor"))
    e, ek = y - f, y - k
    ts = pred["timestamp_mistni"]
    total = {"model": scores(e, y), "kor": scores(ek, y)}
    month = ts.dt.strftime("%Y-%m").to_numpy()
    monthly = [{"mesic": m, "model": scores(e[month == m], y[month == m]),
                "kor": scores(ek[month == m], y[month == m])} for m in sorted(set(month))]
    day = ts.dt.date
    daily = (pd.DataFrame({"d": day, "e": ek})
             .groupby("d")["e"].apply(lambda x: float((x**2).mean() ** 0.5)))
    worst = [{"den": str(d), "rmse": r} for d, r in daily.nlargest(8).items()]

    data = {
        "t": ts.dt.strftime("%Y-%m-%d %H:%M").tolist(),
        "y": [round(v, 1) for v in y],
        "f": [round(v, 1) for v in f],
        "k": [round(v, 1) for v in k],
        "e": [round(v, 1) for v in e],
        "ek": [round(v, 1) for v in ek],
        "vydano": pd.to_datetime(pred["vydano"], utc=True).dt.tz_convert("Europe/Prague")
                    .dt.strftime("%d.%m. %H:%M").tolist(),
        "total": total, "monthly": monthly, "worst": worst,
        "od": str(day.iloc[0]), "do": str(day.iloc[-1]), "dni": int(day.nunique()),
    }
    from plotly.offline import get_plotlyjs
    # plotly.js primo v souboru — graf funguje offline, bez CDN
    html = (Path(__file__).with_name("backtest_template.html").read_text(encoding="utf-8")
            .replace("/*__DATA__*/null", json.dumps(data, separators=(",", ":")))
            .replace("/*__PLOTLY__*/", get_plotlyjs()))
    out = OUT / "predikce_D1.html"
    out.write_text(html, encoding="utf-8")
    print(f"graf: {out}")
    return out


if __name__ == "__main__":
    if "--plot-only" not in sys.argv:
        simulate()
    correct()
    out = plot()
    if "--no-open" not in sys.argv:
        webbrowser.open(out.as_uri())
