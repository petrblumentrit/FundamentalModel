"""Dokumentace modelu v HTML: hlavni myslenka, slozky s obrazky, presnost, navod.

    uv run python explore/dokumentace.py

Obrazky se kresli z aktualniho stavu modelu (provoz/parametry.pkl) a z
vysledku zpetnych simulaci, pokud existuji (simulace/v12*, simulace/intraday*).
Vystup dokumentace/dokumentace.html (mimo repo — obsahuje data portfolia),
offline: plotly.js je vlozene v souboru. `--no-open` neotvira prohlizec.
"""
import json
import sys
import webbrowser
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
import workspace
import etl
import explain
import heating
import operation

OUT = workspace.ROOT / "dokumentace"
BACKTESTS = (("Naměřené počasí (dokonalá předpověď)", "v12"), ("Předpověď počasí z archivu", "v12_predpoved"))
INTRADAY = (("naměřené počasí", "intraday"), ("předpověď počasí", "intraday_predpoved"))

df = etl.load()
state = pd.read_pickle(workspace.ROOT / "provoz" / "parametry.pkl")
params = state["params"]
bp, hp, cp, pp = params
loc = df.index.tz_convert(etl.TZ)
r1 = lambda a: [None if not np.isfinite(v) else round(float(v), 1) for v in np.asarray(a, float)]


def week(col: str) -> dict:
    """Tyden (Po-Ne) posledniho roku s nejvyssim prumerem slozky `col`, rozlozeny na slozky."""
    last = df[df.index >= df.index[-1] - pd.Timedelta(days=365)]
    d = explain.decompose(last, params)
    monday = (last.index.tz_convert(etl.TZ).tz_localize(None).normalize()
              - pd.to_timedelta(last.index.tz_convert(etl.TZ).dayofweek, unit="D"))
    size = d[col].groupby(monday).size()
    best = d[col].groupby(monday).mean()[size >= 7 * 92].idxmax()
    sel = np.asarray(monday == best)
    w = d[sel]
    out = {"t": last.index[sel].tz_convert(etl.TZ).strftime("%Y-%m-%d %H:%M").tolist(),
           "skutecnost": r1(last["baseload"][sel]), "model": r1(w["model"]),
           "kal": r1(w[list(explain.CALENDAR)].sum(axis=1)),
           **{c: r1(w[c]) for c in ("osvetleni", "sero", "topeni", "chlazeni", "fve")},
           "teplota": round(float(last["temp"][sel].mean()), 1)}
    return out


def inertia() -> dict:
    """Tri tydny kolem nejvetsiho ochlazeni posledni zimy: venkovni teplota
    proti teplote, kterou "citi" topeni (osvit, vitr, setrvacnost budov)."""
    ts = pd.Series(heating.t_star(df, hp["a"], hp["b"], hp["w"], hp["tau"]), index=df.index)
    last = df[df.index >= df.index[-1] - pd.Timedelta(days=365)]
    day = last["temp"].groupby(last.index.floor("D")).mean()
    change = (day.shift(-2) - day.shift(2))[day.index.month.isin([11, 12, 1, 2])]
    drop = change.idxmin()        # nejprudsi pokles behem 4 dni v topne sezone
    sel = (df.index >= drop - pd.Timedelta(days=8)) & (df.index < drop + pd.Timedelta(days=13))
    h = df[sel].resample("h").mean(numeric_only=True)
    return {"t": h.index.tz_convert(etl.TZ).strftime("%Y-%m-%d %H:%M").tolist(), "teplota": r1(h["temp"]),
            "pocitova": r1(ts[sel].resample("h").mean()), "topeni": r1(pd.Series(heating.predict(df, hp), index=df.index)[sel].resample("h").mean())}


def backtests() -> list:
    out = []
    for label, name in BACKTESTS:
        path = workspace.ROOT / "simulace" / name / "predikce_D1.csv"
        if not path.exists():
            continue
        d = pd.read_csv(path, sep=";", decimal=",")
        y = d["skutecnost"]
        row = {"rezim": label}
        for key, col in (("model", "predikce"), ("kor", "predikce_kor")):
            e = y - d[col]
            row[key + "_rmse"] = round(float(np.sqrt((e**2).mean())), 1)
            row[key + "_mape"] = round(float((e.abs() / y).mean() * 100), 2)
        out.append(row)
    return out


def horizons() -> dict | None:
    out = {}
    for label, name in INTRADAY:
        path = workspace.ROOT / "simulace" / name / "horizonty.csv"
        if not path.exists():
            return None
        h = pd.read_csv(path, sep=";", decimal=",")
        h = h[h["horizont_h"] <= 12]          # delsi horizonty ma jen cast vydani
        out[label] = {"h": h["horizont_h"].tolist(), "model": r1(h["rmse_model"]), "kor": r1(h["rmse_kor"]),
                      "mape": [round(float(v), 2) for v in h["mape_kor"]]}
    return out


s = operation.end_of_data(df)
yearly = df["baseload"].groupby(loc.year).mean()
data = {
    "zima": week("topeni"), "leto": week("chlazeni"), "setrvacnost": inertia(),
    "curves": explain.curves(params, s), "backtests": backtests(), "horizons": horizons(),
    "params": [{**r, "hodnota": ("0" if abs(r["hodnota"]) < 1e-6 else f"{r['hodnota']:.4g}").replace(".", ",")}
               for r in operation.physical(params)],
    "meta": {
        "od": f"{loc[0]:%d.%m.%Y}", "do": f"{loc[-1]:%d.%m.%Y}", "intervalu": len(df),
        "prumer": round(float(df["baseload"].mean()), 0),
        "min": round(float(df["baseload"].min()), 0), "max": round(float(df["baseload"].max()), 0),
        "n_linearni": int(len(bp["coef"]) + len(hp["k_coef"]) + len(hp["trend_coef"]) + len(cp["k_coef"]) + len(pp["delta"])),
        "n_fyzikalni": len(operation.physical(params)),
        "rmse_trenink": round(float(bp["rmse"]), 1), "r2": round(float(bp["r2"]), 3),
        "stav": f"{state['data_do'].tz_convert(etl.TZ):%d.%m.%Y}", "fve_osvit": explain.PV_REF_SUN,
        "vytvoreno": f"{pd.Timestamp.now():%d.%m.%Y}",
    },
}
from plotly.offline import get_plotlyjs
html = (Path(__file__).with_name("dokumentace_template.html").read_text(encoding="utf-8")
        .replace("/*__DATA__*/null", json.dumps(data, separators=(",", ":"), ensure_ascii=False))
        .replace("/*__PLOTLY__*/", get_plotlyjs()))
OUT.mkdir(exist_ok=True)
path = OUT / "dokumentace.html"
path.write_text(html, encoding="utf-8")
print(f"dokumentace: {path} ({path.stat().st_size / 1e6:.1f} MB); tydny: zima od {data['zima']['t'][0][:10]}, "
      f"leto od {data['leto']['t'][0][:10]}; zpetne simulace: {len(data['backtests'])}, horizonty: {data['horizons'] is not None}")
if "--no-open" not in sys.argv:
    webbrowser.open(path.as_uri())
