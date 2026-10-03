"""Provozni predikce od konce dat (src/operation.py) — spoustet po prichodu mereni.

    uv run python predikce.py

Nacte data, podle stari ulozenych parametru model prefituje (linearni cast
denne, tvar tydne), vyda predikci po 15 min od konce dat, kam saha predpoved
pocasi (nejdal 39 h), a ulozi ji. Volby:

    --prefit=tvar | linearni | ne   vynutit / zakazat prefit (jinak podle stari)
    --korekce=soubor                koeficienty korekce (vychozi z intraday backtestu
                                    simulace/intraday_predpoved/koeficienty.csv); --korekce= vypne
    --konec="2026-07-20 12:00"      prehrani minulosti: data jen pred timto mistnim casem
    --out=provoz                    slozka stavu a vystupu
    --bez-grafu                     nevytvaret vysvetleni.html (~5 MB, plotly.js v souboru)
    --otevrit                       otevrit graf v prohlizeci

Vystupy (mimo repo): parametry.pkl (stav modelu), predikce.csv (posledni
vydani), archiv/predikce_RRRRMMDD_HHMM.csv (vsechna vydani, cas konce dat UTC),
rozklad.csv (slozky modelu za poslednich 14 dni a predikci, spotreba ocistena
o pocasi, vlivy osvitu, vetru a setrvacnosti), vysvetleni.html (offline graf
rozkladu a krivky modelu).
"""
import sys
import time
import webbrowser
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))
import etl
import explain
import operation

OUT = ROOT / "provoz"
COEFS = ROOT / "simulace" / "intraday_predpoved" / "koeficienty.csv"
FORCE, END = None, None
GRAPH = "--bez-grafu" not in sys.argv
OPEN = "--otevrit" in sys.argv
for _a in sys.argv[1:]:
    if _a.startswith("--out="):
        OUT = ROOT / _a.split("=", 1)[1]
    if _a.startswith("--prefit="):
        FORCE = _a.split("=", 1)[1]
    if _a.startswith("--korekce="):
        COEFS = ROOT / _a.split("=", 1)[1] if _a.split("=", 1)[1] else None
    if _a.startswith("--konec="):
        END = pd.Timestamp(_a.split("=", 1)[1]).tz_localize(etl.TZ).tz_convert("UTC")

t0 = time.time()
old = etl.previous()
df = etl.load()
if END is not None:
    df = df[df.index < END]
s = operation.end_of_data(df)
print(f"konec dat: {s.tz_convert(etl.TZ):%Y-%m-%d %H:%M} ({len(df)} intervalu, nacteno za {time.time() - t0:.1f} s)")
# zpetne zmenene hodnoty od minuleho spusteni (zpresneni mereni)
rev = etl.revisions(old, df)
if len(rev):
    loc = rev.index.tz_convert(etl.TZ)
    print(f"zpetne zmeneno {len(rev)} intervalu ({loc[0]:%Y-%m-%d %H:%M} - {loc[-1]:%Y-%m-%d %H:%M}), "
          f"nejvetsi zmena spotreby {rev['baseload'].abs().max():.1f}, teploty {rev['temp'].abs().max():.2f}")

OUT.mkdir(parents=True, exist_ok=True)
path = OUT / "parametry.pkl"
state = pd.read_pickle(path) if path.exists() else None
fitted = state["data_do"] if state else s
state, done = operation.update_params(df, state, FORCE, revised=int((rev.index < fitted).sum()))
if done != "nic":
    pd.to_pickle(state, path)
print(f"parametry: prefit {done} ({time.time() - t0:.0f} s); linearni cast z dat do "
      f"{state['data_do'].tz_convert(etl.TZ):%Y-%m-%d %H:%M}, tvar do {state['tvar_do'].tz_convert(etl.TZ):%Y-%m-%d %H:%M}")

coefs = None
if COEFS is not None and COEFS.exists():
    coefs = pd.read_csv(COEFS, sep=";", decimal=",", index_col=0)
else:
    print("koeficienty korekce nenalezeny — predikce bez korekce (spustit explore/intraday.py --meteo=predpoved)")
t1 = time.time()
f, full = operation.forecast(df, state["params"], operation.weather_forecast(df), coefs)
f.to_csv(OUT / "predikce.csv", sep=";", decimal=",", float_format="%.3f")
(OUT / "archiv").mkdir(exist_ok=True)
f.to_csv(OUT / "archiv" / f"predikce_{s:%Y%m%d_%H%M}.csv", sep=";", decimal=",", float_format="%.3f")
full.to_csv(OUT / "rozklad.csv", sep=";", decimal=",", float_format="%.3f")

loc = f["timestamp_mistni"]
print(f"predikce: {len(f)} intervalu, {loc.iloc[0]:%d.%m. %H:%M} - {loc.iloc[-1]:%d.%m. %H:%M} "
      f"({time.time() - t1:.1f} s, celkem {time.time() - t0:.0f} s)")
show = sorted({k for k in (0, 1, 2, 3, 7, 11, 23, 47, len(f) - 1) if k < len(f)})
cols = ["horizont_h", "baze", "topeni", "chlazeni", "fve", "model", "korekce", "predikce"]
print(f.iloc[show].set_index("timestamp_mistni")[cols].round(1).to_string())
print(f"ulozeno: {OUT / 'predikce.csv'}, rozklad s historii: {OUT / 'rozklad.csv'}")

if GRAPH:
    local = lambda t: t.tz_convert(etl.TZ)
    meta = {
        "konec": f"{local(s):%Y-%m-%d %H:%M}", "konec_text": f"{local(s):%d.%m.%Y %H:%M}",
        "horizont": float(f["horizont_h"].iloc[-1]), "do_text": f"do {loc.iloc[-1]:%d.%m. %H:%M}",
        "data_do": f"{local(state['data_do']):%d.%m.%Y %H:%M}", "tvar_do": f"{local(state['tvar_do']):%d.%m.%Y %H:%M}",
        "fve_osvit": explain.PV_REF_SUN, "params": [],
    }
    out = explain.render(full, explain.curves(state["params"], s), meta, OUT / "vysvetleni.html")
    print(f"graf: {out}")
    if OPEN:
        webbrowser.open(out.as_uri())
