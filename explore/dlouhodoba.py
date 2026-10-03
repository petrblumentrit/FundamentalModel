"""Dlouhodoba predikce od konce dat (src/longterm.py), vychozi rok dopredu.

    uv run python explore/dlouhodoba.py --pocasi=scenare

Volby:
    --pocasi=normal | scenare   pocasi za koncem dat: jedna normalova draha, nebo
                                scenare z historickych let (prumer + pasmo nejistoty)
    --trend=posledni | drzet    pomale slozky (uroven, topna citlivost, kapacita FVE);
                                i vyber `uroven,topeni,fve` a tlumeni `posledni:0.5`
    --dni=365                   delka predikce
    --out=simulace/dlouhodoba   vystupni slozka
    --parametry=provoz          model z provozniho stavu (provoz/parametry.pkl, vcetne
                                pevnych hodnot z config/model.yaml) misto models/*_joint.npz
    --no-open                   neotvirat graf v prohlizeci

Model = ulozeny joint fit ze vsech dat (models/*_joint.npz). Vystupy:
    predikce.csv           15min predikce se slozkami (u scenaru i pasmo)
    predikce_scenare.csv   jednotlive drahy scenaru (jen --pocasi=scenare)
    predikce_mesice.csv    mesicni prumery
    predikce.html          interaktivni offline graf
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
import base
import cooling
import etl
import fit
import heating
import longterm

OUT = workspace.ROOT / "simulace" / "dlouhodoba"
WEATHER, TREND, DAYS = "normal", "posledni", 365
for _a in sys.argv[1:]:
    if _a.startswith("--out="):
        OUT = workspace.ROOT / _a.split("=", 1)[1]
    if _a.startswith("--pocasi="):
        WEATHER = _a.split("=", 1)[1]
    if _a.startswith("--trend="):
        TREND = _a.split("=", 1)[1]
    if _a.startswith("--dni="):
        DAYS = int(_a.split("=", 1)[1])
REF_SHIFT = pd.Timedelta(days=364)   # skutecnost pred rokem, zarovnana na den v tydnu

df = etl.load()
if "--parametry=provoz" in sys.argv:
    params = pd.read_pickle(workspace.ROOT / "provoz" / "parametry.pkl")["params"]
else:
    params = (base.load(base.MODEL_DIR / "base_joint.npz"), heating.load(heating.MODEL_DIR / "heating_joint.npz"),
              fit.load("cooling_joint.npz", cooling.PARAM_NAMES + ["rmse", "r2"]),
              fit.load("pv_joint.npz", ["gamma", "span_days", "knot_days", "rmse", "r2"]))
cutoff = df.index[-1] + longterm.STEP
f = longterm.forecast(df, params, cutoff, DAYS, weather=WEATHER, trend=TREND)
scen = longterm.scenario_columns(f)
has_band = "pasmo_dolni" in f

loc = f.index.tz_convert(etl.TZ)
ref = df["baseload"].reindex(f.index - REF_SHIFT).to_numpy()
mon = longterm.monthly(f)
mon["loni"] = pd.Series(df["baseload"].reindex(f.index - pd.Timedelta(days=365)).to_numpy(), index=f.index).groupby(
    loc.strftime("%Y-%m")).mean()
days_in = pd.Series(1, index=loc.strftime("%Y-%m")).groupby(level=0).size() / 96.0

OUT.mkdir(parents=True, exist_ok=True)
main = f.drop(columns=scen)
main.index.name = "timestamp_utc"
main.insert(0, "timestamp_mistni", loc)
main.to_csv(OUT / "predikce.csv", sep=";", decimal=",", float_format="%.3f")
if scen:
    s = f[scen].copy()
    s.index.name = "timestamp_utc"
    s.to_csv(OUT / "predikce_scenare.csv", sep=";", decimal=",", float_format="%.2f")
mon.assign(dni=days_in).to_csv(OUT / "predikce_mesice.csv", sep=";", decimal=",", float_format="%.2f")

p = f["predikce"]
peak, low = p.idxmax().tz_convert(etl.TZ), p.idxmin().tz_convert(etl.TZ)
last_year = df["baseload"][df.index >= df.index[-1] - pd.Timedelta(days=365)].mean()
print(f"predikce {loc[0]:%d.%m.%Y %H:%M} - {loc[-1]:%d.%m.%Y %H:%M}, pocasi {WEATHER}"
      + (f" ({len(scen)} scenaru)" if scen else "") + f", trend {TREND}")
print(f"prumer {p.mean():.1f} (poslednich 365 dni dat {last_year:.1f}, {100 * (p.mean() / last_year - 1):+.1f} %), "
      f"maximum {p.max():.0f} ({peak:%d.%m.%Y %H:%M}), minimum {p.min():.0f} ({low:%d.%m.%Y %H:%M})")
print(mon.round(1).to_string())


def r1(a) -> list:
    return [None if not np.isfinite(v) else round(float(v), 1) for v in np.asarray(a, float)]


# denni prumery pro siroky zaber grafu (cas = poledne daneho dne)
day_key = loc.strftime("%Y-%m-%d 12:00")
daily = pd.DataFrame({"p": p.to_numpy(), "h": ref}, index=day_key)
if has_band:
    daily["lo"], daily["hi"] = f["pasmo_dolni"].to_numpy(), f["pasmo_horni"].to_numpy()
full_day = daily.groupby(level=0).size() >= 92      # bez neuplneho prvniho a posledniho dne
daily = daily.groupby(level=0).mean()[full_day]

data = {
    "d": {"t": daily.index.tolist(), "p": r1(daily["p"]), "h": r1(daily["h"]),
          "lo": r1(daily["lo"]) if has_band else None, "hi": r1(daily["hi"]) if has_band else None},
    "t": loc.strftime("%Y-%m-%d %H:%M").tolist(),
    "p": r1(p), "h": r1(ref),
    "lo": r1(f["pasmo_dolni"]) if has_band else None,
    "hi": r1(f["pasmo_horni"]) if has_band else None,
    "monthly": [{"mesic": m, "dni": round(float(days_in[m]), 1),
                 **{k: (None if pd.isna(v) else round(float(v), 1)) for k, v in row.items()}}
                for m, row in mon.iterrows()],
    "meta": {"pocasi": WEATHER, "trend": TREND, "scenaru": len(scen), "kvantily": len(scen) >= longterm.MIN_QUANTILE,
             "od": f"{loc[0]:%d.%m.%Y}", "do": f"{loc[-1]:%d.%m.%Y}", "dni": DAYS,
             "prumer": round(float(p.mean()), 1), "loni": round(float(last_year), 1),
             "max": round(float(p.max()), 0), "max_t": f"{peak:%d.%m.%Y %H:%M}",
             "min": round(float(p.min()), 0), "min_t": f"{low:%d.%m.%Y %H:%M}",
             "sirka": round(float((f["pasmo_horni"] - f["pasmo_dolni"]).mean()), 1) if has_band else None,
             "model_sd": longterm.MODEL_SD},
}
from plotly.offline import get_plotlyjs
# plotly.js primo v souboru — graf funguje offline, bez CDN
html = (Path(__file__).with_name("dlouhodoba_template.html").read_text(encoding="utf-8")
        .replace("/*__DATA__*/null", json.dumps(data, separators=(",", ":")))
        .replace("/*__PLOTLY__*/", get_plotlyjs()))
out = OUT / "predikce.html"
out.write_text(html, encoding="utf-8")
print(f"graf: {out}")
if "--no-open" not in sys.argv:
    webbrowser.open(out.as_uri())
