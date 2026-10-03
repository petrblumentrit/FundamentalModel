"""Historie spotreby proti spotrebe po odecteni topeni a chlazeni.

    uv run python explore/ocistena.py

Model z provozniho stavu (provoz/parametry.pkl, fit ze vsech dat) rozlozi
celou historii; od skutecne spotreby se odecte modelovy prispevek topeni a
chlazeni (FVE a ostatni slozky zustavaji). Vystupy do simulace/ocistena/
(mimo repo): ocistena.csv (15 min), ocistena_mesice.csv, ocistena.html
(offline graf). `--out=slozka`, `--no-open`.
"""
import json
import sys
import webbrowser
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
import cooling
import etl
import heating

OUT = ROOT / "simulace" / "ocistena"
for _a in sys.argv[1:]:
    if _a.startswith("--out="):
        OUT = ROOT / _a.split("=", 1)[1]

df = etl.load()
state = pd.read_pickle(ROOT / "provoz" / "parametry.pkl")
_, hp, cp, _ = state["params"]
f = pd.DataFrame({"skutecnost": df["baseload"], "topeni": heating.predict(df, hp),
                  "chlazeni": cooling.predict(df, cp)}, index=df.index)
f["ocistena"] = f["skutecnost"] - f["topeni"] - f["chlazeni"]
loc = f.index.tz_convert(etl.TZ)
cols = ["skutecnost", "ocistena", "topeni", "chlazeni"]


def shares(g: pd.DataFrame) -> pd.DataFrame:
    out = g[cols].mean()
    out["podil_topeni"] = 100 * out["topeni"] / out["skutecnost"]
    out["podil_chlazeni"] = 100 * out["chlazeni"] / out["skutecnost"]
    out["dni"] = len(g) / 96.0
    return out


monthly = f.groupby(loc.strftime("%Y-%m")).apply(shares)
yearly = f.groupby(loc.strftime("%Y")).apply(shares)
total = shares(f)

OUT.mkdir(parents=True, exist_ok=True)
out = f.copy()
out.index.name = "timestamp_utc"
out.insert(0, "timestamp_mistni", loc)
out.to_csv(OUT / "ocistena.csv", sep=";", decimal=",", float_format="%.3f")
monthly.to_csv(OUT / "ocistena_mesice.csv", sep=";", decimal=",", float_format="%.2f")

print(f"historie {loc[0]:%d.%m.%Y} - {loc[-1]:%d.%m.%Y}; model: linearni cast z dat do "
      f"{state['data_do'].tz_convert(etl.TZ):%d.%m.%Y %H:%M}")
print(f"rozpeti dennich prumeru: skutecnost {f['skutecnost'].groupby(loc.date).mean().min():.0f}-"
      f"{f['skutecnost'].groupby(loc.date).mean().max():.0f}, bez topeni a chlazeni "
      f"{f['ocistena'].groupby(loc.date).mean().min():.0f}-{f['ocistena'].groupby(loc.date).mean().max():.0f}")
print(yearly.round(1).to_string())


def r1(a) -> list:
    return [None if not np.isfinite(v) else round(float(v), 1) for v in np.asarray(a, float)]


def rows(t: pd.DataFrame, key: str) -> list:
    return [{key: k, **{c: round(float(v), 2) for c, v in row.items()}} for k, row in t.iterrows()]


# denni prumery pro siroky zaber (cas = poledne), jen uplne dny
day_key = loc.strftime("%Y-%m-%d 12:00")
full_day = f.groupby(day_key).size() >= 92
daily = f.groupby(day_key)[cols].mean()[full_day]
data = {
    "t": loc.strftime("%Y-%m-%d %H:%M").tolist(), **{c: r1(f[c]) for c in cols},
    "d": {"t": daily.index.tolist(), **{c: r1(daily[c]) for c in cols}},
    "yearly": rows(yearly, "obdobi"), "monthly": rows(monthly, "obdobi"),
    "meta": {"od": f"{loc[0]:%d.%m.%Y}", "do": f"{loc[-1]:%d.%m.%Y}", **{c: round(float(v), 2) for c, v in total.items()},
             "t_b": round(float(hp["t_b"]), 1), "t_bc": round(float(cp["t_bc"]), 1)},
}
from plotly.offline import get_plotlyjs
html = (Path(__file__).with_name("ocistena_template.html").read_text(encoding="utf-8")
        .replace("/*__DATA__*/null", json.dumps(data, separators=(",", ":")))
        .replace("/*__PLOTLY__*/", get_plotlyjs()))
path = OUT / "ocistena.html"
path.write_text(html, encoding="utf-8")
print(f"graf: {path}")
if "--no-open" not in sys.argv:
    webbrowser.open(path.as_uri())
