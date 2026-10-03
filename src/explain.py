"""Mezivysledky modelu — rozklad predikce na slozky a krivky modelu.

Model je aditivni, takze rozklad je jen jiny pohled na cleny, ktere se uz
pocitaji; nic se neodhaduje znovu:

    model = uroven + profil + prazdniny + po_volnu + mosty + obdobi   (rozvrh a kalendar)
          + osvetleni + sero + topeni + chlazeni + fve                (fve <= 0)

Je to prirazeni podle modelu, ne mereni: cleny se sdilenym regresorem (osvit
u FVE, sera a solarnich zisku; vecerni spicka u osvetleni a topeni) od sebe
oddeluje jen struktura modelu.

Vliv osvitu, vetru a setrvacnosti na topeni a chlazeni neni samostatny
scitanec — vstupuje dovnitr nelinearni krivky. `weather_effects` ho vycisluje
jako rozdil dvou vypoctu (s vlivem a bez nej); tyto rozdily se navzajem ani
s ostatnimi slozkami nescitaji.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

import base
import cooling
import etl
import heating
import kalendar
import pv

TEMPLATE = Path(__file__).resolve().parent.parent / "explore" / "vysvetleni_template.html"
PV_REF_SUN = 845.0   # [W/m2] osvit, ke kteremu se vztahuje "spicka FVE" (jako backtest.summary)
# sloupec rozkladu -> bloky koeficientu baze (base._slices); obdobi = vsechny obd_*
BASE_PARTS = {
    "uroven": ("level",),
    "profil": ("dt0", "dt1", "dt2", "dt3"),
    "prazdniny": ("praz", "leto_prac", "leto_vol"),
    "po_volnu": ("noc",),
    "mosty": ("most", "most2"),
    "osvetleni": ("svetlo",),
    "sero": ("sero",),
}
CALENDAR = ("uroven", "profil", "prazdniny", "po_volnu", "mosty", "obdobi")   # cast baze nezavisla na pocasi a slunci
WEATHER = ("sero", "topeni", "chlazeni", "fve")                                # cleny zavisle na pocasi
PARTS = (*CALENDAR, "osvetleni", *WEATHER)


def base_parts(df: pd.DataFrame, bp: dict) -> pd.DataFrame:
    """Bazalni slozka rozdelena na bloky (soucet = base.predict)."""
    X = base.design(df, pd.Timestamp(bp["t0"]).date(), bp["span_days"])
    sl = base._slices(bp["span_days"])
    c = bp["coef"]
    out = pd.DataFrame(index=df.index)
    for name, keys in BASE_PARTS.items():
        out[name] = sum(X[:, sl[k]] @ c[sl[k]] for k in keys)
    out["obdobi"] = sum((X[:, sl["obd_" + g]] @ c[sl["obd_" + g]] for g, _ in kalendar.period_groups()),
                        np.zeros(len(df)))
    return out


def decompose(df: pd.DataFrame, params: tuple) -> pd.DataFrame:
    """Aditivni slozky modelu na ose df (sloupce PARTS, `baze` a `model`)."""
    bp, hp, cp, pp = params
    out = base_parts(df, bp)
    out["topeni"] = heating.predict(df, hp)
    out["chlazeni"] = cooling.predict(df, cp)
    out["fve"] = -pv.predict(df, pp)
    out = out[list(PARTS)]
    out["baze"] = out[[*CALENDAR, "osvetleni", "sero"]].sum(axis=1)
    out["model"] = out[list(PARTS)].sum(axis=1)
    return out


def weather_effects(df: pd.DataFrame, params: tuple) -> pd.DataFrame:
    """Vliv osvitu, vetru a setrvacnosti uvnitr topeni a chlazeni jako rozdil
    dvou vypoctu (clen s vlivem minus clen bez nej)."""
    _, hp, cp, _ = params
    h, c = heating.predict(df, hp), cooling.predict(df, cp)
    return pd.DataFrame({
        "topeni_osvit": h - heating.predict(df, {**hp, "a": 0.0}),          # solarni zisky (<= 0)
        "topeni_vitr": h - heating.predict(df, {**hp, "b": 0.0}),           # ochlazovani vetrem (>= 0)
        "topeni_setrvacnost": h - heating.predict(df, {**hp, "w": 1.0}),    # proti okamzite reakci na teplotu
        "chlazeni_osvit": c - cooling.predict(df, {**cp, "a_c": 0.0}),
        "chlazeni_setrvacnost": c - cooling.predict(df, {**cp, "w_c": 1.0}),
    }, index=df.index)


def cleaned(y: np.ndarray, d: pd.DataFrame) -> np.ndarray:
    """Spotreba ocistena o pocasi: skutecnost bez topeni, chlazeni, FVE a sera."""
    return y - d[list(WEATHER)].sum(axis=1).to_numpy()


def curves(params: tuple, now: pd.Timestamp) -> dict:
    """Krivky modelu ve stavu k `now`: odezva na teplotu, citlivost podle
    hodiny dne, denni profily baze a pomale slozky v case."""
    bp, hp, cp, pp = params
    tod = np.arange(0, 24, 0.25)
    day = now.tz_convert(etl.TZ).tz_localize(None).normalize()
    trend_now = heating.trend_profile(now, hp, tod) if "trend_coef" in hp else np.zeros(len(tod))
    k_h = {dt: heating.k_curve(dt, hp, tod) + trend_now for dt in range(4)}
    k_c = {dt: cooling.k_curve(dt, cp, tod) for dt in range(4)}
    # odezva na teplotu v ustalenem stavu (T* = T, bez osvitu a vetru), prumerny pracovni den
    T = np.arange(-15.0, 36.0, 1.0)
    heat = k_h[0].mean() * heating.softplus(hp["t_b"] - T, hp["s"]) / np.maximum(1.0 + hp["alpha"] * T, heating.COP_FLOOR)
    cool = k_c[0].mean() * heating.softplus(T - cp["t_bc"], cp["s_c"]) / np.maximum(
        1.0 - cp["alpha_c"] * (T - cooling.T_REF), cooling.EER_FLOOR)
    weeks = pd.date_range(pd.Timestamp(pp["t0"]), now, freq="7D")
    trend = heating.trend_curve(weeks, hp)
    r = lambda a: [round(float(v), 2) for v in a]
    return {
        "teplota": r(T), "topeni_T": r(heat), "chlazeni_T": r(cool),
        "hodina": r(tod),
        "k_topeni": [r(k_h[dt]) for dt in range(4)], "k_chlazeni": [r(k_c[dt]) for dt in range(4)],
        "profil": [r(base.profile_curve(dt, bp, tod=tod)) for dt in range(4)],
        "osvetleni": r(base.light_curve(bp, tod)), "sero": r(base.gloom_curve(bp, tod)),
        "tyden": [str(d.tz_convert(etl.TZ).date()) for d in weeks],
        "uroven": r(base.level_curve(weeks, bp)),
        "fve_spicka": r(pv.capacity_curve(weeks, pp) * PV_REF_SUN),
        "topna_citlivost": r(heating.k_curve(0, hp, tod).mean() + trend),
        "den": str(day.date()),
    }


def render(full: pd.DataFrame, curves_: dict, meta: dict, path: Path) -> Path:
    """Offline HTML (plotly.js v souboru): rozklad v case + krivky modelu.

    full = rozklad na 15min ose (sloupce PARTS, model, korekce, predikce,
    skutecnost, vlivy pocasi), meta = popisne udaje vcetne tabulky parametru."""
    from plotly.offline import get_plotlyjs
    cols = ["skutecnost", *PARTS, "model", "korekce", "predikce", "ocistena",
            *[c for c in full.columns if c.startswith(("topeni_", "chlazeni_"))]]
    data = {
        "t": full.index.tz_convert(etl.TZ).strftime("%Y-%m-%d %H:%M").tolist(),
        "s": {c: [None if np.isnan(v) else round(float(v), 1) for v in full[c].to_numpy(float)] for c in cols},
        "curves": curves_, "meta": meta,
    }
    html = (TEMPLATE.read_text(encoding="utf-8")
            .replace("/*__DATA__*/null", json.dumps(data, separators=(",", ":"), ensure_ascii=False))
            .replace("/*__PLOTLY__*/", get_plotlyjs()))
    path.write_text(html, encoding="utf-8")
    return path
