"""Provozni predikce: po prichodu mereni vyda predikci od konce dat dal.

Postup pri kazdem spusteni (predikce.py v koreni projektu):

1. data (etl.load) — vstupni CSV se ctou pokazde cela (~1,3 s), takze zpetne
   zpresnena mereni se projevi sama: v rezidualni korekci hned, v parametrech
   pri pristim prefitu (pri zmene REVISION_REFIT a vic intervalu okamzite).
   Konec dat s = prvni interval, kde chybi spotreba nebo namerene pocasi,
2. parametry modelu z ulozeneho stavu; prefit jen kdyz zestarly: linearni
   koeficienty po LINEAR_AGE (~7 s), nelinearni tvar po SHAPE_AGE (~45 s,
   warm start z ulozenych). Bez ulozeneho stavu cely retezec (~100 s),
3. pocasi: do s mereni, od s predpoved z archivu (meteo_forecast) — nejnovejsi,
   ktera v nem je; horizont konci tam, kde konci predpoved, nejdal
   intraday.HORIZON intervalu,
4. model + korekce zavisla na horizontu (src/intraday.py). Koeficienty korekce
   se berou z intraday backtestu (explore/intraday.py --meteo=predpoved ->
   koeficienty.csv); vstupy korekce se pocitaji z rezidui modelu v poslednich
   dnech, takze provoz nepotrebuje historii vydanych predikci.
"""
import numpy as np
import pandas as pd

import yaml

import backtest
import base
import config
import cooling
import etl
import explain
import heating
import intraday
import longterm
import meteo_forecast
import pv

STEP = pd.Timedelta(minutes=15)
_C = config.model()["provoz"]        # hodnoty v config/model.yaml
LINEAR_AGE = pd.Timedelta(days=_C["linearni_prefit_dny"])    # stari dat linearniho prefitu, po kterem se opakuje
SHAPE_AGE = pd.Timedelta(days=_C["tvar_prefit_dny"])     # totez pro nelinearni tvar
REVISION_REFIT = _C["revize_prefit_intervalu"]                  # [intervaly] zpetne zmenenych dat, od kolika se linearni cast prefituje hned
EXPLAIN_DAYS = _C["rozklad_dni"]                    # [dny] historie v rozkladu predikce (rozklad.csv, vysvetleni.html)
FC_HISTORY = pd.Timedelta(days=45)   # predpoved pocasi zpet: klouzavy bias (30 dni) a chyba pocasi v korekci


SECTIONS = {"topeni": (1, heating.PARAM_NAMES), "chlazeni": (2, cooling.PARAM_NAMES), "fve": (3, pv.PARAM_NAMES)}
DIGITS = 6      # platnych cislic v parametry.yaml
HEADER = """# Fyzikalni parametry modelu, se kterymi se prave pocita (generuje predikce.py).
# Hodnotu lze rucne prepsat: pri pristim spusteni se prevezme a linearni cast
# modelu se k ni prefituje. Uprava vydrzi do pristiho prefitu tvaru (tydne);
# trvale se parametr drzi polozkou `pevna` v config/model.yaml.
"""


def _rounded(v: float) -> float:
    return float(f"{float(v):.{DIGITS}g}")


def physical(params: tuple) -> list[dict]:
    """Fyzikalni parametry s popisem z konfigurace: radky {sekce, nazev,
    hodnota, jednotka, popis, puvod} — pro parametry.yaml a tabulku v grafu."""
    rows = []
    for sec, (k, names) in SECTIONS.items():
        sp = config.spec(sec)
        for n in names:
            v, lo, hi = float(params[k][n]), *map(float, sp[n]["meze"])
            if sp[n].get("pevna") is not None:
                src = "pevná hodnota z konfigurace"
            elif abs(v - lo) <= 1e-6 * (hi - lo):
                src = "odhad, na dolní mezi"
            elif abs(v - hi) <= 1e-6 * (hi - lo):
                src = "odhad, na horní mezi"
            else:
                src = "odhad"
            rows.append({"sekce": sec, "nazev": n, "hodnota": _rounded(v), "jednotka": sp[n]["jednotka"],
                         "popis": sp[n]["popis"], "puvod": src})
    return rows


def write_yaml(state: dict, path) -> None:
    """Citelny stav modelu: fyzikalni parametry + odvozene veliciny ke konci dat."""
    bp, hp, cp, pp = state["params"]
    end = pd.DatetimeIndex([state["data_do"]])
    loc = lambda t: f"{t.tz_convert(etl.TZ):%Y-%m-%d %H:%M}"
    out = {"stav": {
        "linearni_cast_z_dat_do": loc(state["data_do"]), "tvar_z_dat_do": loc(state["tvar_do"]),
        "rmse_trenink": round(float(bp["rmse"]), 2),
        "uroven_portfolia": round(float(base.level_curve(end, bp)[0]), 1),
        "topna_citlivost_na_stupen": round(float(heating.k_curve(0, hp).mean() + heating.trend_curve(end, hp)[0]), 2),
        "spicka_fve": round(float(pv.capacity_curve(end, pp)[0] * explain.PV_REF_SUN), 1),
    }}
    for r in physical(state["params"]):
        out.setdefault(r["sekce"], {})[r["nazev"]] = {k: r[k] for k in ("hodnota", "jednotka", "puvod", "popis")}
    path.write_text(HEADER + yaml.safe_dump(out, allow_unicode=True, sort_keys=False), encoding="utf-8")


def read_edits(state: dict, path) -> dict:
    """Rucne prepsane hodnoty v parametry.yaml: {(sekce, nazev): hodnota} tam,
    kde se soubor lisi od ulozeneho stavu (v zapsane presnosti)."""
    if not path.exists():
        return {}
    doc = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    out = {}
    for sec, (k, names) in SECTIONS.items():
        for n in names:
            v = (doc.get(sec) or {}).get(n, {}).get("hodnota")
            if v is not None and float(v) != _rounded(state["params"][k][n]):
                out[(sec, n)] = float(v)
    return out


def apply_edits(state: dict, edits: dict) -> dict:
    """Stav s prevzatymi rucnimi hodnotami (linearni cast je pak nutne prefitovat)."""
    params = [dict(p) for p in state["params"]]
    for (sec, n), v in edits.items():
        params[SECTIONS[sec][0]][n] = v
    return {**state, "params": tuple(params)}


def config_pins(state: dict) -> list[str]:
    """Parametry, jejichz pevna hodnota v config/model.yaml se lisi od stavu."""
    out = []
    for sec, (k, names) in SECTIONS.items():
        fix = config.fixed(sec, names)
        out += [n for n, f in zip(names, fix) if not np.isnan(f) and float(state["params"][k][n]) != f]
    return out


def end_of_data(df: pd.DataFrame) -> pd.Timestamp:
    return df.index[-1] + STEP


def update_params(df: pd.DataFrame, state: dict | None, force: str | None = None,
                  revised: int = 0) -> tuple[dict, str]:
    """Stav {params, data_do, tvar_do} platny pro konec dat df a co se prefitovalo
    ("tvar" | "linearni" | "nic"). force = "tvar" | "linearni" | "ne";
    revised = pocet zpetne zmenenych intervalu v datech, ze kterych se fitovalo.
    Linearni cast se prefituje i tehdy, kdyz se pevna hodnota v konfiguraci
    lisi od stavu (fit.joint ji pri fix_shape prevezme)."""
    s = end_of_data(df)
    if state is not None and force == "ne":
        return state, "nic"
    stale = state is None or state["data_do"] > s     # stav z budoucnosti (prehravani minulosti)
    full = np.ones(len(df), bool)
    if stale or force == "tvar" or s - state["tvar_do"] >= SHAPE_AGE:
        dh = etl.hourly(df)
        shape = backtest.refit(dh, np.ones(len(dh), bool), state["params"] if state else None)
        params = backtest.refit(df, full, shape, fix_shape=True)
        return {"params": params, "data_do": s, "tvar_do": s}, "tvar"
    if force == "linearni" or s - state["data_do"] >= LINEAR_AGE or revised >= REVISION_REFIT or config_pins(state):
        params = backtest.refit(df, full, state["params"], fix_shape=True)
        return {"params": params, "data_do": s, "tvar_do": state["tvar_do"]}, "linearni"
    return state, "nic"


def weather_forecast(df: pd.DataFrame) -> pd.DataFrame:
    """Predpoved pocasi z archivu na 15min ose kolem konce dat; osvit a sero
    kalibrovane z cele historie pred koncem dat."""
    s = end_of_data(df)
    arch = meteo_forecast.load_archive()
    if arch.index.max() < s:
        raise SystemExit(f"archiv predpovedi pocasi konci {arch.index.max()}, pred koncem dat {s}")
    kt = meteo_forecast.calibrate(arch, df, s)
    kt_q = meteo_forecast.calibrate_quantiles(arch, df, s)
    idx = pd.date_range(max(arch.index.min(), s - FC_HISTORY),
                        min(arch.index.max() + 3 * STEP, s + (intraday.HORIZON - 1) * STEP), freq=STEP)
    return meteo_forecast.to_15min(arch, idx, kt, kt_q)


def forecast(df: pd.DataFrame, params: tuple, meteo: pd.DataFrame,
             coefs: pd.DataFrame | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Predikce od konce dat a rozklad na slozky (src/explain.py).

    Vraci (predikce, rozklad): predikce = cile od konce dat po 15 min (slozky,
    model, korekce, vysledek); rozklad = totez i pro poslednich EXPLAIN_DAYS
    dni pred koncem dat (model s namerenym pocasim proti skutecnosti, spotreba
    ocistena o pocasi, vlivy osvitu, vetru a setrvacnosti)."""
    s = end_of_data(df)
    fr = longterm.frame(df, s, days=2).iloc[:longterm.HISTORY_DAYS * intraday.DAY + intraday.HORIZON].copy()
    pos = np.arange(len(fr))
    i = int((fr.index < s).sum())
    y = fr["baseload"].to_numpy(float)
    measured = {c: fr[c].to_numpy(float).copy() for c in intraday.METEO_COLS}
    # od konce dat predpoved; kde chybi, zustane NaN a cil se nevyda
    future = intraday._forecast_columns(fr, df, meteo, s)
    # model s predpovedi pocasi i v minulych dnech — chyba pocasi pro korekci
    past = intraday._forecast_columns(fr, df, meteo, s - pd.Timedelta(days=intraday.W_DAYS))
    for c in intraday.METEO_COLS:
        fr[c] = past[c]
    pw = backtest.components(fr, params)["predikce"].to_numpy()
    for c in intraday.METEO_COLS:
        fr[c] = np.where(pos < i, measured[c], future[c])
    full = explain.decompose(fr, params)
    p = full["model"].to_numpy()
    j = np.flatnonzero((pos >= i) & ~np.isnan(p))
    full["korekce"] = np.nan
    if coefs is not None:
        x = {**intraday.inputs(y - p, i, j),
             **{k + "_w": v for k, v in intraday.inputs(p - pw, i, j).items()}}
        full.iloc[j, full.columns.get_loc("korekce")] = np.asarray(intraday.apply(coefs.iloc[j - i], x), float)
    else:
        full.iloc[j, full.columns.get_loc("korekce")] = 0.0
    full["predikce"] = full["model"] + full["korekce"].fillna(0.0)
    full.insert(0, "skutecnost", y)
    full["ocistena"] = explain.cleaned(y, full)
    full = full.join(explain.weather_effects(fr, params))
    full.insert(0, "horizont_h", np.where(pos >= i, (pos - i + 1) / 4, np.nan))
    full.insert(0, "timestamp_mistni", full.index.tz_convert(etl.TZ))
    full.index.name = "timestamp_utc"
    keep = np.concatenate([np.arange(max(i - EXPLAIN_DAYS * intraday.DAY, 0), i), j])
    out = full.iloc[j].drop(columns=["skutecnost", "ocistena"])
    return out, full.iloc[keep]
