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

import backtest
import config
import etl
import explain
import intraday
import longterm
import meteo_forecast

STEP = pd.Timedelta(minutes=15)
_C = config.model()["provoz"]        # hodnoty v config/model.yaml
LINEAR_AGE = pd.Timedelta(days=_C["linearni_prefit_dny"])    # stari dat linearniho prefitu, po kterem se opakuje
SHAPE_AGE = pd.Timedelta(days=_C["tvar_prefit_dny"])     # totez pro nelinearni tvar
REVISION_REFIT = _C["revize_prefit_intervalu"]                  # [intervaly] zpetne zmenenych dat, od kolika se linearni cast prefituje hned
EXPLAIN_DAYS = _C["rozklad_dni"]                    # [dny] historie v rozkladu predikce (rozklad.csv, vysvetleni.html)
FC_HISTORY = pd.Timedelta(days=45)   # predpoved pocasi zpet: klouzavy bias (30 dni) a chyba pocasi v korekci


def end_of_data(df: pd.DataFrame) -> pd.Timestamp:
    return df.index[-1] + STEP


def update_params(df: pd.DataFrame, state: dict | None, force: str | None = None,
                  revised: int = 0) -> tuple[dict, str]:
    """Stav {params, data_do, tvar_do} platny pro konec dat df a co se prefitovalo
    ("tvar" | "linearni" | "nic"). force = "tvar" | "linearni" | "ne";
    revised = pocet zpetne zmenenych intervalu v datech, ze kterych se fitovalo."""
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
    if force == "linearni" or s - state["data_do"] >= LINEAR_AGE or revised >= REVISION_REFIT:
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
