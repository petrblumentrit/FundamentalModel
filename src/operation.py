"""Provozni predikce: po prichodu mereni vyda predikci od konce dat dal.

Postup pri kazdem spusteni (predikce.py v koreni projektu):

1. data (etl.load) — konec dat s = prvni nezmereny interval,
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
import etl
import intraday
import longterm
import meteo_forecast

STEP = pd.Timedelta(minutes=15)
LINEAR_AGE = pd.Timedelta(days=1)    # stari dat linearniho prefitu, po kterem se opakuje
SHAPE_AGE = pd.Timedelta(days=7)     # totez pro nelinearni tvar
FC_HISTORY = pd.Timedelta(days=45)   # predpoved pocasi zpet: klouzavy bias (30 dni) a chyba pocasi v korekci


def end_of_data(df: pd.DataFrame) -> pd.Timestamp:
    return df.index[-1] + STEP


def update_params(df: pd.DataFrame, state: dict | None, force: str | None = None) -> tuple[dict, str]:
    """Stav {params, data_do, tvar_do} platny pro konec dat df a co se prefitovalo
    ("tvar" | "linearni" | "nic"). force = "tvar" | "linearni" | "ne"."""
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
    if force == "linearni" or s - state["data_do"] >= LINEAR_AGE:
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
             coefs: pd.DataFrame | None = None) -> pd.DataFrame:
    """Predikce od konce dat: slozky modelu, korekce a vysledek po 15 min."""
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
    comp = backtest.components(fr, params)
    p = comp["predikce"].to_numpy()
    j = np.flatnonzero((pos >= i) & ~np.isnan(p))
    out = comp.iloc[j].rename(columns={"predikce": "model"})
    out.insert(0, "horizont_h", (j - i + 1) / 4)
    out.insert(0, "timestamp_mistni", out.index.tz_convert(etl.TZ))
    out["korekce"] = 0.0
    if coefs is not None:
        x = {**intraday.inputs(y - p, i, j),
             **{k + "_w": v for k, v in intraday.inputs(p - pw, i, j).items()}}
        out["korekce"] = np.asarray(intraday.apply(coefs.iloc[j - i], x), float)
    out["predikce"] = out["model"] + out["korekce"]
    out.index.name = "timestamp_utc"
    return out
