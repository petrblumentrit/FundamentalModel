"""Simulace provozni predikce D+1 — klouzavy backtest.

Provozni rezim, ktery simulace napodobuje:

- predikce se vydava v den D v 10:00 mistniho casu na cely den D+1,
- k dispozici je spotreba do D 09:00 (H-1): posledni interval je 08:45-09:00,
  casove znacky jsou zacatky intervalu, takze trenink = radky < D 09:00,
- model se kazdy den **prefituje** z dat do cutoffu: bazalni fit z mirnych dnu
  (levny, linearni) a joint fit s warm startem nelinearnich parametru z
  predchoziho dne; urovnova spline i kapacita FVE koncem treninku koncí,
  za nim drzi (stejne jako v dopredne validaci, krok 7).

Meteo na zbytek dne D a na D+1 se bere **skutecne namerene** — predpoved
pocasi v datech neni, takze simulace meri chybu modelu pri dokonale predpovedi
pocasi. V provozu k ni pribude chyba meteo predpovedi.
"""
import numpy as np
import pandas as pd

import base
import cooling
import etl
import fit
import heating
import pv
import validate

ISSUE_HOUR = 10   # vydani predikce v D [h, mistni cas]
CUTOFF_HOUR = 9   # konec dostupnych dat v D (H-1)


def _wall(day: pd.Timestamp, hour: int) -> pd.Timestamp:
    """Hodina na mistnich hodinach (v den zmeny casu neni D 00:00 + 9 h = 09:00)."""
    return pd.Timestamp(day.year, day.month, day.day, hour).tz_localize(etl.TZ)


def cutoff(issue_day: pd.Timestamp) -> pd.Timestamp:
    return _wall(issue_day, CUTOFF_HOUR)


def train_mask(df: pd.DataFrame, issue_day: pd.Timestamp) -> np.ndarray:
    return (df.index < cutoff(issue_day).tz_convert("UTC"))


def refit(df: pd.DataFrame, mask: np.ndarray, warm: tuple | None) -> tuple:
    """Fit z radku masky; warm = parametry predchoziho dne (jinak cely retezec)."""
    if warm is None:
        return validate.fit_masked(df, mask)
    _, hp, cp, pp = warm
    bp = base.fit(df[mask])
    idx = df.index[mask]
    pp = dict(pp)
    pp["t0"] = str(idx[0])
    pp["span_days"] = float((idx[-1] - idx[0]).total_seconds() / 86400.0) + 1.0
    return fit.joint(df, bp, hp, cp, pp, mask=mask)


def components(df: pd.DataFrame, params: tuple) -> pd.DataFrame:
    """Aditivni slozky predikce na cele ose (setrvacnostni filtry potrebuji historii)."""
    bp, hp, cp, pp = params
    out = pd.DataFrame({
        "baze": base.predict(df, bp),
        "topeni": heating.predict(df, hp),
        "chlazeni": cooling.predict(df, cp),
        "fve": -pv.predict(df, pp),
    }, index=df.index)
    out["predikce"] = out.sum(axis=1)
    return out


def summary(params: tuple, df: pd.DataFrame, mask: np.ndarray) -> dict:
    """Parametry tvaru + uroven a spicka FVE na konci treninku."""
    out = validate.summary_params(params, df)
    end = pd.DatetimeIndex([df.index[mask][-1]])
    out["fve_spicka"] = float(pv.capacity_curve(end, params[3])[0] * 845)
    out["uroven"] = float(base.level_curve(end, params[0])[0])
    out["rmse_trenink"] = float(params[0]["rmse"])
    return out


def run_days(df: pd.DataFrame, issue_days: list[pd.Timestamp],
             log=print) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Postupne pres dny vydani; prvni den cely retezec, dal warm start."""
    loc_date = pd.to_datetime(df["date_local"])
    y = df["baseload"]
    rows, pars, warm = [], [], None
    for d in issue_days:
        mask = train_mask(df, d)
        params = refit(df, mask, warm)
        warm = params
        target = d + pd.Timedelta(days=1)
        sel = (loc_date == target).to_numpy()
        comp = components(df, params)[sel]
        comp.insert(0, "skutecnost", y[sel])
        comp.insert(0, "vydano", _wall(d, ISSUE_HOUR))
        comp.insert(1, "data_do", cutoff(d))
        rows.append(comp)
        s = summary(params, df, mask)
        s["vydano"] = d.date()
        pars.append(s)
        e = comp["skutecnost"] - comp["predikce"]
        log(f"{d.date()} -> {target.date()}  RMSE {np.sqrt((e**2).mean()):6.1f}  "
            f"bias {e.mean():+6.1f}")
    return pd.concat(rows), pd.DataFrame(pars).set_index("vydano")
