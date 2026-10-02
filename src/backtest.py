"""Simulace provozni predikce D+1 — klouzavy backtest.

Provozni rezim, ktery simulace napodobuje:

- predikce se vydava v den D v 10:00 mistniho casu na cely den D+1,
- k dispozici je spotreba do D 09:00 (H-1): posledni interval je 08:45-09:00,
  casove znacky jsou zacatky intervalu, takze trenink = radky < D 09:00,
- model se kazdy den **prefituje** z dat do cutoffu: vsechny linearni
  koeficienty (uroven, profily, osvetleni, kapacita FVE, trend topeni) denne;
  nelinearni parametry tvaru (topna/chladici krivka, setrvacnost; ze dne na
  den se skoro nemeni) jednou tydne pri vydani v SHAPE_REFIT_DOW (a prvni den
  simulace), mezi tim se drzi. Denni prefit tak trva ~2 s misto ~60 s;
  srovnani s dennim nelinearnim prefitem: RMSE 19,95 vs 19,87, predikce se
  lisi v prumeru o 0,7.
- vypocet je dvoufazovy a paralelni (explore/backtest.py): (1) tvar pro
  kazde pondeli, warm start z fitu k prvnimu dni simulace (data pred ni,
  zadny unik z budoucnosti) — pondeli jsou navzajem nezavisla; (2) linearni
  prefit a predikce pro kazdy den s tvarem posledniho pondeli. Vysledek tak
  nezavisi na poctu workeru. Urovnova spline i kapacita FVE koncem treninku konci,
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
SHAPE_REFIT_DOW = 0   # den vydani s plnym (nelinearnim) prefitem, 0 = pondeli
HOLIDAY_CONTEXT = 3   # [dny] okoli svatku, ktere se pri omezenem okne bere cele
PREDICT_HISTORY = 60  # [dny] historie pro vypocet slozek predikce D+1


def _wall(day: pd.Timestamp, hour: int) -> pd.Timestamp:
    """Hodina na mistnich hodinach (v den zmeny casu neni D 00:00 + 9 h = 09:00)."""
    return pd.Timestamp(day.year, day.month, day.day, hour).tz_localize(etl.TZ)


def cutoff(issue_day: pd.Timestamp) -> pd.Timestamp:
    return _wall(issue_day, CUTOFF_HOUR)


def holiday_context(df: pd.DataFrame) -> np.ndarray:
    """Svatky, mosty a Vanoce (22. 12. - 3. 1.) vcetne +-HOLIDAY_CONTEXT dni.

    Okoli je nutne: bez beznych dni kolem by ve starsich letech chybel opora
    pro uroven a efekt svatku by splynul s urovni.
    """
    days = pd.to_datetime(df["date_local"])
    hol = etl._cz_holidays(days.dt.year.min() - 1, days.dt.year.max() + 1)
    special = pd.DatetimeIndex(sorted(hol | etl._bridges(hol)))
    xmas = pd.DatetimeIndex([d for d in pd.date_range(days.min() - pd.Timedelta(days=10),
                                                      days.max() + pd.Timedelta(days=10))
                             if (d.month == 12 and d.day >= 22) or (d.month == 1 and d.day <= 3)])
    anchor = special.union(xmas)
    near = anchor.copy()
    for k in range(1, HOLIDAY_CONTEXT + 1):
        near = near.union(anchor + pd.Timedelta(days=k)).union(anchor - pd.Timedelta(days=k))
    return days.isin(near).to_numpy()


def train_mask(df: pd.DataFrame, issue_day: pd.Timestamp, window: int | None = None) -> np.ndarray:
    """Data do cutoffu; s window [dny] jen poslednich window dni + okoli svatku
    z cele historie (svatky se uci bez omezeni)."""
    end = cutoff(issue_day).tz_convert("UTC")
    mask = df.index < end
    if window is not None:
        recent = df.index >= end - pd.Timedelta(days=window)
        mask &= recent | holiday_context(df)
    return mask


def refit(df: pd.DataFrame, mask: np.ndarray, warm: tuple | None,
          fix_shape: bool = False) -> tuple:
    """Fit z radku masky; warm = parametry predchoziho dne (jinak cely retezec);
    fix_shape = jen linearni prefit s tvarem z warm."""
    if warm is None:
        return validate.fit_masked(df, mask)
    _, hp, cp, pp = warm
    bp = base.fit(df[mask])
    idx = df.index[mask]
    pp = dict(pp)
    pp["t0"] = str(idx[0])
    pp["span_days"] = float((idx[-1] - idx[0]).total_seconds() / 86400.0) + 1.0
    # warm start je blizko — 100 iteraci staci, pojistka proti pomale konvergenci
    return fit.joint(df, bp, hp, cp, pp, mask=mask, fix_shape=fix_shape, max_nfev=100)


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


def shape_days(issue_days: list[pd.Timestamp]) -> list[pd.Timestamp]:
    """Dny vydani s nelinearnim prefitem tvaru: prvni den a kazde SHAPE_REFIT_DOW."""
    return [d for i, d in enumerate(issue_days) if i == 0 or d.dayofweek == SHAPE_REFIT_DOW]


_hourly_cache: dict = {}


def _hourly(df: pd.DataFrame) -> pd.DataFrame:
    """Hodinova agregace df (v ramci procesu jednou)."""
    key = (id(df), len(df))
    if key not in _hourly_cache:
        _hourly_cache.clear()
        _hourly_cache[key] = etl.hourly(df)
    return _hourly_cache[key]


def shape_fit(df: pd.DataFrame, day: pd.Timestamp, warm: tuple | None,
              window: int | None = None) -> tuple:
    """Plny (nelinearni) prefit z dat do cutoffu dne; warm=None = cely retezec.

    Fituje se z hodinove agregace (4x mene radku, ~3-4x rychleji): parametry
    tvaru (bilancni teplota, sirka, setrvacnost) ctvrthodinove rozliseni
    nepotrebuji. Linearni cast a predikce (forecast_day) zustavaji 15min.
    """
    dh = _hourly(df)
    return refit(dh, train_mask(dh, day, window), warm)


def forecast_day(df: pd.DataFrame, day: pd.Timestamp, shape: tuple,
                 shape_from: pd.Timestamp, window: int | None = None,
                 meteo: pd.DataFrame | None = None) -> tuple[pd.DataFrame, dict]:
    """Linearni prefit s danym tvarem a predikce na D+1 (dny jsou nezavisle).

    meteo = predpoved pocasi na 15min ose (src/meteo_forecast.py): od cutoffu
    (zbytek dne D a cele D+1) nahradi namerenou teplotu, osvit a vitr —
    realisticky rezim; trenink zustava na skutecnych datech.
    """
    mask = train_mask(df, day, window)
    params = refit(df, mask, shape, fix_shape=True)
    target = day + pd.Timedelta(days=1)
    # slozky staci pocitat z poslednich PREDICT_HISTORY dni: setrvacnostni
    # filtry (tau ~85 h) zapomenou pocatek za ~2 tydny, zbytek je po radcich
    dates = pd.to_datetime(df["date_local"])
    recent = df[((dates >= target - pd.Timedelta(days=PREDICT_HISTORY)) & (dates <= target)).to_numpy()]
    if meteo is not None:
        recent = recent.copy()
        fut = recent.index >= cutoff(day).tz_convert("UTC")
        fc = meteo.reindex(recent.index[fut])
        for c in ("temp", "sun", "wind"):
            recent.loc[fut, c] = fc[c].fillna(recent.loc[fut, c]).to_numpy()
    sel = (pd.to_datetime(recent["date_local"]) == target).to_numpy()
    comp = components(recent, params)[sel]
    comp.insert(0, "skutecnost", recent["baseload"][sel])
    comp.insert(0, "vydano", _wall(day, ISSUE_HOUR))
    comp.insert(1, "data_do", cutoff(day))
    s = summary(params, df, mask)
    s["vydano"] = day.date()
    s["tvar_z"] = shape_from.date()
    return comp, s
