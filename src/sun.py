"""Poloha slunce a tma — deterministicke, v predikci nepotrebuji meteo.

NOAA aproximace (presnost ~1 min ve vychodu/zapadu) pro stred CR. Pouziti:
clen osvetleni v bazi (spotreba rizena svetlem, ne hodinami — viz
explore/svetlo_test.py: pri zmene casu se ~45 jednotek spotreby posune se
soumrakem, rozvrh zustava na hodinach).
"""
import numpy as np
import pandas as pd

LAT, LON = 49.8, 15.5        # stred CR
DARK_MID = -2.0              # [deg] vyska slunce, kde je "pul tmy"
DARK_WIDTH = 1.5             # [deg] sirka prechodu (soumrak trva ~30-40 min)


def _solar(ts_utc: pd.DatetimeIndex):
    doy = ts_utc.dayofyear.to_numpy()
    hr = ts_utc.hour.to_numpy() + ts_utc.minute.to_numpy() / 60
    g = 2 * np.pi / 365 * (doy - 1 + (hr - 12) / 24)
    eqt = 229.18 * (0.000075 + 0.001868 * np.cos(g) - 0.032077 * np.sin(g)
                    - 0.014615 * np.cos(2 * g) - 0.040849 * np.sin(2 * g))
    decl = (0.006918 - 0.399912 * np.cos(g) + 0.070257 * np.sin(g) - 0.006758 * np.cos(2 * g)
            + 0.000907 * np.sin(2 * g) - 0.002697 * np.cos(3 * g) + 0.00148 * np.sin(3 * g))
    return hr, eqt, decl


def elevation(ts_utc: pd.DatetimeIndex) -> np.ndarray:
    """Vyska slunce [deg] v okamzicich ts_utc."""
    hr, eqt, decl = _solar(ts_utc)
    ha = np.radians((hr * 60 + eqt + 4 * LON) / 4 - 180)
    lat = np.radians(LAT)
    cz = np.sin(lat) * np.sin(decl) + np.cos(lat) * np.cos(decl) * np.cos(ha)
    return 90 - np.degrees(np.arccos(np.clip(cz, -1, 1)))


def sun_times(dates: pd.DatetimeIndex) -> tuple[pd.Series, pd.Series]:
    """Vychod a zapad (UTC) pro kalendarni dny."""
    noon = pd.DatetimeIndex(dates).tz_localize("UTC") + pd.Timedelta(hours=12)
    _, eqt, decl = _solar(noon)
    lat = np.radians(LAT)
    ha0 = np.degrees(np.arccos(np.cos(np.radians(90.833)) / (np.cos(lat) * np.cos(decl))
                               - np.tan(lat) * np.tan(decl)))
    day0 = pd.DatetimeIndex(dates).tz_localize("UTC")
    rise = day0 + pd.to_timedelta(720 - 4 * (LON + ha0) - eqt, unit="min")
    sset = day0 + pd.to_timedelta(720 - 4 * (LON - ha0) - eqt, unit="min")
    return pd.Series(rise, index=dates), pd.Series(sset, index=dates)


def darkness(elev: np.ndarray) -> np.ndarray:
    """Plynula tma: 0 za dne, 1 po obcanskem soumraku."""
    return 1 / (1 + np.exp((elev - DARK_MID) / DARK_WIDTH))


def darkness_15min(index_utc: pd.DatetimeIndex) -> np.ndarray:
    """Prumerna tma v 15min intervalu (znacky = zacatky intervalu)."""
    sub = [darkness(elevation(index_utc + pd.Timedelta(minutes=m))) for m in (2.5, 7.5, 12.5)]
    return np.mean(sub, axis=0)
