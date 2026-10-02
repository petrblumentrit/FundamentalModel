"""Archiv predpovedi pocasi (Analyza/ArchivMeteo.xlsx, mimo repo) pro
realisticky backtest: model se uci ze skutecnych dat, ale na zbytek dne D a
na D+1 dostava predpoved.

Archiv: hodinove, mistni cas (CET/CEST), teplota [C], oblacnost [%],
rychlost vetru [m/s], od 4. 8. 2022. Hodnota k zacatku hodiny (nejlepsi shoda
se skutecnosti bez posunu: teplota MAE 1,2 C).
Predpoved se uklada v D 10:00 na cely nasledujici den D+1 (24 hodnot) — pro
D+1 je to predpoved z D 10:00, pro zbytek dne D z D-1 10:00 (obe dostupne v
okamziku vydani predikce).

Model potrebuje globalni osvit [W/m2], archiv ma oblacnost: osvit = osvit za
jasne oblohy (Haurwitz, z vysky slunce) x pomer jasnosti kt(oblacnost).
kt se kalibruje empiricky z historie (do `until`), protoze vzorec
Kasten-Czeplak vztah u stredni oblacnosti nezachyti (merena stanice: jasno
kt ~0,64, zatazeno ~0,15).
"""
from pathlib import Path

import numpy as np
import pandas as pd

import sun
from etl import TZ

FILE = Path(__file__).resolve().parent.parent / "Analyza" / "ArchivMeteo.xlsx"
KT_BINS = np.arange(0, 101, 10)        # hrany binu oblacnosti [%]


def load_archive() -> pd.DataFrame:
    """Hodinova predpoved na UTC ose (hodiny po rade od mistni pulnoci)."""
    a = pd.read_excel(FILE)
    a.columns = ["ts", "temp", "obl", "wind"]
    a = a.dropna(subset=["temp"])
    ts = pd.to_datetime(a["ts"])
    day = ts.dt.normalize()
    utc = day.dt.tz_localize(TZ).dt.tz_convert("UTC") + pd.to_timedelta(a.groupby(day).cumcount(), unit="h")
    out = a[["temp", "obl", "wind"]].astype(float)
    out.index = pd.DatetimeIndex(utc)
    return out[~out.index.duplicated()]


def clear_sky(index_utc: pd.DatetimeIndex) -> np.ndarray:
    """Osvit za jasne oblohy [W/m2] (Haurwitz) ve stredu intervalu."""
    cz = np.clip(np.sin(np.radians(sun.elevation(index_utc))), 0, None)
    return np.where(cz > 0.01, 1098 * cz * np.exp(-0.057 / np.maximum(cz, 0.01)), 0.0)


def calibrate(arch: pd.DataFrame, actual: pd.DataFrame, until: pd.Timestamp) -> np.ndarray:
    """kt v binech oblacnosti z hodin pred `until` (bez uniku z budoucnosti)."""
    act = actual["sun"].groupby(actual.index.floor("h")).mean()
    j = arch[["obl"]].join(act, how="inner")
    j = j[j.index < until]
    j["cs"] = clear_sky(j.index + pd.Timedelta(minutes=30))
    j = j[j.cs > 50]
    b = np.clip(np.digitize(j.obl, KT_BINS[1:-1]), 0, len(KT_BINS) - 2)
    # pomer souctu (vahy podle osvitu za jasne oblohy), ne prumer pomeru
    return np.array([j.sun[b == k].sum() / j.cs[b == k].sum() for k in range(len(KT_BINS) - 1)])


def to_15min(arch: pd.DataFrame, index_15: pd.DatetimeIndex, kt: np.ndarray) -> pd.DataFrame:
    """Predpoved na 15min ose: teplota, vitr a oblacnost linearne mezi stredy
    hodin, osvit = jasna obloha (stred intervalu) x kt(oblacnost)."""
    mid_h = arch.index + pd.Timedelta(minutes=30)
    mid_q = index_15 + pd.Timedelta(minutes=7.5)
    x, xq = mid_h.asi8.astype(float), mid_q.asi8.astype(float)
    out = pd.DataFrame(index=index_15)
    for c in ("temp", "wind", "obl"):
        out[c] = np.interp(xq, x, arch[c].to_numpy())
    centers = (KT_BINS[:-1] + KT_BINS[1:]) / 2
    out["sun"] = clear_sky(mid_q) * np.interp(out["obl"].to_numpy(), centers, kt)
    return out
