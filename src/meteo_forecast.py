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

Sero pres den (clen baze, src/sun.py) je nelinearni funkce osvitu: predpoved
dava jen stredni osvit pri dane oblacnosti, a sero stredniho osvitu neni
stredni sero (tmave dny by se nikdy nepredpovedely). Proto se pro kazdy bin
oblacnosti kalibruji i **kvantily** kt a sero se pocita jako prumer pres ne —
stredni hodnota nelinearniho clenu pres rozdeleni chyby predpovedi.
"""
from pathlib import Path

import numpy as np
import pandas as pd

import sun
from etl import TZ

FILE = Path(__file__).resolve().parent.parent / "Analyza" / "ArchivMeteo.xlsx"
KT_BINS = np.arange(0, 101, 10)        # hrany binu oblacnosti [%]
KT_QUANTILES = np.arange(0.05, 1.0, 0.1)   # kvantily kt v binu (stredy decilu)


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


def _history(arch: pd.DataFrame, actual: pd.DataFrame, until: pd.Timestamp) -> tuple[pd.DataFrame, np.ndarray]:
    """Hodiny pred `until` se sluncem nad obzorem a jejich bin oblacnosti."""
    act = actual["sun"].groupby(actual.index.floor("h")).mean()
    j = arch[["obl"]].join(act, how="inner")
    j = j[j.index < until]
    j["cs"] = clear_sky(j.index + pd.Timedelta(minutes=30))
    j = j[j.cs > 50]
    return j, np.clip(np.digitize(j.obl, KT_BINS[1:-1]), 0, len(KT_BINS) - 2)


def calibrate(arch: pd.DataFrame, actual: pd.DataFrame, until: pd.Timestamp) -> np.ndarray:
    """kt v binech oblacnosti z hodin pred `until` (bez uniku z budoucnosti)."""
    j, b = _history(arch, actual, until)
    # pomer souctu (vahy podle osvitu za jasne oblohy), ne prumer pomeru
    return np.array([j.sun[b == k].sum() / j.cs[b == k].sum() for k in range(len(KT_BINS) - 1)])


def calibrate_quantiles(arch: pd.DataFrame, actual: pd.DataFrame, until: pd.Timestamp) -> np.ndarray:
    """Kvantily kt v binech oblacnosti (bin x kvantil) — rozdeleni skutecne
    jasnosti pri dane predpovedi oblacnosti."""
    j, b = _history(arch, actual, until)
    ratio = (j.sun / j.cs).to_numpy()
    return np.array([np.quantile(ratio[b == k], KT_QUANTILES) for k in range(len(KT_BINS) - 1)])


def to_15min(arch: pd.DataFrame, index_15: pd.DatetimeIndex, kt: np.ndarray,
             kt_q: np.ndarray | None = None) -> pd.DataFrame:
    """Predpoved na 15min ose: teplota, vitr a oblacnost linearne mezi stredy
    hodin, osvit = jasna obloha (stred intervalu) x kt(oblacnost).

    S kt_q (calibrate_quantiles) navic sloupec `sero_den` = stredni hodnota
    exp(-osvit / I0) pres rozdeleni kt pri dane oblacnosti (bez faktoru tmy)."""
    mid_h = arch.index + pd.Timedelta(minutes=30)
    mid_q = index_15 + pd.Timedelta(minutes=7.5)
    x, xq = mid_h.asi8.astype(float), mid_q.asi8.astype(float)
    out = pd.DataFrame(index=index_15)
    for c in ("temp", "wind", "obl"):
        out[c] = np.interp(xq, x, arch[c].to_numpy())
    centers = (KT_BINS[:-1] + KT_BINS[1:]) / 2
    cs = clear_sky(mid_q)
    out["sun"] = cs * np.interp(out["obl"].to_numpy(), centers, kt)
    if kt_q is not None:
        obl = out["obl"].to_numpy()
        out["sero_den"] = np.mean([np.exp(-cs * np.interp(obl, centers, kt_q[:, q]) / sun.GLOOM_I0)
                                   for q in range(kt_q.shape[1])], axis=0)
    return out
