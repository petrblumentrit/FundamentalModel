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
jasne oblohy (Haurwitz, z vysky slunce) x pomer jasnosti kt(oblacnost, vyska
slunce). kt se kalibruje empiricky z historie (do `until`), protoze vzorec
Kasten-Czeplak vztah u stredni oblacnosti nezachyti (merena stanice: jasno
kt ~0,64, zatazeno ~0,15). Zavislost na vysce slunce je nutna: s kt jen podle
oblacnosti mel odvozeny osvit spatny denni i rocni chod (leto rano +50-70,
v poledne -70 az -100 W/m2; v zime v prumeru dvojnasobek skutecnosti).

Teplota a vitr: predpoved ma proti stanici setrvacny bias (teplota v lete +1
az +1,7 C, autokorelace denni chyby 0,59; vitr +1,7 m/s, coz model cte jako
chladnejsi pocitovou teplotu) — `debias` ho odecita klouzave po hodinach dne
z poslednich BIAS_DAYS dni znamych pri vydani.

Sero pres den (clen baze, src/sun.py) je nelinearni funkce osvitu: predpoved
dava jen stredni osvit pri dane oblacnosti, a sero stredniho osvitu neni
stredni sero (tmave dny by se nikdy nepredpovedely). Proto se pro kazdy bin
oblacnosti kalibruji i **kvantily** kt a sero se pocita jako prumer pres ne —
stredni hodnota nelinearniho clenu pres rozdeleni chyby predpovedi.
"""
from pathlib import Path

import numpy as np
import pandas as pd

from scipy.interpolate import RegularGridInterpolator

import sun
import workspace
from etl import TZ
from etl import _timestamps as etl_timestamps

FILE = workspace.ANALYZA / "ArchivMeteo.xlsx"   # nebo ArchivMeteo.csv vedle nej (ma prednost xlsx)
KT_BINS = np.arange(0, 101, 10)        # hrany binu oblacnosti [%]
EL_BINS = np.array([0.0, 0.1, 0.2, 0.3, 0.45, 0.6, 0.75, 0.9])   # hrany binu sin(vysky slunce)
KT_QUANTILES = np.arange(0.05, 1.0, 0.1)   # kvantily kt v binu (stredy decilu)
MIN_CELL = 30                          # min. hodin v bunce (jinak hodnota jen podle oblacnosti)
MIN_CLEAR_SKY = 10.0                   # [W/m2] hodiny s nizsim osvitem za jasne oblohy se nekalibruji
BIAS_DAYS = 30                         # okno klouzaveho biasu teploty a vetru
BIAS_COLS = ("temp", "wind")


def load_archive() -> pd.DataFrame:
    """Hodinova predpoved na UTC ose (hodiny po rade od mistni pulnoci).

    Zdroj: ArchivMeteo.xlsx, nebo ArchivMeteo.csv ve formatu vstupnich dat
    (oddelovac `;`, desetinna carka, mistni cas) se sloupci cas, teplota,
    oblacnost [%], rychlost vetru v tomto poradi."""
    if FILE.exists():
        a = pd.read_excel(FILE)
        a.columns = ["ts", "temp", "obl", "wind"]
        a = a.dropna(subset=["temp"])
        ts = pd.to_datetime(a["ts"])
    else:
        a = pd.read_csv(FILE.with_suffix(".csv"), sep=";", decimal=",")
        a.columns = ["ts", "temp", "obl", "wind"]
        a = a.dropna(subset=["temp"])
        ts = etl_timestamps(a["ts"])
    day = ts.dt.normalize()
    utc = day.dt.tz_localize(TZ).dt.tz_convert("UTC") + pd.to_timedelta(a.groupby(day).cumcount(), unit="h")
    out = a[["temp", "obl", "wind"]].astype(float)
    out.index = pd.DatetimeIndex(utc)
    return out[~out.index.duplicated()]


def _cos_zenith(index_utc: pd.DatetimeIndex) -> np.ndarray:
    return np.clip(np.sin(np.radians(sun.elevation(index_utc))), 0, None)


def clear_sky(index_utc: pd.DatetimeIndex) -> np.ndarray:
    """Osvit za jasne oblohy [W/m2] (Haurwitz) ve stredu intervalu."""
    cz = _cos_zenith(index_utc)
    return np.where(cz > 0.01, 1098 * cz * np.exp(-0.057 / np.maximum(cz, 0.01)), 0.0)


def _history(arch: pd.DataFrame, actual: pd.DataFrame, until: pd.Timestamp):
    """Hodiny pred `until` se sluncem nad obzorem, jejich bin oblacnosti a
    bin vysky slunce."""
    act = actual["sun"].groupby(actual.index.floor("h")).mean()
    j = arch[["obl"]].join(act, how="inner")
    j = j[j.index < until]
    mid = j.index + pd.Timedelta(minutes=30)
    j["cs"] = clear_sky(mid)
    j["cz"] = _cos_zenith(mid)
    j = j[j.cs > MIN_CLEAR_SKY]
    b = np.clip(np.digitize(j.obl, KT_BINS[1:-1]), 0, len(KT_BINS) - 2)
    e = np.clip(np.digitize(j.cz, EL_BINS[1:-1]), 0, len(EL_BINS) - 2)
    return j, b, e


def _cells(fn, b: np.ndarray, e: np.ndarray) -> np.ndarray:
    """fn(maska) pro kazdou bunku (oblacnost x vyska slunce); ridke bunky se
    doplni interpolaci mezi sousednimi biny oblacnosti pri stejne vysce slunce
    (vyska slunce urcuje kt vic nez oblacnost)."""
    n_o, n_e = len(KT_BINS) - 1, len(EL_BINS) - 1
    cells = [[fn((b == k) & (e == m)) if ((b == k) & (e == m)).sum() >= MIN_CELL else None
              for m in range(n_e)] for k in range(n_o)]
    for m in range(n_e):
        ok = [k for k in range(n_o) if cells[k][m] is not None]
        if not ok:      # vyska slunce bez dat: cely sloupec podle oblacnosti
            for k in range(n_o):
                cells[k][m] = fn(b == k)
            continue
        vals = np.array([cells[k][m] for k in ok], float)
        for k in range(n_o):
            if cells[k][m] is None:
                cells[k][m] = (np.array([np.interp(k, ok, vals[:, q]) for q in range(vals.shape[1])])
                               if vals.ndim == 2 else float(np.interp(k, ok, vals)))
    return np.array(cells, float)


def calibrate(arch: pd.DataFrame, actual: pd.DataFrame, until: pd.Timestamp) -> np.ndarray:
    """kt v bunkach (oblacnost x vyska slunce) z hodin pred `until` (bez
    uniku z budoucnosti)."""
    j, b, e = _history(arch, actual, until)
    sun_a, cs = j.sun.to_numpy(), j.cs.to_numpy()
    # pomer souctu (vahy podle osvitu za jasne oblohy), ne prumer pomeru
    return _cells(lambda m: sun_a[m].sum() / cs[m].sum(), b, e)


def calibrate_quantiles(arch: pd.DataFrame, actual: pd.DataFrame, until: pd.Timestamp) -> np.ndarray:
    """Kvantily kt v bunkach (oblacnost x vyska slunce x kvantil) — rozdeleni
    skutecne jasnosti pri dane predpovedi oblacnosti."""
    j, b, e = _history(arch, actual, until)
    ratio = (j.sun / j.cs).to_numpy()
    return _cells(lambda m: np.quantile(ratio[m], KT_QUANTILES), b, e)


def _lookup(table: np.ndarray, obl: np.ndarray, cz: np.ndarray) -> np.ndarray:
    """Bilinearni interpolace tabulky mezi stredy bunek (za okraji drzi)."""
    c_o = (KT_BINS[:-1] + KT_BINS[1:]) / 2
    c_e = (EL_BINS[:-1] + EL_BINS[1:]) / 2
    f = RegularGridInterpolator((c_o, c_e), table)
    return f(np.column_stack([np.clip(obl, c_o[0], c_o[-1]), np.clip(cz, c_e[0], c_e[-1])]))


def to_15min(arch: pd.DataFrame, index_15: pd.DatetimeIndex, kt: np.ndarray,
             kt_q: np.ndarray | None = None) -> pd.DataFrame:
    """Predpoved na 15min ose: teplota, vitr a oblacnost linearne mezi stredy
    hodin, osvit = jasna obloha (stred intervalu) x kt(oblacnost, vyska slunce).

    S kt_q (calibrate_quantiles) navic sloupec `sero_den` = stredni hodnota
    exp(-osvit / I0) pres rozdeleni kt pri dane oblacnosti (bez faktoru tmy)."""
    mid_h = arch.index + pd.Timedelta(minutes=30)
    mid_q = index_15 + pd.Timedelta(minutes=7.5)
    x, xq = mid_h.asi8.astype(float), mid_q.asi8.astype(float)
    out = pd.DataFrame(index=index_15)
    for c in ("temp", "wind", "obl"):
        out[c] = np.interp(xq, x, arch[c].to_numpy())
    obl, cz, cs = out["obl"].to_numpy(), _cos_zenith(mid_q), clear_sky(mid_q)
    out["sun"] = cs * _lookup(kt, obl, cz)
    if kt_q is not None:
        out["sero_den"] = np.exp(-cs[:, None] * _lookup(kt_q, obl, cz) / sun.GLOOM_I0).mean(axis=1)
    return out


def debias_history(fc: pd.DataFrame, actual: pd.DataFrame,
                   cols: tuple[str, ...] = BIAS_COLS) -> pd.DataFrame:
    """Cela historie predpovedi zbavena klouzaveho biasu — pro kazdy den bias z
    BIAS_DAYS dni koncicich predvcerejskem (zname pri vydani predchozi den).
    Pro uceni modelu na predpovezenem pocasi (backtest --trenink=predpoved-bias)."""
    loc = fc.index.tz_convert(TZ)
    day, hour = pd.DatetimeIndex(loc.date), loc.hour
    out = pd.DataFrame(index=fc.index)
    for c in cols:
        err = (fc[c] - actual[c].reindex(fc.index)).to_numpy()
        tab = pd.Series(err).groupby([day, hour]).mean().unstack()
        bias = tab.rolling(BIAS_DAYS, min_periods=7).mean().shift(2).fillna(0.0)
        b = bias.stack().reindex(pd.MultiIndex.from_arrays([day, hour])).to_numpy()
        out[c] = fc[c].to_numpy() - np.nan_to_num(b)
    if "wind" in out:
        out["wind"] = out["wind"].clip(lower=0.0)
    return out


def debias(fc: pd.DataFrame, actual: pd.DataFrame, cutoff: pd.Timestamp,
           cols: tuple[str, ...] = BIAS_COLS) -> pd.DataFrame:
    """Predpoved teploty a vetru od `cutoff` dal zbavena klouzaveho biasu.

    Bias = prumer (predpoved - skutecnost) po hodinach dne (mistni cas) za
    poslednich BIAS_DAYS dni pred `cutoff` — vse zname v okamziku vydani.
    fc = predpoved na 15min ose (to_15min), actual = namerena data.
    """
    past = fc.index[(fc.index < cutoff) & (fc.index >= cutoff - pd.Timedelta(days=BIAS_DAYS))]
    fut = fc.index[fc.index >= cutoff]
    hour_p, hour_f = past.tz_convert(TZ).hour, fut.tz_convert(TZ).hour
    out = pd.DataFrame(index=fut)
    for c in cols:
        err = fc.loc[past, c] - actual[c].reindex(past)
        bias = err.groupby(hour_p).mean().reindex(range(24)).fillna(0.0)
        out[c] = fc.loc[fut, c] - bias.reindex(hour_f).to_numpy()
    if "wind" in out:
        out["wind"] = out["wind"].clip(lower=0.0)
    return out
