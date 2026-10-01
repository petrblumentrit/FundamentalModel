"""Nacteni vstupnich dat: 15min spotreba portfolia a meteo (CET/CEST -> UTC)."""
from pathlib import Path

import numpy as np
import pandas as pd

DATA_DIR = Path(__file__).resolve().parent.parent / "Data"
TZ = "Europe/Prague"


def _read_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=";", decimal=",")
    ts = pd.to_datetime(df["timestamp"], format="mixed", dayfirst=True)
    # data jsou v lokalnim case: podzimni duplicitni hodina se rozlisi poradim
    ts = ts.dt.tz_localize(TZ, ambiguous="infer")
    df["timestamp"] = ts.dt.tz_convert("UTC")
    return df.set_index("timestamp").sort_index()


def load(with_local_time: bool = True) -> pd.DataFrame:
    """Spotreba + meteo na spolecne UTC ose; volitelne sloupce lokalniho kalendare."""
    bl = _read_csv(DATA_DIR / "baseload.csv").rename(columns={"baseline": "baseload"})
    mt = _read_csv(DATA_DIR / "Meteo15.csv")
    # duplicitni timestampy (artefakty exportu kolem zmen casu) -> prumer
    mt = mt.groupby(level=0).mean()
    # realna data konci tam, kde konci meteo; dal jsou v obou souborech
    # jen prazdne radky sablony (baseline sum kolem nuly, meteo NaN)
    valid_end = mt["temp"].last_valid_index()
    mt = mt.loc[:valid_end]
    # doplneni chybejicich intervalu na plnou 15min osu, kratke diry interpolovat
    full = pd.date_range(mt.index.min(), mt.index.max(), freq="15min", tz="UTC")
    mt = mt.reindex(full).interpolate(limit=4)
    df = bl.join(mt, how="inner").dropna()
    if with_local_time:
        loc = df.index.tz_convert(TZ)
        df["date_local"] = loc.date
        df["tod"] = loc.hour + loc.minute / 60  # cas dne [h]
        df["dow"] = loc.dayofweek  # 0=Po
        # typ dne = role v bloku volna podle toho, zda je volno dnes a zitra:
        #   0 prace->prace (Po-Ct), 1 prace->volno (patek; i ctvrtek pred
        #   patecnim svatkem, mosty), 2 volno->volno (sobota; i svatek pred
        #   vikendem), 3 volno->prace (nedele; i Velikonocni pondeli).
        # Svatky tak prebiraji profily naucene z ~240 vikendu misto nedele pro
        # vsechny; pro bezne tydny se nic nemeni. Volno vcera nese ranni clen
        # po volnu (po_volnu).
        holidays = _cz_holidays(loc.year.min(), loc.year.max() + 1)
        dates = pd.Series(loc.date, index=df.index)
        days = pd.to_datetime(dates)
        off = ((df["dow"] >= 5) | dates.isin(holidays)).to_numpy()
        nxt = days + pd.Timedelta(days=1)
        off_next = ((nxt.dt.dayofweek >= 5) | nxt.dt.date.isin(holidays)).to_numpy()
        df["daytype"] = np.select([~off & ~off_next, ~off & off_next, off & off_next],
                                  [0, 1, 2], default=3)
        # mosty: pracovni den sevreny mezi svatkem a vikendem/svatkem (role 1)
        df["most"] = dates.isin(_bridges(holidays)).to_numpy()
        # letni prazdniny (dominantni skolni volno)
        df["prazdniny"] = loc.month.isin([7, 8])
        # predchozi kalendarni den volny (vikend/svatek) — rezim noci se meni az
        # behem rana, ne o pulnoci (pondelni noc je jeste vikendova)
        prev = pd.Series(pd.to_datetime(dates) - pd.Timedelta(days=1), index=df.index)
        df["po_volnu"] = ((prev.dt.dayofweek >= 5) | prev.dt.date.isin(holidays)).to_numpy()
        # tma (astronomicka, deterministicka) pro clen osvetleni
        import sun
        df["tma"] = sun.darkness_15min(df.index)
    return df


def _bridges(holidays: set) -> set:
    """Pracovni dny sevrene mezi svatkem a vikendem/svatkem."""
    from datetime import timedelta

    out = set()
    for h in holidays:
        for d in (h - timedelta(days=1), h + timedelta(days=1)):
            if d.weekday() >= 5 or d in holidays:
                continue
            prev_off = (d - timedelta(days=1)).weekday() >= 5 or (d - timedelta(days=1)) in holidays
            next_off = (d + timedelta(days=1)).weekday() >= 5 or (d + timedelta(days=1)) in holidays
            if prev_off and next_off:
                out.add(d)
    return out


def _cz_holidays(y0: int, y1: int) -> set:
    """Ceske statni svatky vcetne pohyblivych (Velky patek, Velikonocni pondeli)."""
    from datetime import date, timedelta

    fixed = [(1, 1), (5, 1), (5, 8), (7, 5), (7, 6), (9, 28), (10, 28), (11, 17), (12, 24), (12, 25), (12, 26)]
    out = set()
    for y in range(y0, y1 + 1):
        out.update(date(y, m, d) for m, d in fixed)
        e = _easter(y)
        out.add(e - timedelta(days=2))  # Velky patek
        out.add(e + timedelta(days=1))  # Velikonocni pondeli
    return out


def _easter(year: int):
    """Anonymni gregoriansky algoritmus."""
    from datetime import date

    a, b, c = year % 19, year // 100, year % 100
    d, e = b // 4, b % 4
    f = (b + 8) // 25
    g = (b - f + 1) // 3
    h = (19 * a + b - d - g + 15) % 30
    i, k = c // 4, c % 4
    m = (32 + 2 * e + 2 * i - h - k) % 7
    n = (a + 11 * h + 22 * m) // 451
    month = (h + m - 7 * n + 114) // 31
    day = (h + m - 7 * n + 114) % 31 + 1
    return date(year, month, day)
