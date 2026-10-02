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
        import kalendar
        holidays = kalendar.holidays(loc.year.min() - 1, loc.year.max() + 1)
        dates = pd.Series(loc.date, index=df.index)
        days = pd.to_datetime(dates)
        off = ((df["dow"] >= 5) | dates.isin(holidays)).to_numpy()
        nxt = days + pd.Timedelta(days=1)
        off_next = ((nxt.dt.dayofweek >= 5) | nxt.dt.date.isin(holidays)).to_numpy()
        df["daytype"] = np.select([~off & ~off_next, ~off & off_next, off & off_next],
                                  [0, 1, 2], default=3)
        # mosty: 1-2 pracovni dny sevrene volnem (config/kalendar.yaml)
        runs = kalendar.bridge_runs(holidays)
        df["most"] = dates.isin(set(runs)).to_numpy()
        df["most_delka"] = dates.map(runs).fillna(0).astype(int).to_numpy()
        # letni prazdniny (dominantni skolni volno)
        df["prazdniny"] = kalendar.summer_holidays(dates)
        # zvlastni obdobi (Vanoce ...): priznak pro kazdou skupinu z configu
        for name, flag in kalendar.period_flags(dates, holidays).items():
            df["obd_" + name] = flag
        # predchozi kalendarni den volny (vikend/svatek) — rezim noci se meni az
        # behem rana, ne o pulnoci (pondelni noc je jeste vikendova)
        prev = pd.Series(pd.to_datetime(dates) - pd.Timedelta(days=1), index=df.index)
        df["po_volnu"] = ((prev.dt.dayofweek >= 5) | prev.dt.date.isin(holidays)).to_numpy()
        # tma (astronomicka, deterministicka) pro clen osvetleni
        import sun
        df["tma"] = sun.darkness_15min(df.index)
    return df


def _bridges(holidays: set) -> set:
    """Mosty podle config/kalendar.yaml (obal kvuli starsim volanim)."""
    import kalendar
    return kalendar.bridges(holidays)


def _cz_holidays(y0: int, y1: int) -> set:
    """Statni svatky podle config/kalendar.yaml (obal kvuli starsim volanim)."""
    import kalendar
    return kalendar.holidays(y0, y1)
