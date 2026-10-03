"""Nacteni vstupnich dat: 15min spotreba portfolia a meteo (CET/CEST -> UTC)."""
import time
from pathlib import Path

import numpy as np
import pandas as pd

import workspace

DATA_DIR = workspace.DATA
TZ = "Europe/Prague"
READ_RETRIES = 4   # pokusu o nacteni, kdyz se vstupni soubor pri cteni meni
READ_WAIT = 2.0    # [s] cekani mezi pokusy
SRC_DIR = Path(__file__).resolve().parent
# mezipamet nactenych dat (v Data/, mimo git): plati, dokud se nezmeni vstupni
# soubory, konfigurace kalendare ani kod, ktery je zpracovava
CACHE = DATA_DIR / ".etl_cache.pkl"
_CACHE_DEPS = [DATA_DIR / "baseload.csv", DATA_DIR / "Meteo15.csv",
               workspace.config_file("kalendar.yaml"), workspace.config_file("model.yaml"),
               SRC_DIR / "etl.py", SRC_DIR / "kalendar.py", SRC_DIR / "sun.py"]


def _read_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=";", decimal=",")
    ts = _timestamps(df["timestamp"])
    # data jsou v lokalnim case: podzimni duplicitni hodina se rozlisi poradim
    ts = ts.dt.tz_localize(TZ, ambiguous="infer")
    df["timestamp"] = ts.dt.tz_convert("UTC")
    return df.set_index("timestamp").sort_index()


def _timestamps(s: pd.Series) -> pd.Series:
    """Casove znacky exportu: "01.01.2022" (pulnoc bez casu) nebo
    "01.01.2022 0:15". S pevnym formatem ~0,8 s na soubor; format="mixed"
    parsuje po prvcich (12-17 s) a zustava jen jako zaloha pro jiny zapis."""
    try:
        return pd.to_datetime(s.where(s.str.len() > 10, s + " 0:00"), format="%d.%m.%Y %H:%M")
    except ValueError:
        return pd.to_datetime(s, format="mixed", dayfirst=True)


def step_hours(df: pd.DataFrame) -> float:
    """Casovy krok dat [h] (15 min = 0,25; hodinova agregace = 1)."""
    idx = df.index[:200]
    return float(np.median((idx[1:] - idx[:-1]).total_seconds())) / 3600.0


def hourly(df: pd.DataFrame) -> pd.DataFrame:
    """Hodinova agregace vystupu load(): spojite veliciny prumerem (spotreba,
    meteo, tma), kalendarni sloupce prvni hodnotou. Pro fit tvaru (nelinearni
    parametry nepotrebuji ctvrthodinove rozliseni, 4x mene radku)."""
    num = [c for c in ("baseload", "sun", "wind", "temp", "tma", "sero") if c in df.columns]
    g = df.groupby(df.index.floor("h"))
    out = g[[c for c in df.columns if c not in num]].first()
    out[num] = g[num].mean()
    return out[df.columns]


def _cache_key() -> tuple:
    return tuple((str(p), p.stat().st_mtime_ns, p.stat().st_size) for p in _CACHE_DEPS)


def load(with_local_time: bool = True) -> pd.DataFrame:
    """Spotreba + meteo na spolecne UTC ose; volitelne sloupce lokalniho kalendare.

    Plny vystup (s kalendarem) se uklada do mezipameti — nacteni ~10 s ->
    <1 s; backtest ho nacita v kazdem workeru.
    """
    if with_local_time:
        key = _cache_key()
        try:
            cached_key, cached = pd.read_pickle(CACHE)
            if cached_key == key:
                return cached
        except Exception:
            pass
        df = _load_stable(key)
        key = _cache_key()
        tmp = CACHE.with_suffix(f".{np.random.randint(1 << 30)}.tmp")
        pd.to_pickle((key, df), tmp)
        tmp.replace(CACHE)             # atomicky — workery mohou zapisovat soucasne
        return df
    return _load(False)


def _load_stable(key: tuple) -> pd.DataFrame:
    """_load odolny proti soubehu s aktualizaci vstupnich souboru: kdyz se
    soubor behem cteni zmeni nebo je zapsany jen zcasti, cte se znovu."""
    for attempt in range(READ_RETRIES):
        try:
            df = _load(True)
        except Exception:
            if attempt == READ_RETRIES - 1:
                raise
        else:
            now = _cache_key()
            if now == key:
                return df
            key = now
        time.sleep(READ_WAIT)
    return df


def previous() -> pd.DataFrame | None:
    """Data z minuleho nacteni (mezipamet), i kdyz uz neplati — pro zjisteni,
    co se ve vstupnich souborech zmenilo. Volat pred load()."""
    try:
        return pd.read_pickle(CACHE)[1]
    except Exception:
        return None


def revisions(old: pd.DataFrame | None, new: pd.DataFrame, tol: float = 1e-6) -> pd.DataFrame:
    """Zpetne zmenene hodnoty (zpresneni mereni): radky spolecne obema nactenim,
    kde se spotreba nebo meteo lisi; sloupce = zmena (nova - stara)."""
    cols = ["baseload", "temp", "sun", "wind"]
    if old is None:
        return pd.DataFrame(columns=cols, index=pd.DatetimeIndex([], tz="UTC"))
    idx = old.index.intersection(new.index)
    d = new.loc[idx, cols] - old.loc[idx, cols]
    return d[(d.abs() > tol).any(axis=1)]


def add_calendar(df: pd.DataFrame) -> pd.DataFrame:
    """Doplni sloupce mistniho kalendare (typ dne, svatky, mosty, prazdniny,
    zvlastni obdobi, tma) k ramci s UTC indexem — pro data i pro budouci osu
    dlouhodobe predikce. Meni df na miste a vraci ho."""
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
    # poradi dne v okne prubehu leta (config prazdniny.prubeh), mimo okno -1
    df["leto_den"] = kalendar.summer_position(dates)
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


def _load(with_local_time: bool) -> pd.DataFrame:
    bl = _read_csv(DATA_DIR / "baseload.csv").rename(columns={"baseline": "baseload"})
    mt = _read_csv(DATA_DIR / "Meteo15.csv")
    # duplicitni timestampy (artefakty exportu kolem zmen casu) -> prumer
    mt = mt.groupby(level=0).mean()
    # za koncem realnych dat jsou v obou souborech jen radky sablony: meteo
    # prazdne, baseline zaporny sum kolem nuly (realna spotreba je >= 300).
    # Konec se urcuje pro kazdy soubor zvlast — mereni pocasi a spotreby
    # nemusi prichazet stejne rychle; data konci tim drivejsim z nich.
    valid_end = mt["temp"].last_valid_index()
    mt = mt.loc[:valid_end]
    bl = bl.loc[:bl.index[bl["baseload"] > 0].max()]
    # doplneni chybejicich intervalu na plnou 15min osu, kratke diry interpolovat
    full = pd.date_range(mt.index.min(), mt.index.max(), freq="15min", tz="UTC")
    mt = mt.reindex(full).interpolate(limit=4)
    df = bl.join(mt, how="inner").dropna()
    if with_local_time:
        add_calendar(df)
        import sun
        # sero pres den (tmava obloha za dne) pro clen denniho sviceni
        df["sero"] = sun.gloom(df["tma"], df["sun"])
    return df


def _bridges(holidays: set) -> set:
    """Mosty podle config/kalendar.yaml (obal kvuli starsim volanim)."""
    import kalendar
    return kalendar.bridges(holidays)


def _cz_holidays(y0: int, y1: int) -> set:
    """Statni svatky podle config/kalendar.yaml (obal kvuli starsim volanim)."""
    import kalendar
    return kalendar.holidays(y0, y1)
