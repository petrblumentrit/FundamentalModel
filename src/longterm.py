"""Dlouhodoba predikce — az rok dopredu z modelu odhadnuteho do `cutoff`.

Za horizontem predpovedi pocasi (zatim od konce dat) se pocasi dosazuje dvema
zpusoby (prepinac `weather`):

- **normal**: jedna draha — klimatologicky prumer teploty, osvitu a vetru podle
  dne v roce a casu dne (z historie pocasi). Vysledkem je 15min profil. Hladka
  draha podhodnocuje nelinearni cleny (chlazeni, prechodna obdobi, spicky);
  sero pres den se proto bere jako **prumer skutecneho sera**, ne sero
  prumerneho osvitu.
- **scenare**: model se spusti na skutecnem pocasi kazdeho historickeho roku
  (souvisle okno posunute o cele roky v UTC — poloha slunce sedi, dynamika
  pocasi a setrvacnost zustavaji realne). Aby z kratke historie nebyly jen
  3-4 drahy, bere se kazdy rok i s posunem o OFFSETS dni; teplota se pritom
  prepocte o rozdil normalu (anomalie zdrojoveho dne + normal ciloveho dne).
  Predikce = prumer pres scenare, pasmo = kvantily 10-90 % (pri mene nez
  MIN_QUANTILE scenarich rozpeti). Scenare popisuji jen nejistotu pocasi a
  ty z tehoz roku nejsou nezavisle — samotne pokryvaly v dopredne validaci
  57-69 % skutecnosti misto 80 %. Pasmo se proto rozsiruje o chybu modelu a
  trendu (`model_sd`, z dopredne validace se skutecnym pocasim), scitanou s
  rozptylem pocasi kvadraticky.

Historie pocasi: `Data/meteo_historie.csv` (stejny format jako Meteo15.csv,
libovolne dlouha rada ze stejne stanice), jinak meteo z dat portfolia.

Pomale slozky za koncem treninku (prepinac `trend`):

- **posledni**: pokracuje posledni trend — uroven baze, prirustek topne
  citlivosti a kapacita FVE rostou tempem poslednich TREND_DAYS dni;
- **drzet**: posledni hodnota (jako v predikci D+1);
- vyber slozek carkou (`uroven`, `topeni`, `fve`) a tlumeni cislem za
  dvojteckou, napr. `topeni,fve` nebo `posledni:0.5`.

Kalendarni a behavioralni cast (typy dne, svatky, mosty, Vanoce, prubeh leta,
osvetleni podle slunce) je deterministicka a pocita se pro budouci dny stejne
jako pro minule.
"""
import numpy as np
import pandas as pd

import base
import cooling
import etl
import heating
import pv
import sun

WEATHER_COLS = ("temp", "sun", "wind")
HISTORY_FILE = etl.DATA_DIR / "meteo_historie.csv"
HISTORY_DAYS = 60     # [dny] skutecna data pred cutoffem pro nabeh setrvacnostnich filtru
TREND_DAYS = 365      # [dny] okno "posledniho trendu"
NORMAL_WINDOW = 15    # [dny] vyhlazeni klimatologie pres den v roce
MIN_QUANTILE = 10     # od kolika scenaru ma smysl kvantilove pasmo
OFFSETS = (-14, -7, 0, 7, 14)   # [dny] posuny scenaru kolem celych let
TREND_PARTS = ("uroven", "topeni", "fve")
Z80 = 1.2816          # kvantil normalniho rozdeleni pro pasmo 10-90 %
# chyba modelu a trendu rok dopredu pri znamem pocasi (explore/dlouhodoba_validace.py:
# RMSE 14-17 v hodinovem kroku, RMS mesicniho biasu 5-7 mimo fold s obratem trendu)
MODEL_SD = 15.0
MODEL_SD_MONTH = 8.0
STEP = pd.Timedelta(minutes=15)


def weather_history(df: pd.DataFrame) -> pd.DataFrame:
    """Rada teploty, osvitu a vetru na 15min UTC ose — z HISTORY_FILE, pokud
    existuje (doplnena o novejsi meteo z dat portfolia), jinak z dat portfolia."""
    own = df[list(WEATHER_COLS)]
    if not HISTORY_FILE.exists():
        return own
    ext = etl._read_csv(HISTORY_FILE)[list(WEATHER_COLS)].astype(float)
    ext = ext.groupby(level=0).mean()
    full = pd.date_range(ext.index.min(), ext.index.max(), freq=STEP, tz="UTC")
    ext = ext.reindex(full).interpolate(limit=4).dropna()
    return pd.concat([ext[ext.index < own.index.min()], own])


def frame(df: pd.DataFrame, cutoff: pd.Timestamp, days: int) -> pd.DataFrame:
    """Osa HISTORY_DAYS pred cutoffem + `days` za nim: skutecna data, kde jsou
    (za cutoffem jen pro validaci), kalendar pro vsechny radky."""
    idx = pd.date_range(cutoff - pd.Timedelta(days=HISTORY_DAYS),
                        cutoff + pd.Timedelta(days=days) - STEP, freq=STEP, tz="UTC")
    fr = df[["baseload", *WEATHER_COLS]].reindex(idx)
    past = fr.index < cutoff
    cols = list(WEATHER_COLS)
    fr.loc[past, cols] = fr.loc[past, cols].interpolate(limit_direction="both")
    etl.add_calendar(fr)
    fr["sero"] = sun.gloom(fr["tma"], fr["sun"])
    return fr


def _doy_slot(index: pd.DatetimeIndex) -> tuple[np.ndarray, np.ndarray]:
    """Den v roce (366 -> 365) a ctvrthodina dne v UTC — slunecni cas, zmena
    mistniho casu klimatologii neposouva."""
    return np.minimum(index.dayofyear.to_numpy(), 365), index.hour.to_numpy() * 4 + index.minute.to_numpy() // 15


def climatology(hist: pd.DataFrame) -> dict[str, np.ndarray]:
    """Prumer podle (den v roce, cas dne) vyhlazeny cirkularne pres dny v roce;
    pro kazdou velicinu matice 365 x 96. Krome pocasi i prumerne sero."""
    h = hist.copy()
    h["sero"] = sun.gloom(sun.darkness_15min(h.index), h["sun"])
    doy, slot = _doy_slot(h.index)
    out = {}
    for c in (*WEATHER_COLS, "sero"):
        mat = (pd.Series(h[c].to_numpy(), index=[doy, slot]).groupby(level=[0, 1]).mean()
               .unstack().reindex(index=range(1, 366), columns=range(96)))
        mat = mat.interpolate(limit_direction="both")
        tri = pd.concat([mat, mat, mat], ignore_index=True)
        out[c] = tri.rolling(NORMAL_WINDOW, center=True, min_periods=1).mean().iloc[365:730].to_numpy()
    return out


def fill_normal(fr: pd.DataFrame, cutoff: pd.Timestamp, clim: dict[str, np.ndarray]) -> pd.DataFrame:
    """Radky od cutoffu: normalove pocasi a prumerne sero."""
    out = fr.copy()
    fut = out.index >= cutoff
    doy, slot = _doy_slot(out.index[fut])
    for c in (*WEATHER_COLS, "sero"):
        out.loc[fut, c] = clim[c][doy - 1, slot]
    return out


def scenario_shifts(hist: pd.DataFrame, cutoff: pd.Timestamp, days: int) -> list[int]:
    """Posuny [dny] zpet (cele roky + OFFSETS), pro ktere historie pred
    cutoffem pokryva cele okno predikce."""
    out, k = [], 1
    while cutoff - pd.Timedelta(days=int(round(365.25 * k)) + min(OFFSETS)) >= hist.index.min():
        for off in OFFSETS:
            shift = int(round(365.25 * k)) - off
            start = cutoff - pd.Timedelta(days=shift)
            end = start + pd.Timedelta(days=days) - STEP
            if start < hist.index.min() or end >= cutoff:
                continue
            window = hist[list(WEATHER_COLS)].reindex(pd.date_range(start, end, freq=STEP))
            if window.notna().mean().min() > 0.98:
                out.append(shift)
        k += 1
    return out


def fill_scenario(fr: pd.DataFrame, cutoff: pd.Timestamp, hist: pd.DataFrame, shift: int,
                  clim: dict[str, np.ndarray] | None = None) -> pd.DataFrame:
    """Radky od cutoffu: pocasi z historie posunute o `shift` dni (stejny cas
    UTC). S klimatologii se teplota prepocte o rozdil normalu ciloveho a
    zdrojoveho dne (posun mimo cele roky jinak vnasi sezonni chod)."""
    out = fr.copy()
    fut = out.index[out.index >= cutoff]
    src_idx = fut - pd.Timedelta(days=shift)
    src = hist[list(WEATHER_COLS)].reindex(src_idx).interpolate(limit_direction="both")
    if clim is not None:
        (d_t, s_t), (d_s, s_s) = _doy_slot(fut), _doy_slot(src_idx)
        src["temp"] = src["temp"].to_numpy() + clim["temp"][d_t - 1, s_t] - clim["temp"][d_s - 1, s_s]
    out.loc[fut, list(WEATHER_COLS)] = src.to_numpy()
    out["sero"] = sun.gloom(out["tma"], out["sun"])
    return out


def trend_weights(trend: str) -> dict[str, float]:
    """`posledni` | `drzet` | vyber slozek carkou; tlumeni za dvojteckou."""
    name, _, factor = trend.partition(":")
    f = float(factor) if factor else 1.0
    parts = TREND_PARTS if name == "posledni" else () if name == "drzet" else tuple(name.split(","))
    unknown = set(parts) - set(TREND_PARTS)
    if unknown:
        raise ValueError(f"nezname slozky trendu: {sorted(unknown)}")
    return {p: (f if p in parts else 0.0) for p in TREND_PARTS}


def band(total: np.ndarray, model_sd: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
    """Pasmo 10-90 % z matice scenaru (radky = cas, sloupce = scenare),
    rozsirene o chybu modelu (kvadraticky soucet s rozptylem pocasi)."""
    mean = total.mean(axis=1)
    if total.shape[1] >= MIN_QUANTILE:
        lo, hi = np.quantile(total, [0.1, 0.9], axis=1)
    else:
        lo, hi = total.min(axis=1), total.max(axis=1)
    m = Z80 * model_sd
    return mean - np.sqrt((mean - lo) ** 2 + m**2), mean + np.sqrt((hi - mean) ** 2 + m**2)


def scenario_columns(f: pd.DataFrame) -> list[str]:
    return [c for c in f.columns if c.startswith("scenar_")]


def monthly(f: pd.DataFrame, model_sd: float = MODEL_SD_MONTH) -> pd.DataFrame:
    """Mesicni prumery (mistni cas) slozek a predikce; u scenaru pasmo z
    mesicnich prumeru jednotlivych drah (uzsi nez prumer ctvrthodinoveho pasma)."""
    key = f.index.tz_convert(etl.TZ).strftime("%Y-%m")
    out = f[["baze", "topeni", "chlazeni", "fve", "predikce"]].groupby(key).mean()
    sc = scenario_columns(f)
    if sc:
        out["pasmo_dolni"], out["pasmo_horni"] = band(f[sc].groupby(key).mean().to_numpy(), model_sd)
    return out


def _x_days(fr: pd.DataFrame, t0) -> np.ndarray:
    """Poloha na ose urovnove spline [dny od t0] — stejne jako base.design."""
    return ((pd.to_datetime(fr["date_local"]) - pd.Timestamp(t0)).dt.days.to_numpy(float)
            + fr["tod"].to_numpy() / 24.0)


def _ramp_ref(span_days: float) -> float:
    """Konec posledni plne pozorovane rampy (pv.capacity_basis) [dny od t0]."""
    return float(np.floor(span_days / pv.KNOT_DAYS) * pv.KNOT_DAYS)


def _at(t0: pd.Timestamp, x_days: list[float]) -> pd.DataFrame:
    return pd.DataFrame(index=pd.DatetimeIndex([t0 + pd.Timedelta(days=x) for x in x_days]))


def components(fr: pd.DataFrame, params: tuple, trend: str = "posledni") -> pd.DataFrame:
    """Slozky predikce na ose `fr` s pravidlem pro pomale slozky za koncem treninku."""
    bp, hp, cp, pp = params
    grow = trend_weights(trend)

    # uroven baze: za koncem treninku posledni hodnota + (volitelne) posledni trend
    t0_b, span_b = pd.Timestamp(bp["t0"]).date(), bp["span_days"]
    x = _x_days(fr, t0_b)
    lvl_c = bp["coef"][base._slices(span_b)["level"]]
    lvl = base._spline_basis(x, span_b) @ lvl_c
    l_end, l_prev = base._spline_basis(np.array([span_b, span_b - TREND_DAYS]), span_b) @ lvl_c
    after = np.maximum(x - span_b, 0.0)
    rule = np.where(x > span_b, l_end + grow["uroven"] * (l_end - l_prev) / TREND_DAYS * after, lvl)
    baze = base.predict(fr, bp) - lvl + rule

    # topeni: prirustek citlivosti dk(t, cas dne) pokracuje od konce posledni rampy
    theta = np.array([hp[k] for k in heating.PARAM_NAMES], float)
    shape = heating.shape_term(fr, theta)
    k_tot = heating.k_total(fr, hp)
    if grow["topeni"] and "trend_coef" in hp:
        t0_h, span_h = pd.Timestamp(hp["trend_t0"]), hp["trend_span"]
        ref = _ramp_ref(span_h)
        R = heating.trend_basis(_at(t0_h, [ref, ref - TREND_DAYS]), t0_h, span_h)
        slope = grow["topeni"] * heating._trend_blocks(hp) @ (R[0] - R[1]) / TREND_DAYS   # po slozkach denniho tvaru
        k_tot = k_tot + (heating.trend_tod(fr["tod"].to_numpy()) @ slope) * np.maximum(pv._days(fr, t0_h) - ref, 0.0)
    topeni = k_tot * shape

    # FVE: kapacita roste tempem posledniho roku (neklesa)
    t0_p, span_p = pd.Timestamp(pp["t0"]), pp["span_days"]
    fve = pv.predict(fr, pp)
    if grow["fve"]:
        ref = _ramp_ref(span_p)
        C = pv.capacity_basis(_at(t0_p, [ref, ref - TREND_DAYS]), t0_p, span_p) @ pp["delta"]
        rate = grow["fve"] * max(float(C[0] - C[1]), 0.0) / TREND_DAYS
        w = -pv.design(fr, pp["gamma"], t0_p, span_p)[:, 0]                # osvit x ucinnost
        fve = fve + rate * np.maximum(pv._days(fr, t0_p) - ref, 0.0) * w

    out = pd.DataFrame({"baze": baze, "topeni": topeni, "chlazeni": cooling.predict(fr, cp), "fve": -fve},
                       index=fr.index)
    out["predikce"] = out.sum(axis=1)
    return out


def forecast(df: pd.DataFrame, params: tuple, cutoff: pd.Timestamp, days: int = 365,
             weather: str = "normal", trend: str = "posledni", hourly: bool = False,
             hist: pd.DataFrame | None = None, model_sd: float = MODEL_SD) -> pd.DataFrame:
    """Predikce od `cutoff` na `days` dni.

    weather = normal | scenare | skutecne (skutecne jen pro validaci na datech,
    ktera existuji). hourly = pocitat na hodinove agregaci (parametry z
    hodinoveho fitu). Vraci slozky a `predikce`; u scenaru prumer pres scenare,
    `pasmo_dolni`/`pasmo_horni` a sloupce `scenar_<posun>` jednotlivych drah.
    Historie pocasi se bere jen pred cutoffem (zadny unik z budoucnosti).
    """
    fr = frame(df, cutoff, days)
    if hist is None:
        hist = weather_history(df)
    hist = hist[hist.index < cutoff]

    def run(f: pd.DataFrame) -> pd.DataFrame:
        if hourly:
            f = etl.hourly(f)
        return components(f, params, trend)[lambda c: c.index >= cutoff]

    if weather == "skutecne":
        return run(fr)
    clim = climatology(hist)
    if weather == "normal":
        return run(fill_normal(fr, cutoff, clim))
    if weather != "scenare":
        raise ValueError(f"neznamy rezim pocasi: {weather}")
    shifts = scenario_shifts(hist, cutoff, days)
    if not shifts:
        raise ValueError("historie pocasi nepokryva ani jedno cele okno predikce")
    runs = {s: run(fill_scenario(fr, cutoff, hist, s, clim)) for s in shifts}
    cols = ["baze", "topeni", "chlazeni", "fve", "predikce"]
    out = sum(r[cols] for r in runs.values()) / len(runs)
    total = np.column_stack([r["predikce"].to_numpy() for r in runs.values()])
    out["pasmo_dolni"], out["pasmo_horni"] = band(total, model_sd)
    for s, r in runs.items():
        out[f"scenar_{s}"] = r["predikce"].to_numpy()
    return out
