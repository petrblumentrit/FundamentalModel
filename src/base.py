"""Bazalni slozka B(t) — krok 3 navrhu (NAVRH_MODELU.md).

Aditivni struktura, fit jen z mirnych dnu (mrtve teplotni pasmo, kde
nepusobi topeni ani chlazeni):

    B(t) = level(datum) + profil_daytype(cas dne) + prazdniny_korekce(cas dne)

- level: kubicka B-spline s uzly po KNOT_DAYS dnech a penalizaci druhych
  diferenci koeficientu (P-spline) — pomala sezonnost + trend; penalizace
  drzi hladky prubeh i pres obdobi bez mirnych dnu (zima, vrchol leta)
- profil: Fourierova rada (K_PROFILE harmonickych, perioda 24 h) pro kazdy
  typ dne; profil Po-Ct je bez konstanty (uroven nese spline), ostatni typy
  maji konstantni offset vuci Po-Ct
- prazdniny: konstanta + kratka Fourierova korekce (letni skolni volno)
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import BSpline

from etl import TZ

K_PROFILE = 12  # harmonickych v dennim profilu
K_PRAZ = 4      # harmonickych v prazdninove korekci
KNOT_DAYS = 90  # rozestup uzlu urovnove spline
N_DAYTYPES = 4  # 0 Po-Ct, 1 patek/most, 2 sobota, 3 nedele/svatek
# mrtve pasmo z diagnostiky (09_baze_signal): topeni dozniva ~13-14 C,
# chlazeni nastupuje uz ~18 C (drive nez navrhovy odhad 20-24 C)
MILD_BAND = (14.0, 18.0)   # denni prumer T mirneho dne
PREV_BAND = (12.5, 19.5)   # denni prumer T predchoziho dne (setrvacnost)
SMOOTH = 500.0             # vaha penalizace druhych diferenci spline koeficientu

MODEL_DIR = Path(__file__).resolve().parent.parent / "models"


def _fourier(tod: np.ndarray, k: int) -> np.ndarray:
    w = 2 * np.pi * np.outer(tod, np.arange(1, k + 1)) / 24.0
    return np.hstack([np.cos(w), np.sin(w)])


def _spline_basis(x_days: np.ndarray, span_days: float) -> np.ndarray:
    inner = np.arange(0.0, span_days + KNOT_DAYS, KNOT_DAYS)
    t = np.concatenate([np.repeat(inner[0], 3), inner, np.repeat(inner[-1], 3)])
    x = np.clip(x_days, inner[0], inner[-1] - 1e-9)
    return BSpline.design_matrix(x, t, 3).toarray()


def _n_spline(span_days: float) -> int:
    return len(np.arange(0.0, span_days + KNOT_DAYS, KNOT_DAYS)) + 2


def _slices(span_days: float) -> dict:
    """Indexy bloku koeficientu v celkovem vektoru."""
    out, i = {}, 0
    ns = _n_spline(span_days)
    out["level"] = slice(i, i + ns); i += ns
    out["dt0"] = slice(i, i + 2 * K_PROFILE); i += 2 * K_PROFILE
    for k in range(1, N_DAYTYPES):
        out[f"dt{k}"] = slice(i, i + 1 + 2 * K_PROFILE); i += 1 + 2 * K_PROFILE
    out["praz"] = slice(i, i + 1 + 2 * K_PRAZ); i += 1 + 2 * K_PRAZ
    out["total"] = i
    return out


def design(df: pd.DataFrame, t0, span_days: float) -> np.ndarray:
    x = (pd.to_datetime(df["date_local"]) - pd.Timestamp(t0)).dt.days.to_numpy(float)
    x = x + df["tod"].to_numpy() / 24.0
    tod = df["tod"].to_numpy()
    S = _spline_basis(x, span_days)
    F = _fourier(tod, K_PROFILE)
    dt = df["daytype"].to_numpy()
    blocks = [S]
    for k in range(N_DAYTYPES):
        m = (dt == k).astype(float)[:, None]
        b = F * m
        if k > 0:
            b = np.hstack([m, b])
        blocks.append(b)
    pz = df["prazdniny"].to_numpy(float)[:, None]
    blocks.append(np.hstack([pz, _fourier(tod, K_PRAZ) * pz]))
    return np.hstack(blocks)


def daily_table(df: pd.DataFrame) -> pd.DataFrame:
    d = df.groupby("date_local").agg(
        baseload=("baseload", "mean"), temp=("temp", "mean"),
        n=("baseload", "size"), daytype=("daytype", "max"),
        prazdniny=("prazdniny", "max"),
    )
    d.index = pd.to_datetime(d.index)
    return d


def mild_dates(df: pd.DataFrame) -> pd.DatetimeIndex:
    """Dny v mrtvem pasmu; predchozi den take mirny (tepelna setrvacnost)."""
    d = daily_table(df)
    ok = (d["temp"].between(*MILD_BAND)
          & d["temp"].shift(1).between(*PREV_BAND)
          & (d["n"] >= 92))  # DST dny maji 92/100 intervalu
    return d.index[ok]


def fit(df: pd.DataFrame, smooth: float = SMOOTH) -> dict:
    dates = mild_dates(df)
    sel = df[pd.to_datetime(df["date_local"]).isin(dates)]
    t0 = df["date_local"].min()
    span = float((df["date_local"].max() - t0).days) + 1.0
    X = design(sel, t0, span)
    y = sel["baseload"].to_numpy()
    sl = _slices(span)
    p = sl["total"]
    # normalni rovnice + P-spline penalizace na urovnovem bloku
    XtX = X.T @ X
    ns = sl["level"].stop - sl["level"].start
    D2 = np.diff(np.eye(ns), n=2, axis=0)
    P = np.zeros((p, p))
    P[sl["level"], sl["level"]] = smooth * (D2.T @ D2)
    P += 1e-8 * np.trace(XtX) / p * np.eye(p)  # numericka stabilizace
    coef = np.linalg.solve(XtX + P, X.T @ y)
    resid = y - X @ coef
    return {
        "coef": coef, "t0": str(t0), "span_days": span, "smooth": smooth,
        "k_profile": K_PROFILE, "k_praz": K_PRAZ, "knot_days": KNOT_DAYS,
        "n_mild_days": len(dates), "rmse": float(np.sqrt(np.mean(resid**2))),
        "mae": float(np.mean(np.abs(resid))),
        "r2": float(1 - resid.var() / y.var()),
    }


def predict(df: pd.DataFrame, params: dict) -> np.ndarray:
    return design(df, pd.Timestamp(params["t0"]).date(), params["span_days"]) @ params["coef"]


def level_curve(dates: pd.DatetimeIndex, params: dict) -> np.ndarray:
    """Urovnova slozka (denni prumer baze pro Po-Ct mimo prazdniny)."""
    t0 = pd.Timestamp(params["t0"])
    if dates.tz is not None:  # osa modelu je v lokalnich datech, ne v UTC
        dates = dates.tz_convert(TZ).tz_localize(None)
    x = (dates - t0).days.to_numpy(float) + 0.5
    S = _spline_basis(x, params["span_days"])
    sl = _slices(params["span_days"])
    return S @ params["coef"][sl["level"]]


def profile_curve(daytype: int, params: dict, prazdniny: bool = False,
                  tod: np.ndarray | None = None) -> np.ndarray:
    """Profilova slozka (odchylka od urovne) pro dany typ dne."""
    if tod is None:
        tod = np.arange(0, 24, 0.25)
    sl = _slices(params["span_days"])
    c = params["coef"]
    F = _fourier(tod, K_PROFILE)
    if daytype == 0:
        out = F @ c[sl["dt0"]]
    else:
        b = c[sl[f"dt{daytype}"]]
        out = b[0] + F @ b[1:]
    if prazdniny:
        b = c[sl["praz"]]
        out = out + b[0] + _fourier(tod, K_PRAZ) @ b[1:]
    return out


def save(params: dict, path: Path | None = None) -> Path:
    path = path or MODEL_DIR / "base_params.npz"
    path.parent.mkdir(exist_ok=True)
    np.savez(path, **params)
    return path


def load(path: Path | None = None) -> dict:
    path = path or MODEL_DIR / "base_params.npz"
    raw = np.load(path, allow_pickle=False)
    out = {k: raw[k] for k in raw.files}
    for k in ("t0",):
        out[k] = str(out[k])
    for k in ("span_days", "smooth", "rmse", "mae", "r2"):
        out[k] = float(out[k])
    return out
