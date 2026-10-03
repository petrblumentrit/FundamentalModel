"""FVE — lokalni vyroba, krok 5 navrhu (NAVRH_MODELU.md).

    G_FVE(t) = C(t) * I(t) * eta(T)

- C(t): agregatni kapacita flotily, **monotonne neklesajici** — parametrizovana
  jako pocatecni uroven + soucet nezapornych rampovych prirustku po KNOT_DAYS
  dnech. Nezapornost prirustku (bounds v linearnim reseni) vynuti, ze kapacita
  jen roste; instalace se neodinstalovavaji.
- eta(T) = 1 - gamma*(T - 25): teplotni derating panelu, normalizace eta(25)=1
  (absolutni skala je v C).

Clen je v koeficientech C linearni -> resi se ve vnitrnim linearnim kroku
spolecne s chlazenim (viz `fit.py`; osvit je sdileny regresor).
"""
from pathlib import Path

import numpy as np
import pandas as pd

import config

_C = config.model()["fve"]                   # hodnoty v config/model.yaml
KNOT_DAYS = float(_C["kapacita_uzel_dni"])   # rozestup rampovych prirustku kapacity
T_REF = _C["referencni_teplota"]             # referencni teplota clanku pro eta
PARAM_NAMES = ["gamma"]
LOWER, UPPER, X0 = config.bounds("fve", PARAM_NAMES)

MODEL_DIR = Path(__file__).resolve().parent.parent / "models"


def _days(df: pd.DataFrame, t0: pd.Timestamp) -> np.ndarray:
    return (df.index - t0).total_seconds().to_numpy() / 86400.0


def capacity_basis(df: pd.DataFrame, t0: pd.Timestamp, span_days: float) -> np.ndarray:
    """[1, rampy] — C = basis @ delta je pri delta >= 0 neklesajici.

    Jen rampy, ktere jsou v rozsahu (span_days) pozorovane cele: rampa
    zachycena par dny na konci okna (typicky v prosinci, bez slunce) neni z dat
    identifikovatelna a validace mimo vzorek ji ukazala jako skok +40 jednotek.
    Za koncem rozsahu kapacita drzi posledni hodnotu.
    """
    x = np.minimum(_days(df, t0), span_days)
    cols = [np.ones(len(x))]
    for j in np.arange(KNOT_DAYS, span_days - KNOT_DAYS + 1e-6, KNOT_DAYS):
        cols.append(np.clip((x - j) / KNOT_DAYS, 0.0, 1.0))
    return np.column_stack(cols)


def eta(df: pd.DataFrame, gamma: float) -> np.ndarray:
    return 1.0 - gamma * (df["temp"].to_numpy() - T_REF)


def design(df: pd.DataFrame, gamma: float, t0: pd.Timestamp, span_days: float) -> np.ndarray:
    """Sloupce pro linearni reseni; vyroba se od spotreby odecita -> zaporne."""
    w = df["sun"].to_numpy() * eta(df, gamma)
    return -capacity_basis(df, t0, span_days) * w[:, None]


def predict(df: pd.DataFrame, params: dict) -> np.ndarray:
    """G_FVE >= 0 (vyroba); v aditivnim modelu se odecita."""
    t0 = pd.Timestamp(params["t0"])
    return -design(df, params["gamma"], t0, params["span_days"]) @ params["delta"]


def capacity_curve(dates: pd.DatetimeIndex, params: dict) -> np.ndarray:
    """C(t) — merena v jednotkach spotreby na W/m2 osvitu."""
    t0 = pd.Timestamp(params["t0"])
    fake = pd.DataFrame(index=dates)
    return capacity_basis(fake, t0, params["span_days"]) @ params["delta"]
