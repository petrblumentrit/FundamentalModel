"""Chladici modul — krok 5 navrhu (NAVRH_MODELU.md).

Zrcadlovy protejsek topneho modulu:

    P_chlazeni(t) = k_c(cas dne, typ dne) * g_c(T*_c) / EER(T)

- pocitova teplota chlazeni: T_ef = T + a_c*I  (vetsi role osvitu, vitr se
  neuplatnuje — okna se v horku nezavirala kvuli infiltraci)
- vlastni, **kratsi** setrvacnost: T*_c = w_c*T_ef + (1-w_c)*EMA(T_ef; tau_c),
  tau_c v hodinach (klimatizace reaguje rychleji nez topeni)
- zrcadlovy hladky hinge: g_c(T*) = s_c * ln(1 + exp((T* - T_bc)/s_c))
- EER na okamzitou teplotu (kondenzator je venku, ucinnost klesa s horkem):
  EER(T) = max(1 - alpha_c*(T - 25), EER_FLOOR), normalizace EER(25) = 1

Model je v k_c linearni -> odhad separabilne spolu s FVE (viz `fit.py`).
"""
from pathlib import Path

import numpy as np
import pandas as pd

from heating import softplus

K_TOD = 4          # harmonickych v k_c(cas dne)
N_DAYTYPES = 4
T_REF = 25.0
EER_FLOOR = 0.3

PARAM_NAMES = ["a_c", "t_bc", "s_c", "w_c", "tau_c", "alpha_c"]
LOWER = np.array([0.0, 16.0, 0.5, 0.0, 0.25, 0.0])
UPPER = np.array([0.10, 26.0, 6.0, 1.0, 72.0, 0.06])
X0 = np.array([0.02, 21.0, 2.0, 0.5, 6.0, 0.02])

MODEL_DIR = Path(__file__).resolve().parent.parent / "models"


def t_star(df: pd.DataFrame, a_c: float, w_c: float, tau_c: float) -> np.ndarray:
    te = df["temp"].to_numpy() + a_c * df["sun"].to_numpy()
    alpha_ema = 1.0 - np.exp(-0.25 / tau_c)  # 15min krok
    smooth = pd.Series(te).ewm(alpha=alpha_ema, adjust=False).mean().to_numpy()
    return w_c * te + (1.0 - w_c) * smooth


def k_basis(df: pd.DataFrame) -> np.ndarray:
    """[1, cos/sin harmonicky, dummy dt1..dt3] — sdileny tvar, offsety typu dne."""
    tod = df["tod"].to_numpy()
    wv = 2 * np.pi * np.outer(tod, np.arange(1, K_TOD + 1)) / 24.0
    dt = df["daytype"].to_numpy()
    dums = np.stack([(dt == k).astype(float) for k in range(1, N_DAYTYPES)], axis=1)
    return np.hstack([np.ones((len(df), 1)), np.cos(wv), np.sin(wv), dums])


def shape_term(df: pd.DataFrame, theta: np.ndarray) -> np.ndarray:
    """g_c(T*_c)/EER pro dane nelinearni parametry theta."""
    a_c, t_bc, s_c, w_c, tau_c, alpha_c = theta
    ts = t_star(df, a_c, w_c, tau_c)
    eer = np.maximum(1.0 - alpha_c * (df["temp"].to_numpy() - T_REF), EER_FLOOR)
    return softplus(ts - t_bc, s_c) / eer


def design(df: pd.DataFrame, theta: np.ndarray) -> np.ndarray:
    return k_basis(df) * shape_term(df, theta)[:, None]


def predict(df: pd.DataFrame, params: dict) -> np.ndarray:
    theta = np.array([params[k] for k in PARAM_NAMES], float)
    return design(df, theta) @ params["k_coef"]


def k_curve(daytype: int, params: dict, tod: np.ndarray | None = None) -> np.ndarray:
    if tod is None:
        tod = np.arange(0, 24, 0.25)
    fake = pd.DataFrame({"tod": tod, "daytype": daytype})
    return k_basis(fake) @ params["k_coef"]
