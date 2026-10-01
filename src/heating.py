"""Topny modul — krok 4 navrhu (NAVRH_MODELU.md).

Fituje se na reziduich po odectu baze (y - B) v chladnem a mirnem obdobi:

    P_topeni(t) = k(cas dne, typ dne) * g(T*) / COP(T)

- pocitova teplota:  T_ef = T + a*I - b*v*(T_in - T)+   (T_in = 20 C)
- dvoukanalovy filtr: T* = w*T_ef + (1-w)*EMA(T_ef; tau)
  (rychly kanal = ekviterm, pomaly = vnitrni termostat pres setrvacnost budovy)
- topna krivka:      g(T*) = s * ln(1 + exp((T_b - T*)/s))   (hladky hinge)
- COP na okamzitou teplotu (vyparnik venku): COP(T) = max(1 + alpha*T, 0.2),
  normalizace COP(0) = 1 — absolutni skala je v k
- k(cas dne, typ dne) = Fourier (K_TOD harmonickych) + offsety typu dne;
  pro dane nelinearni parametry je model linearni v k -> separabilni LS
  (vnejsi nelinearni optimalizace jen pres 7 parametru)
- casova zmena citlivosti: k(cas dne, typ dne) + dk(t), kde dk(t) je
  pomaly prirustek od zacatku dat — aditivni clen "pribyle topne zateze"
  (tepelna cerpadla misto plynu po roce 2022, proti tomu zateplovani). Rampy
  po TREND_KNOT_DAYS dnech jako kapacita FVE, bez omezeni znamenka (cisty
  efekt muze jit obema smery), penalizace zmen tempa (instalace mohou
  saturovat); za koncem dat drzi posledni hodnotu, trend se neextrapoluje.
  Prirustek ma **vlastni denni tvar**: dk(t, cas dne) = sum_j phi_j(cas dne)
  * dk_j(t), phi = [1, TREND_K_TOD harmonickych] — nova topna zatez (tepelna
  cerpadla s nocnim utlumem a ekvitermou) nemusi mit denni prubeh jako ta
  puvodni. Kazda slozka ma vlastni rampy se stejnou penalizaci.
  Odhaduje se jen v joint fitu (staged fit ho nema).
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

import pv

K_TOD = 4        # harmonickych v k(cas dne)
TREND_K_TOD = 2  # harmonickych v dennim tvaru prirustku citlivosti
T_IN = 20.0      # vnitrni teplota pro vetrny clen [C]
COP_FLOOR = 0.2
N_DAYTYPES = 4

# poradi nelinearnich parametru a jejich meze
PARAM_NAMES = ["a", "b", "t_b", "s", "w", "tau", "alpha"]
LOWER = np.array([0.0, 0.0, 10.0, 0.5, 0.0, 4.0, 0.0])
UPPER = np.array([0.05, 0.03, 18.0, 6.0, 1.0, 120.0, 0.06])
X0 = np.array([0.01, 0.005, 14.5, 2.5, 0.5, 36.0, 0.02])

MODEL_DIR = Path(__file__).resolve().parent.parent / "models"


def softplus(z: np.ndarray, s: float) -> np.ndarray:
    zs = z / s
    return np.where(zs > 30, z, s * np.log1p(np.exp(np.clip(zs, -700, 30))))


def t_effective(df: pd.DataFrame, a: float, b: float) -> np.ndarray:
    dT = np.maximum(T_IN - df["temp"].to_numpy(), 0.0)
    return df["temp"].to_numpy() + a * df["sun"].to_numpy() - b * df["wind"].to_numpy() * dT


def t_star(df: pd.DataFrame, a: float, b: float, w: float, tau: float) -> np.ndarray:
    """Dvoukanalova efektivni teplota na cele souvisle ose (EMA nelze na vyseku)."""
    te = t_effective(df, a, b)
    alpha_ema = 1.0 - np.exp(-0.25 / tau)  # 15min krok
    smooth = pd.Series(te).ewm(alpha=alpha_ema, adjust=False).mean().to_numpy()
    return w * te + (1.0 - w) * smooth


def k_basis(df: pd.DataFrame) -> np.ndarray:
    """[1, cos/sin harmonicky, dummy dt1..dt3] — sdileny tvar, offsety typu dne."""
    tod = df["tod"].to_numpy()
    wv = 2 * np.pi * np.outer(tod, np.arange(1, K_TOD + 1)) / 24.0
    dt = df["daytype"].to_numpy()
    dums = np.stack([(dt == k).astype(float) for k in range(1, N_DAYTYPES)], axis=1)
    return np.hstack([np.ones((len(df), 1)), np.cos(wv), np.sin(wv), dums])


def shape_term(df: pd.DataFrame, theta: np.ndarray) -> np.ndarray:
    """g(T*)/COP pro dane nelinearni parametry theta (na cele ose df)."""
    a, b, t_b, s, w, tau, alpha = theta
    ts = t_star(df, a, b, w, tau)
    cop = np.maximum(1.0 + alpha * df["temp"].to_numpy(), COP_FLOOR)
    return softplus(t_b - ts, s) / cop


def fit(df: pd.DataFrame, resid: np.ndarray, mask: np.ndarray) -> dict:
    """df = cela souvisla osa; resid = y - B; mask = radky pro fit (bez chlazeni)."""
    Phi = k_basis(df)
    r = resid[mask]

    def inner(theta):
        H = shape_term(df, theta)
        A = Phi[mask] * H[mask, None]
        coef, *_ = np.linalg.lstsq(A, r, rcond=None)
        return coef, r - A @ coef

    def outer(theta):
        return inner(theta)[1]

    sol = least_squares(outer, X0, bounds=(LOWER, UPPER), diff_step=1e-3,
                        x_scale=UPPER - LOWER, verbose=0)
    k_coef, res = inner(sol.x)
    out = dict(zip(PARAM_NAMES, sol.x))
    out.update({
        "k_coef": k_coef, "k_tod": K_TOD,
        "rmse": float(np.sqrt(np.mean(res**2))),
        "r2": float(1 - res.var() / r.var()),
    })
    return out


def trend_basis(df: pd.DataFrame, t0, span_days: float) -> np.ndarray:
    """Rampy prirustku citlivosti (bez interceptu — ten nese k); za koncem
    rozsahu drzi, jen plne pozorovane rampy (stejne jako kapacita FVE)."""
    return pv.capacity_basis(df, pd.Timestamp(t0), span_days)[:, 1:]


def trend_tod(tod: np.ndarray) -> np.ndarray:
    """phi_j(cas dne): [1, cos, sin ...] — denni tvar prirustku."""
    w = 2 * np.pi * np.outer(tod, np.arange(1, TREND_K_TOD + 1)) / 24.0
    return np.hstack([np.ones((len(tod), 1)), np.cos(w), np.sin(w)])


def trend_design(df: pd.DataFrame, t0, span_days: float) -> np.ndarray:
    """Sloupce phi_j(cas dne) * rampa_r(t), poradi j-major (blok na slozku)."""
    R = trend_basis(df, t0, span_days)
    Phi = trend_tod(df["tod"].to_numpy())
    return (Phi[:, :, None] * R[:, None, :]).reshape(len(df), -1)


def k_total(df: pd.DataFrame, params: dict) -> np.ndarray:
    """k(cas dne, typ dne) + dk(t, cas dne)."""
    k_t = k_basis(df) @ params["k_coef"]
    if "trend_coef" in params:
        k_t = k_t + trend_design(df, params["trend_t0"], params["trend_span"]) @ params["trend_coef"]
    return k_t


def predict(df: pd.DataFrame, params: dict) -> np.ndarray:
    theta = np.array([params[k] for k in PARAM_NAMES], float)
    return k_total(df, params) * shape_term(df, theta)


def _trend_blocks(params: dict) -> np.ndarray:
    """Koeficienty trendu jako matice (slozka denniho tvaru x rampa)."""
    c = np.asarray(params["trend_coef"])
    return c.reshape(1 + 2 * TREND_K_TOD, -1)


def trend_curve(dates: pd.DatetimeIndex, params: dict) -> np.ndarray:
    """Prumer dk(t, cas dne) pres den — prirustek topne citlivosti (jednotky k)."""
    if "trend_coef" not in params:
        return np.zeros(len(dates))
    fake = pd.DataFrame(index=dates)
    return trend_basis(fake, params["trend_t0"], params["trend_span"]) @ _trend_blocks(params)[0]


def trend_profile(date, params: dict, tod: np.ndarray | None = None) -> np.ndarray:
    """dk(t, cas dne) v danem dni podle casu dne."""
    if tod is None:
        tod = np.arange(0, 24, 0.25)
    fake = pd.DataFrame(index=pd.DatetimeIndex([pd.Timestamp(date)]))
    r = trend_basis(fake, params["trend_t0"], params["trend_span"])[0]
    return trend_tod(tod) @ (_trend_blocks(params) @ r)


def k_curve(daytype: int, params: dict, tod: np.ndarray | None = None) -> np.ndarray:
    if tod is None:
        tod = np.arange(0, 24, 0.25)
    fake = pd.DataFrame({"tod": tod, "daytype": daytype})
    return k_basis(fake) @ params["k_coef"]


def save(params: dict, path: Path | None = None) -> Path:
    path = path or MODEL_DIR / "heating_params.npz"
    path.parent.mkdir(exist_ok=True)
    np.savez(path, **params)
    return path


def load(path: Path | None = None) -> dict:
    path = path or MODEL_DIR / "heating_params.npz"
    raw = np.load(path, allow_pickle=False)
    out = {k: raw[k] for k in raw.files}
    for k in PARAM_NAMES + ["rmse", "r2"]:
        out[k] = float(out[k])
    if "trend_coef" in out:
        out["trend_t0"] = str(out["trend_t0"])
        out["trend_span"] = float(out["trend_span"])
    return out
