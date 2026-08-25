"""Spolecne odhady pres vice modulu — sdilene regresory.

Osvit je sdileny regresor tri clenu (solarni zisky v topeni, chlazeni, FVE);
navrh (NAVRH_MODELU.md) proto zakazuje odhadovat je nezavisle za sebou.
Chlazeni a FVE se tu fituji **spolecne** jednim separabilnim LS:

- vnejsi nelinearni optimalizace pres parametry tvaru (chlazeni + gamma FVE),
- vnitri linearni krok pres k_c (chlazeni, bez omezeni) a delta (prirustky
  kapacity FVE, omezene na >= 0, aby kapacita jen rostla).

Vnitrni krok se resi pres normalni rovnice (Choleskeho rozklad), takze omezena
uloha ma rozmer poctu koeficientu, ne poctu radku — jinak by kazde vyhodnoceni
vnejsi funkce prochazelo matici 160k x 30.
"""
import numpy as np
import pandas as pd
from scipy.optimize import least_squares, lsq_linear

import base
import cooling
import heating
import pv


def bounded_normal(G: np.ndarray, c: np.ndarray, lb: np.ndarray, ub: np.ndarray,
                   penalty: np.ndarray | None = None, ridge: float = 1e-9) -> np.ndarray:
    """Omezene LS zadane normalnimi rovnicemi (G = A'A, c = A'y)."""
    if penalty is not None:
        G = G + penalty
    scale = np.sqrt(np.diag(G))
    scale[scale <= 0] = 1.0
    G = G / np.outer(scale, scale)
    c = c / scale
    G.flat[:: G.shape[0] + 1] += ridge * np.trace(G) / G.shape[0]
    R = np.linalg.cholesky(G).T
    z = np.linalg.solve(R.T, c)
    sol = lsq_linear(R, z, bounds=(lb * scale, ub * scale))
    return sol.x / scale


def bounded_lstsq(A: np.ndarray, y: np.ndarray, lb: np.ndarray, ub: np.ndarray,
                  penalty: np.ndarray | None = None, ridge: float = 1e-9) -> np.ndarray:
    """min ||Ax - y||^2 + x'Px s omezenim lb <= x <= ub."""
    return bounded_normal(A.T @ A, A.T @ y, lb, ub, penalty, ridge)


PV_SMOOTH = 0.05  # penalizace zmen tempa instalaci FVE (relativne k A'A bloku)


def _pv_penalty(n_c: int, n_p: int, gram_pv: float, lam: float) -> np.ndarray:
    """Prvni diference prirustku kapacity — tempo instalaci se meni hladce."""
    P = np.zeros((n_c + n_p, n_c + n_p))
    D1 = np.diff(np.eye(n_p - 1), axis=0)  # bez interceptu (pocatecni kapacita)
    P[n_c + 1:, n_c + 1:] = lam * gram_pv * (D1.T @ D1)
    return P


def fit_cooling_pv(df: pd.DataFrame, resid: np.ndarray) -> tuple[dict, dict]:
    """resid = y - B - P_topeni na cele ose; vraci (params_chlazeni, params_fve)."""
    t0 = df.index[0]
    span = float((df.index[-1] - t0).total_seconds() / 86400.0) + 1.0
    n_c = cooling.k_basis(df.head(1)).shape[1]
    n_p = pv.capacity_basis(df.head(1), t0, span).shape[1]
    lb = np.concatenate([np.full(n_c, -np.inf), np.zeros(n_p)])
    ub = np.full(n_c + n_p, np.inf)

    def inner(p):
        A_pv = pv.design(df, p[-1], t0, span)
        A = np.hstack([cooling.design(df, p[:-1]), A_pv])
        pen = _pv_penalty(n_c, n_p, float(np.mean((A_pv * A_pv).sum(0))), PV_SMOOTH)
        coef = bounded_lstsq(A, resid, lb, ub, penalty=pen)
        return coef, resid - A @ coef

    lower = np.concatenate([cooling.LOWER, [0.0]])      # gamma >= 0
    upper = np.concatenate([cooling.UPPER, [0.01]])
    x0 = np.concatenate([cooling.X0, [0.004]])
    sol = least_squares(lambda p: inner(p)[1], x0, bounds=(lower, upper),
                        diff_step=1e-3, x_scale=upper - lower, verbose=0)
    coef, res = inner(sol.x)

    c_params = dict(zip(cooling.PARAM_NAMES, sol.x[:-1]))
    c_params.update({"k_coef": coef[:n_c], "k_tod": cooling.K_TOD})
    p_params = {"gamma": float(sol.x[-1]), "delta": coef[n_c:],
                "t0": str(t0), "span_days": span, "knot_days": pv.KNOT_DAYS}
    stats = {"rmse": float(np.sqrt(np.mean(res**2))),
             "r2": float(1 - res.var() / resid.var())}
    c_params.update(stats)
    p_params.update(stats)
    return c_params, p_params


def joint(df: pd.DataFrame, bp: dict, hp: dict, cp: dict, pp: dict,
          max_nfev: int = 400) -> tuple[dict, dict, dict, dict]:
    """Zaverecny spolecny fit vsech modulu (krok 6) s warm startem z kroku 1-5.

    Sekvencni odhad nechava po sdilenych regresorech (hlavne osvitu) systematicke
    chyby — navrh proto joint fit oznacuje za nevynechatelny. Uloha zustava
    separabilni: nelinearni jsou jen parametry tvaru (14), vsechny koeficienty
    urovni a profilu se resi linearne.

    Bazalni blok na nelinearnich parametrech nezavisi, takze se jeho navrh i
    Gramova matice pocitaji jednou; kazde vyhodnoceni pak stavi jen bloky
    topeni, chlazeni a FVE.
    """
    y = df["baseload"].to_numpy()
    t0_b = pd.Timestamp(bp["t0"]).date()
    F = base.design(df, t0_b, bp["span_days"])
    p_F = F.shape[1]
    FtF, Fty = F.T @ F, F.T @ y

    t0_p = pd.Timestamp(pp["t0"])
    span_p = pp["span_days"]
    n_h = heating.k_basis(df.head(1)).shape[1]
    n_c = cooling.k_basis(df.head(1)).shape[1]
    n_p = pv.capacity_basis(df.head(1), t0_p, span_p).shape[1]
    p = p_F + n_h + n_c + n_p

    sl_b = base._slices(bp["span_days"])
    ns = sl_b["level"].stop - sl_b["level"].start
    D2 = np.diff(np.eye(ns), n=2, axis=0)
    pen = np.zeros((p, p))
    # penalizace hladkosti urovne musi rust s poctem radku, jinak je tu proti
    # 161k radkum 7x slabsi nez pri fitu z ~200 mirnych dnu a uroven zacne
    # chytat sezonnost, kterou ma nest topeni (viz varovani v navrhu)
    smooth = bp["smooth"] * len(df) / (float(bp["n_mild_days"]) * 96.0)
    pen[sl_b["level"], sl_b["level"]] = smooth * (D2.T @ D2)
    D1 = np.diff(np.eye(n_p - 1), axis=0)
    pv_lo = p_F + n_h + n_c + 1

    lb = np.full(p, -np.inf)
    lb[p_F + n_h + n_c:] = 0.0   # kapacita FVE neklesajici
    ub = np.full(p, np.inf)

    n_hp = len(heating.PARAM_NAMES)
    n_cp = len(cooling.PARAM_NAMES)

    def blocks(q):
        th_h, th_c, gamma = q[:n_hp], q[n_hp:n_hp + n_cp], q[-1]
        A_pv = pv.design(df, gamma, t0_p, span_p)
        V = np.hstack([heating.k_basis(df) * heating.shape_term(df, th_h)[:, None],
                       cooling.design(df, th_c), A_pv])
        return V, float(np.mean((A_pv * A_pv).sum(0)))

    def inner(q):
        V, gram_pv = blocks(q)
        G = np.empty((p, p))
        G[:p_F, :p_F] = FtF
        FtV = F.T @ V
        G[:p_F, p_F:] = FtV
        G[p_F:, :p_F] = FtV.T
        G[p_F:, p_F:] = V.T @ V
        c = np.concatenate([Fty, V.T @ y])
        P = pen.copy()
        P[pv_lo:, pv_lo:] = PV_SMOOTH * gram_pv * (D1.T @ D1)
        coef = bounded_normal(G, c, lb, ub, penalty=P)
        return coef, y - F @ coef[:p_F] - V @ coef[p_F:]

    x0 = np.concatenate([[hp[k] for k in heating.PARAM_NAMES],
                         [cp[k] for k in cooling.PARAM_NAMES], [pp["gamma"]]])
    lower = np.concatenate([heating.LOWER, cooling.LOWER, [0.0]])
    upper = np.concatenate([heating.UPPER, cooling.UPPER, [0.01]])
    x0 = np.clip(x0, lower + 1e-9, upper - 1e-9)
    sol = least_squares(lambda q: inner(q)[1], x0, bounds=(lower, upper),
                        diff_step=1e-3, x_scale=upper - lower, max_nfev=max_nfev)
    coef, res = inner(sol.x)

    stats = {"rmse": float(np.sqrt(np.mean(res**2))),
             "r2": float(1 - res.var() / y.var())}
    bo = dict(bp); bo["coef"] = coef[:p_F]; bo.update(stats)
    ho = dict(zip(heating.PARAM_NAMES, sol.x[:n_hp]))
    ho.update({"k_coef": coef[p_F:p_F + n_h], "k_tod": heating.K_TOD, **stats})
    co = dict(zip(cooling.PARAM_NAMES, sol.x[n_hp:n_hp + n_cp]))
    co.update({"k_coef": coef[p_F + n_h:p_F + n_h + n_c], "k_tod": cooling.K_TOD, **stats})
    po = {"gamma": float(sol.x[-1]), "delta": coef[p_F + n_h + n_c:],
          "t0": str(t0_p), "span_days": span_p, "knot_days": pv.KNOT_DAYS, **stats}
    return bo, ho, co, po


def save(params: dict, name: str) -> str:
    path = pv.MODEL_DIR / name
    path.parent.mkdir(exist_ok=True)
    np.savez(path, **params)
    return str(path)


def load(name: str, floats: list[str]) -> dict:
    raw = np.load(pv.MODEL_DIR / name, allow_pickle=False)
    out = {k: raw[k] for k in raw.files}
    for k in floats:
        out[k] = float(out[k])
    if "t0" in out:
        out["t0"] = str(out["t0"])
    return out
