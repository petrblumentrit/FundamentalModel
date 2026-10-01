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


# penalizace zmen tempa instalaci FVE (relativne k A'A bloku). Puvodne 0.05;
# validace mimo vzorek ukazala, ze posledni rampa pred koncem okna (v zime bez
# slunce) dostava skok +30 jednotek — 5.0 ho potlaci, in-sample RMSE nemeni
# (18.66 -> 18.68), dopredna chyba roku 2026 klesa 20.0 -> 19.0; vic uz nepomaha
PV_SMOOTH = 5.0
# penalizace zmen tempa prirustku topne citlivosti (relativne k A'A bloku)
HEAT_TREND_SMOOTH = 5.0


def _pv_penalty(n_c: int, n_p: int, gram_pv: float, lam: float) -> np.ndarray:
    """Prvni diference prirustku kapacity — tempo instalaci se meni hladce."""
    P = np.zeros((n_c + n_p, n_c + n_p))
    D1 = np.diff(np.eye(n_p - 1), axis=0)  # bez interceptu (pocatecni kapacita)
    P[n_c + 1:, n_c + 1:] = lam * gram_pv * (D1.T @ D1)
    return P


def fit_cooling_pv(df: pd.DataFrame, resid: np.ndarray,
                   mask: np.ndarray | None = None) -> tuple[dict, dict]:
    """resid = y - B - P_topeni na cele ose; vraci (params_chlazeni, params_fve).

    mask = radky, ze kterych se odhaduje (validace mimo vzorek); regresory se
    stavi na cele ose, protoze setrvacnostni filtry potrebuji souvislou historii.
    Rozsah rampy FVE se bere jen z trenovacich radku — mimo nej kapacita drzi
    posledni hodnotu.
    """
    if mask is None:
        mask = np.ones(len(df), bool)
    idx = df.index[mask]
    t0 = idx[0]
    span = float((idx[-1] - t0).total_seconds() / 86400.0) + 1.0
    r = resid[mask]
    n_c = cooling.k_basis(df.head(1)).shape[1]
    n_p = pv.capacity_basis(df.head(1), t0, span).shape[1]
    lb = np.concatenate([np.full(n_c, -np.inf), np.zeros(n_p)])
    ub = np.full(n_c + n_p, np.inf)

    def inner(p):
        A_pv = pv.design(df, p[-1], t0, span)[mask]
        A = np.hstack([cooling.design(df, p[:-1])[mask], A_pv])
        pen = _pv_penalty(n_c, n_p, float(np.mean((A_pv * A_pv).sum(0))), PV_SMOOTH)
        coef = bounded_lstsq(A, r, lb, ub, penalty=pen)
        return coef, r - A @ coef

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
             "r2": float(1 - res.var() / r.var())}
    c_params.update(stats)
    p_params.update(stats)
    return c_params, p_params


def joint(df: pd.DataFrame, bp: dict, hp: dict, cp: dict, pp: dict,
          max_nfev: int = 400, mask: np.ndarray | None = None,
          fix_shape: bool = False) -> tuple[dict, dict, dict, dict]:
    """Zaverecny spolecny fit vsech modulu (krok 6) s warm startem z kroku 1-5.

    Sekvencni odhad nechava po sdilenych regresorech (hlavne osvitu) systematicke
    chyby — navrh proto joint fit oznacuje za nevynechatelny. Uloha zustava
    separabilni: nelinearni jsou jen parametry tvaru (14), vsechny koeficienty
    urovni a profilu se resi linearne.

    Bazalni blok na nelinearnich parametrech nezavisi, takze se jeho navrh i
    Gramova matice pocitaji jednou; kazde vyhodnoceni pak stavi jen bloky
    topeni, chlazeni a FVE.

    mask = radky pouzite pro odhad (validace mimo vzorek). Regresory vcetne
    setrvacnostnich filtru se stavi na cele ose, do normalnich rovnic vstupuji
    jen radky masky.

    fix_shape = nelinearni parametry tvaru drzet z warm startu a resit jen
    linearni cast (jedno vyhodnoceni misto ~50) — pro denni prefit v provozu,
    kde se tvar meni pomalu a staci ho preodhadnout jednou za cas.
    """
    if mask is None:
        mask = np.ones(len(df), bool)
    y = df["baseload"].to_numpy()[mask]
    t0_b = pd.Timestamp(bp["t0"]).date()
    F = base.design(df, t0_b, bp["span_days"])[mask]
    p_F = F.shape[1]
    FtF, Fty = F.T @ F, F.T @ y

    t0_p = pd.Timestamp(pp["t0"])
    span_p = pp["span_days"]
    n_h = heating.k_basis(df.head(1)).shape[1]
    n_c = cooling.k_basis(df.head(1)).shape[1]
    n_p = pv.capacity_basis(df.head(1), t0_p, span_p).shape[1]
    n_r = heating.trend_basis(df.head(1), t0_p, span_p).shape[1]
    n_j = 1 + 2 * heating.TREND_K_TOD
    n_t = n_j * n_r
    p = p_F + n_h + n_c + n_p + n_t
    sl_pv = slice(p_F + n_h + n_c, p_F + n_h + n_c + n_p)
    sl_tr = slice(sl_pv.stop, p)

    sl_b = base._slices(bp["span_days"])
    ns = sl_b["level"].stop - sl_b["level"].start
    D2 = np.diff(np.eye(ns), n=2, axis=0)
    pen = np.zeros((p, p))
    # penalizace hladkosti urovne musi rust s poctem radku, jinak je tu proti
    # 161k radkum 7x slabsi nez pri fitu z ~200 mirnych dnu a uroven zacne
    # chytat sezonnost, kterou ma nest topeni (viz varovani v navrhu)
    smooth = bp["smooth"] * int(mask.sum()) / (float(bp["n_mild_days"]) * 96.0)
    pen[sl_b["level"], sl_b["level"]] = smooth * (D2.T @ D2)
    D1 = np.diff(np.eye(n_p - 1), axis=0)
    pv_lo = p_F + n_h + n_c + 1
    # zmeny tempa penalizovane zvlast v kazde slozce denniho tvaru prirustku
    D1r = np.diff(np.eye(n_r), axis=0)
    D1t = np.kron(np.eye(n_j), D1r)

    lb = np.full(p, -np.inf)
    lb[sl_pv] = 0.0              # kapacita FVE neklesajici
    lb[sl_b["svetlo"]] = 0.0     # osvetleni za tmy jen pridava spotrebu
    ub = np.full(p, np.inf)

    n_hp = len(heating.PARAM_NAMES)
    n_cp = len(cooling.PARAM_NAMES)

    T_h = heating.trend_design(df, t0_p, span_p)

    def blocks(q):
        th_h, th_c, gamma = q[:n_hp], q[n_hp:n_hp + n_cp], q[-1]
        A_pv = pv.design(df, gamma, t0_p, span_p)[mask]
        H = heating.shape_term(df, th_h)[:, None]
        A_tr = (T_h * H)[mask]
        V = np.hstack([(heating.k_basis(df) * H)[mask],
                       cooling.design(df, th_c)[mask], A_pv, A_tr])
        return V, float(np.mean((A_pv * A_pv).sum(0))), float(np.mean((A_tr * A_tr).sum(0)))

    def inner(q):
        V, gram_pv, gram_tr = blocks(q)
        G = np.empty((p, p))
        G[:p_F, :p_F] = FtF
        FtV = F.T @ V
        G[:p_F, p_F:] = FtV
        G[p_F:, :p_F] = FtV.T
        G[p_F:, p_F:] = V.T @ V
        c = np.concatenate([Fty, V.T @ y])
        P = pen.copy()
        P[pv_lo:sl_pv.stop, pv_lo:sl_pv.stop] = PV_SMOOTH * gram_pv * (D1.T @ D1)
        P[sl_tr, sl_tr] = HEAT_TREND_SMOOTH * gram_tr * (D1t.T @ D1t)
        coef = bounded_normal(G, c, lb, ub, penalty=P)
        return coef, y - F @ coef[:p_F] - V @ coef[p_F:]

    x0 = np.concatenate([[hp[k] for k in heating.PARAM_NAMES],
                         [cp[k] for k in cooling.PARAM_NAMES], [pp["gamma"]]])
    lower = np.concatenate([heating.LOWER, cooling.LOWER, [0.0]])
    upper = np.concatenate([heating.UPPER, cooling.UPPER, [0.01]])
    x0 = np.clip(x0, lower + 1e-9, upper - 1e-9)
    if fix_shape:
        q = x0
    else:
        q = least_squares(lambda q: inner(q)[1], x0, bounds=(lower, upper),
                          diff_step=1e-3, x_scale=upper - lower, max_nfev=max_nfev).x
    coef, res = inner(q)

    stats = {"rmse": float(np.sqrt(np.mean(res**2))),
             "r2": float(1 - res.var() / y.var())}
    bo = dict(bp); bo["coef"] = coef[:p_F]; bo.update(stats)
    ho = dict(zip(heating.PARAM_NAMES, q[:n_hp]))
    ho.update({"k_coef": coef[p_F:p_F + n_h], "k_tod": heating.K_TOD,
               "trend_coef": coef[sl_tr], "trend_t0": str(t0_p), "trend_span": span_p, **stats})
    co = dict(zip(cooling.PARAM_NAMES, q[n_hp:n_hp + n_cp]))
    co.update({"k_coef": coef[p_F + n_h:p_F + n_h + n_c], "k_tod": cooling.K_TOD, **stats})
    po = {"gamma": float(q[-1]), "delta": coef[sl_pv],
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
