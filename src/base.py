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
- prubeh leta: prazdninovy efekt neni po cele leto stejny — pracovni dny jsou
  21. 7. - 9. 8. o ~7 niz (v 7-8 h az -25; dovolene vrcholi), 20.-29. 8. o
  ~6 vys (navrat pred koncem prazdnin), kazdy rok stejne. Hladky clen podle
  dne v lete (kubicke B-spliny v okne z config/kalendar.yaml, na okrajich
  nulove): pro pracovni dny s vlastnim dennim tvarem (konstanta +
  harmonicke), pro volne dny jen uroven. Uroven (uzly po 90 dnech) tak
  rychly prubeh neunese.
- rano po volnu: pracovni den po vikendu/svatku ma do ~NIGHT_END h vlastni
  rezim (noc z nedele na pondeli je jeste vikendova, prechod na pracovni
  rezim probiha az behem rana). Aditivni clen z kubickych B-splin na
  [0, NIGHT_END], v NIGHT_END plynule (hodnota i sklon) odezni na nulu.
  Bez nej je profil Po-Ct kompromisem: pondelni noc model prestreluje o ~10,
  noci Ut-Ct podstreluje o 2-5.
- osvetleni: tma(t) * aktivita(cas dne) — spotreba rizena svetlem, ne
  hodinami (pri zmene casu se ~45 jednotek posune se soumrakem, viz
  explore/svetlo_test.py). tma je astronomicka (src/sun.py), aktivita jsou
  kubicke B-spliny jen v oknech usvitu a soumraku, na okrajich oken plynule
  nulove: mimo okna je tma po cely rok stejna (poledne vzdy svetlo, pozdni
  noc vzdy tma) a clen by tam splyval s profilem.
- sero pres den: sero(t) * aktivita(cas dne) — pri tmave obloze za dne
  (tezka oblacnost, mlha; osvit pod ~100 W/m2) je spotreba vyssi nez za
  jasneho dne se stejnou teplotou, v pracovni i volne dny (sviti se i pres
  den). sero = (1 - tma) * exp(-osvit / I0) (src/sun.py), aktivita kubicke
  B-spliny pres den, >= 0. Na rozdil od osvetleni za tmy zavisi na pocasi.
- most: celodenni aditivni posun (konstanta + K_MOST harmonickych) — mosty
  jsou typovane jako patek, ale byvaji o 20-40 nize (vybirane dovolene);
  zvlast jednodenni a dvoudenni mosty (dvoudenni byvaji slabsi), mimo
  zvlastni obdobi (Vanoce maji vlastni cleny)
- zvlastni obdobi (config/kalendar.yaml, napr. Vanoce): kazda skupina dnu
  aditivni clen s vlastnim dennim tvarem (konstanta + harmonicke z configu)
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import BSpline

import config
import kalendar
import workspace
from etl import TZ, step_hours

_C = config.model()["baze"]              # hodnoty v config/model.yaml
K_PROFILE = _C["harmonicke_profil"]      # harmonickych v dennim profilu
K_PRAZ = _C["harmonicke_prazdniny"]      # harmonickych v prazdninove korekci
KNOT_DAYS = _C["uroven_uzel_dni"]        # rozestup uzlu urovnove spline
N_DAYTYPES = 4  # 0 Po-Ct, 1 patek/most, 2 sobota, 3 nedele/svatek
# mrtve pasmo z diagnostiky (09_baze_signal): topeni dozniva ~13-14 C,
# chlazeni nastupuje uz ~18 C (drive nez navrhovy odhad 20-24 C)
MILD_BAND = tuple(_C["mirny_den"])             # denni prumer T mirneho dne
PREV_BAND = tuple(_C["mirny_predchozi_den"])   # denni prumer T predchoziho dne (setrvacnost)
SMOOTH = _C["uroven_hladkost"]           # vaha penalizace druhych diferenci spline koeficientu
NIGHT_END = _C["rano_po_volnu"]["konec"] # [h] konec rana po volnu (diagnostika: odchylka mizi 8-9 h)
NIGHT_KNOT = _C["rano_po_volnu"]["uzel"] # [h] rozestup uzlu ranniho clenu
K_MOST = _C["harmonicke_most"]           # harmonickych v korekci mostu
# okna, kde se tma behem roku meni (mistni cas, vcetne letniho casu) [h]
# rano jen do 8,5 h: pozdeji je tma jen v prosinci a lednu, tedy hlavne o
# vanocnich prazdninach, a clen chytal vanocni propad (vysel zaporny)
LIGHT_WINDOWS = tuple(tuple(w) for w in _C["osvetleni"]["okna"])
LIGHT_KNOT = _C["osvetleni"]["uzel"]     # [h] rozestup uzlu aktivity osvetleni
LIGHT_RIDGE = _C["osvetleni"]["ukotveni"]  # ridge osvetleni ve fitu z mirnych dnu (relativne)
GLOOM_WINDOW = tuple(_C["sero"]["okno"]) # [h] okno clenu sera pres den (mistni cas)
GLOOM_KNOT = _C["sero"]["uzel"]          # [h] rozestup uzlu aktivity pri seru

SUMMER = kalendar.summer_course()   # okno prubehu leta (None = clen vypnut)

MODEL_DIR = workspace.MODELS


def _fourier(tod: np.ndarray, k: int) -> np.ndarray:
    w = 2 * np.pi * np.outer(tod, np.arange(1, k + 1)) / 24.0
    return np.hstack([np.cos(w), np.sin(w)])


def _spline_basis(x_days: np.ndarray, span_days: float) -> np.ndarray:
    inner = np.arange(0.0, span_days + KNOT_DAYS, KNOT_DAYS)
    t = np.concatenate([np.repeat(inner[0], 3), inner, np.repeat(inner[-1], 3)])
    x = np.clip(x_days, inner[0], inner[-1] - 1e-9)
    return BSpline.design_matrix(x, t, 3).toarray()


def _night_basis(tod: np.ndarray) -> np.ndarray:
    """Kubicke B-spliny na [0, NIGHT_END]; posledni dve vynechane -> v NIGHT_END
    je clen i jeho sklon nulovy, za nim nulovy uplne."""
    inner = np.arange(0.0, NIGHT_END + 1e-9, NIGHT_KNOT)
    t = np.concatenate([np.repeat(inner[0], 3), inner, np.repeat(inner[-1], 3)])
    x = np.clip(tod, 0.0, NIGHT_END - 1e-9)
    B = BSpline.design_matrix(x, t, 3).toarray()[:, :-2]
    return B * (tod < NIGHT_END)[:, None]


# uzlu + 2 kubickych B-splin, bez poslednich dvou
N_NIGHT = len(np.arange(0.0, NIGHT_END + 1e-9, NIGHT_KNOT))


def _window_knots(a: float, b: float, step: float) -> np.ndarray:
    """Rovnomerne uzly pres cele okno, rozestup co nejblize step."""
    return np.linspace(a, b, max(int(round((b - a) / step)), 1) + 1)


def _window_basis(tod: np.ndarray, a: float, b: float, step: float) -> np.ndarray:
    """Kubicke B-spliny na [a, b] bez dvou krajnich na kazde strane — na obou
    okrajich je clen i jeho sklon nulovy, mimo okno nulovy."""
    inner = _window_knots(a, b, step)
    t = np.concatenate([np.repeat(inner[0], 3), inner, np.repeat(inner[-1], 3)])
    x = np.clip(tod, a, b - 1e-9)
    B = BSpline.design_matrix(x, t, 3).toarray()[:, 2:-2]
    return B * ((tod >= a) & (tod < b))[:, None]


def _light_basis(tod: np.ndarray) -> np.ndarray:
    return np.hstack([_window_basis(tod, a, b, LIGHT_KNOT) for a, b in LIGHT_WINDOWS])


# B-splin je (uzlu + 2), bez dvou krajnich na kazde strane
N_LIGHT = sum(len(_window_knots(a, b, LIGHT_KNOT)) - 2 for a, b in LIGHT_WINDOWS)


def _gloom_basis(tod: np.ndarray) -> np.ndarray:
    return _window_basis(tod, *GLOOM_WINDOW, GLOOM_KNOT)


N_GLOOM = len(_window_knots(*GLOOM_WINDOW, GLOOM_KNOT)) - 2


def _summer_basis(pos: np.ndarray) -> np.ndarray:
    """B-spliny podle dne v lete (pos = poradi dne v okne + cast dne), mimo okno nuly."""
    return _window_basis(pos, 0.0, float(SUMMER["dni"]), SUMMER["uzel_dni"])


N_SUMMER = len(_window_knots(0.0, float(SUMMER["dni"]), SUMMER["uzel_dni"])) - 2 if SUMMER else 0
K_SUMMER = SUMMER["harmonicke"] if SUMMER else 0


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
    out["noc"] = slice(i, i + N_NIGHT); i += N_NIGHT
    out["most"] = slice(i, i + 1 + 2 * K_MOST); i += 1 + 2 * K_MOST
    out["most2"] = slice(i, i + 1 + 2 * K_MOST); i += 1 + 2 * K_MOST
    out["svetlo"] = slice(i, i + N_LIGHT); i += N_LIGHT
    out["sero"] = slice(i, i + N_GLOOM); i += N_GLOOM
    out["leto_prac"] = slice(i, i + N_SUMMER * (1 + 2 * K_SUMMER)); i += N_SUMMER * (1 + 2 * K_SUMMER)
    out["leto_vol"] = slice(i, i + N_SUMMER); i += N_SUMMER
    for name, k in kalendar.period_groups():
        out["obd_" + name] = slice(i, i + 1 + 2 * k); i += 1 + 2 * k
    out["total"] = i
    return out


def design(df: pd.DataFrame, t0, span_days: float) -> np.ndarray:
    x = (pd.to_datetime(df["date_local"]) - pd.Timestamp(t0)).dt.days.to_numpy(float)
    x = x + df["tod"].to_numpy() / 24.0
    tod = df["tod"].to_numpy()
    S = _spline_basis(x, span_days)
    F = _fourier(tod, K_PROFILE)

    def fsub(k: int) -> np.ndarray:
        """Prvnich k harmonickych z F (cos/sin se pocitaji jen jednou)."""
        assert k <= K_PROFILE
        return np.hstack([F[:, :k], F[:, K_PROFILE:K_PROFILE + k]])

    dt = df["daytype"].to_numpy()
    blocks = [S]
    for k in range(N_DAYTYPES):
        m = (dt == k).astype(float)[:, None]
        b = F * m
        if k > 0:
            b = np.hstack([m, b])
        blocks.append(b)
    pz = df["prazdniny"].to_numpy(float)[:, None]
    blocks.append(np.hstack([pz, fsub(K_PRAZ) * pz]))
    # rano pracovniho dne po vikendu nebo svatku; most ma vlastni celodenni
    # clen (jinak by jeho celodenni propad stahl ranni clen pondeli)
    most = df["most"].to_numpy()
    after_off = (df["po_volnu"].to_numpy() & (dt <= 1) & ~most).astype(float)[:, None]
    blocks.append(_night_basis(tod) * after_off)
    groups = kalendar.period_groups()
    in_period = np.zeros(len(df), bool)
    for name, _ in groups:
        in_period |= df["obd_" + name].to_numpy()
    for n_run in (1, 2):
        mo = (most & ~in_period & (df["most_delka"].to_numpy() == n_run)).astype(float)[:, None]
        blocks.append(np.hstack([mo, fsub(K_MOST) * mo]))
    blocks.append(_light_basis(tod) * df["tma"].to_numpy()[:, None])
    blocks.append(_gloom_basis(tod) * df["sero"].to_numpy()[:, None])
    if SUMMER:
        ld = df["leto_den"].to_numpy()
        Sb = _summer_basis(np.where(ld >= 0, ld + tod / 24.0, -1.0))
        Phi = np.hstack([np.ones((len(df), 1)), fsub(K_SUMMER)]) if K_SUMMER else np.ones((len(df), 1))
        work = (dt <= 1)[:, None]
        blocks.append((Sb[:, :, None] * Phi[:, None, :]).reshape(len(df), -1) * work)
        blocks.append(Sb * ~work)
    for name, k in groups:
        g = df["obd_" + name].to_numpy(float)[:, None]
        blocks.append(np.hstack([g, fsub(k) * g]) if k else g)
    return np.hstack(blocks)


def daily_table(df: pd.DataFrame) -> pd.DataFrame:
    d = df.groupby("date_local").agg(
        baseload=("baseload", "mean"), temp=("temp", "mean"),
        n=("baseload", "size"), daytype=("daytype", "max"),
        prazdniny=("prazdniny", "max"),
    )
    d.index = pd.to_datetime(d.index)
    return d


def rows_per_day(df: pd.DataFrame) -> float:
    return 24.0 / step_hours(df)


def mild_dates(df: pd.DataFrame) -> pd.DatetimeIndex:
    """Dny v mrtvem pasmu; predchozi den take mirny (tepelna setrvacnost)."""
    d = daily_table(df)
    ok = (d["temp"].between(*MILD_BAND)
          & d["temp"].shift(1).between(*PREV_BAND)
          & (d["n"] >= 23 * rows_per_day(df) / 24))  # den zmeny casu ma 23 hodin
    return d.index[ok]


def fit(df: pd.DataFrame, smooth: float = SMOOTH) -> dict:
    # penalizace je absolutni proti X'X, ktere roste s poctem radku na den —
    # pro hodinova data (24 misto 96 radku) ji umerne zmensit
    rpd = rows_per_day(df)
    smooth = smooth * rpd / 96.0
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
    # mirne dny (jaro, podzim) nepokryvaji hodiny, kdy je tma jen v zime —
    # osvetleni tu stahnout k nule, plne ho odhadne az joint fit
    sl_l = sl["svetlo"]
    P[sl_l, sl_l] += LIGHT_RIDGE * np.trace(XtX) / p * np.eye(sl_l.stop - sl_l.start)
    # sero se z mirnych dnu (malo tmavych dni) odhaduje spatne — stejne ukotveni
    sl_g = sl["sero"]
    P[sl_g, sl_g] += LIGHT_RIDGE * np.trace(XtX) / p * np.eye(sl_g.stop - sl_g.start)
    # prubeh leta: mirnych dni je v lete malo — take ukotvit, odhadne joint fit
    for key in ("leto_prac", "leto_vol"):
        sl_s = sl[key]
        P[sl_s, sl_s] += LIGHT_RIDGE * np.trace(XtX) / p * np.eye(sl_s.stop - sl_s.start)
    coef = np.linalg.solve(XtX + P, X.T @ y)
    resid = y - X @ coef
    return {
        "coef": coef, "t0": str(t0), "span_days": span, "smooth": smooth,
        "k_profile": K_PROFILE, "k_praz": K_PRAZ, "knot_days": KNOT_DAYS,
        "n_mild_days": len(dates), "rows_per_day": rpd,
        "rmse": float(np.sqrt(np.mean(resid**2))),
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


def night_curve(params: dict, tod: np.ndarray | None = None) -> np.ndarray:
    """Ranni clen pracovniho dne po volnu (odchylka od profilu typu dne)."""
    if tod is None:
        tod = np.arange(0, 24, 0.25)
    return _night_basis(tod) @ params["coef"][_slices(params["span_days"])["noc"]]


def light_curve(params: dict, tod: np.ndarray | None = None) -> np.ndarray:
    """Spotreba osvetleni za plne tmy podle casu dne (mimo okna nulova)."""
    if tod is None:
        tod = np.arange(0, 24, 0.25)
    return _light_basis(tod) @ params["coef"][_slices(params["span_days"])["svetlo"]]


def gloom_curve(params: dict, tod: np.ndarray | None = None) -> np.ndarray:
    """Spotreba navic pri plnem seru (nulovy osvit za dne) podle casu dne."""
    if tod is None:
        tod = np.arange(0, 24, 0.25)
    return _gloom_basis(tod) @ params["coef"][_slices(params["span_days"])["sero"]]


def summer_curve(params: dict, pos: np.ndarray, tod: float | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Prubeh leta podle dne v okne: (pracovni dny, volne dny). tod=None =
    denni prumer (jen konstantni slozka), jinak hodnota v danem case dne."""
    sl = _slices(params["span_days"])
    Sb = _summer_basis(np.asarray(pos, float))
    c = params["coef"][sl["leto_prac"]].reshape(N_SUMMER, 1 + 2 * K_SUMMER)
    phi = np.zeros(1 + 2 * K_SUMMER)
    phi[0] = 1.0
    if tod is not None and K_SUMMER:
        phi[1:] = _fourier(np.array([tod]), K_SUMMER)[0]
    return Sb @ (c @ phi), Sb @ params["coef"][sl["leto_vol"]]


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
    for k in ("span_days", "smooth", "rmse", "mae", "r2", "rows_per_day"):
        if k in out:
            out[k] = float(out[k])
    return out
