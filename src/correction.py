"""Korekce predikce D+1 z chyb drivejsich predikci — reziduova vrstva nad modelem.

Chyba predikce je silne setrvacna (autokorelace denniho biasu 0,77 mezi
sousednimi dny): drift urovne a sezonni pohyb baze, ktere fundamentalni model
za cutoffem jen drzi. Fyziku modelu korekce nemeni — je to samostatna aditivni
vrstva nad ulozenymi predikcemi.

Informace v okamziku vydani (D 10:00) pro cilovy den T = D+1:

- x1 = prumerna chyba dne D-1 (predikce vydana v D-2; den je uz cely zmereny),
- x2 = prumerna chyba dne D od 00:00 do 09:00 (predikce vydana v D-1),
- x3 = tvar chyby: odchylka chyby v danem case dne od denniho prumeru,
  prumer pres poslednich SHAPE_DAYS celych dni (D-SHAPE_DAYS .. D-1).

    korekce(T, cas) = a*x1 + b*x2 + c*x3(cas)

Nepravidelne dny (svatky, mosty, vanocni obdobi 21. 12. - 3. 1.) se jako zdroj korekce
nepouzivaji — jejich chyba je kalendarni anomalie, ne drift, a bez vyrazeni
se prenasela do nasledujicich dni (3.-4. 1. 2026 RMSE 66). Misto nich se
bere posledni pravidelny den; do odhadu koeficientu nevstupuji.

Koeficienty se odhaduji **online**: pro kazdy den vydani LS (s ridge ke nule)
jen z cilovych dni, jejichz chyba byla v okamziku vydani znama (T' <= D-1).
Prvnich WARMUP dni bez korekce.
"""
import numpy as np
import pandas as pd

import etl

SHAPE_DAYS = 7
WARMUP = 21        # minimalni pocet znamych cilovych dni pro odhad
RIDGE = 1e-3       # relativne ke stope X'X
FEATURES = ("x1", "x2", "x3")


def irregular_days(days: pd.DatetimeIndex) -> np.ndarray:
    """Svatky, mosty a vanocni obdobi 21. 12. - 3. 1. (model ma chybu 40-90 uz
    od 22. 12.; vlastni vanocni blok v modelu zatim neni)."""
    hol = etl._cz_holidays(days.year.min() - 1, days.year.max() + 1)
    special = hol | etl._bridges(hol)
    xmas = ((days.month == 12) & (days.day >= 21)) | ((days.month == 1) & (days.day <= 3))
    return np.asarray(days.isin(list(special)) | xmas)


def features(bt: pd.DataFrame, cutoff_hour: int = 9) -> pd.DataFrame:
    """bt = vystup backtestu (timestamp_mistni, skutecnost, predikce)."""
    t = pd.to_datetime(bt["timestamp_mistni"], utc=True).dt.tz_convert("Europe/Prague")
    out = pd.DataFrame({
        "den": pd.to_datetime(t.dt.date),
        "slot": t.dt.hour * 4 + t.dt.minute // 15,
        "e": (bt["skutecnost"] - bt["predikce"]).to_numpy(),
    }, index=bt.index)
    day = pd.Timedelta(days=1)
    axis = pd.date_range(out["den"].min(), out["den"].max(), freq="D")
    irr = pd.Series(irregular_days(axis), index=axis)
    em = out.groupby("den")["e"].mean().reindex(axis)
    mm = out[t.dt.hour < cutoff_hour].groupby("den")["e"].mean().reindex(axis)
    dev = out.assign(d=out["e"] - out["den"].map(em)).groupby(["den", "slot"])["d"].mean().unstack()
    dev = dev.reindex(axis)
    # nepravidelne dny ven, misto nich posledni pravidelny den
    em_r, mm_r = em.mask(irr).ffill(), mm.mask(irr).ffill()
    out["x1"] = (out["den"] - 2 * day).map(em_r)
    out["x2"] = (out["den"] - day).map(mm_r)
    # tvar: chyba po slotech minus denni prumer, pres poslednich SHAPE_DAYS
    # pravidelnych dni; pro cil T se bere stav ke dni T-2 (= D-1, posledni cely)
    shape = (dev[~irr.to_numpy()].rolling(SHAPE_DAYS, min_periods=SHAPE_DAYS - 2).mean()
             .reindex(axis).ffill().shift(2))
    out["x3"] = shape.stack().reindex(pd.MultiIndex.from_arrays([out["den"], out["slot"]])).to_numpy()
    out["nepravidelny"] = out["den"].map(irr).to_numpy()
    return out


def online(bt: pd.DataFrame, feats: tuple[str, ...] = FEATURES) -> tuple[np.ndarray, pd.DataFrame]:
    """Korekce pro kazdy radek bt a koeficienty po dnech vydani."""
    f = features(bt)
    cols = list(feats)
    X = f[cols].to_numpy()
    ok = ~np.isnan(X).any(axis=1)
    fit_ok = ok & ~f["nepravidelny"].to_numpy()   # anomalie neodhaduji koeficienty
    corr = np.zeros(len(f))
    coefs = []
    days = np.sort(f["den"].unique())
    den = f["den"].to_numpy()
    for T in days:
        known = fit_ok & (den <= T - np.timedelta64(2, "D"))
        target = (den == T) & ok
        if len(np.unique(den[known])) < WARMUP or not target.any():
            continue
        A, y = X[known], f["e"].to_numpy()[known]
        G = A.T @ A
        G += RIDGE * np.trace(G) / len(cols) * np.eye(len(cols))
        c = np.linalg.solve(G, A.T @ y)
        corr[target] = X[target] @ c
        coefs.append({"cil": pd.Timestamp(T), **dict(zip(cols, c))})
    return corr, pd.DataFrame(coefs).set_index("cil")
