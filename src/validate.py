"""Validace mimo vzorek — cely retezec kroku 3-6 odhadnuty jen z casti dat.

Dve otazky, dva rezimy:

- **vynechany rok** (interpolace): trenink = vse mimo testovany rok. Urovnova
  spline a rampa FVE testovany rok preklenou (penalizace druhych/prvnich
  diferenci), takze se testuje **struktura povetrnostni odezvy** — sedi
  topna/chladici krivka, setrvacnost a FVE i na roce, ktery model nevidel?
- **dopredna** (extrapolace): trenink = vse pred testovanym rokem. Uroven i
  kapacita FVE za koncem treninku **drzi posledni hodnotu** (spline i rampy
  jsou orezane na rozsah treninku), takze chyba obsahuje i neextrapolovany
  trend — to je realna situace predpovedi, bias se proto vykazuje zvlast.

Regresory (vcetne EMA filtru) se vzdy stavi na cele souvisle ose; maska
vybira jen radky, ktere vstupuji do odhadu.
"""
import numpy as np
import pandas as pd

import base
import cooling
import etl
import fit
import heating
import pv

HEAT_MAX_DAYTEMP = 17.0  # stejne kriterium jako explore/heating_fit.py


def year_local(df: pd.DataFrame) -> np.ndarray:
    return df.index.tz_convert(etl.TZ).year.to_numpy()


def fit_masked(df: pd.DataFrame, mask: np.ndarray) -> tuple[dict, dict, dict, dict]:
    """Staged fit (baze -> topeni -> chlazeni+FVE) a joint fit jen z radku mask."""
    bp = base.fit(df[mask])
    y = df["baseload"].to_numpy()
    resid = y - base.predict(df, bp)
    day_temp = df.groupby("date_local")["temp"].transform("mean").to_numpy()
    hp = heating.fit(df, resid, mask & (day_temp <= HEAT_MAX_DAYTEMP))
    resid = resid - heating.predict(df, hp)
    cp, pp = fit.fit_cooling_pv(df, resid, mask)
    return fit.joint(df, bp, hp, cp, pp, mask=mask)


def predict(df: pd.DataFrame, params: tuple[dict, dict, dict, dict]) -> np.ndarray:
    bp, hp, cp, pp = params
    return (base.predict(df, bp) + heating.predict(df, hp)
            + cooling.predict(df, cp) - pv.predict(df, pp))


def score(y: np.ndarray, yhat: np.ndarray) -> dict:
    e = y - yhat
    return {
        "rmse": float(np.sqrt(np.mean(e**2))),
        "bias": float(e.mean()),
        "rmse_bez_biasu": float(e.std()),
        "r2": float(1 - e.var() / y.var()),
    }


def summary_params(params: tuple[dict, dict, dict, dict], df: pd.DataFrame) -> dict:
    """Klicove parametry tvaru + kapacita FVE na konci osy (pro srovnani foldu)."""
    bp, hp, cp, pp = params
    out = {f"top_{k}": float(hp[k]) for k in heating.PARAM_NAMES}
    out.update({f"chl_{k}": float(cp[k]) for k in cooling.PARAM_NAMES})
    out["fve_gamma"] = float(pp["gamma"])
    end = pd.DatetimeIndex([df.index[-1]])
    out["fve_spicka"] = float(pv.capacity_curve(end, pp)[0] * 845)
    return out
