"""Normalove pocasi — dlouhodoby normal pro ocistenou spotrebu (krok 6).

Pro kazdou dvojici (den v roce, cas dne) prumer pres roky, vyhlazeny cirkularne
15dennim oknem. Vysledek je spojita rada na stejne casove ose jako data, takze
do ni lze dosadit setrvacnostni filtry topeni i chlazeni beze zmeny.
"""
import numpy as np
import pandas as pd

import etl

COLS = ("temp", "sun", "wind")
WINDOW = 15  # dnu vyhlazeni sezonniho chodu


def _smooth_circular(mat: pd.DataFrame) -> pd.DataFrame:
    """Klouzavy prumer pres dny v roce, prosinec navazuje na leden."""
    tri = pd.concat([mat, mat, mat], ignore_index=True)
    sm = tri.rolling(WINDOW, center=True, min_periods=1).mean()
    return sm.iloc[len(mat):2 * len(mat)].set_axis(mat.index)


def weather(df: pd.DataFrame) -> pd.DataFrame:
    """Kopie df s meteo sloupci nahrazenymi normalem."""
    loc = df.index.tz_convert(etl.TZ)
    key = pd.MultiIndex.from_arrays([loc.dayofyear, df["tod"].to_numpy()])
    out = df.copy()
    for c in COLS:
        mat = pd.Series(df[c].to_numpy(), index=key).groupby(level=[0, 1]).mean().unstack()
        out[c] = _smooth_circular(mat).stack().reindex(key).to_numpy()
    return out
