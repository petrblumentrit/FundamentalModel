"""Nacitani dat OTE: zbytkove diagramy (ZD), jejich koeficienty (KZD) a
prepoctene typove diagramy dodavek (TDD). Soubory v Analyza/ (mimo repo).

Formaty se lisi podle obdobi: do 6/2024 hodinova data (sloupec "Hodina",
.xls/.xlsx), od 7/2024 ctvrthodinova ("Perioda"). Obe varianty jdou od mistni
pulnoci souvisle (dny zmeny casu maji 23/25 hodin, resp. 92/100 period), takze
UTC = mistni pulnoc + poradi. Vystup je vzdy hodinovy na UTC ose: ZD jako
energie za hodinu (soucet ctvrthodin), KZD a TDD jako prumer (bezrozmerne).

Sloupce: ZD/KZD — 5 = PRE (sit 0031); TDD — 6 = TDD4, 9 = TDD5 Praha,
15 = TDD6, 16 = TDD7, 17 = TDD8 (viz hlavicky souboru).
"""
from pathlib import Path

import pandas as pd

from etl import TZ

DIR = Path(__file__).resolve().parent.parent / "Analyza"


def _read(path: Path, columns: dict[int, str]) -> pd.DataFrame:
    raw = pd.read_excel(path, header=None)
    hdr = next(i for i in range(10) if raw.iloc[i, 0] == "Den")
    step = 15 if str(raw.iloc[hdr, 1]).startswith("Perioda") else 60
    body = raw.iloc[hdr + 1:].dropna(subset=[0])
    body = body[pd.to_datetime(body[0], errors="coerce").notna()]
    midnight = pd.to_datetime(body[0]).dt.tz_localize(TZ).dt.tz_convert("UTC")
    idx = midnight + pd.to_timedelta((body[1].astype(int) - 1) * step, unit="min")
    out = pd.DataFrame({name: body[c].astype(float).to_numpy() for c, name in columns.items()},
                       index=pd.DatetimeIndex(idx, name="utc"))
    out["krok_min"] = step
    return out


def load(prefix: str, columns: dict[int, str], energy: bool) -> pd.DataFrame:
    """Vsechny soubory prefix_*.xls(x), hodinove na UTC ose.

    energy=True: ctvrthodiny se scitaji (ZD), jinak prumeruji (KZD, TDD).
    Prekryvy souboru: plati posledni nacteny (serazeno podle nazvu).
    """
    parts = [_read(f, columns) for f in sorted(DIR.glob(f"{prefix}_*.xls*"))]
    hourly = []
    for p in parts:
        vals = p.drop(columns="krok_min")
        grp = vals.groupby(vals.index.floor("h"))
        hourly.append(grp.sum() if energy and p["krok_min"].iloc[0] == 15 else grp.mean())
    out = pd.concat(hourly)
    return out[~out.index.duplicated(keep="last")].sort_index()


def zd_pre() -> pd.Series:
    return load("ZD", {5: "pre"}, energy=True)["pre"]


def kzd_pre() -> pd.Series:
    return load("KZD", {5: "kzd_pre"}, energy=False)["kzd_pre"]


def tdd() -> pd.DataFrame:
    cols = {6: "TDD4", 9: "TDD5_Praha", 15: "TDD6", 16: "TDD7", 17: "TDD8"}
    return load("Prepoctene_TDD", cols, energy=False)
