"""Nacitani dat OTE: zbytkove diagramy (ZD), jejich koeficienty (KZD),
prepoctene typove diagramy dodavek (TDD) a zatizeni soustavy CEPS. Soubory v
Analyza/ (mimo repo).

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

import workspace

from etl import TZ

DIR = workspace.ANALYZA


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


def ceps(cerpani: bool = False) -> pd.Series:
    """Zatizeni soustavy CEPS [MW], hodinove prumery na UTC ose.

    Zdroj Analyza/CEPS.csv — export z portalu CEPS (hodinova agregace
    prumerem), sloupce "Zatizeni s cerpanim [MW]" a "Zatizeni [MW]" (bez
    cerpani precerpavacich elektraren, vychozi). Soubor je v mistnim case (dny
    zmeny casu 23/25 hodin, podzimni hodina dvakrat); hodiny jdou po rade od
    mistni pulnoci, takze UTC = mistni pulnoc + poradi hodiny ve dni.
    Kontrola: prumer 2025 7 732 MW vs cista spotreba CR dle CSU 6 932 MW
    (zatizeni obsahuje i ztraty).
    """
    c = pd.read_csv(DIR / "CEPS.csv", sep=";", skiprows=2, encoding="utf-8-sig")
    c = c.iloc[:, :3]
    c.columns = ["ts", "s_cerpanim", "zatizeni"]
    ts = pd.to_datetime(c["ts"], format="%d.%m.%Y %H:%M")
    day = ts.dt.normalize()
    midnight = day.dt.tz_localize(TZ).dt.tz_convert("UTC")
    idx = midnight + pd.to_timedelta(c.groupby(day).cumcount(), unit="h")
    col = "s_cerpanim" if cerpani else "zatizeni"
    return pd.Series(c[col].astype(float).to_numpy(), index=pd.DatetimeIndex(idx), name="ceps").dropna()


def entsoe_load() -> pd.DataFrame:
    """ENTSO-E Transparency, BZN|CZ: skutecne zatizeni a denni predpoved CEPS
    [MW], hodinove prumery na UTC ose (Analyza/GUI_TOTAL_LOAD_DAYAHEAD_*.csv).

    MTU je v CET/CEST (podzimni hodina oznacena "(CEST)"/"(CET)"); radky jdou
    po rade od mistni pulnoci, takze UTC = mistni pulnoc + poradi * krok.
    Do 2024 hodinove, od 2025 ctvrthodinove. Duplicitni soubory ("(1)") se
    ignoruji, prazdne hodnoty ("-", "n/e") jsou NaN.
    """
    parts = []
    for f in sorted(DIR.glob("GUI_TOTAL_LOAD_DAYAHEAD_*.csv")):
        if "(" in f.stem:
            continue
        c = pd.read_csv(f, na_values=["-", "n/e"])
        c.columns = ["mtu", "area", "actual", "da"]
        start = c["mtu"].str.slice(0, 16)
        t0 = pd.to_datetime(start, format="%d/%m/%Y %H:%M")
        day = t0.dt.normalize()
        end = pd.to_datetime(c["mtu"].str.extract(r"- (\d\d/\d\d/\d{4} \d\d:\d\d)")[0], format="%d/%m/%Y %H:%M")
        step = int((end.iloc[0] - t0.iloc[0]).total_seconds() // 60)
        midnight = day.dt.tz_localize(TZ).dt.tz_convert("UTC")
        idx = midnight + pd.to_timedelta(c.groupby(day).cumcount() * step, unit="min")
        parts.append(pd.DataFrame({"actual": c["actual"].astype(float).to_numpy(),
                                   "da": c["da"].astype(float).to_numpy()}, index=pd.DatetimeIndex(idx)))
    out = pd.concat(parts)
    out = out.groupby(out.index.floor("h")).mean()
    return out.dropna(how="all").sort_index()
