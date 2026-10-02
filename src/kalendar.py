"""Kalendar z config/kalendar.yaml: svatky, mosty, prazdniny, zvlastni obdobi.

Vse pracuje s kalendarnimi daty (datetime.date), ne s radky dat — funguje
tedy i na podmnozinach (mirne dny, maska treninku).
"""
from datetime import date, timedelta
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

CONFIG = Path(__file__).resolve().parent.parent / "config" / "kalendar.yaml"


@lru_cache(maxsize=1)
def config() -> dict:
    return yaml.safe_load(CONFIG.read_text(encoding="utf-8"))


def _md(s: str) -> tuple[int, int]:
    m, d = s.split("-")
    return int(m), int(d)


def easter(year: int) -> date:
    """Velikonocni nedele (anonymni gregoriansky algoritmus)."""
    a, b, c = year % 19, year // 100, year % 100
    d, e = b // 4, b % 4
    f = (b + 8) // 25
    g = (b - f + 1) // 3
    h = (19 * a + b - d - g + 15) % 30
    i, k = c // 4, c % 4
    m = (32 + 2 * e + 2 * i - h - k) % 7
    n = (a + 11 * h + 22 * m) // 451
    month = (h + m - 7 * n + 114) // 31
    day = (h + m - 7 * n + 114) % 31 + 1
    return date(year, month, day)


def _easter_days(year: int) -> dict[str, date]:
    e = easter(year)
    return {k: e + timedelta(days=v) for k, v in config()["svatky"]["velikonoce"].items()}


def holidays(y0: int, y1: int) -> set[date]:
    """Statni svatky (volno) v letech y0..y1."""
    cfg = config()["svatky"]
    out = set()
    for y in range(y0, y1 + 1):
        out.update(date(y, *_md(s)) for s in cfg["pevne"])
        out.update(_easter_days(y).values())
    return out


def shop_closed(y0: int, y1: int) -> set[date]:
    """Dny se zakazem prodeje ve velkych obchodech."""
    out = set()
    for y in range(y0, y1 + 1):
        ed = _easter_days(y)
        for s in config()["svatky"]["zakaz_prodeje"]:
            out.add(ed[s] if s in ed else date(y, *_md(s)))
    return out


def _off(d: date, hol: set) -> bool:
    return d.weekday() >= 5 or d in hol


def bridges(hol: set) -> set[date]:
    """Pracovni dny v behu nejvys `max_pracovnich_dnu`, sevrenem volnem z obou
    stran, pricemz aspon jedna strana je svatek (ne jen vikend)."""
    return set(bridge_runs(hol))


def bridge_runs(hol: set) -> dict[date, int]:
    """Most -> delka behu pracovnich dni, do ktereho patri (1 nebo 2 ...)."""
    n_max = config()["mosty"]["max_pracovnich_dnu"]
    out = {}
    for h in hol:
        for direction in (1, -1):          # beh pracovnich dni za / pred svatkem
            run, d = [], h + timedelta(days=direction)
            while not _off(d, hol) and len(run) <= n_max:
                run.append(d)
                d += timedelta(days=direction)
            if 0 < len(run) <= n_max and _off(d, hol):
                out.update({x: len(run) for x in run})
    return out


def _in_range(d: date, od: str, do: str) -> bool:
    a, b = _md(od), _md(do)
    md = (d.month, d.day)
    return a <= md <= b if a <= b else (md >= a or md <= b)   # pres Novy rok


def period_groups() -> list[tuple[str, int]]:
    """[(nazev skupiny, harmonickych)] pres vsechna obdobi, v poradi configu."""
    out = []
    for groups in (config().get("obdobi") or {}).values():
        for name, g in groups.items():
            out.append((name, int(g.get("harmonicke", 0))))
    return out


def period_flags(dates: pd.Series, hol: set) -> dict[str, np.ndarray]:
    """Pro kazdou skupinu obdobi priznak po radcich (dates = kalendarni data)."""
    uniq = pd.Series(pd.unique(dates))
    bridge_days = bridges(hol)
    out = {}
    for groups in (config().get("obdobi") or {}).values():
        for name, g in groups.items():
            sel = set()
            for d in uniq:
                if "dny" in g:
                    hit = d.strftime("%m-%d") in g["dny"]
                elif "rozsahy" in g:
                    hit = any(_in_range(d, a, b) for a, b in g["rozsahy"])
                else:
                    hit = _in_range(d, g["od"], g["do"])
                if not hit or (g.get("jen_mosty") and d not in bridge_days):
                    continue
                kind = g.get("dny_v_tydnu", "vse")
                off = _off(d, hol)
                if kind == "vse" or (kind == "pracovni" and not off) or (kind == "vikend" and off):
                    sel.add(d)
            out[name] = dates.isin(sel).to_numpy()
    return out


def summer_holidays(dates: pd.Series) -> np.ndarray:
    p = config()["prazdniny"]["letni"]
    sel = {d for d in pd.unique(dates) if _in_range(d, p["od"], p["do"])}
    return dates.isin(sel).to_numpy()
