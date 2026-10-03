"""Konfigurace modelu z config/model.yaml: rucne zvolene konstanty a fyzikalni
parametry (meze, start, prior, pevna hodnota). Moduly si z ni pri importu
nacitaji sve konstanty; kalendar ma vlastni soubor (src/kalendar.py)."""
from functools import lru_cache
from pathlib import Path

import numpy as np
import yaml

CONFIG = Path(__file__).resolve().parent.parent / "config" / "model.yaml"
SECTIONS = ("topeni", "chlazeni", "fve")   # oddily s fyzikalnimi parametry


@lru_cache(maxsize=1)
def model() -> dict:
    return yaml.safe_load(CONFIG.read_text(encoding="utf-8"))


def spec(section: str) -> dict:
    """Fyzikalni parametry oddilu: nazev -> {popis, jednotka, meze, start, prior, pevna}."""
    return model()[section]["parametry"]


def bounds(section: str, names: list[str]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(dolni mez, horni mez, start) v poradi names."""
    sp = spec(section)
    return (np.array([sp[n]["meze"][0] for n in names], float),
            np.array([sp[n]["meze"][1] for n in names], float),
            np.array([sp[n]["start"] for n in names], float))


def priors(section: str) -> dict[str, tuple[float, float]]:
    """Slabe priory: nazev -> (stred, sm. odchylka)."""
    return {n: (float(p["prior"][0]), float(p["prior"][1])) for n, p in spec(section).items() if p.get("prior")}


def fixed(section: str, names: list[str]) -> np.ndarray:
    """Pevne hodnoty v poradi names; NaN = parametr se odhaduje."""
    sp = spec(section)
    return np.array([np.nan if sp[n].get("pevna") is None else float(sp[n]["pevna"]) for n in names])
