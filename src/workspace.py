"""Pracovni slozka: kde lezi data, konfigurace, stav a vystupy modelu.

Vychozi je koren repozitare (jako dosud). Jinou sadu dat vybere

    uv run python predikce.py --projekt=cesta/ke/slozce

nebo promenna prostredi MODEL_PROJEKT (volba na prikazove radce ma prednost;
relativni cesta se bere od aktualniho adresare). Plati pro vsechny skripty —
zpetne simulace ji predavaji i svym paralelnim procesum.

Struktura pracovni slozky (stejna jako v repozitari):

    Data/baseload.csv, Data/Meteo15.csv      vstupni mereni
    Analyza/ArchivMeteo.xlsx nebo .csv       archiv predpovedi pocasi
    config/model.yaml, config/kalendar.yaml  nepovinne — kdyz chybi, plati konfigurace z repozitare
    provoz/, simulace/, models/, dokumentace/  stav a vystupy (vzniknou samy)

Kod a sablony grafu zustavaji v repozitari.
"""
import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
ENV = "MODEL_PROJEKT"


def _root() -> Path:
    for a in sys.argv[1:]:
        if a.startswith("--projekt="):
            return Path(a.split("=", 1)[1]).resolve()
    return Path(os.environ[ENV]).resolve() if os.environ.get(ENV) else REPO


ROOT = _root()
if ROOT != REPO:
    if not ROOT.is_dir():
        raise SystemExit(f"pracovni slozka neexistuje: {ROOT}")
    os.environ[ENV] = str(ROOT)      # dedi podprocesy (paralelni simulace)
DATA = ROOT / "Data"
ANALYZA = ROOT / "Analyza"
MODELS = ROOT / "models"
PROVOZ = ROOT / "provoz"
SIMULACE = ROOT / "simulace"
DOKUMENTACE = ROOT / "dokumentace"


def config_file(name: str) -> Path:
    """Konfiguracni soubor z pracovni slozky; kdyz v ni neni, z repozitare."""
    own = ROOT / "config" / name
    return own if own.exists() else REPO / "config" / name


def describe() -> str:
    """Jednoradkovy popis pro vypis skriptu (prazdny pro vychozi slozku)."""
    if ROOT == REPO:
        return ""
    cfg = ", ".join(f"{n} {'vlastni' if config_file(n).parent.parent == ROOT else 'z repozitare'}"
                    for n in ("model.yaml", "kalendar.yaml"))
    return f"pracovni slozka: {ROOT} (konfigurace: {cfg})"
