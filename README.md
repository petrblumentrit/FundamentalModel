# Fundamentální model spotřeby elektrické energie

Aditivní fyzikální dekompozice spotřeby portfolia na 15minutových datech:

```
y(t) = B(t) + P_topení(t) + P_chlazení(t) − G_FVE(t) + ε(t)
```

Výstupem je očištěná spotřeba za normálových podmínek. Architektura, průběh
vývoje a interpretace výsledků jsou v **[NAVRH_MODELU.md](NAVRH_MODELU.md)**.

## Prostředí

Projekt používá [uv](https://docs.astral.sh/uv/). Závislosti jsou uzamčené
v `uv.lock`, takže prostředí se vytvoří jedním příkazem:

```bash
uv sync
```

Skripty se pak spouštějí přes `uv run`, bez ruční aktivace venv:

```bash
uv run python explore/joint_fit.py
```

## Data

Vstupy se čekají v `Data/` (mimo repozitář — obsahují data portfolia):

| soubor | obsah |
|---|---|
| `baseload.csv` | 15min spotřeba portfolia, sloupec `timestamp;baseline` |
| `Meteo15.csv` | 15min meteo, sloupce `timestamp;sun;wind;temp` |

Oba v lokálním čase, oddělovač `;`, desetinná čárka. ETL je převádí na UTC,
ořezává na poslední platné meteo a doplňuje kalendář (typy dne, svátky, mosty,
prázdniny).

## Postup

Skripty na sebe navazují — každý čte parametry předchozího z `models/`:

```bash
uv run python explore/explore.py          # 1-2  explorace, ověření mrtvého pásma
uv run python explore/base_fit.py         # 3    bazální složka z mírných dnů
uv run python explore/heating_fit.py      # 4    topný modul
uv run python explore/cooling_pv_fit.py   # 5    chlazení + FVE (společný fit)
uv run python explore/battery_test.py     # 5    detekční test baterií
uv run python explore/joint_fit.py        # 6    joint fit + očištěná spotřeba
uv run python explore/validation.py       # 7    validace mimo vzorek (--plot-only z cache)
uv run python explore/backtest.py         # 8    simulace predikce D+1 za posledni rok (--plot-only)
uv run python explore/svetlo_test.py      # 9    test spotreby rizene svetlem (zmena casu)
uv run python explore/zd_analyza.py       # 14   srovnani s domacnostmi PRE (data OTE v Analyza/)
```

Výsledky simulace predikce (CSV + offline interaktivní HTML graf) jdou do
`simulace/` (mimo repozitář, obsahují data portfolia).

Diagnostické grafy se ukládají do `explore/`, parametry do `models/`
(`*_params.npz` ze sekvenčního odhadu, `*_joint.npz` ze závěrečného
společného fitu — pro produkční použití ty druhé).

## Moduly

| modul | role |
|---|---|
| `src/etl.py` | načtení dat, časová osa, kalendář |
| `src/base.py` | bazální složka: P-spline úroveň + Fourierovy profily |
| `src/heating.py` | topení: softplus + COP + dvoukanálový setrvačnostní filtr |
| `src/cooling.py` | chlazení: zrcadlový softplus + EER |
| `src/pv.py` | FVE: neklesající kapacita × osvit × teplotní derating |
| `src/fit.py` | společné odhady přes moduly (sdílené regresory) a joint fit |
| `src/normal.py` | klimatologie pro normálové podmínky |
| `src/validate.py` | fit z části dat (maska řádků) pro validaci mimo vzorek |
| `src/backtest.py` | simulace provozní predikce D+1 (denní přefit z dat do D 09:00) |
| `src/sun.py` | poloha slunce a tma (astronomicky) pro člen osvětlení |
| `src/ote.py` | načítání dat OTE (zbytkové diagramy, KZD, přepočtené TDD) z `Analyza/` |
| `src/correction.py` | korekce predikce D+1 z chyb dřívějších predikcí (online, reziduová vrstva) |
