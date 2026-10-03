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
uv run python explore/vanoce_analyza.py   # 17   Vanoce 2022-2025 (rezidua po dnech a skupinach)
uv run python explore/ceps_analyza.py     # 19   zatizeni CR (CEPS/ENTSO-E) a predpoved CEPS
uv run python explore/zatop_test.py       # 25   opozdene zatopeni po teplem obdobi (test na reziduich)
uv run python explore/dokumentace.py      # dokumentace modelu v HTML (dokumentace/dokumentace.html)
uv run python explore/ocistena.py          # historie spotreby proti spotrebe bez topeni a chlazeni (graf)
uv run python explore/fve_teplota_test.py  # 39   teplotni koeficient FVE: profil fitu a rezidua za dne
uv run python explore/dlouhodoba_validace.py  # 31  dopredna validace predikce na rok dopredu
uv run python explore/dlouhodoba.py --pocasi=scenare  # 32  predikce na rok dopredu od konce dat
uv run python explore/intraday.py --meteo=predpoved --out=simulace/intraday_predpoved  # 33  intraday backtest (vydani po 15 min)
uv run python predikce.py                 # 34   provozni predikce od konce dat
```

Provozní predikce (`predikce.py`): spouští se po příchodu měření. Parametry
drží ve `provoz/parametry.pkl` a přefituje je podle stáří (lineární část po
dni, tvar po týdnu; `--prefit=tvar|linearni|ne`), počasí bere do konce dat
naměřené a dál z archivu předpovědí, korekci podle horizontu z koeficientů
intraday backtestu (`--korekce=soubor`). Výstup `provoz/predikce.csv` a
`provoz/archiv/`; `--konec="2026-07-20 12:00"` přehraje minulý okamžik.
S uloženými parametry trvá běh ~1 s, po aktualizaci vstupních CSV ~2,5 s
(soubory se čtou celé, takže zpětně zpřesněná měření se projeví sama;
skript vypíše, kolik intervalů se od minulého spuštění změnilo).

Jiná sada dat: `--projekt=cesta` u kteréhokoli skriptu (nebo proměnná
prostředí `MODEL_PROJEKT`) přepne pracovní složku — `Data/`, `Analyza/`,
`provoz/`, `simulace/`, `models/` a volitelně vlastní `config/` se pak berou
z ní (viz `src/workspace.py`). Archiv předpovědí počasí může být i
`Analyza/ArchivMeteo.csv`.

Transparentnost výpočtu:
- `config/model.yaml` — ručně zvolené konstanty a u fyzikálních parametrů
  popis, meze, start, prior a volitelně `pevna` (parametr se drží na zadané
  hodnotě ve všech fitech);
- `provoz/parametry.yaml` — aktuálně odhadnuté fyzikální parametry; ruční
  přepsání hodnoty se při příštím spuštění převezme a lineární část se
  přefituje (vydrží do příštího přefitu tvaru, trvale přes `pevna`);
- `provoz/rozklad.csv` a `provoz/vysvetleni.html` — rozklad predikce na složky
  (rozvrh a kalendář, osvětlení, šero, topení, chlazení, FVE, korekce),
  spotřeba očištěná o počasí a křivky modelu; `--bez-grafu`, `--otevrit`.

Dlouhodobá predikce (`explore/dlouhodoba.py`): `--pocasi=normal` (jedna
normálová dráha počasí, čtvrthodinový profil) nebo `--pocasi=scenare` (počasí
historických let, profil + pásmo nejistoty), `--trend=posledni|drzet`,
`--dni=365`, `--out=slozka`. Delší řadu počasí ze stejné stanice lze dodat
jako `Data/meteo_historie.csv` (formát `Meteo15.csv`) — víc scénářů a lepší
normál.

Volby backtestu: `--out=slozka`, `--meteo=predpoved` (realistický režim s
archivem předpovědí počasí v `Analyza/ArchivMeteo.xlsx`), `--bias=temp,wind`
(klouzavé odstranění biasu předpovědi; `--bias=` vypne), `--okno=dny`,
`--prior=vaha`, `--workers=N`, `--trenink=predpoved-bias` (experiment: učení
na předpovězeném počasí), `--plot-only`, `--no-open`.

```bash
uv run python explore/backtest.py --meteo=predpoved --out=simulace/predpoved_meteo
```

Výsledky simulace predikce (CSV + offline interaktivní HTML graf) jdou do
`simulace/` (mimo repozitář, obsahují data portfolia).

Diagnostické grafy se ukládají do `explore/`, parametry do `models/`
(`*_params.npz` ze sekvenčního odhadu, `*_joint.npz` ze závěrečného
společného fitu — pro produkční použití ty druhé).

Kalendář (svátky, mosty, prázdniny, vánoční skupiny dnů) se nastavuje v
`config/kalendar.yaml`; po změně skupin je potřeba přepočítat fity.

## Moduly

| modul | role |
|---|---|
| `src/etl.py` | načtení dat, časová osa, kalendář |
| `src/base.py` | bazální složka: P-spline úroveň + Fourierovy profily, osvětlení za tmy, šero přes den, průběh léta, zvláštní období |
| `src/heating.py` | topení: softplus + COP + dvoukanálový setrvačnostní filtr, vlastní denní tvar pro volné dny |
| `src/cooling.py` | chlazení: zrcadlový softplus + EER, vlastní denní tvar pro volné dny |
| `src/pv.py` | FVE: neklesající kapacita × osvit × teplotní derating |
| `src/fit.py` | společné odhady přes moduly (sdílené regresory) a joint fit |
| `src/normal.py` | klimatologie pro normálové podmínky |
| `src/validate.py` | fit z části dat (maska řádků) pro validaci mimo vzorek |
| `src/backtest.py` | simulace provozní predikce D+1 (denní přefit z dat do D 09:00) |
| `src/sun.py` | poloha slunce, tma (astronomicky) pro člen osvětlení a šero přes den (z osvitu) |
| `src/meteo_forecast.py` | archiv předpovědí počasí: osvit z oblačnosti a výšky slunce, kvantily jasnosti pro šero, klouzavý bias teploty a větru |
| `src/ote.py` | načítání dat OTE (zbytkové diagramy, KZD, přepočtené TDD), ČEPS a ENTSO-E z `Analyza/` |
| `src/kalendar.py` | kalendář z `config/kalendar.yaml`: svátky, mosty, prázdniny, zvláštní období (Vánoce) |
| `src/correction.py` | korekce predikce D+1 z chyb dřívějších predikcí (online, reziduová vrstva) |
| `src/longterm.py` | dlouhodobá predikce: budoucí osa, počasí normál / scénáře, pravidlo trendu, pásmo nejistoty |
| `src/intraday.py` | intraday predikce: vydání po 15 min, korekce závislá na horizontu (rezidua modelu + chyba předpovědi počasí) |
| `src/workspace.py` | pracovní složka: kde leží data, konfigurace, stav a výstupy (`--projekt=`) |
| `src/config.py` | konfigurace modelu z `config/model.yaml`: konstanty, meze, priory a pevné hodnoty fyzikálních parametrů |
| `src/explain.py` | rozklad predikce na složky, vlivy počasí, křivky modelu, offline graf |
| `src/operation.py` | provozní predikce od konce dat: přefit podle stáří parametrů, počasí, model + korekce |
