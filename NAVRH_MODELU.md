# Fundamentální model spotřeby elektrické energie — návrh

Datum: 2026-08-04. Shrnutí návrhové diskuse před zahájením vývoje.

## Cíl

Matematicko-fyzikální model spotřeby portfolia s **aditivní** dekompozicí (nikoli multiplikativní faktory). Výstupem je očištěná spotřeba nezávislá na venkovních vlivech — spotřeba za "normálových" podmínek. Aditivní tvar umožňuje příspěvky přidávat/ubírat a odpovídá reálné fyzice (přibude tepelné čerpadlo → přičte se člen; přibude FVE → odečte se).

## Data

- Spotřeba: 15min, celé portfolio (ideálně 2+ roky kvůli oběma sezónám).
- Meteo: 10min staniční měření (teplota, globální osvit, rychlost větru), převod na 15 min.
  - Převod časově váženým průměrem přes hranice intervalů, ne nearest-neighbor. Pohlídat konvenci značení intervalů (začátek vs. konec) — posun o interval vyrobí falešné zpoždění, které pak chybně kompenzuje setrvačnostní člen.
  - Osvit = energie za interval; teplota a vítr = průměr okamžitých hodnot — agregují se různě.
  - Portfolio → meteo jako vážený průměr více stanic s vahami podle rozložení odběru (stačí hrubě, např. po krajích). Jedna stanice je největší zdroj šumu v teplotní závislosti.
- Kalendář: české svátky, "mosty", školní prázdniny, přechod letního času, tarifní pásma (NT/VT, HDO).

## Struktura modelu

```
y(t) = B(t) + P_topení(t) + P_TUV(t) + P_chlazení(t) − G_FVE(t) + D_bat(t) + ε(t)
```

### B(t) — bazální složka

Denní profil (96 hodnot, vyhlazený) × typ dne: pracovní den, sobota, neděle/svátek, přechodové dny (pátek, mosty), prázdniny. Svátky mapovat na nedělní profil s vlastní korekcí. Vlastní sezónnost báze (osvětlení, chování) držet úsporně — pomalá spline nebo Kalman — aby "nesnědla" topný signál. Tarifně řízené spotřebiče (akumulačky, tarifní nabíjení baterií) patří sem, ne do povětrnostního bloku.

### Pocitová (efektivní) teplota

Startovní tvar — lineární kompozit:

```
T_ef = T + a·I − b·v·(T_ref − T)
```

- `a`, `b` se odhadují z dat (závisí na charakteru zástavby).
- Větrný člen vázaný na ΔT (infiltrace + konvekce), působí jen v topné sezóně.
- Alternativy pro validaci: sol-air temperature `T + α·I/h(v)`; neparametrická plocha (T, I, v) v GAM jako diagnostika tvaru.

### Tepelná setrvačnost — dvoukanálový filtr

Dvě populace topných jednotek s odlišnou dynamikou, ale stejnou statickou fyzikou:

- **Rychlý kanál**: ekvitermní řízení (čidlo venku) — spotřeba sleduje okamžitou T_ef.
- **Pomalý kanál**: vnitřní termostat — reakce přes tepelnou setrvačnost budovy, exponenciální vyhlazení `T̃_ef` s časovou konstantou τ (typicky 24–72 h).

```
T* = w·T_ef + (1−w)·T̃_ef
```

Jedna statická křivka `g(T*)`, dva parametry navíc (w, τ). Ekvivalent impulzní odezvy špička + exponenciální ocas. Plné dva aditivní členy `k₁·g(T_ef) + k₂·g(T̃_ef)` jen pokud diagnostika ukáže odlišné bilanční teploty populací (riziko kolinearity).

Identifikace w, τ: dynamické epizody (nástup/odeznění studené vlny), denní chod teploty. Pozor na konfundaci s k(čas dne) — identifikovat z rozdílů mezi dny ve stejnou denní hodinu, ne z průměrného denního profilu.

Validace: volná impulzní odezva spotřeby na teplotu (distributed lag 0–96 h s penalizací hladkosti) na reziduích — má vyjít špička + jeden exponenciální ocas. Dvě zřetelné časové škály → přidat třetí kanál.

### P_topení — topná křivka

Hladký hinge (softplus) — portfolio má rozmazanou bilanční teplotu:

```
Q_top = k(čas dne, typ dne) · s · ln(1 + exp((T_b − T*)/s))
P_el  = Q_top / COP(T)
```

- `T_b` ≈ 13–16 °C, `s` = šířka přechodu (rozptyl bilančních teplot v portfoliu).
- `k(čas dne, typ dne)` — časově proměnný koeficient (noční útlumy, pracovní doba, NT pásma); tvar křivky (T_b, s, COP) společný.
- **COP vždy na okamžitou T** (výparník je venku — rychlý i u jednotek řízených zevnitř), tepelná potřeba na filtrovanou T*. COP(T) zhruba lineární, případně Carnot × účinnostní faktor.
- Saturace při extrémním mrazu (konečný instalovaný výkon) — omezení shora.
- Známá položka pro rezidua: odmrazovací cykly (hrbol 0 až +5 °C při vysoké vlhkosti). V první verzi ignorovat.

### P_chlazení

Zrcadlový softplus, vlastní `T_b,chl` ≈ 20–24 °C, vlastní **kratší** setrvačnost (klimatizace reaguje rychleji), větší role osvitu, menší větru. Vynutit `T_b,topení < T_b,chlazení` (mrtvé pásmo ~15–20 °C, kde nepůsobí ani jedno).

### P_TUV

Slabá teplotní závislost (chladnější vstupní voda v zimě) + silný denní/týdenní profil. Malý sezónní člen + kalendářní profil.

### G_FVE — lokální výroba

Člen úměrný osvitu s **časově proměnnou kapacitou** (flotila roste): buď explicitně kapacita × specifická výroba(I, T), nebo regresní koeficient odhadovaný po klouzavých oknech / Kalmanem.

Pozor: osvit je sdílený regresor tří členů — solární zisky (↓ topení, jen v topné sezóně, přes setrvačnost), chlazení (↑, jen léto), FVE (↓ měřenou spotřebu, celoročně jen za světla, úměrně kapacitě). Rozlišují se strukturou; bez ní se slijí.

### D_bat — baterie

Zavádí **mezičasovou vazbu**: spotřeba večer/v noci závisí na osvitu z předchozích hodin. Dva mechanismy se stejnou signaturou: (a) po zataženém dni dobíjení ze sítě v noci; (b) po slunečném dni vybíjení tlumí večerní odběr, po zataženém ne.

- **Cílový model — ekvivalentní baterie se stavem nabití**: agregátní kapacita C(t) (roste v čase), výkonový limit, účinnost η ≈ 85–90 %, dispatch pravidlo (nabíjení z FVE přebytku, vybíjení proti večerní zátěži; volitelně noční dobíjení ze sítě podmíněné denním deficitem osvitu). SOC je spojitý stav — **půlnoc nevnímá**, energie ze slunečného dne se klidně vydává druhý den nad ránem.
- **Detekční test (udělat první)**: distributed lag na reziduích — večerní a ponoční spotřeba proti anomálii osvitu, lagy **přes půlnoc do ~24–36 h**. Záporná závislost večer + růst v čase kopírující instalace = efekt je měřitelný. Tvar odezvy ukáže, jak daleko za půlnoc vliv sahá.
- Fyzikální vazba: co se nabije, to se (× η) později vybije — baterie energii jen přesouvají. Na denní agregaci efekt téměř mizí, ale hranici dne položit na **"solární den"** (před svítáním, ~4:00–5:00, kde je SOC nejspíš vyčerpaný), ne na půlnoc. Zbytkový vícedenní carry-over = známá drobná chyba; kdyby rostl (větší baterie, spotová optimalizace), přejít plně na stavový model i pro fázi odhadu.
- Tarifně/spotově řízené nabíjení není povětrnostní — patří do bazálního bloku.

## Odhad parametrů — staged fitting

1. **Báze** — jen z období s T* v mrtvém pásmu (mírné jarní/podzimní dny = kontrolní skupina bez povětrnostních vlivů).
2. **Topení** — na reziduích v chladném období: softplus + COP + parametry pocitové teploty + w, τ.
3. **Chlazení** — na reziduích v teplém období.
4. **FVE** — osvitový člen s časově proměnnou kapacitou.
5. **Baterie** — detekční regrese, pak stavový model.
6. **Závěrečný joint fit** všech parametrů s warm startem z kroků 1–5 — **nevynechávat** (sdílené regresory, zejména osvit, nechávají po sekvenčním odhadu systematické chyby).

Topný člen identifikovat primárně z den-po-dni odchylek teploty od sezónního normálu (studená vlna vs. mírný týden ve stejném měsíci), ne z ročního chodu — jinak únik mezi sezónností báze a topením.

Nástroje: nelineární LS s omezeními (scipy) jako základ; bayesovská varianta (PyMC) pro priory na fyzikální parametry a nejistoty; GAM jednorázově jako diagnostika tvarů; Kalman / klouzavá okna pro časově proměnné složky (báze, kapacita FVE a baterií).

## Očištěná spotřeba

Dvě definice — vybrat podle použití:

- **Reziduální**: `y_norm = y − P_topení − P_chlazení + G_FVE − D_bat` → zbyde B(t) + ε.
- **Normálové podmínky** (obvyklý význam): `y_norm = y − P(skutečné počasí) + P(normálové počasí)` — povětrnostní členy se vyhodnotí při dlouhodobém normálu / TMY, včetně dosazení do stavové simulace baterií. Díky nulové denní bilanci baterií se jejich normalizace projeví jen v profilu dne, ne v denní energii.

## První kroky vývoje (po nahrání dat)

1. ETL: načtení spotřeby a meteo, převod 10→15 min, kontrola konvencí značení intervalů, vážení stanic. ✅
2. Explorace: scatter spotřeba × teplota po hodinách dne a typech dne; ověření mrtvého pásma; hrubý odhad T_b. ✅
3. Bazální model z mírných dnů. ✅
4. Topný modul (softplus + COP + dvoukanálový filtr), diagnostika volnou impulzní odezvou. ✅
5. Chlazení, FVE, bateriový detekční test. ✅
6. Joint fit, výpočet očištěné spotřeby, analýza reziduí (odmrazování, mosty, anomálie). ✅

## Stav vývoje

### 2026-08-05 — data, ETL, explorace, bazální model (kroky 1–3)

**Data.** Reálný rozsah 1.1.2022 – 4.8.2026 (≈ 4,6 roku, obě sezóny 4–5×) — dostatečné. Oba soubory mají osu až do konce 2028, ale za 4.8.2026 16:00 UTC je jen prázdná šablona (meteo NaN, baseline šum kolem nuly) — ETL ořezává podle posledního platného meteo. Sloupec spotřeby se v dodaných datech jmenuje `baseline` (ETL přejmenovává na `baseload`). Meteo obsahovalo 168 duplicitních + 168 chybějících 15min intervalů (shluky kolem přechodů času 2022–2023, 0,07 %) — duplicity se průměrují, krátké díry (≤ 1 h) interpolují.

**ETL** (`src/etl.py`): lokální čas → UTC (podzimní duplicitní hodina přes `ambiguous="infer"`), společná 15min osa, kalendář: typ dne 0 = Po–Čt, 1 = pátek + mosty (pracovní den sevřený mezi svátkem a vikendem/svátkem), 2 = sobota, 3 = neděle/svátek; příznak letních prázdnin (7–8). Výsledek: 160 964 řádků bez NaN, spotřeba 300–1063.

**Bazální model** (`src/base.py`, krok 3): `B(t) = úroveň(datum) + profil(čas dne, typ dne) + prázdninová korekce`. Úroveň = kubická P-spline (uzly po 90 dnech, penalizace druhých diferencí — hladké překlenutí zim bez mírných dnů); profily = Fourier 12 harmonických na typ dne (Po–Čt bez konstanty, úroveň nese spline; ostatní typy s offsetem); prázdniny = konstanta + 4 harmonické. Fit LS jen z mírných dnů (denní průměr T v mrtvém pásmu, předchozí den též mírný kvůli setrvačnosti). Parametry: `models/base_params.npz`; diagnostika: `explore/base_fit.py` (grafy 06–09).

Výsledky: 225 mírných dnů (137 Po–Čt), RMSE 14,9 (≈ 2,5 % úrovně), R² 0,986. Úroveň 594–620, mělké minimum 2023, bez artefaktů v zimních mezerách. Profily fyzikálně smysluplné: dvouvrcholový pracovní den, pátek s dřívějším odpolednim poklesem, víkend ≈ −15 %, prázdninová korekce −30 až −50 přes den. Rezidua vs. čas dne bez systematiky (krátké špičky ~5:00 — spínání/HDO, Fourier je nerozliší; pro v1 OK).

**Empirické korekce návrhu:**
- **Mrtvé pásmo v denních průměrech je 14–18 °C, ne 15–20** — chladicí větev nastupuje už od ~18 °C denního průměru (sklon nad 18 °C ≈ +9/°C; odpolední hrb s vrcholem ~15:00, v teplých dnech 20–26 °C reziduum odpoledne ~+67 vs. ~+21 v noci — jasná klimatizační signatura). Pozor na souřadnice: den s průměrem 18 °C má odpolední maximum 24–26 °C, takže okamžitá bilanční teplota T_b,chl nejspíš zůstává v návrhovém rozmezí 20–24 °C; „brzký" nástup je z velké části artefakt denního průměrování. Pro výběr mírných dnů (fit báze) ale platí užší pásmo.
- Zbytkový sklon reziduí v pásmu +1,6/°C (≈ 0,3 % báze přes šíři pásma) — pravděpodobně ocas chlazení/osvitu; ponecháno na joint fit.
- Po odečtení báze (`explore/09_baze_signal.png`): čistá topná větev ~+200 při −7 °C, chladicí ~+100 při +30 °C — vstup pro topný modul.

**Další krok:** 4 — topný modul (softplus + COP + dvoukanálový filtr, identifikace w, τ), diagnostika volnou impulzní odezvou.

### 2026-08-05 — topný modul (krok 4)

**Implementace** (`src/heating.py`): `P_topení = k(čas dne, typ dne) · g(T*) / COP(T)`. Pocitová teplota `T_ef = T + a·I − b·v·(20−T)⁺`; dvoukanálový filtr `T* = w·T_ef + (1−w)·EMA(T_ef; τ)` počítaný na celé souvislé ose; `g` = softplus, COP(T) = max(1+α·T, 0.2) na okamžitou teplotu, normalizace COP(0)=1 (škála v k). k = Fourier 4 harmonické + offsety typů dne. Odhad separabilně: vnější nelineární LS (scipy `least_squares`) jen přes 7 parametrů (a, b, T_b, s, w, τ, α), vnitřní lineární řešení k — model je v k lineární. Fit na reziduích po bázi, jen dny s denním průměrem T ≤ 17 °C (74 % řádků, bez chladicí kontaminace). Parametry: `models/heating_params.npz`; diagnostika: `explore/heating_fit.py` (grafy 10–13).

**Odhadnuté parametry:** T_b = 17,2 °C, s = 0,8 °C (ostřejší hinge, než návrh čekal), w = 0,32 / τ = 88 h — **dominuje pomalý kanál** (termostatické populace, τ nad horním okrajem návrhového odhadu 24–72 h), α_COP = 0,020 1/°C. **a = 0,033 °C/(W/m²) je nadhodnocené — osvitový člen zatím absorbuje i FVE** (sluneční den snižuje měřenou spotřebu); po zavedení FVE členu a joint fitu se přerozdělí. **b (vítr) = 0** — větrný efekt z jedné stanice nedetekovatelný; nechat v modelu, ověřit po vážení stanic.

**Výsledky:** chladné dny (průměr < 10 °C) RMSE 110,8 → 25,8; R² topného signálu 0,865. Topná křivka sedí na datech (mírná konvexita od COP), k(čas dne) dvouvrcholový (≈11:30 a 18:00, noc ~½ dne, neděle nejníž), nejchladnější týden (3.–12.1.2026, špičky > 1000) model kopíruje včetně nočních útlumů. Volná impulzní odezva zbytku: malá (±0,2 vs. topný sklon ~10/°C), se špičkami v násobcích 24 h (konfundace s denním profilem, návrh na ni upozorňoval) a mírně zápornou hodnotou v lagu 0 — žádná druhá výrazná časová škála, třetí kanál není potřeba; dořešit v joint fitu.

**Další krok:** 5 — chlazení (zrcadlový softplus, kratší setrvačnost), FVE člen s časově proměnnou kapacitou, bateriový detekční test (distributed lag přes půlnoc).

### 2026-08-05 — chlazení, FVE, baterie (krok 5)

**Chlazení a FVE fitovány společně** (`src/cooling.py`, `src/pv.py`, `src/fit.py:fit_cooling_pv`) — návrh zakazuje odhadovat je za sebou, protože osvit je jejich sdílený regresor. Úloha zůstává separabilní: vnější nelineární optimalizace přes 7 parametrů tvaru, vnitřní lineární krok přes k_c (chlazení) a δ (přírůstky kapacity FVE). Vnitřní krok se řeší přes normální rovnice s Choleskiho rozkladem, takže omezená úloha má rozměr počtu koeficientů (~30), ne počtu řádků (160k) — jinak by každé vyhodnocení procházelo maticí 160k × 30.

- **Chlazení**: T_bc = 19,0 °C, s_c = 1,3, w_c = 0,33, **τ_c = 40 h**, α_EER = 0,0024. Teplé dny (průměr > 20 °C) RMSE 58,1 → 14,9. Pozor: τ_c původně narazilo na horní mez 24 h — mez rozšířena na 72 h, výsledných 40 h je uvnitř. **Chlazení tedy není rychlejší než topení** (τ 84 h), jak návrh předpokládal; setrvačnost budovy je v obou směrech podobná, jen o něco kratší.
- **FVE**: kapacita parametrizována jako počáteční úroveň + součet **nezáporných** rampových přírůstků po 90 dnech (omezení v lineárním kroku), takže kapacita může jen růst — instalace se neodinstalovávají. Bez penalizace vyšla kapacita jako schodiště, proto přidána penalizace prvních diferencí přírůstků (tempo instalací se mění hladce).
- `a_c` (osvit v pocitové teplotě chlazení) i γ (teplotní derating panelů) skončily **na nule** — obojí je degenerované vůči FVE členu ve stejném regresoru; joint fit je jediné místo, kde se to dá rozdělit.

**Baterie — detekční test** (`explore/battery_test.py`): cílen na hodiny 18–05, kdy FVE už nevyrábí, takže cokoliv, co tam osvit vysvětluje, je přesunutá energie. Impulzní odezva přes lagy 0–36 h (přes půlnoc) i roční vývoj večerní citlivosti.

**Výsledek: efekt baterií není měřitelný.** Odezva v lagu 0 je −0,0027 a součet přes lagy 1–36 je −0,003 na W/m²; roční sklony večerní citlivosti (−0,006 / −0,011 / −0,002 / −0,006 / −0,002) jsou všechny v mezích šumu a **bez rostoucího trendu**, který by kopíroval instalace. Ani jedna z podmínek, které návrh stanovil jako důkaz měřitelnosti, není splněna.

Test má prokazatelně sílu: pozitivní kontrola (stejný odhad na denních hodinách a reziduu bez odečtené FVE) najde známý efekt −0,016 na W/m², tedy 6× větší než ta hranice, na které bateriový signál mlčí. **Člen `D_bat` proto ve v1 vypuštěn** — spolu s ním odpadá stavový model SOC i hranice „solárního dne". Test opakovat, až portfolio vyroste o baterie nebo spotovou optimalizaci.

### 2026-08-05 — joint fit a očištěná spotřeba (krok 6)

**Joint fit** (`src/fit.py:joint`): warm start z kroků 1–5, nelineárních parametrů 14 (7 topení + 6 chlazení + γ), všechny koeficienty úrovní a profilů se řeší lineárně. Bazální blok na nelineárních parametrech nezávisí, takže se jeho návrh i Gramova matice počítají jednou a každé vyhodnocení staví jen bloky topení, chlazení a FVE (160k × 44 místo 160k × 174).

**Celkově: RMSE 21,5 → 18,7, R² = 0,984.**

Co joint fit přerozdělil — přesně ta systematická chyba, kvůli které ho návrh označuje za nevynechatelný:

| parametr | staged | joint |
|---|---|---|
| topení `a` (osvit) | 0,0333 | 0,0260 |
| topení T_b / s | 17,2 / 0,79 | 18,0 / 2,31 |
| topení α (COP) | 0,0203 | 0,0107 |
| chlazení T_bc / s_c | 19,0 / 1,34 | 18,2 / 2,30 |
| **špičkový výkon FVE** | **26** | **63** |

**Kapacita FVE se víc než zdvojnásobila** — sekvenční odhad ji podhodnocoval, protože osvitový člen topení odčerpával část solárního signálu. Sklon reziduí na osvit klesl z −0,017 na −0,004 na W/m² (o 75 %). Oba softplusy se rozšířily (s 0,8 → 2,3; s_c 1,3 → 2,3), což odpovídá „rozmazané bilanční teplotě" portfolia, kterou návrh předpokládal — staged fit je držel uměle ostré.

**Nutná oprava při joint fitu:** penalizace hladkosti bazální úrovně musí růst s počtem řádků. Fit z ~200 mírných dnů má 21 600 řádků, joint fit 161 000, takže stejná konstanta je tu 7,5× slabší — úroveň se rozvlnila a začala chytat sezónnost, kterou má nést topení (přesně to, před čím návrh varuje). Po přeškálování je úroveň hladká.

**Očištěná spotřeba** — obě definice z návrhu, normálové počasí z 15denně vyhlazené klimatologie (den v roce × čas dne, `src/normal.py`):

| rok | měřeno | reziduálně | na normál |
|---|---|---|---|
| 2022 | 632,4 | 555,6 | 629,4 |
| 2023 | 616,5 | 544,3 | 615,5 |
| 2024 | 627,9 | 556,9 | 626,1 |
| 2025 | 641,4 | 569,4 | 637,4 |
| 2026 | 658,7 | 582,0 | 646,2 |

Reziduální řada je po odečtení povětrnostních členů plochá (mizí sezónní chod), normálová si sezónnost ponechává — odstraňuje jen **anomálii** počasí, ne roční chod; to je záměr, ne chyba. Rok 2026 je částečný (do 4. 8.), takže jeho průměr není srovnatelný s celými roky.

**Zbývající známé chyby v reziduích:**
- **Prosinec RMSE 39,5** proti 13–20 ve zbytku roku — vánoční týden mezi svátky není modelovaný (svátky se mapují na neděli, ale celý blok 24. 12.–1. 1. má vlastní režim). Největší jednotlivá rezerva modelu.
- Zbytkový sklon reziduí na teplotu +0,03/°C a na osvit −0,004/(W/m²) — malé, ale nenulové.
- `b` (vítr) = 0 a `a_c` = 0, γ = 0 stále na mezích; vítr má smysl zkusit až po vážení více stanic.
- Odmrazovací cykly tepelných čerpadel (hrb 0 až +5 °C) zůstávají neošetřené podle návrhu.

**Další možné kroky:** kalendářní blok pro vánoční období; vážení více meteostanic (mělo by odemknout větrný člen); bayesovská varianta (PyMC) pro nejistoty parametrů; validace mimo vzorek (rok stranou).
