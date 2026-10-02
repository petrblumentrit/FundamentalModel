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

### 2026-09-15 — validace mimo vzorek (krok 7)

**Metoda** (`src/validate.py`, `explore/validation.py`): celý řetězec kroků 3–6 (báze z mírných dnů → topení → chlazení + FVE → joint fit) odhadnutý jen z **masky řádků**; regresory včetně setrvačnostních filtrů se stavějí na celé ose, do normálních rovnic vstupují jen trénovací řádky. Dva režimy, každý odpovídá na jinou otázku:

- **vynechaný rok** (2023, 2024, 2025) — trénink = vše mimo testovaný rok. Úroveň i kapacita FVE mezeru překlenou penalizacemi, takže se testuje **struktura povětrnostní odezvy** (topná/chladicí křivka, setrvačnost, FVE) na roce, který model neviděl.
- **dopředná** (2024, 2025, 2026) — trénink = vše před testovaným rokem. Úroveň i kapacita FVE za koncem tréninku **drží poslední hodnotu** (spline i rampy ořezané na rozsah tréninku), takže chyba obsahuje i neextrapolovaný trend — realistická situace predikce; bias se vykazuje zvlášť.

**Výsledky** (RMSE 15min hodnot; in-sample = joint fit ze všech dat):

| rok | in-sample | vynechaný rok | bias | dopředná | bias | dopředná bez biasu |
|---|---|---|---|---|---|---|
| 2023 | 18,0 | 20,6 | +8,5 | – | – | – |
| 2024 | 18,9 | 19,8 | −3,6 | 22,8 | −9,5 | 20,7 |
| 2025 | 20,3 | 24,0 | −9,7 | 21,8 | 0,0 | 21,8 |
| 2026 (do 4. 8.) | 16,6 | – | – | 19,0 | −2,8 | 18,8 |

R² mimo vzorek 0,978–0,985. Grafy 22 (měsíční RMSE), 23 (stabilita parametrů), 24 (průběh testovaných roků), 25 (úroveň a kapacita FVE po foldech).

**Interpretace:**

- **Struktura se nepřefituje.** Chyba bez biasu je 18,7–22 proti 18–20 in-sample. Parametry tvaru jsou přes foldy stabilní: T_b 17,5–18,0 °C, s 1,8–2,9, w 0,34–0,39, T_bc 17,1–18,7 °C, s_c 2,1–2,4, τ_c 35–43 h. Nejhůř určený je τ topení (63–87 h; dopředné foldy z 2–3 let dávají kratší hodnoty) — pomalý kanál potřebuje víc zim.
- **Bias vynechaného roku je celý z přemostění úrovně**, ne z povětrnostních členů: P-spline s penalizací druhých diferencí přes roční mezeru přestřelí (minimum 2023 na 560 místo 570 → +8,5; letní hrb 2025 na 620 → −9,7, graf 25). Pro odhad chyby predikce je směrodatný dopředný režim.
- **Dopředná predikce:** RMSE 19–23. Bias −9,5 v roce 2024 vzniká držením úrovně ze dna roku 2023 při tréninku jen ze dvou let; v letech 2025 a 2026 je bias 0 až −3. Mimo prosinec a leden je měsíční RMSE 15–19 proti 13–16 in-sample.
- **Nález a oprava — koncová rampa FVE.** Dopředné foldy odhalily neidentifikovanou poslední rampu kapacity: v posledním čtvrtletí před koncem tréninku (říjen–prosinec, málo slunce) dostala skok +30 až +40 jednotek a špička FVE vyšla 97 místo 63 — v plném fitu (konec v srpnu) ji kotví letní slunce, v predikci z konce roku ne. Dvě opravy: (a) `pv.capacity_basis` počítá jen rampy pozorované v rozsahu celé, za koncem rozsahu kapacita drží; (b) penalizace tempa instalací `PV_SMOOTH` 0,05 → 5 (test 0,05 / 0,5 / 5 / 50: in-sample RMSE 18,66 → 18,68, dopředná 2026 20,0 → 19,0, špička 97 → 69; nad 5 už nic nepřidá). Plný fit se prakticky nemění (špička 63, parametry tvaru do ±0,1), kapacita je nyní hladká konkávní křivka místo schodů.

**Zbývající chyby dopředného režimu:** leden 2026 (mrazová vlna: RMSE 29 vs 25 in-sample, špičky podhodnocené o ~30 — topná křivka za nejchladnějším T* v tréninku), prosinec (43–46 shodně ve všech režimech — chybějící vánoční blok nezávisí na tom, kolik dat model viděl), neextrapolovaný trend úrovně.

**Závěr:** in-sample RMSE 18,7 není přefitování; realistická chyba predikce na rok dopředu je **19–22 (R² ≈ 0,98)** plus bias úrovně 0–10 podle vývoje trendu. Rezervy v pořadí velikosti: vánoční blok, extrapolace trendu úrovně (hold vs. lineární pokračování — validovat dopředně), lednové extrémy.

**Další kroky:** vánoční kalendářní blok v bázi; pravidlo pro trend úrovně v predikci; bootstrap nejistot přes bloky týdnů; nová data po 4. 8. 2026 (další léto pro chlazení a FVE).

### 2026-10-01 — simulace provozní predikce D+1 (krok 8)

**Režim** (`src/backtest.py`, `explore/backtest.py`): v den D v 10:00 místního času se vydává predikce na celý den D+1 (96 × 15 min) z dat končících v D 09:00 (H-1; časové značky jsou začátky intervalů, poslední známý interval 08:45–09:00). Model se **každý den přefituje** z dat do cutoffu: bazální fit z mírných dnů + joint fit s warm startem nelineárních parametrů z předchozího dne (první den bloku celý řetězec kroků 3–6). Úroveň i kapacita FVE za cutoffem drží. Cutoff i vydání jsou na místních hodinách, i ve dny změny času. **Meteo je skutečně naměřené** (dokonalá předpověď počasí), v datech předpověď není — v provozu k chybě přibude chyba meteo předpovědi. Klouzavý poslední rok platných dat: cílové dny 4. 8. 2025 – 3. 8. 2026 (365 dní), 10 paralelních bloků, 35 min.

Výstupy (mimo repo, `simulace/`): `predikce_D1.csv` (čas vydání, konec dat, skutečnost, složky báze/topení/chlazení/FVE, predikce), `parametry.csv` (parametry po dnech vydání), `predikce_D1.html` (offline interaktivní graf se zoomem, plotly.js vložený v souboru).

**Výsledky:** RMSE 21,1, MAE 14,7, bias +0,2, MAPE 2,23 %, R² 0,981; denní průměry RMSE 15,5 (MAPE 1,6 %). Bez 23. 12. – 1. 1. RMSE 18,2 — shodné s in-sample. Odpovídá dopředné validaci (19–22).

| měsíc | 08 | 09 | 10 | 11 | 12 | 01 | 02 | 03 | 04 | 05 | 06 | 07 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| RMSE | 14,0 | 13,2 | 18,3 | 23,4 | 42,6 | 26,5 | 17,1 | 14,7 | 18,0 | 16,9 | 14,5 | 14,6 |
| bias | +2,7 | +1,9 | +8,5 | +12,8 | −8,1 | +9,8 | −7,5 | −3,5 | −4,6 | −7,7 | −0,1 | −2,1 |

- **Vánoční blok** dominuje: všech 8 nejhorších dní je 23. 12. – 1. 1. (denní RMSE 53–94; 25. a 31. 12. chyba denního průměru ~80).
- **Bias po měsících ±8–13** i při denním přefitu: podzim podstřeluje (+8,5 / +12,8), jaro přestřeluje (−5 až −8). Hypotéza (neověřeno): konec P-spline úrovně za cutoffem jen drží a nestíhá sezónní pohyb báze; stejně dobře může jít o strukturální chybu přechodových měsíců v topné křivce. Ověřit rozkladem biasu po hodinách dne a T*; kandidát na korekci je aditivní kotva úrovně z reziduí posledních dní před cutoffem.
- **Parametry při denním přefitu stabilní:** T_b 17,3–18,0, w 0,36–0,39, τ 75–85 h, τ_c 34–40 h, špička FVE 56–78 (roste s létem).

### 2026-10-01 — ráno po volnu, mosty, korekce z chyb predikce (kroky 8a, 8b)

**Diagnostika** (rezidua plného fitu po hodinách dne): profil Po–Čt byl kompromis dvou nočních režimů — pondělní noc 0–5 h o 7–10 níž než model, noci Út–Čt o 2–5 výš; v pondělí odchylka doznívá v 8–9 h. Typ dne se láme o půlnoci, ale **režim se mění až během rána** (noc z neděle na pondělí je ještě víkendová). Pátek a sobota bez systematiky. Pracovní dny po svátku jsou o 20–40 níž celý den: mosty (pátek po čtvrtečním svátku, lidé si berou volno) a 27. 12. / 2. 1. (Vánoce).

**8a — báze** (`src/base.py`, `src/etl.py`):
- ETL: příznak `po_volnu` (předchozí kalendářní den víkend/svátek, z kalendáře, ne z řádků dat) a `most`.
- **Ranní člen po volnu**: pracovní den po víkendu/svátku dostane aditivní člen z kubických B-splin na 0–10 h (6 parametrů), v 10 h plynule odezní (hodnota i sklon nulové). Odpovídá myšlence posuzovat noční a denní režim zvlášť, ale úsporně — denní profily zůstávají, přidává se jen přechod. Odhad: −12 až −15 v noci, −7 až −9 v 4–8 h, 0 od 10 h.
- **Most**: celodenní aditivní posun (konstanta + 2 harmonické); bez něj celodenní propad mostů stahoval ranní člen pondělí. Vyšel hluboký (až −73), protože 2 ze 7 mostů jsou vánoční (27. 12. 2024, 2. 1. 2026) — rozdělí ho až vánoční blok.
- In-sample RMSE 18,68 → 18,26; po hodinách je pondělní i úterní–čtvrteční noc ±1.

**8b — korekce z chyb dřívějších predikcí** (`src/correction.py`): chyba predikce je silně setrvačná (autokorelace denního biasu 0,77). Samostatná aditivní vrstva nad uloženými predikcemi, fyziku modelu nemění. V D 10:00 pro cíl D+1: x1 = průměrná chyba dne D−1, x2 = průměrná chyba rána D do 09:00, x3 = tvar chyby podle času dne (odchylka od denního průměru za posledních 7 pravidelných dní); korekce = a·x1 + b·x2 + c·x3. Koeficienty **online** (LS s ridge, jen z cílových dní, jejichž chyba byla v okamžiku vydání známá; prvních 21 dní bez korekce). Ustálené a ≈ 0,54, b ≈ 0,2, c ≈ 0,8.
- **Nepravidelné dny** (svátky, mosty, 21. 12. – 3. 1.) se jako zdroj korekce nepoužívají a neodhadují koeficienty; místo nich poslední pravidelný den. Bez toho se anomálie svátků přenášela do dalších dní (3.–4. 1. 2026 RMSE 66). Užší okno 24. 12. – 1. 1. by v prosinci pomohlo díky setrvačnosti uvnitř Vánoc (Vánoce 42 vs 61), ale mimo Vánoce je horší (14,65 vs 14,40) — Vánoce má řešit vlastní blok v modelu, ne korekce.

**Výsledky backtestu** (RMSE 15min; Vánoce = cílové dny 21. 12. – 3. 1.):

| | celý rok | mimo Vánoce | Vánoce |
|---|---|---|---|
| v1 model | 21,10 | 17,86 | 60,1 |
| v1 + korekce | 18,51 | 14,37 | 61,3 |
| **v2 model** (ranní člen + most) | 20,84 | 17,73 | 58,7 |
| **v2 + korekce** | **18,29** | **14,32** | 59,9 |

- Ranní člen: bias 0–6 h v pondělí −12,4 → −4,5 (stejně jako ostatní dny), RMSE pondělí 0–9 h 17,0 → 14,3. Na celkovém RMSE malý (týká se 1/7 dní, 9 h).
- Korekce: −2,6 RMSE celkem, −3,4 mimo Vánoce; nejvíc pomáhá v listopadu (23,0 → 15,9) a lednu (24,5 → 20,5), kde model drží úroveň za cutoffem.
- Oba kroky se doplňují: korekce tvaru (x3) je průměr přes týden a typ dne nerozliší — na v1 nechávala pondělní noc −6,4 a středeční +5,7; na v2 jsou všechny dny ±2.
- **Zbývá společný noční bias všech dní** (model −3 až −6 v 0–6 h; korekce ho odstraní): sezónní změna tvaru profilu (květen −15, leden +9 na v1), kandidát na krok 3 z plánu — nejdřív ověřit, zda nesedí v nočním poměru topného k(čas dne).

### 2026-10-01 — spotřeba řízená světlem: test a člen osvětlení (krok 9)

**Hypotéza** (uživatel): spotřebu v čase řídí dva procesy — **rozvrh** (práce, školy, obchody; drží se hodin) a **denní světlo** (veřejné i vnitřní osvětlení; drží se slunce). Úsporné zdroje by druhý proces zmenšovaly.

**Test bokem modelu** (`explore/svetlo_test.py`, grafy 26–27; spotřeba očištěná o počasí, Út–Čt):
- **Změna času jako přirozený experiment** — přes noc se hodiny posunou o hodinu vůči slunci, počasí, sezóna i rozvrh zůstanou. Rozdíl profilu 14 dní po a před změnou minus stejný rozdíl o 14 dní dřív (placebo bez změny ±6). Na jaře spotřeba klesne až o ~55 přesně v hodině mezi starým a novým západem slunce (~18:15–19:25), na podzim vzroste až o ~50 (~16:40–17:45), ráno zrcadlově menší efekt (podzim −18 v 6:45–7:40). Zbytek dne se nehne. **Potvrzeno:** soumraková složka průměrně ~43 jednotek (~7 % zátěže v té hodině), jednotlivé změny 21–62.
- **Měsíční profily**: celkový tvar dne se kryje mnohem líp na ose hodin než na ose slunce — dominuje rozvrh; v říjnu až březnu na něm sedí večerní hrb, který na ose slunce leží ve všech měsících 0–1 h po západu.
- **Úsporné zdroje**: pokles soumrakové složky 2022–2026 neprokázán (rozptyl mezi změnami větší než jakýkoli trend). Výkyv jaro 2024 nejspíš kvůli Velikonocům ve stejný den jako změna času (neověřeno).
- Regrese „tma × aktivita“ přes celý den se neosvědčila (v poledne tma nenastává, křivka neurčená) — vyřazena.

**Člen osvětlení v bázi** (`src/sun.py`, `src/base.py`): `L(t) = tma(t) · aktivita(čas dne)`. Tma je astronomická (výška slunce, NOAA, střed ČR; logistický přechod kolem −2°, průměr přes 15min interval) — deterministická, v predikci nepotřebuje meteo. Aktivita = kubické B-spliny jen v oknech, kde se tma během roku mění (3–8,5 h a 15–23 h), na okrajích plynule nulové (12 parametrů). Úpravy po prvním fitu:
- ranní okno zkráceno z 9,5 na 8,5 h — v 8–9 h je tma jen v prosinci a lednu, tedy o vánočních prázdninách, a člen vyšel záporný (−49, chytal Vánoce);
- v joint fitu omezení ≥ 0 (osvětlení jen přidává);
- ve fitu báze z mírných dnů ridge na osvětlení — jaro a podzim nepokrývají hodiny, kdy je tma jen v zimě; bez něj neurčené koeficienty rozbily postupný fit topení (RMSE 94).

**Výsledky:** osvětlení za plné tmy ~50 v 16–19 h, ~40 ve 20–21 h, ~15 v 5–7 h. Soumraková složka v reziduích nového modelu při změnách času průměrně −0,3 (dřív ~45). **Topný k(čas dne) ztratil večerní špičku v 18 h** (je teď přes den plochý) — topení dřív neslo osvětlení; T_b 17,8 → 16,1 °C, špička FVE 64 → 71. In-sample RMSE 18,26 → 17,27; sezónní vzorec reziduí (listopad 15–18 h +20 → +10, říjen 18–21 h +19 → +5, duben–červen večer −10 → ~0) zmizel z velké části.

| backtest RMSE | celý rok | mimo Vánoce | 15–21 h | Vánoce |
|---|---|---|---|---|
| v2 model | 20,84 | 17,73 | 23,2 | 58,7 |
| **v3 model** (+ osvětlení) | 20,26 | 16,99 | 20,6 | 58,8 |
| v2 + korekce | 18,29 | 14,32 | – | 59,9 |
| **v3 + korekce** | 18,11 | 14,09 | – | 59,8 |

Korekce z v3 získá méně, protože část soumrakového vzorce dřív chytala sama (tvar chyby za posledních 7 dní). Noční bias modelu (−5 v 0–6 h) zůstává — to je drift tvaru přes roky (návrh B), ne osvětlení.

**Další kroky:** drift tvaru přes roky — nejdřív časová změna topné složky (tepelná čerpadla místo plynu po roce 2022), pak obecný drift tvaru; vánoční blok (uživatel dodá starší data); validace mimo vzorek (krok 7) pro aktuální verzi modelu.

### 2026-10-01 — trend topné citlivosti, role dnů v bloku volna, zrychlení backtestu (kroky 10–12)

**10 — časová změna topné citlivosti** (`src/heating.py`, `src/fit.py`): `P_top = [k(čas dne, typ dne) + Δk(t)] · g(T*)/COP(T)`. Δk(t) je aditivní přírůstek od začátku dat („přibylá topná zátěž“ — tepelná čerpadla místo plynu po roce 2022; proti tomu zateplování, které model zachytí implicitně v čistém efektu). Rampy po 90 dnech jako kapacita FVE, **bez omezení znaménka**, penalizace změn tempa (instalace mohou saturovat); za koncem dat drží, trend se neextrapoluje. Odhaduje se jen v joint fitu (lineární krok).
- Výsledek: **topná citlivost vzrostla od zimy 2021/22 do 2025/26 o ~25 %** (zimy +5,8 / +12,9 / +18,4 / +24,6 %), robustní na penalizaci (0,5 / 5 / 50 → 26,8 / 24,7 / 23,6 %), téměř lineárně ~6 %/rok, saturace zatím nevidět (graf 28). Úroveň báze za stejnou dobu vzrostla jen o ~6 %, takže nejde jen o růst portfolia.
- Noční drift tvaru (noc rok od roku relativně níž) zůstal a je **stejný i v létě** — nesouvisí s topením, samostatný jev k řešení.
- In-sample 17,27 → 17,10; backtest model mimo Vánoce 16,99 → 15,78 (leden 25,2 → 23,2, únor 16,5 → 14,2, listopad 21,2 → 18,5).

**11 — typ dne = role v bloku volna** (`src/etl.py`): role podle toho, zda je volno dnes a zítra — 0 práce→práce (Po–Čt), 1 práce→volno (pátek; i čtvrtek před pátečním svátkem, mosty), 2 volno→volno (sobota; i svátek před víkendem, Velký pátek až Velikonoční neděle), 3 volno→práce (neděle; i Velikonoční pondělí, svátek uprostřed týdne). Volno včera nese ranní člen po volnu. Svátky tak přebírají profily naučené z ~240 víkendů, bez nových parametrů; běžné týdny beze změny.
- Motivace (uživatel): čtvrtek 30. 4. 2026 před svátkem měl odpoledne jiný tvar. Diagnostika: všech 23 pracovních dní před svátkem (mimo Vánoce) mělo odpoledne a večer −15 až −25; po změně ~0 (zbývá −11 ve 12–15 h), 30. 4. večer −49 → −18.
- Svátky v pracovní den jsou dopoledne stále o +12 až +16 výš než profil soboty (1. 5. 2026 +39) → příště aditivní korekce konkrétních svátků se smršťováním (hierarchicky: svátek obecně + konkrétní svátek + den v týdnu) a příznak zákazu prodeje (zavřené obchody jen o některých svátcích).
- Vánoce (22. 12. – 3. 1.) se budou modelovat samostatně: pracovní dny mezi svátky, čisté svátky 25.–26. 12. a specifický večer 24. 12. (příprava večeře → špička, večeře → pokles, televize → nárůst).
- In-sample 17,10 → 16,71. Backtest (13 dní kolem svátků mimo Vánoce): model 23,3 → 18,6, s korekcí 22,3 → 17,2; ostatní dny 15,43 → 15,26.

**12 — zrychlení backtestu** (`src/backtest.py`, `explore/backtest.py`): nelineární parametry tvaru se přefitují **jednou týdně** (pondělí) a první den simulace, lineární část denně (`fit.joint(..., fix_shape=True)`, ~2 s místo ~60 s). Ověřeno na modelu v4: RMSE 19,95 vs 19,87 při denním nelineárním přefitu, predikce se liší v průměru o 0,7. Výpočet je třífázový a paralelní: (0) celý řetězec k prvnímu dni simulace v hlavním procesu, (1) tvar pro každé pondělí — warm start z fáze 0 (data před simulací, žádný únik z budoucnosti), pondělí nezávislá, nejdřív nejvzdálenější, 24 workerů, (2) lineární přefit a predikce pro každý den s tvarem posledního pondělí. Výsledek nezávisí na počtu workerů. **84 min → 10,3 min** (fáze 0: 2 min, fáze 1: 6,3 min, fáze 2: ~2 min); úzké hrdlo je fáze 1 (fity tvaru 23–258 s, sdílená propustnost paměti).

**Souhrn backtestu** (RMSE 15min, poslední rok, Vánoce = 21. 12. – 3. 1.):

| verze | model | mimo Vánoce | + korekce | + korekce mimo Vánoce |
|---|---|---|---|---|
| v1 (výchozí) | 21,10 | 17,86 | 18,51 | 14,37 |
| v2 ráno po volnu, mosty | 20,84 | 17,73 | 18,29 | 14,32 |
| v3 osvětlení | 20,26 | 16,99 | 18,11 | 14,09 |
| v4 trend topení | 19,87 | 15,78 | 18,39 | 13,84 |
| **v5 role dnů** | **19,47** | **15,36** | **18,00** | **13,21** |

**Další kroky:** noční drift tvaru přes roky (není z topení — je i v létě); korekce konkrétních svátků (metoda 2); vánoční blok (uživatel dodá starší data); validace mimo vzorek (krok 7) pro aktuální verzi.

### 2026-10-01 — délka trénovacího okna, data OTE (PRE), trend topení s denním tvarem (kroky 13–15)

**13 — délka trénovacího okna** (`src/backtest.py`, `--okno=dny`): trénink jen z posledních N dní, svátky, mosty a Vánoce (22. 12. – 3. 1.) i s ±3 dny okolí z celé historie (bez okolí by ve starších letech chyběla opora pro úroveň). Backtest (RMSE / MAPE modelu; s korekcí mimo Vánoce): celá historie 19,47 / 1,99 % (13,21), 3 roky 19,37 / 1,98 % (13,39), 2 roky 19,67 / 2,01 % (13,72), 1 rok 20,20 / 2,12 % (14,40). **Okno nezkracovat** — zima potřebuje dlouhou paměť (leden 20,2 vs 24,0 při 1 roce), ani léto kratší okno nezlepší. Noční bias zůstal i při okně 1 rok → nevzniká průměrováním starých let.

**Rozklad chyby backtestu**: chyba modelu fitovaného ze všech dat je na stejných dnech prakticky stejná jako chyba predikce → chyba je **strukturální**, ne predikční (výjimka podzim: +4 až +6 zpoždění nástupu topné sezóny). Fit ze všech dat a fit k cutoffu se liší v rozdělení báze vs topení až o ~20 i v létě (topná složka má letní chvost ~14) — na predikci vliv nemá, zkresluje očištěnou spotřebu; řešit zvlášť.

**14 — data OTE pro PRE** (`src/ote.py`, `explore/zd_analyza.py`, grafy 29–30; data v `Analyza/` mimo repo): zbytkový diagram (ZD, domácnosti + veřejné osvětlení), jeho koeficienty KZD (poměr ZD k očekávání podle TDD přepočtených na skutečnou teplotu) a přepočtené TDD, 1. 1. 2022 – 31. 5. 2026 (do 6/2024 hodinově, pak 15min; vše převedeno na hodiny v UTC).
- Tvar dne portfolia koreluje s domácnostmi PRE 0,76–0,85; topná citlivost domácností je ~2× vyšší.
- **Dlouhodobý drift tvaru je malý**: portfolio očištěné o počasí mění podíl noci a večera o ≤ 1,5 p. b. (hlavně 2022 → 2024, od 2024 stabilní); domácnosti posouvají noc → večer hlavně na jaře.
- **Topná citlivost** (sklon, index první zima = 100) 2021/22 → 2025/26: domácnosti PRE 100 → 128, trend v modelu 100 → 125 (sedí); hrubý sklon portfolia 179 je zkreslený mrazivým lednem 2026 (nelinearita, COP).
- **KZD roste**: domácnosti spotřebovávají čím dál víc nad TDD (2022 ~0,96 → 2024–25 ~1,03 → I–III/2026 1,05–1,11), nejvíc v zimě — souhlasí s růstem tepelných čerpadel.
- **Energetická krize 2022**: KZD na podzim 2022 0,92 — domácnosti šetřily; vysvětluje mělké minimum úrovně 2023. Část odhadovaného růstu topné citlivosti od 2022 může být návrat po krizovém šetření (pro predikci nevadí, pro výklad a saturaci ano). Uživatel zmiňuje i přechod na přímotopy kvůli cenám plynu — data ho výrazně neukazují (KZD v zimě 2022/23 nevzrostl), mohl být menší nebo krátký.
- TDD7 (domácnosti přímotop/TČ): citlivost na teplotu v noci ~0,7 denní — stejný poměr jako topný k(čas dne) v modelu. Noční útlum TČ (nižší bilanční teplota v noci) z TDD neověřitelný (OTE přepočítává lineárně); v reziduích je jen mírný podstřelený ranní náběh (+3 až +7 v 6–9 h v topných dnech).
- KZD pro predikci nepoužitelné (publikováno zpětně), vhodné pro sledování vývoje.

**15 — trend topení s vlastním denním tvarem** (`src/heating.py`, `src/fit.py`): Δk(t, čas dne) = Σ_j φ_j(čas dne) · Δk_j(t), φ = [1, 2 harmonické], každá složka vlastní rampy se stejnou penalizací změn tempa (80 koeficientů, lineární krok joint fitu). **Přírůstek topné zátěže je soustředěný odpoledne a večer**: na konci dat noc +17 až +19 %, den +25 %, 15–21 h +30 až +35 % průměrného k. In-sample 16,71 → 16,58.

| backtest | RMSE | MAPE | mimo Vánoce | + korekce | + korekce MAPE | + korekce mimo Vánoce |
|---|---|---|---|---|---|---|
| v5 role dnů | 19,47 | 1,99 % | 15,36 | 18,00 | 1,74 % | 13,21 |
| **v6 trend s tvarem** | **19,31** | **1,92 %** | **14,61** | **17,88** | **1,70 %** | **12,89** |

Zimní bias má vyrovnaný tvar (dřív noc −0,3 / odpoledne +15,5, nyní +3 až +9 přes celý den) — zbývá úroveň: **Vánoce (20. 12. – 6. 1.) jsou in-sample o −23,6 pod modelem a stahují zbytek zimy (+2,8)**, protože model vyrovnává zimu jako celek. Podzimní odpoledne +11 (část zpoždění nástupu sezóny), jarní noc −4,3.

**Další kroky:** vánoční blok (pracovní dny 22. 12. – 3. 1., svátky 25.–26. 12., Štědrý večer: příprava jídla → špička, večeře → pokles, televize → nárůst) — pomůže i zbytku zimy; letní chvost topení vs úroveň; validace mimo vzorek pro aktuální verzi.

### 2026-10-02 — kalendář v konfiguraci, vánoční blok, prior na tvar, data ČEPS (kroky 16–19)

**16 — konfigurace kalendáře** (`config/kalendar.yaml`, `src/kalendar.py`): svátky (pevné + Velikonoce posunem od neděle), dny se zákazem prodeje (zatím jen informativně), mosty, prázdniny a zvláštní období se skupinami dnů. ETL i báze berou kalendář z konfigurace; změna skupin mění strukturu báze (přepočítat fity).
- **Mosty rozšířeny na 1–2 pracovní dny** sevřené volnem (aspoň z jedné strany svátek) — pokrývá např. Po–Út před středečním Štědrým dnem (22.–23. 12. 2025) nebo Čt–Pá po středečním Novém roce (2.–3. 1. 2025). Jednodenní a dvoudenní mosty mají oddělené celodenní členy, mosty uvnitř zvláštních období se z nich vyjímají.

**17 — vánoční blok** (analýza `explore/vanoce_analyza.py`, grafy 31–32): Vánoce 2022–2025 jsou mezi roky velmi konzistentní — 24. 12. −7,4 % ±0,5 (propad v 19–20 h až −155 = večeře, návrat do půlnoci = televize), 25. 12. −11 %, 26. 12. −8 %, 27.–30. 12. pracovní ~−10 % (víkend poloviční), 31. 12. −7 % ±4,4 (závisí na dni v týdnu), 1. 1. −8 % ±0,5 (v noci +40 oslavy, dopoledne −120). Okrajové dny (21.–23. 12., 2.–3. 1.) jsou nízko jen jako mosty. Skupiny v konfiguraci (každá aditivní člen s vlastním denním tvarem): `stedry_den` (8 harmonických), `bozi_hod`, `sv_stepan`, `mezi_svatky_prac`, `mezi_svatky_vikend`, `silvestr` (přičítá se k mezi svátky), `novy_rok` (8 harm.) — 75 parametrů. Vánoce dřív stahovaly celou zimu (in-sample −23,6 o Vánocích, zbytek zimy +2,8); po zavedení odhad růstu topné citlivosti +23 → +34 %.
- In-sample RMSE 16,58 → 13,14 (prosinec 38 → 17, leden 18,4 → 13,1).

**18 — prior na parametry tvaru a paměť workerů** (`src/fit.py`, `heating.PRIOR`, `cooling.PRIOR`, `explore/backtest.py`): s vánočními členy a trendem topení s denním tvarem se tvar topné křivky v joint fitu posouval po téměř plochém údolí (šířka s → horní mez 6, T_b dolů, COP → 0) a fity tvaru trvaly až 86 min (i samostatně > 28 min). Slabý prior (T_b 17 ± 1,5 °C, s 2 ± 1, τ 80 ± 20 h, α 0,012 ± 0,008; T_bc 18 ± 1,5, s_c 2 ± 1 — z validace mimo vzorek): 1 sm. odch. stojí ~1 % součtu čtverců, data rozhodují, prior ukotví ploché směry → fit 75–81 s, 12 iterací. Pojistka max. 100 iterací při warm startu. Worker backtestu zabere ~1,3 GB — při 24 workerech na 32 GB systém odkládal na disk (~20× zpomalení); počet workerů nyní podle volné paměti (~1,4 GB/worker, `--workers=N`). Backtest 14,3 min.

| backtest | RMSE | MAPE | mimo Vánoce | Vánoce | + korekce | + korekce MAPE | + korekce mimo Vánoce |
|---|---|---|---|---|---|---|---|
| v6 trend s tvarem | 19,31 | 1,92 % | 14,61 | 66,1 | 17,88 | 1,70 % | 12,89 |
| **v7 Vánoce + prior** | **14,48** | **1,67 %** | **13,85** | **25,6** | **13,09** | **1,49 %** | **12,52** |

Vánoce 2025/26 (naučené jen z 2022–2024), RMSE dne v6 → v7: 24. 12. 84 → 24, 25. 12. 112 → 33, 29.–30. 12. 77 → 8–10, 31. 12. 86 → 15, 1. 1. 75 → 21; leden 23,0 → 14,2. Zbývá: dvoudenní most 22.–23. 12. (~36), most 2. 1. 2026 se zhoršil (19 → 32), duben 12,5 → 13,9 (ověřit vliv prioru).

**19 — zatížení ČR (ČEPS, ENTSO-E)** (`src/ote.py`: `ceps()` z CSV portálu ČEPS, `entsoe_load()` z ENTSO-E Transparency; `explore/ceps_analyza.py`): CSV ČEPS je v CET/CEST, hodinově, MW (průměr 2025 7 732 MW vs čistá spotřeba ČSÚ 6 932 MW — zatížení obsahuje ztráty); ENTSO-E (skutečné zatížení + denní předpověď ČEPS, do 2024 hodinově, od 2025 15min) ~0,92× CSV ČEPS, korelace 0,996. Denní předpověď ČEPS má hodinové MAPE 1,2–1,5 %.
- Model na zatížení ČR (jedna meteostanice — topení jen orientačně, parametry na mezích): **krize 2022 celostátně −5,4 %** (portfolio −2,7 %), ČR se dosud nevrátila, portfolio od 2025 nad 2022; **růst topné citlivosti je specifický pro portfolio** (ČR bez jasného trendu, portfolio +24/+39/+50 % noc/poledne/18 h); portfolio je citlivější na osvětlení (9,8 vs 5,8 % průměru v 17 h); FVE v ENTSO-E zatížení nevidět (zahrnuje spotřebu krytou decentrální výrobou).
- Nevysvětlené odchylky portfolia a ČR korelují 0,32–0,33. „Překvapení ČEPS“ (předpověď ČEPS / model ČR z dat před backtestem − 1) jako vstup korekce: RMSE 13,09 → 12,85, ale jen díky Vánocům, mimo Vánoce 12,52 → 12,60, MAPE beze změny → **zatím nepřidávat**. Backtest ale používá skutečné meteo, předpověď ČEPS předpovězené — přínos jeho předpovědi počasí tu nejde změřit; vrátit se s archivem meteo předpovědí a ověřit čas zveřejnění.

**Další kroky:** mosty kolem Vánoc (22.–23. 12., 2. 1.); vliv prioru na jaro; letní chvost topení vs úroveň; archiv meteo předpovědí (realistický backtest, přínos ČEPS); validace mimo vzorek pro aktuální verzi.

### 2026-10-02 — zrychlení backtestu 14,3 → 5,0 min (krok 20)

Profil: fit tvaru vyhodnocuje model ~150–260× (numerická Jacobiho matice, 13 parametrů); vyhodnocení ~1,2 s, z toho ~0,4 s opakované skládání bloků nezávislých na parametrech. Načtení dat 11–32 s v každém procesu. Denní krok ~2 s (matice báze 3×, predikce na celé historii).

Úpravy (A–D beze změny výsledků, ověřeno porovnáním — predikce shodné do 4·10⁻⁷):
- **A** mezipaměť `etl.load()` v `Data/.etl_cache.pkl` (klíč: čas změny a velikost vstupních dat, `config/kalendar.yaml` a zpracovávajícího kódu) — 0,03 s místo 11–24 s;
- **B** v `fit.joint` bloky nezávislé na nelineárních parametrech (denní tvary topení/chlazení, rampy FVE, trend topení) spočítané jednou, matice V předalokovaná — sestavení bloků 0,39 → 0,11 s na vyhodnocení;
- **C** v `base.design` Fourierovy funkce jednou (nižší řády jako výřezy) — sestavení báze ~2×;
- **D** složky predikce D+1 jen z posledních 60 dní (`backtest.PREDICT_HISTORY`; setrvačnostní filtry τ ~85 h zapomenou počátek za ~2 týdny) — denní krok 2,0 → 1,3 s;
- **E** fit tvaru (fáze 0 a 1) z **hodinové agregace** (`etl.hourly`, 4× méně řádků); lineární část a predikce zůstávají 15min. Setrvačnostní filtry počítají s krokem dat (`etl.step_hours`), penalizace hladkosti úrovně a výběr mírných dnů s počtem řádků na den. Ověřeno: parametry tvaru shodné na setiny (T_b 16,70/16,67 °C, τ 81,6/81,6 h), přesnost predikce 7 dní shodná (12,60/12,60; 8,71/8,71), plný fit 68 → 13 s.

Backtest: fáze 0 136 → 57 s, fáze 1 531 → 112 s, fáze 2 ~190 → ~130 s; RMSE / MAPE beze změny (14,48 / 1,67 %, s korekcí 13,09 / 1,49 %), predikce se liší v průměru o 0,13. Zbývající rezerva hlavně ve fázi 2 (matice báze se skládá pro každý den znovu).

### 2026-10-02 — mosty kolem Vánoc (krok 21)

Diagnostika (in-sample): mosty v roce model zvládá (reziduum ~0, obecné členy jednodenní −39 / dvoudenní −13), mosty kolem Vánoc jsou hlubší — reziduum −8 až −23 navíc, celkový efekt ~−35 až −57 na den (mezi běžným mostem a režimem mezi svátky −70; lidé mají volno celý týden). Nová skupina v `config/kalendar.yaml`: `most_vanoce` = mosty ve dnech 20.–23. 12. a 2.–6. 1. (nové volby skupin `rozsahy` a `jen_mosty`), konstanta + 2 harmonické; tyto dny se vyjímají z obecného mostového členu. Pokrývá 23. 12. 2024, 2.–3. 1. 2025, 22.–23. 12. 2025, 2. 1. 2026.

Backtest v8: RMSE 14,30 / MAPE 1,66 % (v7 14,48 / 1,67 %), s korekcí 12,93 / 1,48 %; Vánoce 25,7 → 23,1, mimo Vánoce beze změny. 22. 12. 2025 37 → 16, 23. 12. 35 → 12. **2. 1. 2026 beze změny (32)** — model přestřeluje celý den (748 vs 724, ráno −51) i následující víkend 3.–4. 1. (17–21), takže jde spíš o nižší úroveň po Novém roce (školní prázdniny do 2. 1., čtyřdenní volno od čtvrtka) než o samotný most; vánoční mostový člen se zatím učí jen ze tří dnů 2024/25 (mělčích). Zpřesní se s dalšími Vánoci.

### 2026-10-02 — ověření prioru, graf bez WebGL (krok 22)

**Prior a duben** (backtest v8 bez prioru, `--prior=0`): celkově s priorem mírně lépe (RMSE 14,30 vs 14,39, MAPE 1,66 vs 1,67 %, s korekcí 12,93 vs 12,98); duben s priorem 14,0 vs 13,6 bez něj — prior vysvětluje jen ~0,4 ze zhoršení dubna proti v6 (12,5), zbytek přinesl vánoční blok (přerozdělení jara mezi topení 78 → 50 a bázi 554 → 581; v dubnu model podstřeluje den o +10 až +14 v 9–15 h a odečítá víc FVE). Bez prioru parametry v průběhu roku ujíždějí (s 1,6–6,0 °C, T_b 14,0–16,4 °C, α až 0, T_bc 16,1–18,7 °C), s priorem stabilní (T_b 16,6–16,8, s 2,0–2,1); fáze 1 backtestu 2 vs 9,5 min. **Prior ponechán.** Dubnové podstřelení dne — samostatný úkol (dělení FVE vs topení na jaře, denní tvar trendu topení).

**Graf backtestu**: čáry přes SVG místo WebGL (`scattergl`) — v prohlížeči bez dostupného WebGL (vypnutá akcelerace, ovladač) se čáry nevykreslily, fungoval jen tooltip.

### 2026-10-02 — realistický backtest s předpovědí počasí, ČEPS znovu, slunce a klimatizace (krok 23)

**Archiv předpovědí počasí** (`Analyza/ArchivMeteo.xlsx`, mimo repo; `src/meteo_forecast.py`): hodinově, CET/CEST, teplota, **oblačnost [%]**, vítr, od 4. 8. 2022 (včetně následujících dní → předpověď). Shoda se skutečností: teplota MAE 1,2 °C (bias +0,55 v období backtestu), vítr korelace 0,76 (bias +1,7 m/s — model na vítr prakticky necitlivý). Osvit z oblačnosti: osvit za jasné oblohy (Haurwitz, z výšky slunce) × poměr jasnosti kt(oblačnost) **kalibrovaný empiricky** z dat před backtestem (Kasten–Czeplak vztah u střední oblačnosti nezachytí; měřená stanice: jasno kt 0,67, zataženo 0,31 v hodinových binech). Odvozený osvit: hodinová korelace 0,84, denní 0,94, denní chyba ~28 %. **Ověřeno u uživatele:** archiv se ukládá v D 10:00 na celý následující den D+1 (24 hodnot) — backtest pro D+1 používá předpověď vydanou v D 10:00, pro zbytek dne D předpověď z D−1 10:00; obě jsou v okamžiku vydání dostupné, režim tedy není optimistický.

**Realistický režim backtestu** (`--meteo=predpoved`): model se učí ze skutečných dat, od D 09:00 do konce D+1 dostane předpovězenou teplotu, osvit a vítr. Výsledek (v8): RMSE 14,30 → **19,42**, MAPE 1,66 → **2,17 %**; s korekcí 12,93 → **18,01**, MAPE 1,48 → **2,01 %** — realistická přesnost v provozu MAPE ~2,0 %. Nejvíc trpí jaro (11,2 → 17,9) a zima (14,2 → 19,6); do chyby se promítá hlavně topení (rozdíl složek RMSE 7,8), pak chlazení 6,8 a FVE 6,1.

**Předpověď ČEPS v realistickém režimu** (`explore/ceps_analyza.py --bt=simulace/predpoved_meteo --meteo=predpoved`; model ČR pro „překvapení ČEPS“ také s předpovězeným počasím, jinak by unikla informace o skutečném počasí): s korekcí 18,01 → 17,85 jen díky Vánocům, mimo Vánoce 17,69 → 17,69, MAPE beze změny → **nezapojovat** ani realisticky.

**Slunce a klimatizace** (hypotéza uživatele: v létě slunce zvyšuje chlazení víc, než v zimě snižuje topení): model má pro chlazení vlastní parametr a_c, odhad je 0 (topení a = 0,021). Večerní test (19–23 h, FVE nevyrábí, teplé pracovní dny, kontrola teploty): rozdíl po slunečném a zataženém dni +0,3 ± 2,3 — nulový. Profil věrohodnosti (a_c napevno, fit z hodinových dat): letní RMSE 12,96 (a_c 0) → 13,16 / 13,35 / 13,78 / 14,32 (0,005 / 0,01 / 0,02 / 0,03), odhad špičky FVE 64 → 77 / 93 / 133 / 148 (vynucené chlazení se kompenzuje nerealistickou FVE), s rostoucím a_c klesá i topné a (0,021 → 0,005) — osvit je sdílený regresor. **Vliv slunce na chlazení nad rámec teploty v portfoliu neměřitelný**; případný malý efekt nelze bez nezávislé kapacity FVE oddělit.

### 2026-10-02 — rozklad zbývající chyby, oprava korekce (krok 24)

**Rozklad chyby backtestu v8** (skutečné počasí, s korekcí, RMSE 12,93): úroveň dne (denní průměr chyby) 8,6 = 45 % rozptylu, tvar uvnitř dne pomalejší než 2 h 8,5 = 43 %, čtvrthodinový šum 3,6 = 8 % (šum je i v samotné spotřebě: 3,4–3,9 po odečtení průměrného profilu → tvrdé dno MAPE ~0,4 %). Autokorelace denní chyby po korekci: 1 den +0,36 (v okamžiku vydání známá jen do 09:00), 2 dny +0,03 — z vlastní historie chyb už úroveň zlepšit nejde. Chyba není soustředěná: nejhorších 10 % dní = 35 % součtu čtverců, medián denního RMSE 10,4. Po hodinách: noc 8–9, den 14–17 (maximum ve 14–15 h); sobota 13,5 a neděle 14,4 proti pátku 10,6. Bias přechodných měsíců: říjen–listopad +4 až +4,7, březen–duben +2,6 až +3,5. **Předpověď počasí tvoří zhruba polovinu rozptylu provozní chyby** (12,9 → 18,0; ve 12–14 h RMSE 15–16 → 26; srpen MAPE 1,75 → 2,77 %).

Výhrada: kroky 8–23 se ladily a hodnotily na stejném roce backtestu (8/2025 – 8/2026) → čísla mohou být mírně optimistická; čistý test dají data po 4. 8. 2026.

**Oprava korekce** (`src/correction.py`): `irregular_days` porovnávalo `DatetimeIndex` s množinou `datetime.date`, takže svátky a mosty mimo Vánoce se od zavedení kalendáře nevyřazovaly (vyřazovalo se jen 21. 12. – 3. 1.). Po opravě (10 dní v backtestu): s korekcí 12,93 → **12,84**, MAPE 1,482 → **1,473 %**; s předpovědí počasí 18,01 → 17,94.

### 2026-10-02 — opožděné zatopení po teplém období (krok 25)

**Pozorování z praxe** (uživatel): po teplém víkendu nebo týdnu se v administrativních budovách po ochlazení netopí hned — zatopí se během dne nebo až další den. Test na reziduích plného fitu (`explore/zatop_test.py`, mimo Vánoce):
- **První chladný pracovní den po teplém období** (průměr ≤ 12 °C po 2–5 dnech ≥ 15 °C; 7–10 událostí za 4,5 roku): reziduum v pracovní době −2 ± 3 (podle definice −7 až 0), v noci −3 až −9, **den poté −4 až −8**. Směr odpovídá pozorování, velikost je do 1 % zátěže a na hranici šumu.
- **Stav ústředního topení podle pravidla vyhlášky** (zahájení při denním průměru < 13 °C dva dny po sobě, přerušení při > 13 °C dva dny po sobě): chladné dny ve stavu „neběží“ (56 dní, hlavně září, začátek října, květen) mají reziduum −2,0 ± 1,4; pracovní dny celkově 0, volné dny −6,4 ± 2,6, podzimní pracovní doba −5,2 ± 2,5 (n = 14).
- **Závěr: samostatný člen nezavádět.** Většinu zpoždění už nese pomalý kanál teplotního filtru (váha 64 %, τ ≈ 82 h); zbytek je ~5 jednotek ve ~12 dnech ročně — na ročním RMSE pod 0,01 — a z tak mála událostí by se člen odhadoval hůř, než kolik by přinesl. Část chytá korekce z ranní chyby dne vydání (x2). Vrátit se, až bude událostí víc, nebo pokud by se dal použít skutečný údaj o zahájení topné sezóny.
- Vedlejší nález: 24. 9. 2025 (nejhorší den backtestu, −32) není zatopení — spotřeba je o 35–40 níž rovnoměrně od 1 h do 16 h včetně noci (jednorázová událost v portfoliu nebo v datech); následující den naopak +25 až +40 v pracovní době.

### 2026-10-02 — šero přes den: spotřeba při tmavé obloze (krok 26)

**Diagnostika polední chyby** (rezidua plného fitu, 9–15 h, mimo Vánoce): reziduum závisí na osvitu nelineárně — při osvitu pod ~10 W/m² za dne +7 až +14, mezi 100 a 400 W/m² záporné (zima −10 až −12, přechod −3, léto −1), nad 400 opět kladné. V zimě koreluje polední reziduum s jasností oblohy −0,55, v létě 0. Rozdíl „tmavá minus světlá čtvrthodina“ je +6 až +12 v 8–17 h, **v pracovní i volné dny a ve všech sezónách**; největší kolem poledne.

**Zamítnutá hypotéza — sklon panelů FVE**: model násobí kapacitu osvitem na vodorovnou plochu, skloněné panely při nízkém zimním slunci dostávají víc. Přepočet osvitu do roviny panelů (rozklad na přímou a difuzní složku podle Erbse, sklon 20/35/50° k jihu) in-sample nepomohl (hodinový fit 12,62 → 12,60 až 12,66, zimní korelace s jasností −0,50 až −0,53); teplota článku místo teploty vzduchu 12,62 → 12,59. Do modelu nezavedeno.

**Člen „šero přes den“ v bázi** (`src/sun.py`, `src/base.py`, `src/etl.py`): `šero(t) = (1 − tma) · exp(−osvit / I0)`, `I0 = 40 W/m²`; příspěvek = šero × aktivita(čas dne), aktivita z kubických B-splin v 5–21 h (6 parametrů), ≥ 0, společná pro všechny typy dnů. Nejpravděpodobnější výklad je svícení přes den při těžké oblačnosti a mlze (přímo neověřeno — může v tom být i část solárních zisků, které lineární člen v topení nevystihne). Na rozdíl od osvětlení za tmy závisí na počasí; ETL ho počítá jako sloupec `sero`.
- Hledání I0 (hodinový fit celého řetězce): bez členu 12,62; I0 = 15 / 40 / 80 / 150 → 12,20 / 12,22 / 12,27 / 12,36. Zvoleno 40 (plošší, méně citlivé na chybu předpovědi).
- Při plném šeru je spotřeba vyšší o ~4 v 6 h, ~19 v 10 h, **~30–33 ve 12–14 h**, ~21 v 16 h, ~6 ve 20 h. Zimní poledne RMSE 15,8 → 13,7, korelace rezidua s jasností −0,55 → −0,32.
- **Solární koeficient topení a klesl z 0,021 na 0,0135 °C/(W/m²)** — část efektu tmavé oblohy dřív neslo topení jako „chybějící solární zisky“.

| backtest | RMSE | MAPE | + korekce | + korekce MAPE |
|---|---|---|---|---|
| v8 (po opravě korekce) | 14,30 | 1,66 % | 12,84 | 1,473 % |
| **v9 šero přes den** | **13,85** | **1,63 %** | **12,42** | **1,438 %** |

**Předpověď počasí a nelineární člen** (`src/meteo_forecast.py`, `src/backtest.py`): šero spočítané ze *středního* předpovězeného osvitu realistický backtest nezlepšilo (s korekcí 17,94 → 17,99, model sám 19,42 → 19,93) — předpověď z oblačnosti dává střední osvit a šero středního osvitu není střední šero, tmavé dny se tak nikdy nepředpoví. Proto se pro každý bin oblačnosti kalibrují i **kvantily** poměru jasnosti a šero se v predikci počítá jako průměr přes ně (střední hodnota nelineárního členu přes rozdělení chyby předpovědi): s korekcí 17,94 → **17,77**, MAPE 2,006 → **1,992 %**. Volba `--sero=stredni-osvit` vrací původní výpočet.

Očištěná spotřeba: šero je součástí báze, při přepočtu na normálové počasí zůstává skutečné (rozdíl proti normálu je v ročním průměru malý).

Plný fit na 15min datech (uložené parametry, `models/*_joint.npz`): RMSE 13,08 → **12,71**; T_b 16,8 °C, s 2,0, τ 82 h beze změny.

### 2026-10-02 — kalibrace předpovědi počasí: osvit podle výšky slunce, bias teploty a větru (krok 27)

**Chyby předpovědi z archivu proti stanici** (období před backtestem i v něm):
- **Osvit z oblačnosti měl špatný denní i roční chod**: poměr jasnosti kt závisel jen na oblačnosti, ale při stejné oblačnosti silně roste s výškou slunce (zataženo: 0,01 u obzoru → 0,43 v létě v poledne; jasno 0,05 → 0,78). Léto: ráno +50 až +70 W/m², v poledne −70 až −100; zima: v průměru dvojnásobek skutečnosti (+14 při průměru 12–33 W/m²).
- **Teplota**: předpověď je teplejší než stanice, hlavně v létě odpoledne a v noci (+1 až +1,7 °C, každý rok), přes den málo; denní chyba je setrvačná (autokorelace 0,59 / 0,41 / 0,36 pro 1 / 2 / 3 dny).
- **Vítr**: předpověď o ~1,7 m/s vyšší než stanice — model to čte jako chladnější pocitovou teplotu (topení v lednu +5,9, rovnoměrně přes den).

**Úpravy** (`src/meteo_forecast.py`, `src/backtest.py`; vše jen z dat známých při vydání):
- kt v buňkách **oblačnost × sin(výšky slunce)** (10 × 7, bilineární interpolace, řídké buňky doplněné ze sousedních oblačností), stejně i kvantily pro šero. Bias osvitu: léto v poledne −98 → −11 W/m², zima +14 → +0,5; korelace předpovězeného a skutečného šera 0,84 → 0,90.
- **Klouzavé odstranění biasu** teploty a větru (`meteo_forecast.debias`): průměr (předpověď − skutečnost) po hodinách dne za posledních 30 dní před cutoffem. Chyba předpovědi teploty RMSE 1,57 → 1,37 °C (okna 7 / 14 / 30 / 60 dní: 1,44 / 1,39 / 1,37 / 1,40), bias +0,55 → 0. Volba `--bias=temp,wind` (výchozí), `--bias=` vypne.

| realistický backtest (předpověď počasí) | RMSE | MAPE | + korekce | + korekce MAPE |
|---|---|---|---|---|
| v9 (šero, kvantily) | 19,43 | 2,18 % | 17,77 | 1,992 % |
| + osvit podle výšky slunce, bias teploty | 18,32 | 2,05 % | 17,31 | 1,911 % |
| **+ bias větru** | **18,27** | **2,04 %** | **17,25** | **1,909 %** |

Po měsících (model): červenec 22,0 → 18,2, srpen 21,2 → 16,1, květen 17,6 → 14,2; bias červenec −10,7 → −3,5, srpen −12,0 → −2,7. Samotný osvit s teplotou zhoršil leden a únor (bias −3,6 → −8,6): nadhodnocený zimní osvit dřív náhodou kompenzoval nadhodnocený vítr; s opravou větru je zimní bias stejný jako se skutečným počasím (leden −3,8 vs −3,7). Červen se zhoršil (14,3 → 17,3) — bias +8 odpovídá běhu se skutečným počasím (+5,6), dřív ho kryla chyba předpovědi.

Rozdíl proti běhu se skutečným počasím (12,42) zůstává 4,8 RMSE: po složkách topení 7,1, chlazení 4,8, FVE 4,7, báze (šero) 3,2 (v zimě v poledne šero z předpovědi v průměru o 3–4 jednotky nižší než skutečné) — převážně náhodná chyba předpovědi, kterou kalibrace neodstraní; další zlepšení už jen lepším vstupem (osvit přímo z numerického modelu, víc zdrojů předpovědi).

### 2026-10-02 — vlastní denní tvar topení a chlazení pro volné dny (krok 28)

**Diagnostika** (rezidua plného fitu po hodinách, typu dne a režimu): denní tvar topné citlivosti k(čas dne) i chladicí k_c byl společný pro všechny typy dnů (typ dne jen posouval úroveň). O víkendech a svátcích se ale ranní náběh posouvá: v topných dnech (topná složka > 100) reziduum volných dnů **−12 až −16 v 7–9 h a +9 v 11–12 h**, pracovní dny naopak +8 v 7–8 h (společný tvar byl kompromis); u chlazení zrcadlově — volné dny +16 až +17 v 6–8 h a −7 až −9 v 11–12 h, pracovní −12 v 7 h. Sobota a neděle se chovají stejně. Proto byly víkendy v backtestu horší (sobota 13,5, neděle 14,4 proti pátku 10,6).

**Úprava** (`src/heating.py`, `src/cooling.py`): k i k_c dostaly Fourierovu odchylku pro volné dny (typ dne 2 a 3; 4 harmonické, +8 parametrů v každém členu, lineární krok fitu).

| backtest | RMSE | MAPE | + korekce | + korekce MAPE |
|---|---|---|---|---|
| v9 šero přes den | 13,85 | 1,63 % | 12,42 | 1,438 % |
| **v10 tvar topení/chlazení pro volné dny** | **13,41** | **1,58 %** | **11,94** | **1,381 %** |
| realisticky (předpověď počasí) v9 + krok 27 | 18,27 | 2,04 % | 17,25 | 1,909 % |
| **realisticky v10** | **17,95** | **2,00 %** | **16,93** | **1,868 %** |

Hodinový fit celého řetězce 12,22 → 11,77. Leden 12,8 → 11,6 (s korekcí), červenec 13,8 → 12,8, listopad 12,7 → 11,9.

**Korekční vrstva — tvar chyby podle třídy dne** (zkoušeno, nezavedeno): tvar chyby x3 počítaný zvlášť z posledních pracovních a zvlášť z posledních volných dnů dával na v9 12,42 → 12,20; po opravě v modelu už nic (11,94 beze změny vs 12,03, realisticky 16,93 vs 17,00). Chyba byla strukturální, ne drift — korekce zůstává beze změny. Nepomohl ani pomalý bias (průměr chyby posledních 7 dní) a chyba stejného dne minulého týdne.

Plný fit na 15min datech: RMSE 12,71 → **12,26**.
