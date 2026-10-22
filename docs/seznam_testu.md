# Seznam testů Causality Lab

Vygenerováno příkazem `python run.py report` (2026-10-18 11:17). Neupravuj ručně: hlavičky testů jsou v `phases/`, plán v `docs/planovane_testy.csv`, historie v `results/historie_3x/prehled.csv`.

## Implementované testy

| ID | Fáze | Testovaná otázka | Úroveň závěru | Měřidla | Závislosti | Stav (poslední výsledek) |
| --- | --- | --- | --- | --- | --- | --- |
| K.2 (E47) | K Měřidla a kalibrace | Měří čelo kauzální budoucnosti zdroje zpoždění a ohyb, který OTR předpovídá pro metriku H1? | kalibrace | — | — | d3: -; d4: - |
| K.3 (E45c) | K Měřidla a kalibrace | Přečte délka nejdelšího řetězce malých diamantů profil Π, který byl do hustoty vložen? | kalibrace | — | — | d3: -; d4: - |
| K.4 (E45b) | K Měřidla a kalibrace | Počítá kvadratura očekávanou hodnotu rozmazaného BD operátoru správně, včetně lokální limity a konvergence? | kalibrace | — | — | quadrature: -; shape_d3: -; shape_d4: - |
| K.5 (E45b) | K Měřidla a kalibrace | Je diskrétní BD operátor jednoho sprinklingu v průměru roven očekávanému operátoru? | kalibrace | K.4 | — | d3_e0.2: -; d3_e0.05: -; d4_e0.2: - |
| K.6 (E45a) | K Měřidla a kalibrace | Zachovává routing tok přesně, s kapacitou i bez ní? | kalibrace | — | — | d3: -; d4: - |
| K.7 (E46) | K Měřidla a kalibrace | Mění čtení posunu každého funkcionálu znaménko přesně při zrcadlení bodové množiny (bez pole)? | kalibrace | — | — | d3: -; d4: - |
| K.8 (E55) | K Měřidla a kalibrace | Měří značená sonda propustnost T_out/T_in bez vlivu své velikosti a bez posunu pozadí? | kalibrace | K.6 | — | d3: - |
| K.9 (E55b) | K Měřidla a kalibrace | Dělají nulové modely přepojení to, co deklarují (platnost, stupně, promíchání)? | kalibrace | — | — | d3: - |
| F1.2 (E45d) | F1 Transport a pole | Blíží se agregovaný BD operátor jednoho uspořádání aplikovaný na pole zdroje očekávanému operátoru, s chybou klesající jako (L/ℓ)^(−p) a exponentem p z nezávislé kalibrace? | mechanismus | K.4, K.5 | F1.1 | d3: -; d4: - |

## Plánované testy

| ID | Fáze | Testovaná otázka | Úroveň závěru | Závislosti | Poznámka |
| --- | --- | --- | --- | --- | --- |
| K.1 | K Měřidla a kalibrace | Od jaké velikosti je geometrie čitelná z řádu ve 3+1? | kalibrace | — | výsledek z E54 platí (od N* = 32 000), měřidlo se převede, až ho bude potřebovat F0 nebo F6 |
| F0.4 | F0 Substrát | Je nelokálnost vlastnost substrátu, nebo měřidla? | mechanismus | K.4, K.5 | škálování ℓ_k a nelokálnosti s N, kalibrace na plochých strukturách |
| F0.5 | F0 Substrát | Operační definice tuhosti I | kalibrace | K.3 | hustota diamantů ve světotrubici, nezávislá na Π a na klasifikaci oblasti |
| F1.3 | F1 Transport a pole | Jaká je asymptotika pole (BD-3)? | mechanismus (formální) | K.4 | Π − Π_vac ∝ −1/r pro r ≫ ℓ_k, odlišit 1/r, 1/r², 1/r³, sklon −1,04 (3+1) a −0,12 (2+1) na [5, 20] ℓ_k už spočten |
| F1.4 | F1 Transport a pole | Sčítají se pole dvou zdrojů? | mechanismus | F1.2 | superpozice v toleranci |
| F1.5 | F1 Transport a pole | Škáluje pole se silou zdroje? | mechanismus | F1.2 | Π − 1 ∝ M |
| F2.5 | F2 Pohyb v poli | Přenos geodetického chování do 3+1 a závislost na N | mechanismus | K.1, K.7 | vysvětlit přebytek 10 % u Σ Π |
| F2.6 | F2 Pohyb v poli | Netriviální pole (dva zdroje, nesouměrná konfigurace) | mechanismus | F2.5 |  |
| F2.7 | F2 Pohyb v poli | Uzavírací test zdroj → pole → pohyb | řetězec | F1.2, F2.5 | Π vytvořené zdrojem, ne vložené |
| F2.8 | F2 Pohyb v poli | Newtonova limita pro zkušební tělesa | řetězec | F1.3, F2.7 | a ∝ M/r², nezávislé na tělese, přechod F1/F2 |
| F3.3 | F3 Vznik struktury | Vytváří zpětná vazba lokalizovanou trvalou strukturu? | mechanismus | F3.2, K.8, K.9 | Π_sat, exponent m a směr efektu určí autor před během, zapnuto proti vypnuto na stejných grafech (dříve E56) |
| F3.4 | F3 Vznik struktury | Vzniká struktura bez ladění? | mechanismus | F3.3 | sken předem daného širokého rozsahu |
| F3.5 | F3 Vznik struktury | Vlastnosti struktury | mechanismus | F0.5, K.8, F3.3 | velikost, stabilita, ω, I a Π měřené nezávisle, Z = I/Π až z nich |
| F4.1 | F4 Setrvačnost | Klade struktura odpor změně pohybu? | mechanismus | F3.3 nebo vložená struktura |  |
| F4.2 | F4 Setrvačnost | Je setrvačná hmotnost úměrná ω? | mechanismus | F4.1, F3.5 |  |
| F4.3 | F4 Setrvačnost | Zpomalují se vnitřní hodiny pohybující se struktury? | mechanismus | F4.1 | ω′ = ω √(1 − v²/c²) |
| F4.4 | F4 Setrvačnost | Platí E = mc²? | mechanismus | F4.2 | samostatné kritérium |
| F4.5 | F4 Setrvačnost | Je dráha nejmenšího odporu táž jako dráha maximálního řetězce? | řetězec | F2.5, F4.1 | R = Σ Z Δτ proti F2 |
| F5.1 | F5 Gravitace a ekvivalence | Je struktura zdrojem pole? | řetězec | F1.2, F3.3 | profil a síla podle A5′ se silou z ω |
| F5.2 | F5 Gravitace a ekvivalence | Přitahují se dvě struktury? | řetězec | F5.1, F4.1 |  |
| F5.3 | F5 Gravitace a ekvivalence | Shodují se setrvačná, pasivní a aktivní hmotnost? | řetězec | F4.2, F5.1 | struktury různé velikosti a složení |
| F5.4 | F5 Gravitace a ekvivalence | Newtonova limita ze struktur | řetězec | F5.1–F5.3 | bez vloženého zdroje pole |
| F6.4 | F6 Dynamika řádu | Má soubor uspořádání s globální akcí fázi podobnou varietě? | mechanismus | K.1, K.4 | Benincasa–Dowkerova akce |
| F6.5 | F6 Dynamika řádu | Vzniká geometrie, nebo je vynucená výběrem? | mechanismus | F6.4 |  |
| F6.6 | F6 Dynamika řádu | Plyne z D_eff(χ) MONDovská relace v⁴ = GMa₀? | mechanismus | F6.4 | bez ladění podle dat |
| F6.7 | F6 Dynamika řádu | Odvození A1, A2 a a₀ = cH₀/2π | mechanismus (formální) | — | motivováno, neodvozeno |
| F7.1 | F7 Světlo a vlny | Vyplývá zúžení kuželů c_eff = c·Π²? | mechanismus | F5.1, K.2 | parametr γ = 1, bez c_eff jen poloviční ohyb (R86) |
| F7.2 | F7 Světlo a vlny | Shapirovo zpoždění | mechanismus | F7.1 |  |
| F7.3 | F7 Světlo a vlny | Tenzorové vlny (H2) | mechanismus | F5.2, F6 |  |
| F8.1 | F8 Konfrontace s daty | Sluneční soustava: γ (Cassini), β (perihelium Merkuru), η (MICROSCOPE) | predikce | F5, F7 | F6 není podmínkou |
| F8.2 | F8 Konfrontace s daty | Rotační křivky, RAR, a₀(z) ∝ H(z) | predikce | F5, F6.6, F6.7 |  |
| F8.3 | F8 Konfrontace s daty | Bílí trpaslíci a práh χ_rig | predikce | F3, F5 | okno 38–40 až po odvození mapy χ_K ↔ χ_obs, NICER samostatně |
| F8.4 | F8 Konfrontace s daty | Kupy galaxií | predikce | F5 | včetně Bullet Cluster |
| F8.5 | F8 Konfrontace s daty | CMB | predikce | kosmologické pozadí a poruchy, případně F6 |  |
| F8.6 | F8 Konfrontace s daty | Termodynamika horizontu | predikce (formální) | dle odvození | entropie ∝ plocha, Hawkingova teplota |

## Výsledky GTS Lab 2.x a 3.x (historie)

Testy samotné jsou v archivu; přežívají jejich výsledky. Uložené záznamy 3.x patří do `results/historie_3x/`.

| Nové ID | Původní označení | Testovaná otázka | Výsledek | Úroveň závěru |
| --- | --- | --- | --- | --- |
| K.1 | E54 | Od jaké velikosti je geometrie čitelná z řádu ve 3+1? | kalibrováno od N* = 32 000, replikováno | kalibrace |
| K.2 | E47 | Funguje kauzální kužel jako měřidlo ohybu? | kalibrováno: 0,913 a 0,919 hodnoty OTR | kalibrace |
| K.3 | E45c | Lze Π zpětně přečíst z hustoty řetězců? | kalibrováno v testovaném rozsahu, replikováno | kalibrace |
| K.5 | E45b BD-1 | Je diskrétní BD operátor nestranný? | ověřeno: +0,004 ± 0,008 (klastrováno po sprinklingu) | kalibrace |
| K.6 | E45a (identity) | Zachovává routing tok a je okno ustálené? | ověřeno: reziduum ~10⁻¹⁵, doznívání rychlostí ~1 | kalibrace |
| K.7 | E46 (brána) | Jsou funkcionály řetězců zrcadlově souměrné? | ověřeno: přesně 0 | kalibrace |
| K.8 | E55 (sonda) | Měří sonda propustnost T_out/T_in? | ověřeno: bilance 0, citlivost 9·10⁻⁶ | kalibrace |
| K.9 | E55b (audit) | Jak věrné jsou nulové modely topologie? | DAG null 90 % nekrycích hran, užší chvost toku, přísný null promíchá 18 % hran – vliv krycích relací nerozhodnut | kalibrace |
| K.10 | E49 | Jsou verdikty stabilní při jiném seedu? | nerobustní: E40, E41-A, E42b | kalibrace |
| F0.1 | E43 | Nese řád signaturu geometrie? | potvrzeno: překryv minulostí 1,42–1,66 proti 1,00 | mechanismus |
| F0.2 | E44a | Je hloubka dána nejdelším řetězcem? | potvrzeno: nejkratší 2–5, nejdelší 7–66 | mechanismus |
| F0.3 | E42b | Lze tuhost I číst ze stupně události? | nedrží: stupeň roste s N, shlukování 0 | mechanismus |
| F1.1 | E45a | Platí Gaussův zákon pro čistý tok trvalého zdroje? | potvrzeno 2+1 i 3+1: J/S_w = 1,003–1,004, null 0,77–0,88 | mechanismus |
| F1.6 | E45b | Splňuje pole A5′ bod po bodu? | nelze: při dosažitelném N rozhoduje šum (INSTRUMENT_FAILS) | mechanismus |
| F1.7 | E44d | Jaký profil má pole kolem zdroje? | nerozhodnuto: měřena zátěž, ne Π | mechanismus |
| F1.8 | E40, E41-A | Omezuje kapacita tok, ne počet vazeb? | podpořeno, nerobustní vůči seedu | mechanismus |
| F2.1 | E44b | Ohýbají se maximální řetězce v poli Π jako geodetiky (1+1)? | potvrzeno | mechanismus |
| F2.2 | E46, E46b-A | Geodetické ohýbání v konformní reprezentaci 2+1 | potvrzeno: ohyb sleduje hustotu substrátu, každý maximální řetězec, replikováno | mechanismus |
| F2.3 | E46, E46b-A | Specifičnost axiomu 6 v konformní reprezentaci | NOT_SPECIFIC: bez prokázané specifické predikční hodnoty | mechanismus |
| F2.4 | E46b-B | Predikční hodnota akce ve vážené reprezentaci 2+1 | potvrzeno: Σ Π² → Π⁴η, 0,0770 proti 0,0776, replikováno | mechanismus |
| F2.9 | E46b-G | Sleduje geodetiku lokální hladový krok? | ne: 0,46 proti 0,043 | mechanismus |
| F3.1 | E55a | Šíří se saturace od přetíženého zdroje do okolí? | potvrzeno stínění: k/k_pred = 1,0054, vně 0 z 79 921 saturováno | mechanismus |
| F3.2 | E55b | Vznikne bez zdroje a bez zpětné vazby lokální struktura? | potvrzeno: Δ_pers +0,006 ± 0,006, E_loc +0,004 ± 0,004 (mez 0,05) | mechanismus |
| F5.0 | E44c | Je setrvačná a pasivní hmotnost totéž (vložené těleso, 1+1)? | potvrzeno podmíněně | mechanismus |
| F6.1 | E35–E41-B, E39, E42d, E43 | Vyroste z lokálního pravidla varieta? | nedrží: 0 ze 7 pravidel | mechanismus |
| F6.2 | E04b, E31 | Vznikne z dynamiky trojrozměrný prostor? | nedrží | mechanismus |
| F6.3 | E37d, E41-B | Roste šířka přítomnosti (2. zákon)? | nedrží | mechanismus |
