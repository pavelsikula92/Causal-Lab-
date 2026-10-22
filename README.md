# Causality Lab

Laboratoř Kauzální fyziky. Testuje, zda gravitace a setrvačnost jsou dva projevy jedné
dynamiky kauzálního odporu, na sprinklovaných kauzálních uspořádáních. Testy jsou seřazené podle fází
roadmapy a každá fáze patří k jedné části teorie (`docs/roadmapa.md`).

## Instalace

Python 3.10 nebo novější, knihovny numpy a scipy:

    pip install numpy scipy
    python run.py core          rychlé kontroly jádra (asi 10 s)

## Spuštění

    python run.py list                     fáze, testy a jejich poslední stav
    python run.py show F1.2                otázka, předpověď, kritérium, parametry testu
    python run.py F1.2                     test, všechny jeho části
    python run.py F1.2 --part d3           jen jedna část
    python run.py F1.2 --quick             zmenšený běh: ověří, že test běží (nic nerozhoduje)
    python run.py F1.2 --seed 7            nezávislá replikace jiným seedem
    python run.py F1.2 --workers 6         počet procesů (jinak CAUSALITY_WORKERS nebo všechna jádra)
    python run.py results F1.2             uložené výsledky testu
    python run.py report                   přegeneruje docs/seznam_testu.md


## Jak číst výsledek

Každý test končí jedním stavem: **PASSED** (prošel; u měřidla kalibrováno v deklarovaném rozsahu),
**FAILED** (neprošel), **UNDECIDED** (rozlišení nestačí, nebo kritérium splní i nulový model),
**INSTRUMENT_FAILED** (selhala kontrola měřidla, bez závěru), **RECORDED** (jen zaznamenaná čísla).

Úroveň závěru (kalibrace, mechanismus, řetězec, predikce) říká, co výsledek smí tvrdit. Výsledek se uloží
do `results/<ID>/` jako JSON a text, s otiskem testu i jádra, parametry, seedem a verzemi knihoven;
nic se nepřepisuje. Zmenšené běhy jdou do `results/<ID>/quick/` a nepočítají se.

## Struktura

    core/          jádro: řád, tok, operátory, řetězce, měřidla, statistika, otisky, rámec testů
    phases/        testy po fázích (K_meridla, F0_substrat, F1_pole, …); test = jeden soubor
    results/       výsledky; historie_3x/ pro výsledky GTS Lab 2.x a 3.x
    docs/          roadmapa, list nálezů, plánované testy, generovaný seznam testů
    core_tests/    rychlé kontroly jádra

Hlavička každého testu (`TEST = {...}`) obsahuje fázi, část teorie, otázku, metodu, předpověď, kritérium,
nulový model, úroveň závěru, použitá měřidla a závislosti. Kritérium se po běhu nemění; změna kritéria je
nová revize testu s novým ID (např. F1.2b).


Dokumentace česky, kód anglicky. Licence MIT.
