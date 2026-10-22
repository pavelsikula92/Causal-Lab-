# Roadmapa Kauzální fyziky

Stav 18. 10. 2026. Úplný návrh s tabulkami testů je v `Causality_Lab_navrh.pdf`; aktuální stav všech testů
generuje `python run.py report` do `docs/seznam_testu.md`.

## Cíl

Ověřit, zda gravitace a setrvačnost jsou dva projevy jedné dynamiky kauzálního odporu
(*Physics of Causality and GTS*, únor 2026; matematická báze *Kauzální fyzika, Díl I*):
R = Σ (I_v/Π_v) Δτ = Σ Z_v Δτ, hmota jako trvalá struktura vysoké impedance, setrvačnost i gravitace jako
dráhy nejmenšího kumulativního odporu. Laboratoř teorii testuje, nezachraňuje; do výsledku nesmí vložit
geometrii ani hledaný zákon.

Řetězec: **kauzální struktura → transport a rezistence → efektivní geometrie → gravitační a setrvačné chování →
nezávislé pozorovatelné predikce.**

## Tři úrovně závěru

| Úroveň | Smí tvrdit | Sama nedokazuje |
| --- | --- | --- |
| mechanismus | mechanismus funguje v deklarovaném modelu a rozsahu | že mechanismy tvoří uzavřenou teorii |
| řetězec | navazující mechanismy společně reprodukují předpovězené chování | že řetězec odpovídá přírodě |
| predikce | předem předpovězený nezávislý výsledek uspěl podle kritéria | ostatní důsledky teorie |

Měřidla (fáze K) mají vlastní kategorii **kalibrace** bez fyzikálního závěru. Fáze F8 je predikcí jen tehdy,
když nic nebylo kalibrováno na testovacích datech.

## Fáze a brány

| Fáze | Část teorie | Brána |
| --- | --- | --- |
| K Měřidla a kalibrace | průřezová | každé měřidlo použité ve fázi má záznam kalibrace v rozsahu, ve kterém se použije |
| F0 Substrát | axiomy 1, 2, 8 | měřidla substrátu ověřena v deklarovaném rozsahu; nelokálnost otevřená, dokud ji nerozliší F0.4 |
| F1 Transport a pole | axiom 3, 4, 5 a A5′ | tři oddělené úrovně: integrální A5′ proti operátoru (F1.2), asymptotika BD-3 (F1.3, formálně), numerika r ≫ ℓ_k jen kde je dosažitelná; superpozice a škálování samostatně |
| F2 Pohyb v poli | axiom 6 | geodetické chování ověřeno ve 2+1; predikční hodnota akce jen ve vážené reprezentaci 2+1; otevřeno 3+1 a uzávěr zdroj → pole → pohyb (F2.7) |
| F3 Vznik struktury | axiom 4, nelineární člen A5′ | struktura vzniká v předem daném širokém rozsahu, zpětná vazba zapnuto proti vypnuto; I, Π, ω měřeny nezávisle; jinak hmota zůstává postulátem |
| F4 Setrvačnost | R = Σ Z Δτ, m = ω | každý cíl (odezva, m ∝ ω, dilatace hodin, E = mc²) s vlastní definicí, tolerancí a nulovým modelem |
| F5 Gravitace a ekvivalence | ústřední tvrzení | uzavřený gravitační řetězec na syntetických datech (5.1–5.4), Newtonova limita bez vloženého zdroje; vložená struktura = výsledek podmíněný |
| F6 Dynamika řádu | axiom 7, D_eff ≈ 3 | vlastní brána; selhání nezneplatňuje F1–F5 na pevném uspořádání; úspěch ve 2D není důkaz ve 3+1 |
| F7 Světlo a vlny | H1, H2 | do odvození a numerického ověření zůstávají hypotézami |
| F8 Konfrontace s daty | problémy 6, 7, 10, sluneční soustava | odvozené, kalibrované a předpovězené odděleno; data a kritéria dána před vyhodnocením |

## Závislosti

```
K ─▶ F0 ─▶ F1 ─▶ F2 ─┐
       └─▶ F3 ─▶ F4 ─┴─▶ F5 ─▶ F7 ─▶ F8
           F3 ──────────▶ F5 ─────────▶ F8
F6 (souběžná větev) ┄┄▶ F7 (částečné vstupy), F8.2 a F8.5
```

F8.1 (sluneční soustava) potřebuje F5 a F7, F6 není podmínkou.

## Pravidla

- Kritérium je v hlavičce testu a po běhu se nemění. Změna kritéria = nová revize testu (F1.2b) s důvodem v hlavičce;
  starý výsledek zůstává a platí podle původního kritéria.
- Verdikt platí až po replikaci jiným seedem (`--seed K`).
- Úroveň závěru je povinná; výstup nesmí tvrdit víc, než úroveň dovoluje.
- Měřidlo se používá jen v rozsahu své kalibrace.
- Nová hypotéza může vzniknout po neúspěchu, ale nesmí zpětně změnit kritéria, podle kterých neuspěl starý test.
- Negativní výsledky se zachovávají.

## Nejbližší rozhodující testy

F1.2 (integrální A5′) uzavírá pole, F3.3 (zpětná vazba) otevírá hmotu. F2.7 je první test na úrovni „řetězec“.
