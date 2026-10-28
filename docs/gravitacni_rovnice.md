# Gravitační rovnice: současná formulace teorie

Sepsáno 27. 10. 2026 po replikaci F7.1c. Je to analytický souhrn, ne výsledek nového testu. Text shrnuje, co rovnice
říká, co z ní je odvozené, co laboratoř testovala a co je zatím vstup. Rozhodnutí, na kterých stojí: 1 (rovnice pole,
čtvrté kolo), 6 (zdroj) a 12 (gravitace mění kauzální řád), viz `docs/rozhodnuti_autora.md`.

Konvence: signatura (−, +, +, +), c = 1, ve 3+1. Π je rychlost hodin (dτ = Π dt podél statické světočáry) a L = ln Π.
U hmoty je L < 0. M značí GM/c² (délka). Δ a ∇ jsou ploché souřadnicové operátory izotropních souřadnic.

## 1. Rovnice

**Diskrétně (kauzální množina):**

    B_g ln Π = κ ρ_s

B_g je Benincasův–Dowkerův operátor na uspořádání s geometrií g. ρ_s je hustota nasycených událostí (rozhodnutí 6:
S_v = 1, když ΔQ_v ≥ ΔQ_crit).

**Metrika (H1), kterou rovnice předpokládá:**

    ds² = −Π² dt² + Π⁻² dx² = −e^{2L} dt² + e^{−2L} dx²

**Spojitá limita:** B_g → □_g − R/2. V metrice H1 platí přesně:

    R = 2 e^{2L} (ΔL − |∇L|²)
    (□_g − R/2) L = e^{2L} [ (1 − L) ΔL + L |∇L|² ] = κ ρ_s

**Slabé pole:** ΔL = κρ_s, tedy L = −U (newtonovský potenciál). Kalibrace je κρ_s ↔ 4πGρ.

**Vakuum:** (1 − L) ΔL + L |∇L|² = 0.

## 2. Odkud je metrika: princip zachování kauzální propustnosti

Podrobně v `proc_gamma1.md` (oddíly 4 a 7). Ve zkratce:

- Metrika se z uspořádání čte třemi měřidly: hodiny (počet událostí v diamantu) → φ; hranice relace (rychlost
  kuželu) → φ + ψ; souřadnicová hustota → φ − d_s ψ.
- Princip: hustota událostí krát rychlost kuželu je konstantní (ρ_c · c = konst; kniha 5.3, 23, 26.1). Ve slabém
  poli odtud plyne γ = 2/(d_s − 1), ve 3+1 tedy γ = 1.
- S hodinami Π je H1 dané přesně: A·B = 1, c_eff = Π², ρ_c = Π⁻².
- A6 (S = Σ Π Δτ) dává efektivní metriku Π² g_u. H1 nastane právě tehdy, když má g_u rychlost kuželu Π².
- Schwarzschild splňuje A·B = 1 jen do prvního řádu: A·B = (1 − M²/4r²)².
- Dvě alternativy dávají jiné γ: pevný řád dává γ = −1 (světlo se neohne) a samotné zúžení kuželů γ = 1/d_s.

## 3. Kontroly, které rovnice prošla (analyticky, 27. 10.)

Všechny body byly přepočítány v kontejneru (sympy a scipy).

**(a) Křivost a tvar rovnice.** R i levá strana (□_g − R/2)L byly spočteny z metriky (Christoffelovy symboly →
Ricciho tenzor) pro obecné L(r). Rozdíl proti vzorcům z oddílu 1 je identicky 0.

**(b) První integrál a uzavřené vakuové řešení (nové).**

- Se sférickou symetrií a w = r²L′ zní vakuová rovnice (1 − L) w′ + L L′ w = 0, tedy d ln w / dL = 1 − 1/(1 − L).
- Odtud plyne přesný první integrál:

      r² L′ = M e^{L} (1 − L)

  M je asymptotický náboj (r²L′ → M pro r → ∞).
- Druhá integrace vede na exponenciální integrál:

      Ei(1 − L) = Ei(1) + e M / r

- Numerická integrace ODR od r = 200 M se s uzavřeným tvarem shoduje na 5·10⁻⁸.
- Rozvoj: L = −M/r + M³/(6r³) − M⁴/(12r⁴) + … Koeficient u M²/r² je nulový, to je a = 0 z dřívějška.

**(c) PPN.**

- g_tt = −e^{2L} = −(1 − 2M/r + 2M²/r² − …), tedy β = 1.
- g_rr = e^{−2L} = 1 + 2M/r + …, tedy γ = 1.
- γ teď plyne z principu propustnosti (oddíl 2), ne jen z ansatzu.
- β = 1 + a, kde a je koeficient u M²/r² v L. Vakuová rovnice vynutí a = 0.
- β ale není nezávislá předpověď, protože f = ln Π bylo vybráno podle perihelia (rozhodnutí 1, čtvrté kolo).
  Z rodiny (Π^k − 1)/k dává β = 1 jen ln Π.

**(d) Druhý post-newtonovský řád (předpověď, zatím netestovaná).**

| izotropně, x = M/r | naše rovnice | obecná relativita |
|---|---|---|
| −g_tt | 1 − 2x + 2x² − **1**·x³ | 1 − 2x + 2x² − **3/2**·x³ |
| g_rr | 1 + 2x + **2**·x² + **1**·x³ | 1 + 2x + **3/2**·x² + **1/2**·x³ |

První odchylka je v g_rr ve druhém řádu. Je to přesně neplatnost A·B = 1 pro Schwarzschilda. Měřitelná by byla ve
druhém řádu ohybu světla nebo Shapirova zpoždění; tak přesné měření zatím neexistuje.

## 4. Silné pole (vakuové řešení bodové hmoty)

> **Podmíněné (doplněno 27. 10., `tenzorova_vazba.md`):** celý tento oddíl platí, jen pokud je skalární rovnice
> nelineárně přesná (rozcestí B). Pokud je princip propustnosti zákonem lineárního řádu a nelineární doplnění je
> obecná relativita (rozcestí A, doporučení laboratoře), hrdlo, chybějící horizont, +0,91 % stínu i 2PN koeficienty
> z bodu 3d nejsou předpověďmi teorie.

Hodnoty plynou z uzavřeného řešení v bodě 3b.

| veličina | naše rovnice | obecná relativita |
|---|---|---|
| horizont | žádný: g_tt = e^{2L} > 0 pro každé r > 0 | r_areal = 2M |
| fotonová sféra | rL′ = 1/2: r_iso = 1,798 M, L = −0,535, r_areal = 2M(1 − L) = **3,070 M** | 3 M |
| kritický parametr dopadu | b_crit = min r e^{−2L} = **5,2435 M** | 3√3 M = 5,1962 M (+0,91 %) |
| hrdlo | rL′ = 1: r_iso = 0,612 M, L = −1,341, g_tt = 0,068, r_areal = M(1 − L) = **2,341 M** | – |

**Čtení.**

- **Stín.** Je o 0,91 % větší než v obecné relativitě. Na současnou přesnost snímků stínu (EHT, chyby řádu
  10 %) to nerozliší, ale je to konkrétní číslo.
- **Hrdlo.**
  - Pod 2,341 M plocha sfér znovu roste.
  - g_tt jde k 0 až pro r → 0, a to v nekonečné vlastní vzdálenosti (∫ e^{−L} dr diverguje).
  - Čistě vakuová bodová hmota je tedy hrdlo jako u exponenciální metriky, ne černá díra s horizontem.
  - Pro srovnání: exponenciální metrika L = −M/r má hrdlo v e·M = 2,718 M, fotonovou sféru 3,297 M a b_crit 5,437 M.
    Naše řešení leží mezi ní a obecnou relativitou.
- **Fyzikální relevance hrdla.**
  - Závisí na tom, zda může těleso s povrchem pod fotonovou sférou nebo pod hrdlem vůbec existovat staticky.
  - F1.7 N našel pro gaussovský zdroj mez statických řešení (ln Π: u_lin(0) mezi −0,3 a −0,5).
  - Převod na kompaktnost GM/(Rc²) ale udělaný není a zdroj nemá tlak. Srovnání s neutronovými hvězdami a se
    zhroucením proto zatím nejde.
- **Ohraničení.** Všechny body tohoto oddílu jsou spojitá limita. Diskrétní řešení v silném poli laboratoř
  nemá.

## 5. Strana zdroje (otevřená)

> **Doplněno 27. 10. (`tenzorova_vazba.md`):** v lineárním řádu je zdrojem T_ab toku, ne skalár ρ_s. Ve 3+1 je
> princip propustnosti totéž co čistá tenzorová vazba (h̄_ab = 4U u_a u_b). Skalár jako úplný dynamický zákon je
> vyloučený: kvadrupólové vyzařování kruhové dvojhvězdy je 1/6 hodnoty obecné relativity a strhávání chybí.
> Body níže popisují skalární formulaci.

- ρ_s je skalár, hustota nasycených událostí. Rovnice je skalární v pevném ansatzu metriky H1, není tenzorová.
- Světlo, tlak a pohyb zdroj nenesou, dokud se nedefinuje, jak přispívají k nasycení.
  - Ohyb světla hmotou (γ) to neovlivní.
  - Ovlivní to ale to, zda světlo samo gravituje.
  - Ovlivní to aktivní hmotu horkých a tlakových zdrojů. F1.7 N: nelineární pole přispívá k vnější hmotě
    (1,02–1,51 × lineární hodnota).
- K.15 má měřidlo neseného toku T^ab. Vazba T^ab → ρ_s ale definovaná není.
- κ ↔ G je kalibrace (κρ_s ↔ 4πGρ), ne odvození.
- **Dynamika a vlny chybí.** Rovnice je zatím statická.
  - Skalární teorie gravitace obvykle vyzařují monopólově nebo dipólově. Pozorování dvojitých pulsarů (úbytek
    periody souhlasí s kvadrupólovým vyzařováním obecné relativity) jsou proto pro rovnici hlavní budoucí test.
  - Bez dynamické formulace se k nim nedá nic říct.

## 6. Co je odvozené, co testované a co je vstup

| | položka | stav |
|---|---|---|
| **vstup** | princip zachování propustnosti ρ_c · c = konst | kniha 26.1 jako podmínka Einsteinova režimu; proč ji substrát drží, je otevřené |
| **vstup** | hodiny Π, akce S = Σ Π Δτ (A6) | axiom |
| **vstup** | BD operátor (□ − R/2 v kontinuu) | převzato z teorie kauzálních množin; K.11, K.13a |
| **vstup** | f = ln Π | vybráno podle perihelia, β = 1 je proto vstup |
| **vstup** | ρ_s = ρ★ S_v (nasycené události), κ ↔ G | rozhodnutí 6; kalibrace |
| **odvozené** | H1 z principu a hodin; γ = 2/(d_s − 1), ve 3+1 γ = 1 | `proc_gamma1.md` |
| **odvozené** | R a spojitý tvar rovnice v H1 | sympy, identita (3a) |
| **odvozené** | první integrál r²L′ = M e^L(1 − L), řešení Ei(1 − L) = Ei(1) + eM/r | 3b, numericky ověřeno |
| **odvozené** | a = 0 (β = 1 při dané volbě ln Π), 2PN koeficienty | 3c, 3d |
| **odvozené** | žádný horizont, fotonová sféra 3,070 M, b_crit 5,2435 M, hrdlo 2,341 M | oddíl 4 |
| **testované** | měřidlo γ z uspořádání: 2+1 (K.16b) a 3+1 (K.16c) | PROŠLO |
| **testované** | samotné zúžení kuželů dává γ = 1/d_s | F7.1a 2+1: 0,47 / 0,43 (SELHAL podle předpovědi); F7.1c 3+1: 0,349 / 0,328 a 0,318 / 0,310 (předpověď 1/3) |
| **testované** | propustnost dává γ = 2 (2+1) | F7.1b 2×NEROZHODNUTO: 1,81 / 1,94 a 1,77 / 1,96; systematický posun A dolů a B nahoru |
| **testované** | propustnost dává γ = 1 (3+1) | F7.1c 2×NEROZHODNUTO: 1,000 / 0,972 a 0,948 / 0,979; cesta A projde v obou bězích, cesta B neprojde kvůli šumu srovnání s kontrolou (SE ~0,1) |
| **testované** | BD operátor v geometrii H1: zrušení lineárního členu f = Π | K.13a |
| **testované** | diskrétní odezva mimo zdroj c ≈ 0 (β = 1) | F1.9a MĚŘIDLO SELHALO (kvadratura σ = 4); popisně c souhlasí se spojitou hodnotou a odchylka klesá jako (ℓ_k/σ)² |
| **testované** | mez statických řešení silného pole | F1.7 N (spojité, gaussovský zdroj) |

Shrnutí: γ = 1 se z ansatzu posunulo na princip, který laboratoř na surových hodnotách potvrzuje ve 2+1 i 3+1.
Formální verdikt chybí jen kvůli statistické slabosti srovnání. β = 1 zůstává vstupem. Všechno v silném poli je
odvozená předpověď spojité limity, bez diskrétního testu.

## 7. Co by teorii posunulo nejvíc

V pořadí podle toho, co nejvíc rozhoduje o pravdě:

1. **Dynamika.** Časově závislá forma rovnice a otázka vyzařování (monopól a dipól proti dvojitým pulsarům). Tady
   může skalární formulace padnout. Je to nejdůležitější nevyřešený bod.
2. **Strana zdroje.** Jak světlo, tlak a pohyb přispívají k ρ_s (přes T^ab z K.15). Bez toho nejde srovnat
   kompaktní tělesa ani kosmologii.
3. **Proč substrát drží propustnost.** Odvození principu z axiomů místo převzetí z knihy 26.1.
4. **Formální verdikt F7.1c.** Jen pokud o něj autor stojí: nový test s předem danou statistikou (víc sprinklingů
   pro kontrolu i kandidáta). Kritéria F7.1c se nemění (R132).
