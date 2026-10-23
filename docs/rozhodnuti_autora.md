# Rozhodnutí autora k otevřeným otázkám (22. 10. 2026)

Záznam rozhodnutí a jejich přímých důsledků pro laboratoř. Rozhodnutí neměníme zpětně podle výsledků (R132).

## 1. Rovnice pole a ohyb světla ve 3+1

**Rozhodnutí:** γ = 1 zůstává fyzikálním cílem; rovnice pole se přepracuje. Současná rovnice (a) se nezachraňuje
dodatečným laděním a člen −R/2 se nepřidává ani neodebírá jen kvůli číslu.

- H1 se odděluje na dvě věci: **zákon šíření** c_eff = cΠ² (pracovní postulát pro kužely, F6.3b) a **efektivní
  gravitační metriku**, ze které plyne časové chování, prostorová křivost i pohyb světla. Že H1 plyne z konformní
  metriky Π²η, tvrdit nelze (konformní metrika má kužele ploché).
- Slabopolní cíl: ds² = −(1 + 2φ/c²)c²dt² + (1 − 2ψ/c²)dx², ψ = φ (γ = 1).
- **Proč (a) dává ve 3+1 při γ = 1 nulovou odezvu (odvozeno 22. 10., sympy, slabé pole Π = 1 + εu):** v metrice
  −Π²dt² + Π^(−2γ)dx² s d_s prostorovými rozměry (D = d_s + 1 celkem)
  - □_g Π → +1·∇²u (nezávisle na γ a d_s),
  - −(R/2) Π → (1 − γ(d_s − 1))·∇²u,
  - součet λ₁ = 2 − γ(d_s − 1).

  d_s = 3, γ = 1 ⇒ λ₁ = 0: −R/2 dává přesně −∇²u a ruší □_g, newtonovská odezva zmizí (přesné algebraické vyrušení
  lineárního členu, ne numerický artefakt). d_s = 2, γ = 1 ⇒ λ₁ = 1: příspěvek −R/2 je nulový, proto ve 2+1 (a)
  s H1 funguje. Člen −R/2 není libovolný: je to spojitá limita Benincasova–Dowkerova operátoru
  (□ − R/2, K.11), takže jeho odstranění by znamenalo opustit BD operátor.
- **Rozhodnutí (druhé kolo, 22. 10.): BD operátor zůstává; mění se, na co působí: B_g f(Π) = ρ_s.**
  V lineárním řádu □_g f(Π) → f′(1)∇²u a −(R/2) f(Π) → f(1)(1 − γ(d_s − 1))∇²u, ve 3+1 a γ = 1 tedy
  λ₁ = f′(1) − f(1): f = Π dává 0, f = Π^k dává k − 1, a každá funkce s f(1) = 0, f′(1) = 1 (ln Π, Π − 1,
  (Π^k − 1)/k) dává 1. Lineární odezva je tedy nerozliší.
  - primární kandidát f = ln Π (motiv: rychlosti hodin se skládají násobně, ln je převede na aditivní veličinu;
    motiv, ne odvození – musí se ukázat, že definice Π v teorii takové skládání má); kandidát, ne nový axiom;
  - povinná lineární kontrola f = Π − 1; nelineární kontroly (Π^k − 1)/k pro několik k;
  - rozdíl se hledá mimo lineární řád: B_g f(Π) = f′(Π) □_g Π + f″(Π)(∇Π)² − (R/2) f(Π),
    (∇Π)² = g^{μν}∂_μΠ ∂_νΠ; ln Π nese navíc −(∇Π)²/Π², Π − 1 nic.
  - Benchmark (statický sféricky symetrický zdroj ve 3+1, stejné ρ_s(r), okrajové podmínky a normalizace slabého
    pole pro všechny kandidáty): (1) slabé pole – stejná lineární odezva; (2) nelineární řešení – regulární profil
    vracející se k Π → 1; (3) druhý potenciál ψ vede k γ = 1; (4) přijatelné řešení při změně velikosti a síly
    zdroje bez dodatečného ladění f. ψ a γ musí vyplynout z celého modelu, ne z volby pro jednotlivý test.
  - Poznámka laboratoře k bodu (3): pro f s f(1) = 0 člen −R/2 v lineárním řádu vypadne úplně (R je O(ε), násobí
    ho f(1) = 0), takže ve slabém poli je rovnice □ f(Π) = ρ_s bez ohledu na geometrii; rozdíl BD proti □_g se
    projeví až nelineárně. Druhý potenciál: pokud je metrika zadaná ansatzem H1 (−Π²dt² + Π⁻²dx²), je ψ = φ
    a γ = 1 vstupem, ne výsledkem, a test (3) by byl tautologický. Netriviální test γ vyžaduje číst oba potenciály
    nezávisle – časový z hodin (Π, délky řetězců), prostorový z prostorových vzdáleností nebo z šířky kuželů
    v dynamickém uspořádání (F6.3b) – a teprve pak porovnat.

- **Rozhodnutí (třetí kolo, 22. 10.): K.13a zůstává hlavním diskrétním testem, K.13b v původní podobě vypadá,
  poznámka o ln Π je samostatné analytické zdůvodnění.**
  - Důvod: BD operátor je lineární v poli. Pro f(1) = 0 je E_g[f(Π)] = ε f′(1) B_0[u] + O(ε²) a rozdíl
    ln Π − (Π − 1) ve druhém řádu je −½ ε² B_0[u²] (plochý operátor; ∇²(u²) až v kontinuu). Odezva operátoru tedy
    f nerozliší. Rozliší ho jen princip a celá nelineární rovnice (F1.7 N).
  - K.13a měří operátor sám: E_g[B_g w] = B_0[w] + ε δB_u[w] + O(ε²) pro pevná pole w nezávislá na ε. Pole w = 1
    dává −R/2 = −ε[γ(d_s − 1) − 1]∇²u, gaussovské w dává opravu □_g a změna ℓ_k/σ konvergenci; systematika konečné
    škály se uvádí zvlášť. Povinné kontroly: Π = 1 (plochý operátor) a γ = −1 (K.11).
  - Princip pro f: skládání poměrů hodin (docs/derivations/ln_pi_potencial.md). Je to motivace, ne odvození.
    Jednoznačnost by vyžadovala nový axiom aditivity proměnné pole. Akce a zachování toku jsou dvě nezávislé
    kontrolní cesty.
  - Poznámka laboratoře: v H1 je první řád operátoru lineární v γ, δB = (1 − γd_s) P + (1 + γ) Q. P je hustotní
    část (ploché kužele, přístroj K.11), Q je část ze změny rychlosti kuželů. K.13a měří obě části zvlášť.

## 2. Tok v dynamickém uspořádání

**Rozhodnutí:** globální normalizace podle počtu budoucích potomků není základní fyzikální pravidlo; v F6.3b
zůstává jen jako referenční výpočet. Fyzikální kandidát: lokální rezervoárový model (zásoba Q_v, při vzniku vazby
se přidělí část dostupné zásoby, zbytek zůstává), bilance Σ Φ_in + s_v = Σ Φ_out + ΔQ_v v každém kroku.

Podmínky: (1) žádná znalost budoucích potomků, (2) bilance v každém kroku, (3) výsledky nezávislé na pořadí
výpočetních aktualizací.

**Napětí, které je třeba vyřešit dřív, než se pravidlo zvolí (poznámka laboratoře):** potomci události vznikají
v pořadí, které je lineárním rozšířením uspořádání, a různé rozšíření (např. jiná inerciální soustava) je seřadí
jinak. Pravidlo typu „každá nová vazba dostane podíl α zbývající zásoby“ dá potomkům podíly závislé na tomto pořadí,
tedy poruší podmínku 3. Pevný podíl z původní zásoby na vazbu je na pořadí nezávislý jen do vyčerpání zásoby; pak
rozhoduje, kteří potomci přišli dřív. Podmínky 1–3 dohromady tedy pravidlo silně omezují; nezávislost na pořadí
se musí měřit (bitová shoda v různých lineárních rozšířeních, jak to laboratoř dělá od F3.3b).

## 3. Pozadí toku

**Rozhodnutí:** univerzální tok pozadí se nezavádí. Základní test je při s_bg = 0. Režimy: s_bg = 0 (samostatná
kapka), s_bg > 0 (kapka udržovaná okolním tokem), sken s_bg (kritická hodnota nebo robustní oblast). Kapka, která
existuje jen při s_bg > 0, je „struktura udržovaná vnějším tokem“, ne samovázaný objekt. Kosmologický výklad
pozadí až po vyloučení numerických příčin a jen z dalšího pravidla. (F6.3b A má s_bg = 0 v mřížce.)

## 4. Kopírovaný pulz

**Rozhodnutí:** pulz není axiom, je to hypotéza laboratoře; fyzikální význam (frekvence, fáze) se mu zatím
nepřisuzuje. Rozlišují se: stav pulzu p_v (přenášená značka), tok Φ (bilancovaný), zásoba Q_v, odezva hodin Π_v.
Kopírování informačního stavu nesmí kopírovat konzervovanou veličinu.

**Důsledek pro dosavadní modely (poznámka laboratoře):** v F3.3c–F3.3h a F6.3/F6.3b pulz působí jako vnitřní zdroj,
s_eff = s + P_in, takže každá událost s pulzem přidá do toku ω. Bilance v kódu sedí jen proto, že je pulz započten
jako zdroj; fyzikálně jde o tok vytvořený kopírováním stavu, což toto rozhodnutí vylučuje. Samoudržující fáze
F3.3f po vypnutí zdroje i kapka F6.3b jsou poháněné právě tímto vstřikem. Výsledky těchto testů platí pro model
„pulz jako zdroj“; model podle rozhodnutí 4 (pulz mění odezvu, např. Π nebo kapacitu, ale tok nepřidává) je nový
a zatím netestovaný.

## 5. Parametry Π_sat, m, ΔQ_crit

**Rozhodnutí:** pracovní parametry, ne konstanty. F3.3e (sken ΔQ_crit, výběr předem daným skóre, ověření na nových
seedech) je hotový: optimum 0,2 na nových datech NEPROŠLO. F3.4: citlivost na Π_sat a m včetně interakce s prahem;
kombinace se nevybírá a neověřuje na stejných datech. Hledá se stabilní oblast, ne ostrý vrchol; varianta, která na
nových datech neobstojí, se odmítne nebo přeformuluje.

## 6. Saturovaná hmota jako zdroj gravitace

**Rozhodnutí:** nejjednodušší mapování S_v = 1, když ΔQ_v ≥ ΔQ_crit, jinak 0; ρ_s(v) = ρ★ S_v (ρ★ jediná
normalizace, neodvozená; první test sleduje tvar pole, ne absolutní sílu). Smyčka: zatížení → S → ρ_s → Π →
transport → nové zatížení, s lokální bilancí po každém kroku, určeným pořadím aktualizací a podmínkou konvergence.
Při neúspěchu nejdřív určit, zda selhává mapování, rovnice pole, nebo zpětná vazba.

Poznámka laboratoře: řešič pole (K.11c, F2.7, K.12) je statický a radiální; pro kapku v ose je radiální
přiblížení možné (zdroj zprůměrovaný po sprinklinzích a v čase), obecný nestatický zdroj vyžaduje nový řešič.
Platí rozhodnutí 1: dokud se nerozhodne rovnice pole ve 3+1, smyčka se testuje ve 2+1.

## 7. Růst uspořádání

**Rozhodnutí:** ano, ale ne v první verzi. F6.3b: dynamika relací (pevné body). Další fáze: dynamika zrození
událostí s lokálním pravidlem; nová událost nesmí vzniknout v minulosti existujících; předpovědi nesmí záviset na
libovolném pořadí narození (předloha: modely sekvenčního růstu kauzálních množin). Růst N není sám o sobě
expanze: je nutné rekonstruovat prostorovou geometrii a ukázat změnu fyzikálních vzdáleností.

## 8. Délka ℓ_d a most k SI

**Rozhodnutí:** G není předpovězená konstanta. G = κc³ℓ_d²/(4πλ₁ħ) zůstává podmíněným mostem; ℓ_d = ℓ_k se
nezavádí bez odvození. Cesty: odvodit ℓ_d z mechanismu (např. stabilní defekt), nebo ji určit z jednoho
pozorování a pak předpovídat jiné veličiny bez dalšího ladění (G pak není předpověď).

## 9. Gravitace jako termodynamika kauzální struktury (22. 10.)

**Rozhodnutí:** zkusit Jacobsonovu cestu (Jacobson 1995): entropie lokálních horizontů úměrná ploše + Clausiova
rovnost + Unruhova teplota ⇒ Einsteinovy rovnice jako stavová rovnice. Entropii dodá kauzální množina sama jako
počet molekul horizontu. Rozbor předpokladů: `docs/derivations/jacobson_cesta.md`; první měřidlo K.14.

Poznámky laboratoře:
- Volba substrátu (audit axiomů, oddíl 0) se zatím formálně nezapisuje. Autor ji podmínil tím, že se najde
  podstata gravitace. Jacobsonova cesta stojí na Lorentzovsky invariantním substrátu (lokální Rindlerovy horizonty
  v každém směru), takže její úspěch by volbu A podpořil.
- Cesta rovnice pole B_g f(Π) (rozhodnutí 1, K.13a, F1.7) zůstává v platnosti jako srovnání. Pokud Jacobsonova
  cesta projde, rovnice pole se nevolí, ale odvozuje, a γ = 1 přestává být vstupem.
- Chybějící kus je hmota s tenzorem energie a hybnosti (předpoklad 5). Bez hybnosti nelze hmotu teorie vložit
  do Clausiovy rovnosti.
