# List nálezů a objevů

Stav 21. 10. 2026: výsledky GTS Lab 2.x a 3.x a první běhy Causality Lab (kalibrace K.2–K.9, K.11, K.11b, F1.2, F1.2b, F1.3, F1.4, F1.5 předběžně, F2.5 d4 a d3, F2.7, K.12 s16, K.11c, F3.3 všechny tři kanály, F3.3b q90, F3.3c main_q99, F3.3d celý sken, F3.3e search a confirm, F3.3f capacity a source, F3.3g A); vývojová diagnóza K.11b a odvozený rozpor (a) s H1 ve 3+1.
Úplný seznam testů a jejich posledních stavů je v `docs/seznam_testu.md`.

Nejsilnější výsledek: dráhy v poli Π sledují geodetiku, tok zdroje splňuje Gaussův zákon a pole zdroje má daleko
newtonovský tvar se silou podle Gaussova zákona (F1.3); jedno uspořádání pole reprodukuje ve 2+1 (F1.2b).
Největší mezera: žádná struktura podobná hmotě zatím nevznikla z dynamiky sama.

## Drží

- **Pole lokalizovaného zdroje je daleko Greenova funkce se sílou podle Gaussova zákona** (F1.3, mechanismus,
  očekávaný nelokální operátor, okno 10–40 l_k): 3+1 g = −1,02 (1/r), síla 1,0045 předpovědi Q/4π; 2+1 g = −0,08
  (ln r), síla 1,055 předpovědi Q/2π. Dvě rozlišení shodná. Relativní odchylka od spojitého řešení roste k okraji
  (3+1: 2 %, 7 %, 20 % v r = 10, 20, 40), protože spojité u → 0 u stěny a pole nese konstantní posun od okraje;
  tvar ani síla tím zatížené nejsou.
- **Integrální A5′: zdroj vytváří pole podle předem vypočteného operátoru** (F1.2b, mechanismus, N = 32 000, jeden
  sprinkling, hladký gaussovský zdroj σ = 0,05). Platí ve 2+1; ve 3+1 je výsledek NEROZHODNUTO (viz níže). Zprůměrovaný operátor aplikovaný na pole
  vrací v průměru zdroj na všech šesti škálách průměrování (max |z| 1,32 ve 2+1, 1,16 ve 3+1) a chyba klesá
  stejným exponentem jako u pole bez zdroje (p_test / p_cal = 0,63 / 0,48 ve 2+1, 0,76 / 0,81 ve 3+1).
  Nulový model s permutovanými vahami jádra selhal (z do 19,7 ve 2+1). Ve 3+1 seed 0 prošel jen těsně
  (nulový model selhal na jediné škále, z = 3,23); replikace seedem 1 skončila NEROZHODNUTO: pole A i B splnilo
  (max |z| 1,51; p_test 0,51 ± 0,21 proti p_cal 0,44 ± 0,15), ale nulový model splnil A také (max |z| 2,64).
  Ve 3+1 při N = 32 000 (L/ℓ = 1,6–7,1) tedy jedno uspořádání pole od permutovaných vah spolehlivě nerozliší;
  je to nedostatek rozlišení, ne rozpor s A5′. Předchozí verze s ostrým zdrojem (F1.2) skončila jako selhání měřidla, viz níže.
- **Gaussův zákon pro tok trvalého zdroje** (F1.1, E45a): 2+1 i 3+1, N do 256 000, J/S_w = 1,003–1,004;
  model bez zachování toku 0,77–0,88.
- **Dráhy sledují geodetiku metriky Π²η** (F2.2, F2.4; E46, E46b, replikováno). V konformní reprezentaci ji sleduje
  každý maximální řetězec; ve vážené reprezentaci rozhoduje tvar akce (Σ Π² → Π⁴η: 0,0770 proti 0,0776).
- **Ve 3+1 rozhoduje o dráze tvar akce stejně jako ve 2+1** (F2.5 d4, mechanismus, vážená reprezentace, vložené
  pole měkkého bodového zdroje, 256 seedů na velikost, PROŠEL). Pravidlo E46b části B při N0 = 64 000: SPECIFIC.
  H1 (Σ Π) = 0,0524 ± 0,0026 proti P1 = 0,0432 (geodetika Π²η, pásmo ± 0,0108); H2 (Σ Π²) = 0,0799 ± 0,0025 proti
  P2 = 0,0776 (geodetika Π⁴η); H2 − H1 = 0,0275 ± 0,0023, tedy funkcionály se zřetelně liší. Při 32 000 také
  SPECIFIC: H1 0,0466 ± 0,0026, H2 0,0760 ± 0,0026. Měřidla prošla (zrcadlová identita přesně 0, váhy přesně,
  kontrola 1+1 E44b 0,0678 ± 0,0102 proti 0,0584, žádný neplatný seed). Tři výhrady: (1) podmínka „chyba H1
  neroste“ prošla jen díky statistické rezervě: bodově chyba vzrostla z 0,0035 na 0,0093 (dovoleno do 0,0109);
  (2) H1 leží při 64 000 jen 0,0016 pod horním okrajem pásma; (3) u velké části seedů je H1 = H2 nebo H1 = 0, protože
  ve 3+1 jsou řetězce krátké a nejtěžší řetězec bývá tentýž pro obě váhy, případně přímý. Geodetický zákon Π²η ve 3+1
  tím platí podle předem daného pravidla, ale přesnost čtení H1 ve 3+1 je na hraně a s N se nezlepšila.
- **Π je zpětně čitelná z hustoty řetězců** (K.3, E45c, replikováno); geometrie je čitelná z řádu ve 3+1 od N = 32 000 (K.1, E54).
- **Kauzální kužel jako měřidlo** (K.2, kalibrace, 2+1 i 3+1): zpoždění 0,982 a 0,986 hodnoty OTR, ohyb 0,874
  a 0,844 (v E47 0,913 a 0,919). Zpoždění je přesné; ohyb čte soustavně o 10–15 % níž, protože lokální
  lineární fit sklonu vyhlazuje strmou část profilu 1/y. Měřidlo rozliší γ = 1 od γ = 0,5, menší rozdíly ne.
- **Kalibrační měřidla Causality Lab prošla**: očekávaný operátor (K.4), nestrannost operátoru (K.5), identity
  toku (K.6), zrcadlová identita řetězců (K.7), sonda (K.8), nulové modely topologie (K.9).
- **Kauzální množina ve vzniklé geometrii Π²η dává ve spojité limitě □_g − ½R** (K.11 local, kalibrace): poměr
  operátoru k lokálnímu tvaru rovnice (a) 0,77 → 0,89 → 0,95 (2+1) a 0,87 → 0,954 → 0,9998 (3+1) pro šířky pole
  4, 8, 16 l_k; □_g bez −½R by dal 0,41 a 0,29. Vlastní objem intervalů (průměr Π^D z přesné projekce diamantu)
  souhlasí s Monte Carlem (K.11 diamond, max |z| 3,22 a 2,41; nejhrubší úroveň kvadratury neprošla).
  Konvergence po řádcích prošla (K.11 conv). Ve 2+1 je zakřivený operátor zkalibrovaný i pro vzdálené pole
  (K.11b far_d3: objemově vážený zbytek Q_far/Q_src = +0,015, −0,0001, +0,002 ve třech úrovních).
- **Zakřivený operátor je ve 3+1 zkalibrovaný i pro vzdálené pole, s hladkou (C²) geometrií** (K.11c, kalibrace,
  20. 10.). Kritérium K.11b beze změny (|Q_far(L3) − Q_far(L2)| <= 0,01 Q_src), jen geometrie mezi uzly je C² splajn
  místo lomené čáry. 3+1, síť krok 1 l_k / poměr 1,08 (far_d4_g1), PROŠEL: Q_far/Q_src = −0,0107, −0,0132, −0,0101
  (L1, L2, L3), |L3 − L2| = 0,0031; hladká reference při L2 −0,0106 (odchylka splajnu 0,0026); stačí už L1. Stejné pole
  s lomenou geometrií dává +2,26, tedy příčina selhání K.11b je potvrzená. 2+1 prošlo na obou sítích (far_d3_g1,
  far_d3_g2). 3+1 na hrubší síti krok 2 l_k / poměr 1,12 (far_d4_g2): MĚŘIDLO SELHALO podle kritéria
  (|L3 − L2| = 0,0122 a odchylka splajnu od reference při L2 0,0121, obojí nad 0,01). Popisně: L1 +0,0322, L3 +0,0301 a
  reference při L2 +0,0300 se shodují; odchyluje se jen splajn při L2 (+0,0179), téměř celé v jediném řádku r = 43 l_k
  (3,8e-5 proti 6,0e-5 v ostatních). Splajn sám je tedy v pořádku, nespolehlivá je úroveň L2 na řídké síti; i na síti g1
  je L2 nejdál od ostatních. Pravidlo: pole ve 3+1 počítat se splajnem na síti g1 (kalibrovaná); ve 2+1 stačí obě sítě.
  Tím je odblokované 3+1 pro F1.4 a F1.5 (jejich 3+1 části zatím neexistují).
- **Superpozice v rovnici (a) selhává ve druhém řádu, jak (a) předpovídá** (F1.4, 2+1, mechanismus): jádro
  (σ = 8 l_k) a slupka (r0 = 30, šířka 4 l_k); deficit Δ = u12 − u1 − u2 je kladný všude (společná studna je
  mělčí než součet: zdroj ve studni druhého slábne), ve středu +0,046 proti +0,057 spojité předpovědi (a);
  RMS(Δ_cs − Δ_a) = 0,0101 při toleranci 0,0111, tedy těsně. Π(0): jádro 0,7099 (spojitě 0,7099), slupka 0,915
  (0,901), obojí 0,671 (0,668). Dvě výhrady: (1) kvadratura F1.4 [160, 8, 48, 6, 24, 1] leží podle pozdějšího
  K.11b o 1,3 % zdroje mimo jeho toleranci; (2) slupka o šířce 4 l_k je na nelokálnost citlivá (K.11 local
  při šířce 4: 0,77), její pole vyšlo asi o 14 % slabší než spojité a deficit o 18 % menší; srovnání je s lokální
  spojitou rovnicí, kritérium se nemění. Replikace s kvadraturou podle K.11b by výhradu (1) odstranila.
- **PŘEDBĚŽNĚ: pole v rovnici (a) má vazbu κ/D a nelinearitu podle (a)** (F1.5, 2+1, mechanismus). Běh starší
  verze souboru F1.5 (kvadratura [160, 8, 48, 6, 24, 1], krok uzlů 1 l_k, 48 uzlů, 94 řádků), ne aktuální verze
  s kvadraturou podle K.11b a krokem 2 l_k. P1: k = 0,3363 proti 1/D = 0,3333 (□_g bez −½R by dal 0,99).
  P2: síla vzdáleného pole na jednotku zdroje s(q) = 0,826, 0,709, 0,560, 0,402 proti (a) 0,825, 0,709, 0,562,
  0,408 pro Π(0) = 0,89 až 0,58; RMS rozdílu 0,0034 při toleranci 0,10. Π(0) sedí se spojitou (a) na 0,003.
  Předběžné, protože kvadratura leží podle K.11b o 1,3 % zdroje mimo jeho toleranci a potvrzující běh aktuální
  verze ještě neproběhl; nelze to vydávat za nezávislé potvrzení. Vzdálené pole plochého operátoru je o 3 %
  silnější než spojité (A_flat 65,9 proti 63,9; F1.3 při σ = 1 l_k: 5,5 %).
- **Setrvačná a pasivní hmotnost jsou totéž** při faktorizaci tělesa (F5.0, E44c, 1+1, podmíněně).
- **Překryv minulostí je signatura geometrie**: 1,42–1,66 proti 1,00 (F0.1, E43).
- **Stínění saturací** (F3.1, E55a): zdroj na kapacitě omezí vlastní odtok, vnější tok klesne přesně v poměru
  C/F0_tube (k/k_pred = 1,0054), vně se nesaturuje nic. Bez zpětné vazby se oblast odporu nešíří.
- **Bez zpětné vazby nevzniká struktura** (F3.2, E55b): Δ_pers = +0,006 ± 0,006, E_loc = +0,004 ± 0,004,
  ekvivalentní nule v mezi 0,05; základní stav pro F3.3.
- Kapacita omezuje tok, ne počet vazeb (F1.8, E40, E41-A), verdikt ale není robustní vůči seedu (K.10, E49).

## Nedrží

- **Zpětná vazba přes kapacitu (C_v = C·Π_v) strukturu nevytvoří** (F3.3 capacity, mechanismus, 2+1, 128 sprinklingů,
  Π_sat = 0,5, m = 2, NEPROŠEL). Zpětná vazba saturuje 54 % událostí (bez ní 0,9 %), ale rovnoměrně: párové rozdíly
  proti vypnuté zpětné vazbě na stejných grafech d_pers = −0,007 ± 0,008, d_loc = −0,010 ± 0,004 (obojí EKVIVALENTNÍ
  nule, mez 0,05); s ní E_loc = −0,004 ± 0,001 (rovnoměrnější než náhoda), Δ_pers = +0,001 ± 0,002. Po vypnutí zdroje
  saturace zmizí během jedné vrstvy (0,54 → 0,006 → 0). Větev bez zpětné vazby zopakovala F3.2 s novými seedy
  (Δ_pers +0,008 ± 0,008, E_loc +0,006 ± 0,004). Výklad (popisný, ne test): v jednom kauzálním průchodu je stav
  události dán jen jejím vlastním zatížením, takže kanál A jen posune práh saturace (pro m = 2 na F/C ≈ 0,54) a
  saturace kopíruje obálku toku; saturovaná událost navíc pošle dál méně, čímž zatížení následníků sníží. Událost
  nepředává svůj stav potomkům (zásoba ΔQ se nepřenáší), takže trvalá struktura nemá nosič.
- **Zpětná vazba přes zdroj (s_v·Π_v) saturaci odstraní** (F3.3 source, mechanismus, 2+1, 128 sprinklingů, Π_sat = 0,5,
  m = 2, NEPROŠEL podle pravidla „saturace zmizí“): saturovaný podíl 0,9 % bez zpětné vazby, se zpětnou vazbou
  1 × 10⁻⁶ (saturované události jen ve 2 ze 128 sprinklingů). Záporná zpětná vazba: blízko kapacity hodiny zpomalí
  a zdroj klesne až na polovinu, takže zatížení kapacitu nepřekročí. Bez saturovaných událostí není co lokalizovat.
- **Směrování k nízké impedanci (váhy Π_k/I_k, I = 1) saturaci odstraní** (F3.3 routing, mechanismus, 2+1,
  128 sprinklingů, Π_sat = 0,5, m = 2, stacionární stav, NEPROŠEL podle pravidla „saturace zmizí“). Saturovaný podíl
  0,93 % bez zpětné vazby, se zpětnou vazbou 0 ve všech 128 sprinklinzích; nejmenší Π 0,546 > Π_sat. Stacionární
  stav zkonvergoval (21 průchodů), měřidla prošla. Výklad (popisný): událost posílá tok přednostně potomkům
  s rychlejšími hodinami, tedy méně zatíženým, takže zatížení se vyrovná pod kapacitu (záporná zpětná vazba, jak
  hlavička testu předpokládala). Všechny tři kanály F3.3 tedy NEPROŠLY: A zaplaví, B a C saturaci odstraní; vítěz
  pro krok D (kombinace) není.
- **Paměť přes přenášenou zásobu strukturu nevytvoří (při C = Q_0,9)** (F3.3b q90, mechanismus, 2+1, 128 sprinklingů,
  NEPROŠEL). Zásoba Q se přenáší po linkách k potomkům rovnoměrně a znovu nabízí; nic jiného se nemění. Paměť je
  měřitelná: saturovaný podíl 8,0 % → 9,9 %, průměrná držená zásoba +32 %. Struktura ale ne: párově d_pers =
  +0,0006 ± 0,0016, d_loc = −0,0015 ± 0,0014 (obojí EKVIVALENTNÍ nule s velkou rezervou, mez 0,05); absolutně
  Δ_pers +0,001 ± 0,002, E_loc −0,008 ± 0,002 (saturované události jsou mírně rovnoměrnější než náhoda, i bez
  paměti −0,007). Po vypnutí zdroje saturace i zásoba zmizí během jedné vrstvy 0,25 v obou větvích (10,5 % → 0,1 %).
  Měřidla prošla včetně bitové shody při vyhodnocení v pořadí boostované soustavy. Výklad (popisný): rovnoměrné
  dělení mezi potomky rozprostře zásobu do budoucího kuželu stejně jako tok, takže paměť se ředí a nevytvoří
  nosič, který by zůstal soustředěný podél časupodobného směru.
- **Kopírovaný pulz (pravidlo od hlavního rodiče) zaplaví oblast, ale lokalizovanou strukturu nevytvoří** (F3.3c
  main_q99, mechanismus, 2+1, 128 sprinklingů, NEPROŠEL). Pulz ω = 1 se nedělí, každý potomek dostane pulz rodiče,
  od kterého přitéká nejvíc toku, žije jen v saturovaných událostech a působí jako zdroj. Saturace 0,85 % → 30,5 %,
  ale rovnoměrně: párově d_pers = 0,000 ± 0,008, d_loc = −0,004 ± 0,004 (obojí EKVIVALENTNÍ nule); s pulzem
  Δ_pers +0,005 ± 0,002, E_loc +0,003 ± 0,002. Nové proti F3.3–F3.3b: struktura přežívá vypnutí zdroje. Saturovaný
  podíl po vypnutí 44 % → 8,1 % → 1,6 % → 0,05 % → 0 ve vrstvách po 0,25 (bez pulzu zmizí hned) a vzor po vypnutí se
  podobá vzoru před ním (r_s +0,66, +0,41, +0,12). Pulz je tedy paměť, která se neředí a nese i polohu, ale za
  buzení tvoří plošné zaplavení, ne lokalizovaný uzel. Měřidla prošla včetně bitové shody v pořadí boostované soustavy.
- **Práh materializace ΔQ_crit = 0,25 dá slabou, ale skutečnou lokalizaci, pod mezí 0,05** (F3.3d crit025,
  mechanismus, 2+1, 128 sprinklingů, NEPROŠEL). Model F3.3c main_q99 s jedinou změnou: pulz žije jen při
  ΔQ ≥ 0,25, jinak hned zanikne. Zaplavení zmizí: saturace 0,92 % → 4,4 % (F3.3c 30,5 %), pulz nese 2,3 % událostí.
  Párově d_pers = +0,010 ± 0,008 a d_loc = +0,019 ± 0,004 (obojí EKVIVALENTNÍ, mez 0,05). Lokalizace je přitom
  4,7 SE nad nulou: první nenulový d_loc v řadě F3 (s pulzem E_loc +0,020 ± 0,003, F3.2 +0,004 ± 0,004), zhruba
  2,5× pod mezí. Cena je paměť: po vypnutí zdroje saturace 6,5 % → 0,27 % → 0,01 % → 0 a pulz 3,6 % → 0,11 % → 0
  (F3.3c 44 % → 8,1 % → 1,6 %); r_s s posledním úsekem před vypnutím +0,13 v první vrstvě (málo buněk, jen hlášeno).
  Výklad (popisný, ne test): práh vyměňuje vytrvalost za soustředění; pulz přežije jen v nejsilněji zatížených
  místech, ta jsou shlukovaná, ale bez buzení zatížení pod prahem nevydrží. Měřidla prošla včetně bitové shody.
- **Silnější práh ΔQ_crit = 0,5 lokalizaci nezvýší, ale sníží** (F3.3d crit050, mechanismus, 2+1, 128 sprinklingů,
  NEPROŠEL). Pulz nese 0,16 % událostí, saturace 0,92 % → 1,2 %. Párově d_pers = +0,002 ± 0,004, d_loc = +0,007 ± 0,002
  (obojí EKVIVALENTNÍ, mez 0,05; d_loc 3,0 SE nad nulou, asi třetina hodnoty při 0,25). Paměť po vypnutí zmizí během
  jedné vrstvy (saturace 1,3 % → 0,01 % → 0, pulz 0,19 % → 0). Spolu s crit025: nad 0,25 slábne s prahem lokalizace
  i paměť, okno nad 0,25 v předem daném skenu není.
- **Slabší práh ΔQ_crit = 0,1 vrací zaplavení a paměť, lokalizace klesne** (F3.3d crit010, mechanismus, 2+1,
  128 sprinklingů, NEPROŠEL). Saturace 0,92 % → 17,3 %, pulz nese 14,8 % událostí. Párově d_pers = +0,007 ± 0,007,
  d_loc = +0,008 ± 0,004 (obojí EKVIVALENTNÍ, mez 0,05). Po vypnutí saturace 27,7 % → 3,6 % → 0,41 % → 0, pulz
  24,3 % → 2,8 % → 0,28 %, vzor podobný vzoru před vypnutím (r_s +0,58, +0,25, +0,16): mezi F3.3c a crit025.
- **Krajní práh ΔQ_crit = 1,0 pulz prakticky vypne** (F3.3d crit100, NEPROŠEL): pulz nese < 0,01 % událostí,
  saturace 0,92 % → 0,93 %, d_pers = +0,001 ± 0,001, d_loc = +0,000 ± 0,001 (EKVIVALENTNÍ); model = F3.2.
- **Předem daný sken prahu materializace NEPROŠEL v celém rozsahu** (F3.3d, ΔQ_crit 0,1 / 0,25 / 0,5 / 1,0, s F3.3c
  jako ΔQ_crit = 0). Párové d_loc: −0,004 / +0,008 / +0,019 / +0,007 / +0,000; saturace 30,5 / 17,3 / 4,4 / 1,2 /
  0,9 %; saturace v první vrstvě po vypnutí 8,1 / 3,6 / 0,27 / 0,01 / 0 %. Lokalizace má jediné maximum kolem 0,25
  (souhlasí s odhadem autora 0,15–0,25), velikostí asi 2,5× pod mezí 0,05; paměť s prahem monotónně klesá. Uzel,
  který by se soustředil a zároveň vydržel bez buzení, v rozsahu skenu nevzniká. Maximum hledá a ověřuje F3.3e.
- **Ani optimum prahu nedosáhne meze** (F3.3e search, průzkum, ZAZNAMENÁNO, 2+1, 64 nových sprinklingů, mřížka
  ΔQ_crit 0,05–0,40 po 0,025 na stejných kauzálních množinách). Hledaná veličina min(d − 2 SE) má široký, mělký
  vrchol kolem 0,15–0,2 (souhlasí s odhadem autora): nejlépe ΔQ_crit = 0,2, d_pers = +0,034 ± 0,009,
  d_loc = +0,015 ± 0,006, min(d − 2 SE) = +0,002; průchod by potřeboval > 0,05. d_pers je kladné 3–3,6 SE
  v celém rozsahu 0,1–0,275 (v F3.3d při 0,1 a 0,25 v rámci chyby nula; rozdíl ~1,3 σ). Paměť po vypnutí klesá
  s prahem monotónně (saturace v 1. vrstvě 5,3 % při 0,05 → 0,62 % při 0,2 → 0,03 % při 0,4). Měřidla prošla.
  Záznam se při prvním běhu neuložil (zastaralé core\provenance.py u autora); čísla z výpisu, běh se opakuje.
- **Práh materializace neprojde ani v optimu na nezávislých datech** (F3.3e confirm, mechanismus, 2+1, 128 nových
  sprinklingů, ΔQ_crit = 0,2 vybraná z dat, NEPROŠEL). Párově d_pers = −0,009 ± 0,007, d_loc = +0,010 ± 0,004 (obojí
  EKVIVALENTNÍ, mez 0,05). Saturace 0,92 % → 6,1 %, pulz nese 3,9 % událostí; po vypnutí saturace 10,4 % → 0,57 % →
  0,01 % → 0. Kladné d_pers ze search (+0,034 ± 0,009) se nezopakovalo; d_loc se drží slabě kladné (2,2 SE). Tím je
  model prahu (kopírovaný pulz od hlavního rodiče, C = Q_0,99) na pevném kauzálním uspořádání uzavřen v celém rozsahu
  i v optimu: lokalizace existuje, ale je ~5× pod mezí, a stálost nad úroveň bez pulzu nevzniká.
- **Práh + zpětná vazba přes kapacitu vytvoří samoudržující, ale plošnou saturovanou fázi, ne uzel** (F3.3f capacity,
  mechanismus, 2+1, 128 sprinklingů, ΔQ_crit = 0,2 z dat F3.3e, C·Π, Π_sat = 0,5, m = 2, NEPROŠEL). Saturace: bez
  pulzu 1,0 %, jen práh 7,6 %, práh + kanál 99,1 %; všechny saturované události mají Π = Π_sat = 0,5. Párově proti
  větvi bez pulzu d_pers = +0,004 ± 0,009, d_loc = +0,036 ± 0,005 (obojí EKVIVALENTNÍ, mez 0,05). d_loc je největší
  v řadě F3, ale při 99 % saturovaných událostí měří E_loc (Moranovo I je symetrické vůči χ ↔ 1 − χ) shlukování
  zbývajícího 1 % nesaturovaných míst, tedy „děr“, ne uzlu. Nové: stav přežije vypnutí zdroje – saturace 99 % →
  51 % → 34 % → 17 % → 6,8 % → 3,4 % ve vrstvách po 0,25 (samotný práh 12,7 % → 0,9 % → 0,04 % → 0), pulz nese
  totéž, vzor slabě podobný vzoru před vypnutím (r_s ≈ +0,1). Výklad (popisný): kladná zpětná vazba C·Π s pulzem
  jako zdrojem dává bistabilní médium – oblast přejde celá do saturované fáze s pomalými hodinami a ta se sama drží
  a jen pomalu odtéká (pravděpodobně od okraje pohlcování). Při homogenním buzení zaplaví celou oblast; zda může
  existovat ohraničená „kapka“ této fáze, tento test neměří. Měřidla prošla včetně bitové shody.
- **Práh + zpětná vazba přes zdroj saturaci odstraní** (F3.3f source, mechanismus, 2+1, 128 sprinklingů, s_eff =
  (s + P_in)·Π, NEPROŠEL podle pravidla „saturace zmizí“). Ve větvi s kanálem 0 saturovaných událostí ve všech
  128 sprinklinzích (nejmenší Π 0,519 > Π_sat); větve bez pulzu (1,0 %) a jen s prahem (7,6 %) shodné s částí
  capacity (stejné sprinklingy). Pomalé hodiny ztlumí zdroj i pulz pod kapacitu dřív, než práh pulz zapálí –
  záporná zpětná vazba jako v F3.3 source. Krok D je tím uzavřen: kladný kanál (C·Π) dá samoudržující plošnou fázi,
  záporný (s·Π) nic.
- **Bistabilní saturovaná fáze nemá hranici: kapka buď vyprchá, nebo zaplaví oblast** (F3.3g A, mechanismus, 2+1,
  128 sprinklingů, model F3.3f capacity bez pole a bez zpětné vazby na směrování, zapálení disku r0 = 0,15 (s_ign = 5,
  t ∈ [3; 4)) na pozadí s_bg = 0–0,8, NEPROŠEL pro všechny hodnoty). Při s_bg ≤ 0,2 kapka zanikne (při 0,2 asi za
  2,5 času); při 0,25 vyprchává asi 4,5 času (93 % zaniklo do pozdního okna), při 0,3 a víc vyplní celou oblast
  (R_rms → 0,35–0,40, prstenec u okraje 18–94 % saturovaný). Přechod leží mezi 0,25 a 0,3 a stojící kapka u něj
  nevzniká. Podstatné je, jak se to děje: hned po vypnutí zapálení se kapka rozplyne do řídkého oblaku – R_rms skočí
  z 0,15 na ~0,31 během 0,25–0,5 času (zhruba rychlostí světla), kompaktnost B klesne z ~0,8 na ~0,45 a saturace
  uvnitř C_drop z ~0,5 na 0,13–0,25 nezávisle na pozadí. Saturovaná „fáze“ tedy netvoří ostrou doménu, ale řídký
  oblak, jehož hustota pak podle pozadí klesá (vyprchá) nebo roste (zaplaví). Výklad (popisný): událost rozdělí tok
  rovnoměrně mezi potomky v celém budoucím kuželu, takže vliv každé saturované oblasti se šíří rychlostí světla a
  hranici nic nedrží; omezit to může jen zúžení kuželů. Část B (poruchy) se nespouští: její předpoklad (kapka
  z části A) nenastal. Měřidla prošla.
- **Poznámka k měřidlu Δ_pers (otevřené, ne verdikt):** větev bez pulzu, tedy týž model F3.2 na různých nezávislých
  sadách sprinklingů, dává Δ_pers +0,006 ± 0,006 (F3.2), +0,008 ± 0,008 (F3.3), +0,000 ± 0,007 (F3.3d),
  −0,009 ± 0,008 (F3.3e search) a +0,021 ± 0,006 (F3.3e confirm). Rozptyl je větší, než odpovídá uvedeným SE
  (χ² ≈ 11 kolem váženého průměru +0,008 při 4 stupních volnosti, p ≈ 0,03); bootstrapová SE může Δ_pers podhodnocovat zhruba 1,5×. Verdikty
  F3 to nemění (všechny jsou EKVIVALENTNÍ s rezervou), ale tvrzení typu „x SE nad nulou“ u d_pers a d_loc je třeba
  brát opatrně, dokud to nerozhodne kalibrace rozptylu na nezávislých sadách.
- Lokální růst nevytvoří prostoročas: 0 ze 7 pravidel (F6.1).
- Trojrozměrný prostor z dynamiky nevznikl (F6.2).
- Mez počtu vazeb (stupeň ≤ K₀) neplatí.
- Původní místní zákon Π = f(F/C) dává zrychlení ~1/r³; nahradila ho A5′.
- Bodová A5′ je při dosažitelném N neprokázatelná, rozhoduje šum (F1.6, E45b).
- Lokální hladový krok geodetiku nesleduje (F2.9).
- Saturace se nešíří podle nekorigovaného profilu toku (F0(r_sat) = C).
- **Čtení Π z hustoty řetězců ve 3+1 při nejnižší hustotě** (K.3, část d4): při ρ = 4·10⁶ chyba 0,114 proti
  limitu 0,10. Ve 3+1 je čtení Π použitelné až od ρ = 1,6·10⁷; ve 2+1 prošlo v celém rozsahu.

## Selhání měřidel (bez závěru o teorii)

- **K.12 s16: článek zdroj → pole ve 2+1 není při σ = 16 ℓ_k kalibrován pro test typu F2.7, ale odchylka klesá s σ**
  (kalibrace, sobě podobná rodina σ = 16, q = 0,0025, R = 200, NEPROŠEL). Geodetika v poli kauzální množiny
  P_cs = 3,771 proti spojitému P_a = 3,354, poměr 1,124 (tolerance P_a/16, tj. 1,0625); při σ = 8 (pole F2.7) 1,206.
  Přebytek tedy klesl z 0,206 na 0,124 při zdvojnásobení σ (formálně ~σ^(−0,73); zda klesá k nule, rozhodne s32).
  Kvadratura zkonvergovala (kontrolní a hlavní úroveň se liší o 0,0001 proti P_a/64 = 0,052), Π(0) 0,7202 proti 0,7099.
  Poměr sklonů cs/spojité v r = 1,25σ / 2,5σ / 3,75σ: 1,158 / 1,097 / 1,049 (při σ = 8: 1,245 / 1,129 / 1,056).
  Gaussův deficit toku 1 − tok(R/2)/tok(3σ) zůstává ~0,10 při obou σ (spojité −0,014): v sobě podobné rodině je
  R = 12,5σ, takže okrajový jev by škáloval s σ a nezmenšoval se; rozliší s08_R200.
- **F1.2, integrální A5′ s ostrým zdrojem**: reziduum očekávaného pole 0,0875 q proti toleranci 0,05.
  Ostrá hrana zdroje není na škále nelokality ℓ_k zobrazitelná. Toleranci jsem nastavil bez měření na plné
  mřížce (moje chyba). Revize F1.2b změnila jen tvar zdroje na gaussovský, kritéria zůstala stejná.
  Popisně: A i B byly splněné, nulový model selhal (z do 21), ale bez platného měřidla se to nepočítá.

- **Zakřivený operátor ve 3+1 nekonverguje ve vzdáleném poli** (K.11b far_d4): objemově vážený zbytek
  vzdálených řádků Q_far/Q_src = 3,5, 2,3, 1,9 ve třech úrovních kvadratury (tolerance 0,01). Falešný zdroj je
  tedy větší než skutečný a s rozlišením klesá jen pomalu; F1.4 a F1.5 ve 3+1 s tímto měřidlem počítat nejde.
  Příčina nalezená při vývoji K.11c (19. 10., vývojový běh, ne test): zbytek nevzniká z dlouhých intervalů
  samých, ale z po částech lineární interpolace geometrie Π mezi uzly. Operátor reaguje na zakřivení geometrie
  (člen −R/2, tedy ∇²Π^D) a lineární interpolace má ∇²Π v uzlech jako delta funkce; dlouhé téměř světelné
  intervaly jsou tenké jehly, které je sbírají, ve 3+1 s vahou ~r. Řádek r = 40 l_k: lineární geometrie
  +3,4e-3 a při zdvojení libovolného parametru kvadratury 1,3e-3 až 4,8e-3; hladká geometrie −4,15e-5, stabilní
  na 0,2 %; C² splajn přes tytéž uzlové hodnoty −4,152e-5; C¹ splajn (PCHIP) +4,4e-4. Oprava měřidla je
  v K.11c (geometrie jako C² splajn); K.11b zůstává NEPROŠEL podle svého kritéria. Ve 2+1 se lineární
  a hladká geometrie liší ve vzdálených řádcích asi o 1e-4 (r = 30, síť F1.5); K.11c far_d3 změří, kolik to
  dělá v objemově váženém zbytku, tedy jak velkou výhradu to přináší F1.4 a F1.5.

## Objevy

- **Ohyb dráhy je vlastnost substrátu**: ohýbají se i nejkratší řetězce a průměr přes všechny maximální řetězce,
  ten je navíc asi 5× méně šumný.
- **Stínění saturací** (výše): bez zpětné vazby se saturace za oblast buzení nešíří.
- **Těžký chvost toku je vlastnost topologie linků, ne stupňů**: graf se stejnými stupni a prohozenými vazbami má
  užší rozdělení (99. percentil 4,45 proti 5,06) a skoro nesaturuje.
- **Graf linků kauzální množiny je tuhý**: výměny vazeb, které zachovají platné uspořádání, projdou jen v 0–0,8 % pokusů.
- **Odezva toku na impulz doznívá exponenciálně** rychlostí asi 1 za jednotku času, téměř nezávisle na N.
- **Plný ohyb světla vyžaduje c_eff = c·Π² u hmoty**: váha Π sama dá poloviční ohyb, konformní hustota nulový.
- **Součet Π leží nad spojitou předpovědí, součet Π² na ní, ve 2+1 i ve 3+1** – zatím nevysvětleno. Ve 2+1 asi o 10 %;
  ve 3+1 (F2.5 d4) o 8 % při N0 = 32 000 a o 22 % při 64 000 (rozdíl −0,0058 ± 0,0037, nevýznamný), zatímco H2
  sedí na 3 %. Stejný vzorec ve dvou dimenzích ukazuje spíš na vlastnost čtení Σ Π (malé rozdíly vah, shody řetězců)
  než na náhodu. F2.5 d3 (2+1, 768 seedů na velikost, ZAZNAMENÁNO) ukázal, že s N nemizí: přebytek H1 − P1 =
  +0,0044 ± 0,0016, +0,0029 ± 0,0015, +0,0044 ± 0,0013 při N0 = 8 000, 16 000, 32 000 (+10 %, +7 %, +10 %; trend
  e(8 000) − e(32 000) = 0,0000 ± 0,0021 → PERSISTENT), zatímco H2 − P2 = +0,0004 až +0,0011 (pod 1,5 %). Pravidlo B
  je SPECIFIC při každé velikosti. Přebytek je tedy soustavný a v tomto rozsahu na N nezávislý. Teorie tvrdí, že
  maximální řetězce Σ Π konvergují ke geodetice Π²η (spojitá limita akce); pro Σ Π² to data potvrzují, pro Σ Π
  zůstává posun asi 10 %, který s N neklesá. Buď konvergence Σ Π je velmi pomalá, nebo limita Σ Π není přesně
  geodetika Π²η; nerozhodnuto.
- **Pohyb sleduje geodetiku pole, které kauzální množina sama vytvořila; to pole je ale strmější než spojité (a)**
  (F2.7, řetězec, 2+1, 512 seedů, NEROZHODNUTO). Zdroj σ = 8, q = 0,01; samokonzistentní pole (a) na kauzální
  množině zkonvergovalo (13 iterací, reziduum 3,7e−5), Π(0) = 0,7113 proti spojitému 0,7099. Pohyb (b = 10, T = 40):
  H1 = +2,027 ± 0,071 proti P_a = +1,677 (geodetika spojitého řešení (a)), poměr 1,21; |H1 − P_a| ± 2 SE =
  0,207 až 0,493 proti pásmu d_alt/4 = 0,419 → podle kritéria nerozhodnuto (ani nulové pole, ani □_g bez −R/2
  s P = 3,49 to ale nejsou). Hlášeno mimo kritérium: geodetika v poli kauzální množiny P_cs = +2,023, tedy
  H1 − P_cs = +0,004 ± 0,071 – pohyb leží přesně na geodetice pole, ve kterém se řetězce vážily. Celý rozdíl proti
  P_a vzniká v kroku zdroj → pole: Π_cs podél dráhy r = 10…30 je 0,757 → 0,883, spojité (a) 0,746 → 0,856, tj.
  vzrůst o ~15 % strmější (P_cs/P_a = 1,21). H2 = +3,740 ± 0,068 proti P2_a = 3,322 (poměr 1,13, hlášeno; P2 v poli
  kauzální množiny nespočteno). Napětí s F2.5 d3: tam Σ Π leželo ~10 % nad geodetikou daného pole, zde na ní
  (jiné nastavení; nerozhodnuto, zda jde o týž jev). K.11c kalibroval ve 2+1 dalekou bilanci toku operátoru,
  ne tvar samokonzistentního profilu v r = 10–30.
  Rozbor uloženého profilu (22. 10., bez nového běhu): ve 2+1 má (a) tvar ∇²Π = (κ/3)Π²ρ_s, takže mimo zdroj
  musí být tok r·dΠ/dr konstantní (Gaussův zákon; spojité řešení 0,1197 od r ≈ 30 až ke stěně). Pole kauzální
  množiny tuto vlastnost nemá: r·dΠ/dr vrcholí 0,129 u r ≈ 22 a pak ke stěně klesá (0,124 při r = 36, 0,111 při 56,
  0,073 při 79, 0,006 při 99); Π_cs dosáhne 0,9993 už u r = 89 (spojité 0,986). Poměr sklonů cs/a klesá plynule
  z 1,28 (r = 4) přes 1,24 (10), 1,13 (20), 1,06 (30) k 1,0 u r ≈ 43 a pod 1 dál. Hloubka jámy přitom souhlasí
  (1 − Π(0): 0,2887 proti 0,2901). Nejde tedy o jinou konstantní vazbu (ta by škálovala sklon všude stejně), ale
  o profil, který mimo zdroj není harmonický a „cítí“ stěnu R = 100: jako kandidáti zbývají okrajový jev řádků
  u stěny (uříznutá minulost v zakřivené geometrii, kterou odečet vakua E[B_flat]1 neruší přesně) nebo stínění
  nelokálním operátorem; rozliší je sken R.
- **Bistabilní saturovaná fáze** (F3.3f capacity, pozorování k NEPROŠEL výsledku): kopírovaný pulz s prahem a se
  zpětnou vazbou kapacity C·Π vytvoří první stav v řadě F3, který po vypnutí zdroje sám vydrží několik vrstev
  (polovina saturace po ~1 vrstvě 0,25, stopa ještě po 1,25). Je plošný, ne lokalizovaný.
- **Šum operátoru je prostorově korelovaný** (F1.2b, pozorování, ne test): chyba průměrovaného pole klesá
  s exponentem p ≈ 0,5–0,8, hluboko pod d/2 pro nekorelovaný šum (1,5 ve 2+1, 2 ve 3+1). Přesné čtení pole
  proto potřebuje velké L/ℓ, tedy hustší sprinkling, ne víc řádků operátoru.

## Odvozené důsledky (bez testu)

- **Rovnice pole (a) a hypotéza H1 nemohou ve 3+1 platit současně** (odvozeno 19. 10. z lokální limity
  □_g − R/2, kterou K.11 numericky potvrdil jen pro konformní geometrii Π²η; pro jiné geometrie se opírá o
  Belenchia et al. 2016). Pro statickou metriku −Π²dt² + Π^(−2γ)dx² (γ je parametr PPN; Π²η má γ = −1,
  H1 má γ = 1) dává symbolický výpočet lokální tvar (□_g − R/2)Π:
  2+1: γ = −1: 3Π⁻²∇²Π; γ = 0: 2∇²Π + Π⁻¹|∇Π|²; γ = 1: Π²∇²Π + 2Π|∇Π|²;
  3+1: γ = −1: 4Π⁻²∇²Π + 2Π⁻³|∇Π|²; γ = 0: 2∇²Π + Π⁻¹|∇Π|²; γ = 1: 2Π|∇Π|² (bez členu s ∇²Π).
  Lineární vazba je κ/c(γ), c(γ) = 2 − γ(n − 1), n počet prostorových rozměrů. Ve 3+1 tedy (a) v geometrii
  Π²η (vazba κ/D, ve 2+1 předběžně F1.5) dává γ = −1, tj. žádný ohyb světla ani Shapirovo zpoždění (Cassini: γ = 1);
  v geometrii H1 dává rovnici 2Π|∇Π|² = κρ_s, takže mimo zdroj je Π konstantní a newtonovské pole nevzniká.
  Ve 2+1 se H1 s (a) nevylučuje (vazba κ). Důsledek: F7.1 ve smyslu „H1 při rovnici (a)“ nemá ve 3+1 řešení;
  buď se změní tvar rovnice pole (na co operátor působí, nebo zda nese −R/2), nebo H1. Rozhodnutí je na autorovi.
  Odvození vazby κ/D → κ = 16π (kapitola 4.7) platí jen pro konformní geometrii.

## Otevřené otázky

- Rozliší jedno uspořádání ve 3+1 pole od nulového modelu při vyšším N (větší L/ℓ)?
- Je kladné vychýlení A5′ na velkých škálách (v obou dimenzích, nevýznamné) náhoda, nebo systematika?
- Vznikne se zpětnou vazbou trvalá lokalizovaná struktura bez ladění (F3.3)?
- Má taková struktura setrvačnost a gravitační pole (F4, F5)?
- Proč maximální řetězce Σ Π leží asi o 10 % nad geodetikou Π²η, když Σ Π² leží přesně na geodetice Π⁴η, a proč
  se to s N nemění (F2.5 d3)? Je spojitá limita akce W = Σ Π opravdu geodetika Π²η?
- Proč samokonzistentní pole (a) na kauzální množině ve 2+1 nemá mimo zdroj konstantní tok r·dΠ/dr a ke stěně
  R zplošťuje, zatímco u zdroje je strmější (F2.7)? Je to okrajový jev měřidla (posune se se stěnou při větším R),
  nebo nelokální stínění (zůstane na stejných r)?
- Je nelokálnost vlastnost přírody, nebo měřidla (F0.4)?
- Udrží saturovanou oblast pohromadě zúžení budoucích kuželů v místě pomalých hodin (c_eff = cΠ² z F7), když
  rovnoměrné dělení toku do celého kuželu ji rozplyne rychlostí světla (F3.3g A)?
- Jaká je dynamika řádu: lokální pravidlo, nebo vážený součet přes uspořádání (F6)?
