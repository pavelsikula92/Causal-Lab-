# List nálezů a objevů

Stav 20. 10. 2026: výsledky GTS Lab 2.x a 3.x a první běhy Causality Lab (kalibrace K.2–K.9, K.11, K.11b, F1.2, F1.2b, F1.3, F1.4, F1.5 předběžně, F2.5 d4); vývojová diagnóza K.11b a odvozený rozpor (a) s H1 ve 3+1.
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
  než na náhodu; zda s N mizí, změří F2.5 d3 (2+1, 8 000 až 32 000).
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
- Mizí přebytek H1 s rostoucím N (F2.5 d3), nebo jde o soustavnou vlastnost čtení Σ Π?
- Je nelokálnost vlastnost přírody, nebo měřidla (F0.4)?
- Jaká je dynamika řádu: lokální pravidlo, nebo vážený součet přes uspořádání (F6)?
