# List nálezů a objevů

Stav 21. 10. 2026: výsledky GTS Lab 2.x a 3.x a první běhy Causality Lab (kalibrace K.2–K.9, K.11, K.11b, F1.2, F1.2b, F1.3, F1.4, F1.5 předběžně, F2.5 d4 a d3, F2.7, K.12 s16, s08_R200 a s32, F6.3 V1 a V2, F6.3b V1, V2 a A, K.14 všech 5 částí, F3.6 rest a move, F3.10 unit, attraction, momentum, universality a law, F3.10b law, K.13a všechny tři části, F1.7 V, W a N, F0.6 fixed, F3.7 search, F3.8 search, F3.9 search, K.11c, F3.3 všechny tři kanály, F3.3b q90, F3.3c main_q99, F3.3d celý sken, F3.3e search a confirm, F3.3f capacity a source, F3.3g A); vývojová diagnóza K.11b a odvozený rozpor (a) s H1 ve 3+1.
Úplný seznam testů a jejich posledních stavů je v `docs/seznam_testu.md`.

Nejsilnější výsledek: dráhy v poli Π sledují geodetiku, tok zdroje splňuje Gaussův zákon a pole zdroje má daleko
newtonovský tvar se silou podle Gaussova zákona (F1.3); jedno uspořádání pole reprodukuje ve 2+1 (F1.2b).
Největší mezera: žádná struktura podobná hmotě zatím nevznikla z dynamiky sama.

## Drží

- **F3.10b law: pravidlo předává pole partnera jako sílu se správným průběhem na vzdálenosti** (mechanismus, 2+1,
  d = 8 / 12 / 16 / 20, 64 sprinklingů na d, okno lineární odezvy [1,25 d; 2,5 d], PROŠEL 23. 10., seed 2;
  replikace chybí). Přijatá hybnost dvojice (B − A) proti předpovědi prvního řádu z větve A (pole skutečného
  partnera, svazky osamoceně, žádná data z B): exponent k_meas = 0,845, k_pred = 0,946, Δk = −0,101 ± 0,051,
  |Δk| + 2 SE = 0,203 ≤ 1/4. Spojitá statická úsečka ve stejných oknech dává 0,927 (jádro samo 1; rozdíl dělá
  retardované nabíhání pole v okně). Poměry naměřené / předpovězené: 0,931 ± 0,018 / 0,908 ± 0,030 /
  0,964 ± 0,038 / 1,021 ± 0,051; posun separace v okně 0,4–2,1, hluboko pod d/4. Měřidlo na všech 256
  sprinklinzích: C = A, zrcadlo i soustava bitově, 0 remíz. Čtení je omezené:
  - jádro 1/(2πτ) a znaménko vazby jsou vstup; test ukazuje, že diskrétní dynamika pole převádí na sílu věrně
    (do ~10 %), ne že zákon 1/r plyne z řádu;
  - při d = 8 a 12 je síla o 7–9 % menší než předpověď (2–4 SE), při d = 16 a 20 sedí; příčina neznámá
    (odhad gradientu z rodičů proti centrální diferenci, druhý řád posunu);
  - vysvětluje trend staré části F3.10 law (0,61 → 0,98): ten dělala hlavně vadná kvazistatická předpověď.
  Část force (ρ_d = 1 pro každé d) má vlastní verdikt a se stejným seedem počítá tatáž data.
- **F3.10 attraction: dva oddělené svazky se přes svá retardovaná pole přitahují** (mechanismus, 2+1, 32 sprinklingů
  ~41 500 událostí, d = 12, PROŠEL, seed 0; replikace seedem 8 PROŠLA). Relativní zrychlení po odečtení driftu osamocených
  svazků a = −0,0080 ± 0,0011. Předpověď z retardovaného pole statické úsečky (ε = 0,1 zamčené z F3.8/F3.9,
  λ a μ z kalibrace jednoho svazku) je −0,0111, poměr 0,72. Měřidlo: C = A volba po volbě (vazba jen polem),
  přepočet v soustavě β' = 0,6 i zrcadlo x → −x volí identicky. Žádný nový mechanismus: dynamika je pole z F3.8.
  Čtení je omezené:
  - přitahování plyne z toho, že pravidlo ohýbá nesený směr k pomalejším hodinám a pole partnera hodiny zpomaluje;
    jádro 1/(2πτ) a znaménko vazby jsou vstup (F3.8), ne výstup;
  - samo přitahování ještě neříká nic o zákonu síly, hybnosti ani univerzalitě (samostatné části níže a v plánu).
  Replikace seedem 8 (23. 10., PROŠEL; nezávislé sprinklingy i kalibrace): a = −0,0092 ± 0,0014 (a + 2 SE = −0,0064),
  zrychlení z hybnosti +0,0038 a −0,0041 (seed 0: +0,0037 a −0,0036), Ṗ_L + Ṗ_R = −0,051 ± 0,042; měřidlo OK, 0 remíz.
  Nezávislá kalibrace dala φ_ref 26,9 místo 24,7, tj. λ o 9 % menší (0,00371): šum kalibrace λ mezi seedy je
  ~9 %. Pro znaménko přitahování ani pro F3.10b (stejné λ v simulaci i předpovědi) nevadí; poměry ke kvazistatické
  aproximaci (tady 0,92) jsou proto i z tohoto důvodu jen orientační.

  Část momentum (d = 12, toky 1 a 2, okno lineární odezvy [15, 30], 32 sprinklingů, PROŠEL, seed 0): svazky si
  předávají opačné a stejně velké hybnosti. Ṗ_L = +1,552 ± 0,091, Ṗ_R = −1,545 ± 0,120, součet +0,007 ± 0,158
  proti toleranci 0,387; rozlišení bilance je tedy ~20 % jednotlivé změny. Zrychlení z hybnosti +0,0094 a −0,0047:
  dvakrát těžší svazek zrychluje poloviční rychlostí (to plyne z Ṗ_L = −Ṗ_R a M ∝ tok, nic navíc). Čtení je omezené:
  - vazba je symetrická (stejné λ, stejné jádro, každý svazek je zdrojem i příjemcem úměrně svému toku), takže
    akce a reakce je v pravidle zabudovaná pro statické zdroje; test ukazuje, že ji retardace, diskrétnost a
    ztráty při slučování v lineární odezvě nerozbijí. Pole hybnost nenese, při rychlejším pohybu by retardace
    bilanci porušit mohla; to tu rozlišit nejde;
  - separace se v okně zkrátí o 1,49 ± 0,18 (lineární odezva drží); kvazistatická aproximace tu dává poměr 1,06
    (těžiště) a 0,84 (hybnost L), bez verdiktu.
  Replikace seedem 8 (23. 10.): NEROZHODNUTO, tedy **nepotvrzeno**. Ṗ_L = +1,467 ± 0,056, Ṗ_R = −1,543 ± 0,152,
  součet −0,076 ± 0,179; |součet| + 2 SE = 0,434 nad tolerancí 0,376, |součet| − 2 SE pod ní. Bilance je v obou
  bězích slučitelná s nulou (+0,007 ± 0,158 a −0,076 ± 0,179), test na to ale nemá sílu: chyba Ṗ_R je sama skoro
  polovina tolerance. Seed 0 prošel těsně (0,323 proti 0,387). Rozhodnutí by vyžadovalo nový běh s předem
  stanoveným větším počtem sprinklingů; kritéria se nemění (R132).

  Část universality (d = 12, vnější zdroj R s tokem 2, který reaguje jen na své pole, L s tokem 1 a 2 na týchž
  sprinklinzích, okno lineární odezvy, 32 sprinklingů, PROŠEL, seed 0): zrychlení z nesené hybnosti
  a_P,L = +0,00919 ± 0,00046 (tok 1) a +0,00904 ± 0,00062 (tok 2), párový rozdíl +0,00014 ± 0,00066 proti
  toleranci 0,00152. Těžiště souhlasí (+0,0092 ± 0,0026 a +0,0105 ± 0,0022). Přijatá hybnost se s tokem zdvojí
  (1,55 → 3,02), zrychlení zůstane. Vnější pole je v obou konfiguracích bitově stejné. Čtení je omezené:
  ohyb nezávisí na toku nosiče, takže univerzalita odezvy je v pravidle zabudovaná; test ukazuje, že ji vlastní
  pole (u toku 2 silnější), slučování a diskrétnost nerozbijí. Nezávislé potvrzení principu ekvivalence to není.
  Kvazistatická aproximace (počítá s pohyblivým R, tady R stojí, poměry proto nejsou srovnatelné) dává 0,82 a
  0,78, bez verdiktu; s částí momentum (R reaguje) souhlasí a_P,L 0,0094 v rámci chyby.
  Replikace seedem 8 (23. 10.): NEROZHODNUTO, tedy **nepotvrzeno**. a_P,L = +0,00908 ± 0,00029 (tok 1) a
  +0,00980 ± 0,00052 (tok 2), párový rozdíl −0,00072 ± 0,00052 proti toleranci 0,00157; |rozdíl| + 2 SE = 0,00176.
  V seedu 0 byl rozdíl +0,00014 ± 0,00066 (opačné znaménko), dohromady tedy bez systematického rozdílu, ale
  replikace kritérium nesplnila. Tok 1 dal v obou bězích stejné zrychlení (0,00919 a 0,00908).

  Část law (d = 8 / 12 / 16 / 20, 24 sprinklingů na d) skončila NEROZHODNUTO. Všechna d přitahují; exponent
  těžiště je k = 0,77 ± 0,57 proti k_pred = 1,11 (stejný fit předpovědi), šum těžiště při d = 16 a 20 je 36–45 %.
  Diagnostika bez verdiktu: zrychlení z hybnosti Ṗ/M je zhruba 5× přesnější. Jeho poměr k bodové předpovědi
  roste se vzdáleností (0,61 / 0,76 / 0,91 / 0,98, průměr L a R) a exponent je ~0,58 místo 1,11. Bilance hybnosti
  u symetrické dvojice: Ṗ_L + Ṗ_R je v mezích ~10 % jednotlivé změny při všech d.
  Běh seedem 8 (23. 10., nahrazená část, jen pro záznam): NEROZHODNUTO, k = 1,51 ± 0,49 proti aproximaci 1,12;
  těžiště při d = 16 nerozlišené (−0,0011 ± 0,0010). Zrychlení z hybnosti 0,0052 / 0,0042 / 0,0032 / 0,0029
  (seed 0: 0,0058 / 0,0042 / 0,0037 / 0,0034) – průběh se reprodukuje, plošší než aproximace; rozhodne F3.10b.

  QA review autora (23. 10., po běhu attraction a law) našel vadu předpovědi, ne dynamiky. Kvazistatická
  předpověď dosazuje současnou vzdálenost do vzorce statického zdroje a svazek bere jako úsečku. S kalibrovaným λ
  a μ bodové úsečky při d = 8 projdou v okně skrz sebe (nejmenší r = −4,7) a při d = 12 se separace zkrátí z 12
  na 5,5. Poměry a/a_pred i trend diagnostiky jsou proto zkreslené samotnou předpovědí; vysvětlení „rozlehlým
  zdrojem“ je jen jedna z možností. Verdikt law zůstává NEROZHODNUTO (R132), law se nereplikuje a nahrazuje ji
  F3.10b: předpověď prvního řádu ze skutečného diskrétního pole skutečného partnera v okně lineární odezvy.
  Další nálezy QA:
  - model zachovává skalární tok, ne čtyřhybnost: při slučování různých směrů je |P| > F a přebytek se zahodí;
    změřeno ~1–3·10⁻⁴ toku, ~0,03–0,05 % x-hybnosti na událost, proti efektům zanedbatelné. Pole hybnost nenese;
  - pořadí zpracování je lineární rozšíření kauzálního řádu (teď ověřené v kódu); boost řád zachová, takže
    kontrola soustavy porovnává stejný výpočet;
  - kontroly C, zrcadla a soustavy teď porovnávají i tok a směry, rozdíl začínající numerickou remízou se hlásí
    zvlášť; přibyla deterministická část unit (λ = 0, bilance toku, znaménko kopnutí). unit PROŠEL 23. 10.
    (4 sprinklingy, vše OK, 0 remíz; 93–96 % kopnutých událostí míří k partnerovi; přebytek čtyřhybnosti při
    slučování 1·10⁻⁵ až 3·10⁻⁴ toku); znovu na nezávislých sprinklinzích (seed 8) PROŠEL 23. 10. (vše OK, 0 remíz,
    87–97 % kopnutých událostí k partnerovi, přebytek 8·10⁻⁵ až 1,7·10⁻⁴).
- **F3.6 (rest i move): tok, který nese svůj časupodobný směr jako stav, zůstane u světočáry zdroje a pohybuje se
  s ním (setrvačnost); směr převzatý z posledního kroku ne** (mechanismus, 32 sprinklingů 2+1, ~32 800 událostí, zdroj 204 událostí, PROŠEL). Rychlost rozšiřování
  R_rms v okně [5, 11] v násobcích c: rovnoměrné dělení (řada F3) 0,770 ± 0,003, směr z posledního kroku
  0,393 ± 0,007, nesený směr **0,054 ± 0,006** (tolerance lokalizace 1/4). Nesený svazek drží celý tok (~300 v každé
  vrstvě), šířka roste pomalu z 3,0 na 4,0 za 24 ℓ (zhruba difuzně). Nosiče se slévají z 73 na 27 v jedné vrstvě:
  svazek se stahuje k menšímu počtu světočar. Přepočet ve vztažné soustavě β' = 0,6 volil v každé události stejně
  (pravidlo nevybírá soustavu).

  Kopnutí jedné vazby je ~0,45–0,52 rapidity. Hybnost proto nemůže být posunutí mezi sousedními událostmi (taková
  paměť se rozpadne a tok se šíří rychlostí ~0,4), ale musí být vnitřní stav toku. Čtení je omezené:
  - pravidlo používá geometrii na škále vazby, ne jen řád a počty (kritérium (a) auditu), a volba potomka potřebuje
    budoucnost (b), stejně jako celá řada F3;
  - lokalizace neznamená vázaný objekt: svazek drží, protože tok se zachovává a nese směr. Jeho šířku nic neváže.

  Část move (zdroj β = 0,5, ~49 200 událostí, PROŠEL):
  - nesený směr je lokalizovaný (v_s 0,050 ± 0,006) a v soustavě zdroje stojí: rychlost těžiště +0,0075 ± 0,0043
    proti toleranci β/4 = 0,125 a alternativě „zůstane v klidu v laboratoři“ −0,5;
  - šířka se vyvíjí stejně jako u zdroje v klidu (2,94 → 3,74), jak to vyžaduje princip relativity;
  - směr z posledního kroku se šíří (0,349 ± 0,012), rovnoměrné dělení také (0,782).

  Hypotéza „lokalizace v invariantním substrátu vyžaduje paměť směru (hybnost)“ tím má oporu a nesený směr dává
  i setrvačnost. Replikace rest seedem 1 (23. 10., PROŠEL): nesený směr 0,0534 ± 0,0073, krok 0,367 ± 0,012,
  rovnoměrné dělení 0,774 ± 0,004; šířka 2,94 → 3,92, nosiče 73 → 26 v jedné vrstvě, soustava β' = 0,6 volí
  identicky – souhlasí se seedem 0 ve všech číslech. Replikace move seedem 1 (23. 10., PROŠEL): nesený směr
  lokalizovaný (v_s 0,0528 ± 0,0056), v soustavě zdroje stojí (v_c +0,0045 ± 0,0061 proti toleranci 0,125 a
  alternativě −0,5), šířka 2,89 → 3,77; krok 0,369, rovnoměrné dělení 0,777. F3.6 je tím uzavřený: první
  výsledek laboratoře, který nebyl vložen do pravidla, potvrzený na nezávislých datech.
- **F1.7 V a W: řešič rovnice B_g f(Π) = κρ_s (3+1, H1, kontinuum) je ověřený a ve slabém poli dávají všichni
  kandidáti s f(1) = 0 newtonovskou odezvu, rovnice (a) žádnou** (mechanismus, PROŠEL 23. 10.). V: symbolické tvary
  reprodukují C5 a λ₁(γ) = 2 − γ(d_s − 1); řešič proti analytickému lineárnímu řešení 3,3·10⁻⁹; ln Π při σ = 4 a 16
  dává stejný profil Π(r/σ) do 1,1·10⁻¹² (Π(0) = 0,6253 při u_lin(0) = −0,3). W (u_lin(0) = −10⁻⁴): ln Π, Π − 1,
  (Π⁻¹ − 1)/(−1), (Π² − 1)/2, (Π³ − 1)/3 se od lineárního řešení liší o 2·10⁻⁵ až 1,8·10⁻⁴ (mez 10⁻²); negativní
  kontrola f = Π má lineární koeficient přesně 0. Kandidáty to podle očekávání nerozliší; rozhodne N. γ = 1 je tu
  vstup ansatzu, ne výsledek.
- **F1.7 N: silné pole – každý kandidát má regulární statické řešení jen do konečné síly zdroje; ln Π je druhý
  nejlepší z pěti** (mechanismus, ZAZNAMENÁNO 23. 10.; 3+1, H1, kontinuum, gaussovský zdroj, síla jako u_lin(0)).
  Poslední síla s regulárním řešením / první bez něj:
  - (Π⁻¹ − 1)/(−1): −0,5 / −0,7; ln Π: −0,3 / −0,5; Π − 1 a (Π² − 1)/2: −0,2 / −0,3; (Π³ − 1)/3: −0,1 / −0,2;
  - pořadí je monotónní v k rodiny (Π^k − 1)/k (ln Π je k → 0): menší k, širší rozsah;
  - u slabopolní větve se objevuje druhá, silnopolní větev (Π(0) ≈ 0,09–0,40); obě se k mezní síle přibližují
    (u k = −1 při −0,5: 0,41 a 0,28), což vypadá jako záhyb (dvě řešení splynou a zaniknou);
  - vnější náboj roste nad lineární hodnotu (ln Π: 1,02 / 1,08 / 1,22 / 1,51 při −0,03 / −0,1 / −0,2 / −0,3):
    nelineární pole samo přispívá k hmotnosti viděné zvenku.
  Čtení: konečná mez sama o sobě vadou není (i obecná relativita má mez kompaktnosti statických těles), ale zdroj
  tu nemá tlak a převod u_lin(0) na kompaktnost GM/(Rc²) není udělaný, takže srovnání s neutronovými hvězdami
  zatím nejde. Kandidáty rozliší až druhý řád vnějšího pole (Odvozené důsledky: β). Ohraničení: γ = 1 je vstup.
- **K.15 (všechny čtyři části): tenzor neseného toku T^ab = (1/ΔV₄) Σ F u^a u^b τ_vazby měří prach, protiběžné proudy
  a izotropní plyn správně a s pravidlem sloučení F = |P| se 4-hybnost zachovává přesně** (beam: kalibrace,
  3+1, 64 sprinklingů ~28 700 událostí, PROŠEL 23. 10.). V klidové soustavě svazku (Landauova soustava) p = 0 a π = 0,
  v laboratoři E/ρ = γ² = 1,5625, |q|/E = v = 0,6, S_xx/E = v² = 0,36 a p/E = 0,12 v každé buňce 2 i 3 bez rozptylu;
  kontrola soustavy (tenzor ze zboostovaných souřadnic) do 1e-9. Obě pravidla sloučení dávají totéž, protože stejné
  směry se slévají beze ztráty. Ověřuje algebru rozkladu (E, q, p, π).
  counter (dvě protiběžné populace ±v, PROŠEL 23. 10.): q_x/E = −0,001 ± 0,009 (0), S_yy = 0 přesně, p/E = 0,1127 ±
  0,0006, S_xx/E = 0,338, π_xx/E = 0,225 – vše o ~6 % pod hodnotou emitovaných směrů (0,12 / 0,36 / 0,24), v toleranci.
  Deficit vysvětluje směšování: 9,1 % toku v okně už nemá rychlost v (pravidlo F = |P|; součtové 7,5 %), protože
  se protiběžné řetězce slily. Nepružné sloučení přemění usměrněný pohyb na klidovou hmotu (při F = |P| se energie
  zachová, tlak klesne): populace se chladí. Buňky 2 a 3 dávají p/E 0,109 a 0,113 s rozptylem 0,025 a 0,017.
  iso (izotropní plyn, PROŠEL 23. 10.): q_x, q_y, q_z / E = +0,003 / +0,000 / +0,009 (± 0,005), π_xx, π_yy, π_xy / E
  pod 0,0012 ± 0,0017 (izotropie), p/E = 0,1119 ± 0,0005 proti 0,12 (o 7 % níž; směšování 18,7 % toku, sloučení 7,9 %).
  conserve (izotropní plyn, PROŠEL 23. 10.): s pravidlem sloučení F = |P| se 4-hybnost na řezech t = 2,5–4,5 zachovává
  na 6·10⁻¹⁷ (přesně); dnešní pravidlo (tok = součet toků) ztrácí energii průměrně 0,67 % (nejvýš 1,9 %) za 2 jednotky
  času při 7 % sloučeného toku. Žádný tok v okně nekončí na stěně. Důsledek: dnešní pravidlo přenosu (řada F3) porušuje
  zachování 4-hybnosti, kdykoli se slévají různé směry; odvozené pravidlo F = |P| ho splňuje přesně. Zda ho přijmout
  jako pravidlo laboratoře, je rozhodnutí autora.
  Replikace seedem 2 (23. 10., všechny čtyři části PROŠLY): beam přesně; counter p/E 0,1132 ± 0,0006 (seed 0: 0,1127),
  q_x/E −0,015 ± 0,009, směšování 8,2 %; iso p/E 0,1114 ± 0,0005 (0,1119), složky q a π do 0,004 ± 0,005,
  směšování 19,8 %; conserve F = |P| 7,6·10⁻¹⁷, součtové pravidlo ztrácí 0,81 % energie (nejvýš 2,0 %).
  K.15 je uzavřený.
- **K.13a: diskrétní BD operátor v geometrii H1 (γ = 1) konverguje v prvním řádu ke spojité limitě □_g − R/2,
  hustotní i kuželová část zvlášť** (kalibrace, 3+1, všechny tři části PROŠLY 23. 10.; deterministická kvadratura,
  replikace seedem nemá smysl). lin, poměr k limitě při σ = 4 / 8 / 16 / 32 ℓ_k:
  - w = 1, r = 0: 0,743 / 0,894 / 0,961 / 0,987; r = σ: 0,806 → 0,992; r = 2σ: 0,537 → 0,977;
  - gaussovské w: r = 0: 0,596 → 0,977; r = σ: 0,689 → 0,986; r = 2σ: 0,700 → 0,986;
  - odchylka klesá ostře ve všech šesti řádcích, místní exponent log2 1,56–1,65 (pomaleji než (ℓ_k/σ)²); při σ = 32
    je 0,8–2,4 % proti mezi d_alt/4 (12,5–57,5 %); kvadratura hlavní proti kontrolní ≤ 0,02 (mez d_alt/16);
  - hustotní část P i kuželová část Q (odezva na změnu rychlosti kuželů, tj. na prostorovou část metriky) konvergují
    každá zvlášť; koeficient u ∇²u pro w = 1 jde k −1 (člen −R/2 = −∇²u je reprodukovaný).
  Substrát tedy nese i kuželovou část, na které závisí γ: geometrie s ohybem světla je na kauzální množině
  vyjádřitelná. Důsledek pro rovnici (a) viz Odvozené důsledky. geom a valid (σ = 4):
  - optický čas proti optimalizované dráze konverguje kvadraticky (log2-poměry 1,97–2,00 při γ = 1 i 0, mez 1,75);
  - objemy intervalů se zakřivenými kužely proti Monte Carlu (20 diamantů, 400 000 vzorků): max |z| 2,73;
  - tenké intervaly proti mpmath: 40/40, nejhorší relativní chyba 2,1·10⁻¹⁰.
  valid:
  - Π = 1 proti K.4 na relativní úrovni ~10⁻¹² (řádky r = 0, 4, 8);
  - γ = −1 proti centrální diferenci přístroje K.11: rozdíly 2·10⁻⁷ až 2,4·10⁻⁵ proti tolerancím ~10⁻³–10⁻².
  Při σ = 4 dává diskrétní operátor 0,74 (w = 1) a 0,66 (gaussovské w) spojité hodnoty na ose: to je systematika
  ℓ_k/σ, jejíž pokles s σ testuje lin.
- **K.14: počet molekul horizontu měří plochu horizontu, v zakřivené geometrii vlastní plochu** (kalibrace,
  všech 5 částí PROŠLO; měřidlo entropie pro Jacobsonovu cestu; převzatá konstrukce Barton et al. 2019, ne vlastní
  tvrzení teorie). Plochý prostor proti konstantě a^(d) z nezávislé kvadratury, na jednotku plochy v jednotkách ℓ^(d−2):
  - 1+1: 0,3325 ± 0,0043 a 0,3241 ± 0,0043 při ρ = 1 a 4 (a = 1/3; 20 000 sprinklingů);
  - 2+1: 0,2182 ± 0,0089 a 0,3612 ± 0,0127 (a ρ^(1/3) = 0,2189 a 0,3474);
  - 3+1: 0,1716 ± 0,0014 a 0,3407 ± 0,0027 (a ρ^(1/2) = 0,1732 a 0,3464; a^(4) = √3/10).

  Všechna |z| <= 2,2. Dvě z = −2,17 při ρ = 4 kontrolní běh v kontejneru nepotvrdil jako soustavné (1+1, ρ = 4,
  60 000 sprinklingů: 0,3329 ± 0,0025, z = −0,17). Konformně plochý prostor Π(y) = 0,75 + 0,25 cos(2πy/40)
  (Π 0,5–1): exponent vlastní plochy Π^k
  - 2+1: k = 1,060 ± 0,065 (cíl 1), C/a = 0,992, model bez volných parametrů χ² 18,7 / 16 pásů, souřadnicová plocha
    1151, objem 778;
  - 3+1: k = 1,964 ± 0,055 (cíl 2), C/a = 0,959, χ² 26,0 / 16, souřadnicová plocha 11 746, objem 1 261.

  První předpoklad Jacobsonova argumentu (entropie úměrná vlastní ploše s univerzální konstantou) tak kauzální
  množina dodává aspoň tam, kde se mění hustota a kužele zůstávají ploché. Zbývají zakřivené kužele, hmota
  s hybností a Clausiova rovnost (docs/derivations/jacobson_cesta.md). V 3+1 je amplituda o 4 % pod hlavním řádem;
  může jít o konečnou škálu ℓ/λ, rozhodl by jiný poměr λ.
- **F6.3b V1 a V2: přesná konstrukce dynamického uspořádání je ověřená ve vakuu i v limitách** (validace, PROŠEL obě).
  V1: přesné vakuum = vazby sprinklingu ve 4 sprinklinzích N = 36 000 (žádná vazba u stěny se tentokrát nevyskytla),
  hrubá síla s proměnným Π bez rozdílu v 6 sprinklinzích. V2: při Π_sat = 1 končí řešení po 1. iteraci a je prvek po
  prvku rovné pevnému uspořádání ve všech 4 sprinklinzích; okna 1,5 a 3,0 dávají na řešení s hmotou totožné vazby
  (~1,95–1,97 mil., vakuum ~1,52 mil.: hmota přidá ~28 % vazeb). Podíl událostí, jimž se mění Π mezi iteracemi,
  klesá z ~0,5 na 0,26 za 9 iterací (první sprinkling). Část A smí běžet.
- **F6.3 V1: uspořádání z kuželů závislých na hodinách je ve vakuu přesně kauzální množina** (validace, PROŠEL).
  Při Π = 1 dává přesná konstrukce (tranzitivní uzávěr místní relace, Hasseho diagram) vazby prvek po prvku totožné
  s vazbami sprinklingu ve 4 sprinklinzích N = 36 000 (~1,52 mil. vazeb každý); s proměnným Π (jáma hodin do 0,5)
  se shoduje s hrubou silou (úplný uzávěr) v 6 sprinklinzích po 1 200 událostech bez jediného rozdílu (~20 000
  vazeb každý); všechny vazby jdou dopředu v čase (acykličnost).
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

- **F3.6c (eta03 i eta06): slučování s přesnou 4-hybností (F = |P|) svazek z horkého zdroje neudrží; rozbíhání
  jen zpomalí zhruba na polovinu** (mechanismus, 2+1, 32 sprinklingů, ~61 000 událostí, 202 událostí zdroje,
  obě části NEPROŠEL; měřidlo ve všech podmínkách v pořádku). Část eta03 (η = 0,3). Rychlost rozšiřování
  v okně [6, 18] (těžišťová soustava, R_rms kolem těžiště):
  - sym: B 0,112 ± 0,005, A 0,114 ± 0,005, volný pohyb C 0,220 ± 0,002 (práh vázanosti C/4 = 0,055);
  - asym: B 0,117 ± 0,007, A 0,118 ± 0,007, C 0,214 ± 0,004;
  - cold (kontrola): A = B = C bit po bitu, 0,049 ± 0,003 (difuzní chvění polohy, stejné jako F3.6).

  Bilance B 4·10⁻¹⁵ a 1·10⁻¹⁵, B − A −0,0011 ± 0,0007 a −0,0010 ± 0,0016 (bez měřitelného rozdílu). Slučování
  (~136 sloučení, průměrná ekvivalentní rapidita 0,24) proudy zchladí: vnitřní rozptyl klesne z 0,30 na 0,13 v t = 5
  a dál pomalu na 0,088 v t = 29. Počet proudů klesá z 78 na 35. Přírůstek klidové hmotnosti v B je 3,9 %, staré
  pravidlo A ztrácí 3,7 % energie. Šířka B roste lineárně (3,05 → 5,56 za 24), tedy balisticky, ne difuzně.

  Čtení podle pravidel autora: pravidlo F = |P| funguje účetně (4-hybnost přesně), šířku neřeší a v tomto režimu se
  od starého pravidla měřitelně neliší. Slučování je jen chlazení: proudy, které se jednou rozejdou, se už nepotkají
  a nic je nevrací. Pro Jacobsonovu cestu je F = |P| dál nutné (A ztrácí energii), na lokalizaci nic nemění.

  Část eta06 (η = 0,6, NEPROŠEL):
  - sym: B 0,237 ± 0,007, A 0,249 ± 0,007, C 0,465 ± 0,003 (práh C/4 = 0,116);
  - asym: B 0,243 ± 0,007, A 0,251 ± 0,006, C 0,440 ± 0,003.

  Bilance B 4·10⁻¹⁵ a 1·10⁻¹⁵. Tady se B od A měřitelně liší: B − A = −0,0115 ± 0,0023 a −0,0079 ± 0,0021, tedy B je
  o 3–5 % užší. A ztrácí 12 % energie, B ji převádí na klidovou hmotnost (+14 %). Vnitřní rozptyl klesne z 0,60 na
  0,27 v t = 5 a pak se skoro zastaví (0,22 v t = 29). Počet proudů se ustálí kolem 50: slučování se rychle vyčerpá.

  Souhrn obou částí: poměr v_s(B)/v_s(C) je 0,51–0,55 při obou η. Slučování zchladí jen proudy, které se potkají
  hned u zdroje, pak už se rozbíhající proudy nepotkají. F = |P| je správné účetnictví a při vyšším η svazek
  nepatrně zúží (víc váhy těžším, už sloučeným proudům), vazbu ale nedává. Co drží šířku, zůstává otevřené.

- **F3.9 search: pole s dosahem spolu s konečnou kapacitou kanálu vázaný stav nedá; kapacita jen tlumí odraz a krátce
  drží plošinu, pak se svazek pomalu hroutí** (mechanismus, průzkum, 32 sprinklingů 2+1, ZAZNAMENÁNO; kandidát podle
  pravidla z hlavičky: žádný). F_ref = 7,0, φ_ref = 43,5. Sklon šířky v W1 / W2:
  - free: +0,052 / +0,031;
  - replikace F3.8: ε 0,1 −0,010 / −0,031 (kolaps), ε 0,5 −0,030 / +0,235 (kolaps a odraz) — souhlasí s F3.8;
  - **cap ε 0,1 κ 1: +0,0018 ± 0,0030 (vázaný v W1) / −0,027 ± 0,004 (kolaps v W2)**: šířka drží ~3,4 od t ≈ 13 do 27,
    pak klesá k 2,95;
  - cap ε 0,1 κ ¼: +0,018 / −0,019 (vystoupá na 3,9 a vrací se);
  - cap ε 0,5: kolaps a odraz, kapacita odraz tlumí (W2 +0,205 / +0,102 / +0,030 pro κ 4 / 1 / ¼ proti +0,235 bez
    kapacity);
  - cap ε 0,01: jako free.

  Kontrola soustavy prošla ve všech větvích. Čtení (popisné): kapacita v této podobě posílá přebytek méně souběžnými
  cestami, ale se stejným neseným směrem. Rozšiřuje polohy, nenese však hybnost ven, takže to není tlak v pravém
  smyslu. Slévání nosičů navíc průměruje směry, tedy ztrácí relativní pohyb (nepružné, „chladí“). Odpovídá tomu
  pomalé smršťování po plošině, podobné kontrakci chladnoucí samogravitující látky.
- **F3.7 search: ohyb k pomalejším hodinám s hodinami z místní zátěže svazek nesváže, ale zahřeje; kovariantní
  zúžení kuželů ho nesváže také** (mechanismus, průzkum, 32 sprinklingů 2+1, ~49 100 událostí, ZAZNAMENÁNO; kandidát
  pro confirm podle pravidla z hlavičky: žádný). Hodiny Π = 1 − 0,5 min(F/C, 1)², C = κ F_ref, F_ref = 6,0. Sklon
  šířky v oknech W1 [12, 28] a W2 [28, 44]:
  - free (F3.6): 0,037 ± 0,007 a 0,038 ± 0,007;
  - lens κ = 4 / 1 / ¼: 0,12 / 0,47 / 0,59 v W1, 0,43 / 0,43 / 0,33 v W2. Všechny rostou rychleji než volný svazek.
    Při κ = ¼ a 1 se svazek nejdřív smrští (R 3,0 → 2,5 do t ≈ 9), pak exploduje až ke stěně. Ohyb tedy míří
    dovnitř, ale bez tlumení dráhy prolétnou středem a šum odhadu gradientu z rodičů svazek zahřívá (střední ohyb
    0,05–0,30 rapidity na krok, Π klesá až na 0,5);
  - narrow κ = 4 / 1: šíří se 0,17 / 0,37 až ke stěně. Při κ = ¼ je kužel v jádře tak úzký, že platí záložní volba
    a jádro se chová jako free (0,044 ± 0,007), zatímco okraj se rozlévá (nosičů ~2 000).

  Kontrola soustavy (β′ = 0,6) prošla ve všech větvích. Čtení: s hodinami přímo z místní zátěže je pole drsné
  a krátkodosahové (jen uvnitř svazku), jeho gradient je hlavně šum diskrétní zátěže. Vázání by potřebovalo hladké
  pole s dosahem, tedy pole buzené tokem přes rovnici pole (A5′, retardovaná Greenova funkce kauzální množiny), nebo
  tlumení (výměnu hybnosti). Obojí je nová hypotéza po neúspěchu, s vlastními kritérii.
- **F0.6 fixed: zúžené kužele kreslené v souřadnicích (F3.3h) vybírají preferovanou soustavu** (mechanismus,
  brána relativity, 128 sprinklingů N = 72 000, NEPROŠEL; všech 6 rozhodujících srovnání PŘEKROČENÍ). Tytéž
  události, tentýž řád, buzení i odečet v klidové soustavě zdroje; jen pravidlo kreslí kužele v soustavě β:
  - s_bg = 0: v klidové soustavě přežije 86 % sprinklingů (N pozdě 8,6), při β = 0,3 i 0,6 zaniknou všechny.
    Rámcový efekt −8,56 ± 0,48 je stejně velký jako celý efekt mechanismu (+8,56; tolerance 2,14).
  - s_bg = 0,2: N pozdě 116 → 78,9 → 18,7 pro β = 0 / 0,3 / 0,6 (rámcový efekt −37 ± 1 a −97 ± 1 proti toleranci
    29). Útvar ujíždí podél osy pohybu: těžiště +0,224 ± 0,002 a +0,339 ± 0,002 proti toleranci r0/4 = 0,0375;
    R_late roste o 0,08 a 0,19.
  - Změna soustavy změní počet povolených potomků u 93–96 % událostí v saturované oblasti (5 % při s_bg = 0).

  Pozitivní kontrola (Π_sat = 1, kužel = světelný kužel) byla v obou soustavách bitově shodná a bilance sedí:
  rozdíl dělá pravidlo, ne transformace. Důsledek: výsledky postavené na tomto pravidle (F3.3h, kompaktní doména
  F6.3b A, jejíž lokální relace kreslí kužele stejně) nejsou fyzikou. Myšlenka zúžení kuželů tím nepadá, padá jen tato
  souřadnicová podoba. Kovariantní verze by zužovala kužel vůči nesenému směru toku (klidové soustavě hmoty, F3.6).
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

- **F1.9a: β diskrétního operátoru mimo zdroj, MĚŘIDLO SELHALO jen v kontrole kvadratury při σ = 4**
  (mechanismus, 25. 10., 252 úloh, 5,8 h). Jediná nesplněná podmínka: |c_hlavní − c_kontrolní| v řádku 6σ při
  σ = 4 je 0,090 (H1) a 0,076 (γ = −1) proti toleranci 1/16. Při σ = 8 a 16 je to 0,023 a 0,004. Stěna (rozdíl
  ≤ 0,002), kontrola γ = −1 (odchylka 0,014 a 0,012 při σ = 16) a brána přístroje (K.13a lin PASSED se stejným
  core) jsou v pořádku. Podle protokolu se o β nerozhoduje. Popisně, bez verdiktu:
  - H1, řádky rozhodnutí: c_disc(σ = 4, 8, 16) ve 4σ −0,114, −0,029, −0,002 proti spojité hodnotě směsi +0,008;
    v 6σ +0,016, +0,047, +0,055 proti +0,058. Odchylka klesá s faktorem 3,3–3,9 na zdvojení σ (zhruba (ℓ_k/σ)²),
    takže by jinak splnila PROŠEL (≤ 1/4 při σ = 16, neroste). β_disc při σ = 16: 1,001 a 0,973. Odchylka od 1
    v 6σ je vlastnost směsi (spojitá hodnota 0,971, |∇L| klesá s exponentem −2,09 místo −2), ne diskrétnosti.
  - Kontrola γ = −1 jde ke 2 stejně: 4σ +2,081, +2,045, +2,022 (spojitě 2,008).
  - Řádek 3σ (jen hlášen, ocas zdroje): H1 −0,419, −0,281, −0,224 proti −0,195, odchylka klesá.

  Čtení: selhal nejhrubší bod žebříku (pole šíře 4 ℓ_k s kontrolní kvadraturou [80, 12, 8, 3, 8, 4]). Podle R132 se
  kritérium zpětně nemění. Formální verdikt o β by vyžadoval nový test s předem daným protokolem (např. žebřík
  σ = 8, 16, 32 se stejnými pásmy). Hodnota takového testu je potvrzovací, pásmo Δc = 1 je od naměřených čísel daleko.

- **F6.3b A: v uspořádání ze zúžených kuželů vzniká kompaktní, pomalu rostoucí doména saturované fáze, ale jen
  s vnějším tokem; samostatná kapka nevzniká** (mechanismus, 32 sprinklingů, N = 36 000, MĚŘIDLO SELHALO: soubor
  řešení není při s_bg = 0,1 stacionární, C_late se mezi polovinami průměrovacího okna liší o −0,117 ± 0,014 proti
  toleranci 0,05). Verdikty po pozadích (kritérium F3.3h A):
  - s_bg = 0: dynamické uspořádání ZANIKLO (100 %), ale později než pevné: doména drží do t ≈ 6–7 (pevné do ≈ 4,5).
  - s_bg = 0,1 a 0,2: NEROZHODNUTO. Doména přežila na všech sprinklinzích a je doménou fáze (C_drop pozdě 0,69 ± 0,01
    a 0,84 ± 0,004, kompaktnost B 0,87–0,95), ale stále roste: růst R 0,035 ± 0,002 a 0,039 ± 0,001 je těsně kolem
    r0/4 = 0,0375 a poloměr se na konci okna dál plíživě zvětšuje.
  - s_bg = 0,3 a 0,4: NEPROŠEL (růst nad r0/4 a zapálení pozadí, prstenec 0,08 a 0,31).

  Pevné uspořádání se stejným buzením se rozplyne v řídký oblak (jako F3.3g); vazeb dynamického uspořádání je o 7–43 %
  víc než ve vakuu. Popisně: zúžení kuželů mění rozplynutí rychlostí světla na kompaktní doménu, která roste
  rychlostí ~0,01. Čtení je omezené (audit axiomů 22. 10.): kužele jsou kreslené v souřadnicích (F0.6, rychlé běhy:
  na soustavě závislé), pulz působí jako zdroj toku (v rozporu s rozhodnutím 4) a řešení je samosouhlasné přes celý
  prostoročas. Platné je proto jen negativní čtení: ani s těmito výhodami nevznikne kapka bez vnějšího toku
  (s_bg = 0). Doména při s_bg > 0 je podle rozhodnutí 3 struktura udržovaná vnějším tokem. Opakování s opraveným
  měřidlem stacionarity se nedoporučuje, dokud neprojde F0.6.
- **F6.3 V2: limita Π_sat → 1 se lišila o jednu vazbu kvůli vadě starého měřidla u stěny** (validace, NEPROŠEL;
  revize F6.3b). Při Π_sat = 1 se dynamické řešení (vakuové uspořádání přesné konstrukce) v 1 ze 4 sprinklingů lišilo
  od pevného uspořádání z core.order.links o 1 vazbu z 1 966 510. Rozbor: dvojice (4725, 5857) leží u stěny válce
  (r = 0,492 a 0,498), polovina jejího Alexandrovova intervalu je mimo nasypanou oblast, efektivní ρV = 15,4
  (nominálně 30,2), přímý počet 0 – je to skutečná vazba. core.order.links vybírá kandidáty podle nominálního
  objemu intervalu v neomezené oblasti (λ = 25), a takové vazby u stěny proto vynechá; přesná konstrukce je najde.
  F6.3 V1 prošel jen proto, že jeho 4 sprinklingy takovou vazbu neměly. Vada se týká jen vnější stěny (za poloměrem
  pohlcování 0,4 v testech F3, kde tok odchází), takže dřívější výsledky prakticky neovlivní; core.order.links
  zůstává beze změny. Část (b) prošla: okna 1,5 a 3,0 dávají na řešení s hmotou totožné vazby (~1,96 mil.) ve všech
  4 sprinklinzích. Revize F6.3b používá ve všech větvích přesnou konstrukci a V1 připouští jen prázdné vazby u stěny.
- **K.12 s16: článek zdroj → pole ve 2+1 není při σ = 16 ℓ_k kalibrován pro test typu F2.7, ale odchylka klesá s σ**
  (kalibrace, sobě podobná rodina σ = 16, q = 0,0025, R = 200, NEPROŠEL). Geodetika v poli kauzální množiny
  P_cs = 3,771 proti spojitému P_a = 3,354, poměr 1,124 (tolerance P_a/16, tj. 1,0625); při σ = 8 (pole F2.7) 1,206.
  Přebytek tedy klesl z 0,206 na 0,124 při zdvojnásobení σ (formálně ~σ^(−0,73); zda klesá k nule, rozhodne s32).
  Kvadratura zkonvergovala (kontrolní a hlavní úroveň se liší o 0,0001 proti P_a/64 = 0,052), Π(0) 0,7202 proti 0,7099.
  Poměr sklonů cs/spojité v r = 1,25σ / 2,5σ / 3,75σ: 1,158 / 1,097 / 1,049 (při σ = 8: 1,245 / 1,129 / 1,056).
  Gaussův deficit toku 1 − tok(R/2)/tok(3σ) zůstává ~0,10 při obou σ (spojité −0,014): v sobě podobné rodině je
  R = 12,5σ, takže okrajový jev by škáloval s σ a nezmenšoval se; rozliší s08_R200.
- **K.12 s08_R200 a s32: přebytek pole kauzální množiny v oblasti pohybu není efekt stěny a podle předem daného
  pravidla má dno ~10 % (s32 zatím bez certifikace)** (kalibrace, 2+1, rovnice (a), 24. 10.). s08_R200 NEPROŠEL:
  P_cs/P_a = 1,227 při stěně R = 200 proti 1,206 při R = 100. Čtení „stěna“ = nepodpořena: ve stejné oblasti
  r = 48 Gaussův deficit 0,088 → 0,074, přebytek 0,206 → 0,227 (obojí změna pod čtvrtinu). Blíž ke stěně deficit
  na stěně závisí (r = 64: 0,190 → 0,110; r = 80: 0,414 → 0,143), tok kauzální množiny ale klesá s r i při R = 200
  (0,116 v r = 24 → 0,099 v r = 80, spojité 0,106 konstantní): pole mimo zdroj není harmonické.
  s32: P_cs/P_a = 1,105 (kvadratura shodná na 0,0003 proti toleranci 0,105), ale MĚŘIDLO SELHALO, protože se
  kontrolní řešení z vakua nedotáhlo do konvergence ve 30 iteracích. Čtení „konvergence“ (předběžné): přebytky
  0,206 / 0,124 / 0,105 pro σ = 8 / 16 / 32, Aitkenovo dno 0,099 > 1/16 → **dno**. Hloubka jádra Π(0) se od spojité
  hodnoty 0,7099 s rostoucím σ vzdaluje (0,7113 / 0,7202 / 0,7256). Gaussův deficit v R/2: 0,098 / 0,095 / 0,135.
  Čtení podle předem dané formulace: pokud se dno potvrdí, diskrétní operátor v tomto režimu (2+1, rovnice (a),
  samokonzistentně, silné pole Π(0) ≈ 0,71) nekonverguje k předpokládanému spojitému; stěna v oblasti pohybu
  vyloučena; příčina otevřená (okrajové podmínky, normalizace zdroje, škálování rodiny, nebo operátor). Rozdíl
  nezávislý na σ v sobě podobné rodině ukazuje spíš na jinou spojitou rovnici než na lokální opravu řádu ℓ_k/σ.
  Dál (deklarováno 24. 10. před běhy): s32 znovu s max_iter = 60 (certifikace) a slabopolní rodina w08/w16/w32
  (q desetkrát menší, Π(0) ≈ 0,97): je dno vlastností lineárního operátoru, nebo nelineární samokonzistentní vazby?
- **Rozklad „dna“ K.12: nelokální stínění operátoru je vyloučené, dno je součet složek přístroje a otevřeného
  nelineárního zbytku** (analýza a lineární model, 24. 10., `docs/derivations/k12_dno.md`, skripty
  `docs/derivations/k12_diag/`; statusy a čtení K.12 beze změny, R132). Přesný symbol očekávaného BD operátoru
  2+1 na statickém poli: B(k) = −k²(1 − 1,155 ℓ_k k + …), M₀ = M₂ = 0, r(0) = 1. Vzdálené pole vidí celý náboj
  a operátor u zdroje antistíní ∝ ℓ_k/σ (bez sítě a stěny: tok v 1,25σ 1,159 / 1,077 / 1,038 pro σ = 8 / 16 / 32).
  Lineární model schématu K.12 (běhy autora, jemná síť 1 ℓ_k): sklon u zdroje 1,161 / 1,077 / 1,037, tedy
  operátor na 0,002 a stěna oblast pohybu nezasahuje. Na sobě podobné síti K.12 (h = σ/4) přibývá +0,007 /
  +0,021 / +0,051: chyba sítě roste se σ, odezva operátoru klesá, proto se přebytek zastaví (Aitken 0,099 je
  přechod, ne asymptota). Gaussův deficit v R/2 ≈ 0,10 nezávisle na σ i na jemné síti, i bez PL báze; mizí
  s lokální odezvou geometrie. Původ: nelokální operátor vidí lom pole nulově pokračovaného za stěnou (ve 2+1
  nese celý tok), referencí je přitom lokální Dirichletova úloha. Daleké pole nelokálních variant není
  zkonvergované v kroku sítě u stěny (σ = 8: tok 0,9R 0,494 → 0,562 při kroku 1 → 0,5). Stejný mechanismus
  vysvětluje hloubku jádra, která se se σ vzdaluje (model 0,978 / 0,924 / 0,895, K.12 0,995 / 0,965 / 0,946).
  Definice zdroje je v pořádku (varianta s lokální geometrií konverguje: u(0) 1,0006 při σ = 32). Otevřený
  zbytek: sklon v 1,25σ v K.12 je o 0,077 / 0,060 / 0,055 nad lineárním modelem na téže síti (silné pole,
  nelinearita). Výhrada: nelokální odezvu geometrie model zastupuje plochým operátorem na splajnu. Závěr: K.12
  není falsifikace spojité limity Path A, ani její potvrzení. Rozhodne silné a slabé pole s krokem u zdroje
  ≤ 1 ℓ_k.
- **K.12 h16: krok u zdroje 1 ℓ_k odstranil složku sítě, jak rozklad předpověděl; pole stále neprošlo** (kalibrace,
  2+1, σ = 16, q = 0,0025, R = 200, 76 uzlů, 25. 10., NEPROŠEL). P_cs/P_a = 1,111 (s16 se sítí σ/4: 1,124; tolerance
  1,0625). Sklon v 1,25σ 1,141 proti předpovědi 1,137 (operátor +0,077 a nelineární zbytek +0,060 z rozboru).
  Čtení 'grid' (deklarováno 24. 10.) = **podpořeno**: rozdíl sítí +0,018 proti lineárnímu modelu +0,021. Kvadratura
  shodná (|L2 − L1| = 0,0000), řešení zkonvergovalo (6 iterací). Gaussův deficit v R/2 0,077 (s16 0,095), Π(0) 0,7193.
  Po odečtení sítě zbývá při σ = 16 přebytek nad tolerancí. Tvoří ho operátor (∝ ℓ_k/σ) a nelineární zbytek,
  o kterém rozhodnou h32 a slabopolní hw16/hw32.
- **K.12 uzavřen 25. 10. bez h32 a slabopolní rodiny** (rozhodnutí po rozboru, běh h32 přerušil pád počítače).
  Pro teorii z něj platí: diskrétní operátor konverguje ke spojitému a nestíní (analyticky, r(0) = 1), dno bylo z velké
  části chybou přístroje (síť, stěna) a předpověď rozkladu se potvrdila na h16. Otevřený zůstává nelineární zbytek
  ~0,06 ve sklonu v 1,25σ v silném poli 2+1 se starou rovnicí (a). Statusy částí beze změny. Poučení pro další pole:
  krok u zdroje pevný v ℓ_k, stěna daleko nebo vnější pokračování pole, reference se symbolem operátoru.
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

- **F3.8 search: retardované pole s dosahem svazek poprvé stahuje, vázaný stav ale nevznikl: ε = 0,1 se pomalu
  hroutí, ε = 0,5 se zhroutí a odrazí** (mechanismus, průzkum, 32 sprinklingů 2+1, ~49 200 událostí, ZAZNAMENÁNO;
  kandidát podle pravidla z hlavičky: žádný). Pole φ = Σ F/(2πτ) přes tok v minulosti (τ >= 1), ln Π = −λφ,
  λ = ε/φ_ref, φ_ref = 44,15. Sklon šířky v W1 [12, 28] a W2 [28, 44] (tolerance vázanosti s_free/4 ≈ 0,010):
  - free: +0,042 ± 0,005 a +0,040 ± 0,006 (roste);
  - ε = 0,01: +0,037 a +0,027 (roste, o trochu pomaleji);
  - **ε = 0,1: −0,004 ± 0,005 (nerozlišeno, |s| + 2 SE = 0,014 těsně nad tolerancí) a −0,036 ± 0,005 (kolaps)**;
    šířka se zastaví na ~3,35 a pak klesá k 2,7;
  - ε = 0,5: −0,042 ± 0,009 (kolaps, R 3,0 → 2,1 do t ≈ 23), pak +0,278 ± 0,017 (odraz a rozpínání ke 7,6).

  Kontrola soustavy prošla ve všech větvích. Je to první mechanismus v laboratoři, u kterého vlastní pole svazku
  přemůže difuzi diskrétnosti a svazek přitahuje k sobě. Průběh odpovídá studenému prachu bez tlaku: v gravitaci
  se hroutí (ε = 0,1) a bez tlumení se po zhroucení rozletí (ε = 0,5). Rovnovážný (vázaný, stacionární) stav
  potřebuje něco, co vzdoruje stlačení: tlak, rozptyl rychlostí, nebo konečnou kapacitu (saturaci) teorie. To je nová
  hypotéza po neúspěchu, s vlastními kritérii. Výhrada předem: ve 2+1 je potenciál logaritmický, fyzikální čtení
  vyžaduje 3+1.
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

- **Horizontová cesta k a₀ (Deser–Levin, Jacobson, Milgrom) přímé odvození pozorovaného a₀ nedává** (odvozeno
  24. 10., revize autora, `docs/derivations/a0_horizont.md`). Důsledně podle Jacobsona (κ pomocného horizontu) se
  teplota vykrátí a vzniká obecná relativita; modifikace vyžaduje dodatečně ztotožnit κ s fyzikálním g. Plná
  Deserova–Levinova teplota jako G_eff dává hluboký režim MOND neexistující (pro g_N < cH bez řešení). Rozdíl teplot
  (Milgrom 1999) jako G_eff dává g_N = cH_dS(√(1 + (g/cH_dS)²) − 1): hluboký MOND s a₀ = 2cH_dS (9–11× víc než
  pozorované) a konstantní radiální člen ~cH_dS ≈ 5–7·10⁻¹⁰ m/s² i ve vysokém zrychlení, silně v napětí s planetárními
  rezidui (Iorio: 10⁻¹⁵ Mars, 10⁻¹⁴ m/s² Saturn; Carloni a kol.: perihelia vylučují univerzální zrychlení velikosti
  Pioneerovy anomálie, radarový limit bez výhrad 3·10⁻⁹). Faktor 1/(2π) z cesty nevzniká. Nevyvrací to Milgromovu
  modifikovanou setrvačnost ani všechny horizontové mechanismy. K.14-dS jen jako samostatný geometrický test.
- **Princip fokusace (každý světelný směr k se fokusuje úměrně toku energie hmoty T_kk) vynucuje ve slabém
  statickém poli γ = 1 a Newtonův zákon, bez skalární rovnice pole a bez ansatzu H1** (odvozeno 23. 10., sympy,
  `docs/derivations/fokusace_gamma.py`). S nezávislými potenciály φ (hodiny) a ψ (prostor) je
  R_kk = ∇²φ + ∇²ψ + nⁱnʲ∂ᵢ∂ⱼ(ψ − φ); stejná fokusace ve všech směrech pro hmotu v klidu dává ψ = φ a
  2∇²φ = 8πGρ. Pro pohybující se hmotu vyžaduje tenzorové pole Σ F u u (gravitomagnetismus). Prostorová část
  metriky je v kauzální množině tvar kuželů, tedy deformace řádu kolem hmoty; na pevném sprinklingu chybí.
  Konzistence principu vyžaduje, aby hmota zachovávala 4-hybnost: odvozená oprava pravidla sloučení F = |P|.
  Hybnost pole nese deformace řádu, která se šíří. Plán: `docs/derivations/prostorova_cast.md`.
- **Jacobsonova podmínka (Clausius ⇔ R_kk = 8πG T_kk) v lineárním řádu vyžaduje, aby pole neslo 4-rychlost
  zdrojů: H1 jako skalární Π v pevné soustavě ji pro pohybující se hmotu porušuje, H1 zapsané s vlastní
  4-rychlostí zdroje ji splňuje přesně** (odvozeno 23. 10., sympy, `docs/derivations/jacobson_linear_check.py`;
  slabé pole h_ab = −2U(η_ab + 2 n_a n_b), prach rychlostí v = 0,6, čtyři světelné směry: rozdíly −8,2 / +32,9 /
  +2,4 / −6,8 proti přesně 0). Kovariantní tvar je linearizovaná obecná relativita pro prach. Pro mnoho zdrojů
  je pole tenzorové, h̄_ab ∝ Σ F u_a u_b G_ret, s gravitačními vlnami s polarizacemi + a × a s nesouměrnou
  deformací kuželů. Hybnost hmoty, kterou dodal F3.6, je tedy nutnou podmínkou Jacobsonovy cesty (krok 3)
  i gravitačních vln. Platí v lineárním řádu a ve 3+1; ve 2+1 obecná relativita místní gravitaci nemá.
  Podrobnosti: `docs/derivations/jacobson_cesta.md`, oddíl 6.
- **Status kotvy a₀ = cH₀/2π je v textu teorie rozporný a mechanismus komprese stojí na vyloučeném K₀**
  (rozbor F6.7 na papíře, 23. 10., `docs/derivations/a0_kotva.md`). Text dává a₀ roli konvence vykazování
  (nemění dynamiku) i fyzikálního prahu přechodu D_eff 3 → 2. Testovat lze jen práh: F6.7 pak znamená spočítat
  a*_jádro = 1/(2π) v jednotkách cH₀ z dynamiky. BTFR v⁴ = GMa₀ vyjde jen tehdy, když se dimenze mění při
  pevném zrychlení (při pevné délce vyjde v² ∝ M). Plato D_eff ≈ 3 je v textu dané saturací u K₀, které testy
  vyloučily; komprese potřebuje nový základ (F6.4, nebo efektivní dimenze toku). Varianta a₀(z) ∝ H(z), ke které
  vede motivace „kauzální hloubky“, je v napětí s rotačními křivkami při z ≈ 2 pro každý mechanismus s přechodem
  při a₀(z); konstantní varianta (H₀ nebo Λ) s daty souhlasí. Rozhodnutí je na autorovi.
- **Perihelium vybírá f = ln Π: v rovnici B_g f(Π) = κρ_s s H1 (γ = 1) je parametr PPN β = (1 − f″(1))/2, pro
  rodinu (Π^k − 1)/k tedy β = 1 − k/2** (odvozeno 23. 10. po F1.7 N, sympy: vnější vakuové řešení
  Π = 1 − U + cU², U = M/r, dosazené do radiální rovnice do druhého řádu). Výsledky: ln Π: c = 1/2, β = 1
  (Π = e^(−U) do druhého řádu, g₀₀ = 1 − 2U + 2U², jako obecná relativita); Π − 1: β = 1/2; (Π⁻¹ − 1)/(−1): β = 3/2;
  (Π² − 1)/2: β = 0; (Π³ − 1)/3: β = −1/2. Stáčení perihelia je úměrné (2 + 2γ − β)/3; měření (planetární
  efemeridy) dávají |β − 1| řádu 10⁻⁴, tedy |k| ≲ 2·10⁻⁴ a f″(1) = −1. Důsledky:
  - lineární kontrola Π − 1 (a s ní stažený návrh B_g Π − Π·B_g 1, který se s ní shoduje v lineárním řádu) by
    dala o 1/6 větší stáčení perihelia než obecná relativita (β = 1/2) a je pozorováním vyloučená;
  - ln Π už není jen motiv (skládání poměrů hodin): je to jediný člen rodiny slučitelný s periheliem. Pozorování
    ale fixuje jen f″(1); jiné f s f(1) = 0, f′(1) = 1, f″(1) = −1 dají v prvních dvou řádech totéž;
  - ln Π v N nemá nejširší rozsah řešení (k = −1 má širší), ale k = −1 má β = 3/2.
  Platí za předpokladu H1 (γ = 1 vstupem) a spojité limity operátoru; numerická kontrola z vnějších profilů N
  (fit 1 − U + cU²) je levná a zatím neudělaná.
  Autor 23. 10. přijal ln Π jako kanonickou skalarizaci (docs/rozhodnuti_autora.md, bod 1, čtvrté kolo).
- **Pravidlo ohybu řady F3 (F3.8–F3.10) je pohyb po geodetikách konformní geometrie Π²η (PPN γ = −1), světlo
  tedy neohýbá** (odvozeno 23. 10. po F3.10 universality). Řetězec ohýbá nesený směr k −∇ ln Π rychlostí |∇ ln Π|
  na jednotku vlastního času. Pro řetězec s rychlostí v kolmo ke gradientu statického pole: kopnutí jedné vazby
  Δu_x = |∇ ln Π| τ, laboratorní čas vazby γτ, takže zrychlení je |∇ ln Π| / γ² = |∇ ln Π| (1 − v²). Přesně to dává
  Eulerova–Lagrangeova rovnice pro L = −Π √(1 − v²), tj. geodetiky metriky Π²η. Obecná relativita dává
  (1 + v²), pro světlo dvojnásobek newtonovského ohybu; Π²η dává pro světlo nulu (nulové geodetiky konformní
  geometrie jsou ploché). Důsledky:
  - univerzalita volného pádu (F3.10 universality) je u geodetického pravidla automatická, nezávislý test
    principu ekvivalence vyžaduje tělesa s jiným složením hmoty (vazebná energie, Nordtvedtův jev), tedy vázané
    struktury, které zatím nemáme (F3.7–F3.9);
  - dynamika F3 realizuje právě tu geometrii, kterou výše uvedený důsledek rovnice (a) vylučuje pozorováním
    (Cassini γ = 1). Ohyb světla vyžaduje i prostorovou část metriky (H1: optický čas L − (1 + γ)εU), kterou
    pravidlo F3.8 nepoužívá. K.13a to nerozhodne: ověřuje, že diskrétní BD operátor na geometrii H1 (γ = 1)
    konverguje ke spojité limitě; s ověřeným přístrojem se teprve rozhodne, jaká rovnice pole je s γ = 1
    slučitelná (rozpor H1 a (a) ve 3+1 výše). K.13a PROŠEL 23. 10.: překážkou ohybu světla není substrát, ale rovnice (a).
  Ověřitelné levně: rychlá sonda kolmo k poli, ohyb jako funkce v (pravidlo 1 − v², GR 1 + v²).
- **Rovnice pole (a) a hypotéza H1 nemohou ve 3+1 platit současně** (odvozeno 19. 10. z lokální limity
  □_g − R/2, kterou K.11 numericky potvrdil jen pro konformní geometrii Π²η; pro jiné geometrie se opírá o
  Belenchia et al. 2016). Pro statickou metriku −Π²dt² + Π^(−2γ)dx² (γ je parametr PPN; Π²η má γ = −1,
  H1 má γ = 1) dává symbolický výpočet lokální tvar (□_g − R/2)Π:
  2+1: γ = −1: 3Π⁻²∇²Π; γ = 0: 2∇²Π + Π⁻¹|∇Π|²; γ = 1: Π²∇²Π + 2Π|∇Π|²;
  3+1: γ = −1: 4Π⁻²∇²Π + 2Π⁻³|∇Π|²; γ = 0: 2∇²Π + Π⁻¹|∇Π|²; γ = 1: 2Π|∇Π|² (bez členu s ∇²Π).
  Lineární vazba je κ/c(γ), c(γ) = 2 − γ(d_s − 1), d_s počet prostorových rozměrů (D = d_s + 1). Ve 3+1 tedy (a) v geometrii
  Π²η (vazba κ/D, ve 2+1 předběžně F1.5) dává γ = −1, tj. žádný ohyb světla ani Shapirovo zpoždění (Cassini: γ = 1);
  v geometrii H1 dává rovnici 2Π|∇Π|² = κρ_s, takže mimo zdroj je Π konstantní a newtonovské pole nevzniká.
  Ve 2+1 se H1 s (a) nevylučuje (vazba κ). Důsledek: F7.1 ve smyslu „H1 při rovnici (a)“ nemá ve 3+1 řešení;
  buď se změní tvar rovnice pole (na co operátor působí, nebo zda nese −R/2), nebo H1. Rozhodnutí je na autorovi.
  Odvození vazby κ/D → κ = 16π (kapitola 4.7) platí jen pro konformní geometrii.
  Potvrzeno vlastním přístrojem (K.13a, 23. 10.): v geometrii H1 ve 3+1 jde koeficient u ∇²u v B_g 1 k −1, takže
  v B_g Π = B_0[εu] + ε δB_u[1] se lineární členy ε∇²u − ε∇²u ruší (c(1) = 0). Rozpor tedy už nestojí jen na
  Belenchia et al. 2016; platí v prvním řádu a na σ ≥ 4 ℓ_k s odchylkou pod 2,4 % při σ = 32.
- **Odezva BD operátoru nerozliší argument f rovnice B_g f(Π) = ρ_s** (odvozeno 22. 10., vede ke zrušení plánovaného
  K.13b). Operátor je lineární v poli. Pro f(1) = 0 je E_g[f(Π)] = ε f′(1) B_0[u] + O(ε²), takže geometrie (H1 i jiná)
  se v lineárním řádu neprojeví. Rozdíl ln Π − (Π − 1) ve druhém řádu je −½ ε² B_0[u²], tedy plochý operátor na u²
  (∇²(u²) až v kontinuu); geometrie vstupuje až v ε³. Kandidáty f rozliší jen princip
  (docs/derivations/ln_pi_potencial.md: aditivita poměrů hodin, c² ln N je přesný potenciál zrychlení statických
  pozorovatelů) a celá nelineární rovnice (F1.7 N).
- **První řád operátoru v rodině −Π²dt² + Π^(−2γ)dx² je lineární v γ** (sympy C6, C7, 22. 10.):
  δ(□_g − R/2)w = 2γ u∇²w + (1 − γ(d_s − 2)) ∇u·∇w + (1 − γ(d_s − 1)) w∇²u = (1 − γd_s) P + (1 + γ) Q.
  P = (−2u∇²w + (d_s − 1)∇u·∇w + d_s w∇²u)/D je odezva na změnu hustoty √−g při plochých kuželech (konformní
  přístroj K.11). Q = (2u∇²w + 2∇u·∇w + w∇²u)/D je odezva na změnu rychlosti kuželů Π^(1+γ). Diskrétní přístroj
  K.13a počítá obě části zvlášť.

- **Text teorie a laboratoř stojí na dvou různých substrátech; testy volbu fakticky učinily** (audit axiomů 22. 10.,
  `docs/audit_axiomu.md`, opraveno po připomínce autora). *Physics of Causality* předpokládá konečný stupeň K₀,
  lokalitu a dimenzi z grafové vzdálenosti, čas jako počet aktualizací a Lorentzův faktor z hodin v klidu vůči
  substrátu (klidová soustava). Laboratoř používá Lorentzovsky invariantní sprinkling s nekonečně mnoha vazbami na
  událost. Podle věty Bombelliho, Hensona a Sorkina se konečný stupeň a invariance substrátu vylučují. Testy už
  rozhodly ve prospěch invariantního substrátu: K₀ neplatí (F0.3), kapacita omezuje tok, ne vazby (F1.8), lokální
  růst typu grafu varietu nevytvořil (F6.1, F6.2). Nevysloveno zůstává: text teorie nese zbytky substrátu s klidovou
  soustavou a pravidla laboratoře také (zúžené kužele v souřadnicích, dělení toku podle počtu potomků; srv. F0.6).
- **Funkcionál odporu R = Σ (I/Π) Δτ zatím není funkcionál pole a variační princip sám nic nevybírá** (audit 24. 10.,
  `docs/derivations/variacni_audit.md`). Při pevném grafu a I nezávislém na Π je R v u = ln Π přísně monotónní
  (δR = −Σ I e^(−u) η Δτ ≠ 0), takže nemá vnitřní stacionární bod. Pokud Δτ = Π Δt, Π z R úplně vypadne
  (R = Σ I Δt). Otevřená rozhodnutí: A1 (co je Δτ), A2 (operační definice I). Každá jednosložková rovnice
  A∇²L + B|∇L|² = f je po integračním faktoru μ′/μ = (2B − A′)/A variační, takže existence akce nevybere ani β,
  ani ln Π. Statická druhá variace rovnice (a) je kladná (∫ |∇η|² + (2κ/3) ρ Π η²) a o stabilitě struktury nic
  neříká. Test hypotézy D_eff ≈ 3 patří do fáze C (stabilita kauzální množiny vůči diskrétní akci typu
  Benincasa–Dowker, předem daný protokol). Literární souvislosti (Kleitman–Rothschild, fázové chování) jsou
  neověřená rešerše.

- **β v diskrétním poli: ve druhém řádu ho určují jen dvě už kalibrované části a konvence zdroje ho neovlivní**
  (odvozeno 25. 10., bez testu; rozhodnutí o F1.9 ponechal autor na laboratoři). Pro L = ln Π = εL₁ + ε²L₂
  rovnice B_g L = κρ_s dává B₀[L₁] = κρ_s a B₀[L₂] = −dB_{L₁}[L₁] + (zdrojový člen ∝ ρ_s). β se čte mimo zdroj, kde
  je ρ_s = 0. Jakákoli konvence zdroje ve druhém řádu (vlastní nebo souřadnicová hustota, faktor √−g = e^{(1−3γ)L})
  proto mění jen aktivní hmotu (harmonická odezva 1/r), ne β. Mimo zdroj platí ∇²L₁ = 0, z odezvy prvního řádu
  δ(□_g − R/2)w zbývá (1 − γ(d_s − 2))∇u·∇w = (1 − γ)|∇L₁|² (d_s = 3). Pro S = c|∇L₁|² je L₂ = −(c/2)L₁² a
  β = 1 − c/2, pro H1 (γ = 1) tedy c = 0 a β = 1. Kandidáti f = (Π^k − 1)/k přidávají k|∇L₁|² (β = 1 − k/2, shodně
  s F1.7). Diskrétní β ve druhém řádu tak závisí jen na odezvě prvního řádu v H1 (K.13a: konverguje ke spojitému tvaru,
  0,8–2,4 % při σ = 32) a na plochém operátoru (K.4, F1.3). **Neověřeno:** K.13a měřil s gaussovskými poli u zdroje
  (r ≤ 2σ), ne s polem 1/r mimo zdroj, kde se β čte. Nelokální ocasy tam ve 3+1 klesají rychle, přímé měření
  c(r) ale chybí. Možný test (F1.9a): c(r) = dB_{L₁}[L₁]/|∇L₁|² v řádcích mimo zdroj, předpověď 0, kontrola stejnou
  kvadraturou γ = −1 s předpovědí 2. Konvence zdroje: laboratoř dál používá ρ_s jako pevnou souřadnicovou funkci bez
  √−g (F1.7, K.12). Pro definici hmoty (κ ↔ G) se musí stanovit samostatně.

## Otevřené otázky

- Přijme autor výslovně Lorentzovsky invariantní substrát, ke kterému testy vedou (F0.3, F1.8, F6.1, F6.2), a přepíše
  podle něj axiomy (lokalita, dimenze, čas, odvození Lorentzova faktoru)? (audit axiomů, oddíl 0)

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
