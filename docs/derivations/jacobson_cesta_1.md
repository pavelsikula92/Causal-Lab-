# Jacobsonova cesta: gravitace jako stavová rovnice kauzální struktury

Sepsáno 22. 10. 2026 na pokyn autora („tak to zkusme“). **Status: plán a rozbor předpokladů, nic z toho zatím
není výsledek laboratoře.** Měřidlo prvního kroku je K.14.

## 1. Jacobsonův argument (1995)

T. Jacobson, *Thermodynamics of Spacetime: The Einstein Equation of State*, PRL 75, 1260 (1995). Jednotky
c = k_B = 1.

Předpoklady:

1. **Entropie horizontu je úměrná ploše.** Každý lokální kauzální horizont (Rindlerův horizont zrychleného
   pozorovatele kolem libovolného bodu, v libovolném směru) nese entropii S = η·A s **univerzální** konstantou η.
2. **Clausiova rovnost** δQ = T δS pro teplo, které horizontem proteče, v lokální rovnováze.
3. **Unruhova teplota** T = ħκ/2π, kde κ je zrychlení pozorovatele. Pochází z kvantové teorie pole.
4. **Geometrie nulových svazků** (Raychaudhuriho rovnice). Plocha řezu horizontu se mění podle toho, jak ji
   zakřivení fokusuje: δA = −∫ λ R_ab k^a k^b dλ dA.
5. **Tok energie** přes horizont: δQ = ∫ λ T_ab k^a k^b dλ dA (T_ab je tenzor energie a hybnosti hmoty).

Výsledek: požadavek 2 pro všechny body a všechny nulové směry k vyžaduje R_ab k^a k^b = (2π/ħη) T_ab k^a k^b.
Odtud plynou **Einsteinovy rovnice**

    R_ab − ½ R g_ab + Λ g_ab = 8πG T_ab,     G = 1/(4ħη),

s kosmologickou konstantou Λ jako integrační konstantou (není určena).

Einsteinovy rovnice tu nejsou základní zákon. Jsou to stavové rovnice, podobně jako rovnice ideálního plynu
vyplývá ze statistiky molekul, které jednotlivě neznáme.

## 2. Co z toho pro teorii plyne

Kdyby předpoklady platily pro substrát kauzální fyziky, byla by gravitace **odvozená**, ne vložená:

- žádná rovnice pole, kterou bychom volili (B_g f(Π), volba f, ln Π proti Π − 1);
- žádný ansatz H1: γ = 1, Shapirovo zpoždění, stáčení perihelia i gravitační vlny by vyplynuly z obecné relativity;
- obě poloviny gravitace (hodiny i prostor) by vznikly současně, takže Nordströmův problém (γ = −1 pro skalár)
  nenastane.

**Most k SI (rozhodnutí 8).** Je-li entropie počet molekul horizontu krát entropie s₀ jedné molekuly, pak
η = s₀ a^(4)/ℓ² a podmínka G = 1/(4ħη) dává

    ℓ² = 4 s₀ a^(4) ℓ_P² = (2√3/5) s₀ ℓ_P²,    ℓ ≈ 1,075 √s₀ ℓ_P

(a^(4) = √3/10, ℓ_P Planckova délka). Délka diskrétnosti už není volný parametr, ale je svázaná s G až na s₀.
Předpovědí se to stane teprve tehdy, až se ℓ nebo s₀ změří nezávisle (např. ze Sorkinova odhadu fluktuací Λ ∝ 1/√V,
nebo z nelokálnosti operátoru). Do té doby je to vztah konzistence.

**Souvislost s pomalými hodinami.** Ve statickém gravitačním poli platí Tolmanův–Ehrenfestův vztah: lokální
teplota T_lok = T_∞/N, N = √(−g₀₀). V metrice laboratoře je N = Π. Kde běží hodiny pomaleji, je lokálně tepleji.
Výchozí intuice teorie (gravitace = zpomalení hodin) má v termodynamickém čtení přesný význam: Π je faktor
rudého posuvu teploty.

## 3. Které předpoklady laboratoř umí dodat a které ne

| Předpoklad | Co v kauzální množině | Stav |
|---|---|---|
| 1. S ∝ A, univerzální η | počet molekul horizontu (Dou–Sorkin; 1-molekula Barton et al. 2019) | **K.14** (plochý a konformně plochý prostor); dál zakřivené kužele |
| 2. Clausiova rovnost | požadavek na dynamiku řádu v lokální rovnováze | vstup (princip), zatím netestovatelný: laboratoř nemá dynamiku |
| 3. Unruhova teplota | kvantové pole na kauzální množině (Sorkin–Johnston) | vstup z kvantové teorie pole |
| 4. Raychaudhuri | platí ve spojité limitě sprinklingu | převzato |
| 5. T_ab hmoty | tok laboratoře **a hybnost**, kterou laboratoř nemá (audit, oddíl 3) | **chybí**: bez hybnosti nelze hmotu vložit do rovnice |

Nejcitlivější body jsou dva. Univerzalita η (předpoklad 1) je to, co může selhat numericky: kdyby hustota molekul
na jednotku vlastní plochy závisela na zakřivení už v hlavním řádu, argument padá. Předpoklad 5 je to, co chybí
koncepčně: hmota teorie musí mít tenzor energie a hybnosti, tedy i hybnost.

## 4. Kroky

1. **K.14:** počet molekul = a^(d) × plocha v plochém prostoru (2, 3, 4 rozměry, dvě hustoty) a úměrnost
   **vlastní** ploše v konformně plochém prostoru (alternativy: souřadnicová plocha, objemová hustota).
2. **Zakřivené kužele:** horizont v geometrii, kde se mění i kužele, ne jen hustota (slabé pole Schwarzschilda /
   H1). Potřebuje sprinkling s nekonformní kauzální relací. Pro začátek stačí první řád, jako v K.13a.
3. **Hmota s hybností:** definovat tok energie přes horizont z vlastní hmoty teorie (tok + hybnost). Tady se
   Jacobsonova cesta potkává s řadou F3.
4. **Clausius na kauzální množině:** v geometrii s tokem hmoty přes horizont porovnat změnu počtu molekul
   s δQ/T. Pro geometrii, která už řeší Einsteinovy rovnice, je to při platnosti kroku 1 automatické. Netriviální
   je to až s hmotou z kroku 3 a s dynamikou řádu.

## 5. Výhrady

- Jacobsonův argument předpokládá lokální rovnováhu. Mimo ni dává opravy (viskozita, vyšší členy křivosti;
  Eling, Guedens, Jacobson 2006).
- Ztotožnění entropie s počtem molekul je hypotéza. Jinou možností je entropie provázanosti kvantového pole na
  kauzální množině (Sorkin, Yazdi), která potřebuje useknutí spektra.
- Argument říká, jak musí gravitace vypadat v rovnováze. Neříká, jaká mikroskopická dynamika Clausiovu rovnost
  vynucuje. Růst řádu (audit, oddíl 5) zůstává otevřený.
- Λ zůstává neurčená. To je ale v souladu s kauzálními množinami, kde se Λ objevuje jako fluktuace (Sorkin).
