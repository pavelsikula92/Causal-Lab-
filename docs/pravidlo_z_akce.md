# Pravidlo substrátu z akce: odvození místo vložení

Sepsáno 28. 10. 2026 na pokyn autora. Je to analytický rozbor; výpočty byly ověřeny v sympy. Navazuje na
`tenzorova_vazba.md`, F7.4a a rozhodnutí 13.

## 1. Co se má odvodit

Pravidlo, které F7.4a ověřil na kauzální množině, má tři části:

- **(a) vlastní soustava:** každý tok zužuje kužel ve své klidové soustavě;
- **(b) sčítání:** deficity kapacity od různých toků se sčítají;
- **(c) poměr:** zúžení kuželu a nasycení (přidání událostí) jsou v poměru, který dává γ = 1 (princip propustnosti).

V F7.4a bylo pravidlo vstupem. Níže je ukázáno, že (a) a (b) plynou z axiomu A6 a z linearity bez dalšího
předpokladu o řádu. (c) je jediné číslo a určí ho akce samotného řádu.

Konvence: statické slabé pole ds² = −(1 + 2Φ) dt² + (1 − 2Ψ) dx², d prostorových rozměrů, G = c = 1. Newtonovsky
Φ = −U. Platí γ = Ψ/Φ.

## 2. Krok 1: strana hmoty (A6) dává vlastní soustavu a sčítání

A6: akce toku je vlastní čas podél jeho světočáry, S_m = −Σ m τ.

- Na kauzální množině se vlastní čas čte počtem: délkou nejdelšího řetězce nebo objemem diamantu (K.16 hodiny).
- Pro malou poruchu g = η + h je δτ = −½ h_ab u^a u^b τ, kde u je 4-rychlost toku.
- Variace akce hmoty podle h je tedy ½ Σ m u_a u_b δ(světočára) = ½ T_ab.

Plyne z toho:

- **(a)** Řád „vidí“ tok jen přes u_a u_b. Vlastní soustava toku je proto soustava, ve které se měří jeho vlastní čas,
  a není potřeba ji předpokládat.
- **(b)** Akce je součet přes toky, takže v lineárním řádu se příspěvky sčítají.

Pro tuto část není potřeba nic vědět o dynamice řádu. Stačí A6 a to, že řád odpovídá lineárně.

## 3. Krok 2: Lorentzova invariance nechává jediné číslo

Lineární, lokální a Lorentzovsky invariantní odezva h_ab = G ∗ T má dva kanály: T_ab a η_ab T. Pro statický prach
to vede na jednoparametrickou rodinu (`tenzorova_vazba.md`):

    h̄_ab = U [ (1 − γ) η_ab + 2(1 + γ) u_a u_b ]

Struktura pravidla F7.4a je tím vynucená až na jedno číslo γ:

- γ = −1: čistá stopová vazba (Nordström, pevný řád F3);
- γ = 1: čistá vazba na T_ab.

## 4. Krok 3: γ určuje akce řádu

Obecná statická kvadratická akce řádu s vazbou na prach je

    L = c_ΦΦ |∇Φ|² + c_ΦΨ ∇Φ·∇Ψ + c_ΨΨ |∇Ψ|² − κ' ρ Φ

(vazba jen přes Φ, protože prach v klidu má jen T_00). Variace podle Ψ dává c_ΦΨ ΔΦ + 2 c_ΨΨ ΔΨ = 0, tedy

    γ = Ψ/Φ = −c_ΦΨ / (2 c_ΨΨ)

nezávisle na c_ΦΦ. Variace podle Φ určí jen amplitudu, tedy vztah κ ↔ G.

**Einsteinova–Hilbertova akce** (sympy, √−g R do druhého řádu, d = 2, 3, 4):

    c_ΦΦ = 0,   c_ΦΨ = −2(d − 1),   c_ΨΨ = (d − 1)(d − 2)   ⇒   γ = 1/(d − 2)

Ve 3+1 je γ = 1. Ve 2+1 je c_ΨΨ = 0 a Φ nemá zdroj, takže tam není newtonovské přitahování.

**Princip propustnosti** vyžaduje γ = 2/(d − 1), tedy c_ΨΨ / c_ΦΨ = −(d − 1)/4. S Einsteinovou–Hilbertovou akcí se to
shoduje jen pro d = 3.

**Fyzikální čtení:**

- **Hodiny nemají vlastní tuhost (c_ΦΦ = 0).** Časový potenciál Φ vystupuje jen jako vazba a čas se určí z prostoru.
- **Propustnost je důsledek, ne postulát.** Prach nemá napětí a řád má Einsteinovu–Hilbertovu tuhost, proto se
  prostorová trace-reversed část nedeformuje (h̄_ij = 0). Ve 3+1 to je totéž co ρ_c·c = konst.
- **Jen hustota jako stupeň volnosti dává Nordströma.** Kdyby řád měl jako dynamický stupeň volnosti jen hustotu
  (konformní faktor), vyšel by γ = −1, tedy přesně pevný řád laboratoře.

## 5. Kandidát akce řádu: Benincasova–Dowkerova(–Glaserova) akce

- **Co je známé.** Očekávaná hodnota BD(G) akce na sprinklinzích jde v kontinuální limitě k Einsteinově–Hilbertově
  akci plus hraničním členům. Benincasa a Dowker to ukázali pro BD akci (2010). Pro malé kauzální diamanty to
  prověřuje arXiv:2007.13192: d = 2 je odvozeno explicitně, d = 3, 4, 5 souhlasí do prvního řádu v křivosti.
- **Laboratoř stejný operátor už používá:** B_g z K.11 a K.13a.
- **Výsledek řetězce:**

      A6 (vlastní čas toku počítaný řádem) + BD akce řádu  ⇒  pravidlo F7.4a v lineárním řádu (γ = 1)

- **Co z toho plyne:** pravidlo substrátu je lineární odezva řádu s BD tuhostí na akci toků. Část (c), tedy
  propustnost, je důsledek toho, že BD akce má Einsteinův–Hilbertův limit.

## 6. Co je nové a co ne

- **Známé (učebnice):** v kontinuu je řetězec A6 + Einstein–Hilbert ⇒ linearizovaná obecná relativita. Je známé
  i to, že BD akce má Einsteinův–Hilbertův limit.
- **Nové pro teorii:** pravidlo laboratoře a princip propustnosti nejsou samostatné postuláty. Jsou to lineární
  odezva řádu s BD tuhostí na akci A6. Z postulátů se tak stává odvození, ale v kontinuální limitě.
- **Nové a publikovatelné by bylo:** ukázat, že diskrétní BD akce na kauzální množině sama vybere γ = 1, tedy změřit
  c_ΦΦ, c_ΦΨ, c_ΨΨ přímo z řádu. Nejostřejší předpověď je **c_ΦΦ = 0**: statická porucha samotných hodin nemá
  kvadratickou akci. Nezávisí na normalizaci a pevný řád ani skalární akce ji nesplní.
- **Hlavní překážka je šum.**
  - Signál je řádu a² ∫|∇Φ|² d³x · T ~ 0,1–1 pro a = 0,1 v boxu laboratoře (ℓ = 1).
  - Fluktuace lokální BD akce rostou s počtem prvků, odhadem √N · O(1–10) ~ 10²–10³ pro N ~ 5·10⁴.
  - Přímo je to tedy o několik řádů mimo. Pomoci může rozmazaná (nelokální) BD akce a korelované vzorkování (stejné
    body pro různá Φ, Ψ). Je to odhad řádu velikosti a musí ho potvrdit test proveditelnosti.
- **Mimo dosah laboratoře:** skutečná dynamika (dráhový integrál s BD akcí a nemanifoldní fáze) je otevřený problém
  teorie kauzálních množin. Řetězec v oddíle 5 předpokládá varietní fázi, jako F6.4b.

## 7. Návrhy (nic neběží)

1. **K.17, test proveditelnosti:** rozptyl rozmazané BD akce na sprinklingu 3+1 jako funkce N a ε a rozptyl rozdílu
   při korelovaném vzorkování dvou blízkých metrik. Rozhodne, zda je měření c_ΦΦ, c_ΦΨ, c_ΨΨ vůbec možné.
2. **F6.8 (jen pokud K.17 projde):** kvadratická forma diskrétní BD akce ve statickém sektoru (Φ, Ψ). Předpovědi:
   c_ΦΦ = 0 a γ = −c_ΦΨ/(2c_ΨΨ) = 1 ve 3+1. Alternativy: skalární akce hustoty γ = −1, princip propustnosti mimo 3+1.
3. **Rešerše novosti:** vazba toků přes délku řetězce (A6) a odezva BD akce na hmotu ve slabém poli. Zatím nalezeno
   jen Sverdlov–Bombelli 2009 (skalární pole, ne částice) a práce o limitě BD(G) akce.
