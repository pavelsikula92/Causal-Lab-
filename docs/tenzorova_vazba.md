# Co pole nese: skalár nasycení, nebo tenzor toku

Sepsáno 27. 10. 2026 (fundamentální výzkum, na pokyn autora). Je to analytický rozbor, ne výsledek laboratoře.
Výpočty byly přepočítány v kontejneru (sympy). Navazuje na `gravitacni_rovnice.md`, `proc_gamma1.md`
a `jacobson_cesta.md` (oddíl 6).

Konvence: signatura (−, +, +, +), c = 1, g = η + h. Statické slabé pole má h_00 = 2U a h_ij = 2γU δ_ij, U > 0
je newtonovský potenciál (U = −L). Trace-reversed pole je h̄_ab = h_ab − ½ η_ab h. u je 4-rychlost zdroje.

## 1. Otázka

Rovnice B_g ln Π = κρ_s je skalární. Má jedno pole Π, jeden zdroj ρ_s a metriku H1 zapsanou v jedné soustavě. Je to
úplný dynamický zákon gravitace, nebo jen statický řez něčeho většího? Rozhodují o tom tři pozorování, která skalár
v pevné soustavě nemůže splnit:

- **Strhávání soustav** (g_0i od pohybující se nebo rotující hmoty). Gravity Probe B naměřila −37,2 ± 7,2 mas/rok
  proti −39,2 v obecné relativitě. Skalár Π v pevné soustavě má g_0i = 0. Laboratoř to už ukázala lineárně:
  `jacobson_linear_check.py`, pro pohybující se zdroj podmínka neplatí.
- **Úbytek periody dvojitých pulsarů.** Hulseův–Taylorův pulsar souhlasí s kvadrupólovým vzorcem obecné
  relativity na ~0,2 %, dvojitý pulsar J0737−3039 na ~0,01 %.
- **Polarizace gravitačních vln.** Měření se třemi detektory (GW170814 a další) silně preferuje čistě tenzorové
  polarizace před čistě skalárními.

## 2. Co by dal skalár jako úplný zákon

Skalární pole normované tak, aby staticky dalo Newtona (ΔL = 4πGρ), vyzařuje kvadrupólově výkonem

    P_skalár = (G/60) [ (tr I⃛)² + 2 I⃛_ij I⃛_ij ]

I_ij = ∫ ρ x_i x_j. Obecná relativita má P = (G/5) Q⃛_ij Q⃛_ij s bezstopým Q. Pro kruhovou dvojhvězdu vychází
P_skalár / P_OTR = **1/6** (sympy, přesně).

K tomu přistupuje monopólové a dipólové vyzařování podle toho, jak ρ_s závisí na rychlosti a na vazebné energii
těles. Pro kompaktní tělesa (neutronové hvězdy) bývá dipól ještě silnější než kvadrupól.

Výhrada: rovnice laboratoře nemá stanovenou akci pole, takže normalizace toku energie je převzatá ze standardního
skalárního pole. Faktor řádově 1/6 se ale normalizací, která dává Newtona, opravit nedá. Opravit ho může jen
jiný tvar pole.

**Závěr:** skalární rovnice jako úplný dynamický zákon je pozorováním vyloučená. Platí jako statický řez v klidové
soustavě zdroje.

## 3. Hlavní výsledek: γ je poměr tenzorové a stopové vazby

Statické slabé pole s parametrem γ v klidové soustavě zdroje, d_s prostorových rozměrů:

    h̄_00 = U (1 + γ d_s),     h̄_ij = U (1 + 2γ − γ d_s) δ_ij

Ve 3+1 (d_s = 3):

    h̄_ab = U [ (1 − γ) η_ab + 2(1 + γ) u_a u_b ]

- **γ = −1 (pevný řád, A6 na sprinklingu):** h̄_ab = 2U η_ab. Pole se váže jen na stopu tenzoru energie a hybnosti
  (spin 0, Nordströmova teorie). Světlo se neohýbá. Světlo samo negravituje (elektromagnetické pole má nulovou stopu).
- **γ = 1 (princip propustnosti):** h̄_ab = 4U u_a u_b. Pole se váže na celý T_ab prachu bez stopové části. To je
  přesně linearizovaná obecná relativita (□h̄_ab = −16πG T_ab, T = ρ u u).
- **γ = 1/3 (jen užší kužele):** směs obou.

Princip propustnosti √−g · c = 1 dává γ = 2/(d_s − 1). Podmínka čisté tenzorové vazby (h̄_ij = 0 pro prach) dává
γ = 1/(d_s − 2). Obojí se shoduje **jen pro d_s = 3**.

Tři důsledky:

1. **Výsledek F7.1c má fyzikální obsah, který přesahuje čtení γ.** Substrát, který drží propustnost, se ve 3+1
   v lineárním řádu váže na zdroj stejně jako obecná relativita: na T_ab, ne na jeho stopu.
   - Princip propustnosti je totéž co h̄_ij = 0: deformace řádu leží jen podél proudu (u_a u_b) a napříč proudu
     je nulová.
2. **F7.1b (2+1) testuje princip propustnosti, ne obecnou relativitu.** Obecná relativita ve 2+1 nemá newtonovské
   přitahování (jacobson_cesta, oddíl 6). Že se princip a obecná relativita potkají právě ve 3+1, je buď náhoda,
   nebo náznak, proč má svět tři prostorové rozměry. Odvozené to není.
3. **Otázka strany zdroje je v lineárním řádu zodpovězená.** Zdrojem není skalár ρ_s, ale T_ab toku. Z toho plyne:
   - světlo gravituje;
   - tlak gravituje jako ρ + 3p (aktivní hmota horké a tlakové hmoty);
   - pohyb zdroje vytváří g_0i (strhávání).

   Nasycené události jsou v tomto čtení T_ab u^a u^b v klidové soustavě toku. Lokalizovaný tok nese svou 4-rychlost
   jako stav (F3.6), K.15 měří T^ab a pravidlo sloučení F = |P| zachovává 4-hybnost přesně. Zachování ∂_a T^ab = 0
   je přesně podmínka, bez které tenzorová vazba není konzistentní.

## 4. Tenzorové čtení v lineárním řádu: co předpovídá

Pole je součet přes toky, h̄_ab(x) = Σ 4 U_s u_a u_b (retardovaně), tedy zobecnění F3.8 Σ F/(2πτ) na tenzor.
Předpovědi jsou pak předpovědi linearizované obecné relativity:

- strhávání soustav s hodnotou obecné relativity;
- kvadrupólový vzorec vyzařování přesně (při ∂T = 0);
- dvě příčné polarizace + a × šířící se rychlostí c;
- ohyb světla s γ = 1 a gravitace světla i tlaku.

Kniha (36.2) říká totéž slovy: v rovnovážném režimu (plateau) dynamika impedance reprodukuje tenzorové vlnové
průběhy obecné relativity.

## 5. Nelineární řád: rozcestí

Takzvaný bootstrap (Gupta 1954, Kraichnan 1955, Feynman, Weinberg 1965, Deser 1970) tvrdí: nehmotné pole spinu 2,
které se konzistentně váže na celkový zachovávaný T_ab včetně energie pole samotného, vede jednoznačně na obecnou
relativitu. Předpoklady jsou Lorentzova invariance, lokalita a nejvýše dvě derivace. Ve 4D to podporuje
i Lovelockova věta.

Obecná relativita přitom princip propustnosti v přesném tvaru A·B = 1 nesplňuje. Schwarzschild má
A·B = (1 − M²/4r²)², takže odchylka je ve druhém řádu. Z toho plyne rozcestí:

**A. Propustnost je zákon plateau (lineární řád) a nelineární doplnění je obecná relativita.**
- Odpovídá to Jacobsonově cestě (rozhodnutí 9) i knize 36.2.
- Skalární rovnice B_g ln Π = κρ_s zůstává správná do 1PN (β = 1 souhlasí).
- Její silnopolní a 2PN předpovědi (hrdlo 2,341 M, žádný horizont, stín +0,91 %, g_rr s 2x²) by ale nebyly
  předpověďmi teorie. Byly by to artefakty prodloužení slabopolního principu za jeho platnost.
- Teorie by pak měla odvodit obecnou relativitu (Clausius na kauzální množině), ne ji volit.

**B. Propustnost platí přesně i nelineárně.**
- Teorie se pak od obecné relativity liší od 2PN dál a musí porušovat některý předpoklad bootstrapu: nelokalitu
  BD operátoru na škále ℓ, nebo preferovanou soustavu.
- Zároveň musí ukázat, že vyzařování a strhávání přesto souhlasí s pozorováním.
- Je to výrazně riskantnější cesta a dnes pro ni neexistuje mechanismus.

Doporučení laboratoře: **A**. Stojí na hotových výsledcích (K.14, K.15, F3.6, lineární Jacobsonova kontrola,
F7.1c), nekoliduje s větami o jednoznačnosti a souhlasí s knihou. Rozhodnutí je na autorovi a do té doby platí
silnopolní oddíl `gravitacni_rovnice.md` jen podmíněně.

## 6. Další fundamentální kroky (návrhy, nic neběží)

1. **Diskrétní test tenzorové vazby (F7.2): gravitují protiběžné proudy podle T_ab?**
   - Dva protiběžné svazky (konfigurace K.15) s rychlostí ±v mají v klidové soustavě součtu aktivní hmotu
     T_00 + T_ii.
   - Tři kandidáti dávají různé předpovědi:
     - tenzorová vazba: 2ργ²(1 + v²);
     - klidová hustota nasycených událostí (Lorentzův skalár): 2ρ;
     - hustota v laboratorní soustavě: 2ργ.
   - Měřidlem je K.16c (φ z hodin).
   - Předpoklad: substrátové pravidlo musí rozhodnout, jak se skládají nasycení s různými u. Právě tím se kandidáti
     liší a to je ta fundamentální otázka.
2. **Jacobsonova cesta, krok 4 nelineárně:** Clausiova rovnost na kauzální množině s tokem přes horizont. Rozhoduje
   mezi A a B: pokud substrát Clausiovu rovnost drží, vychází obecná relativita a varianta A.
3. **Proč substrát drží propustnost:** ve tvaru h̄_ij = 0 je to otázka, proč nasycení deformuje řád jen podél proudu.
   Kandidátem je zachování toku (axiom IV): tok bez příčného napětí (prach) nepředává hybnost napříč.
