# Gravitační demonstrace: gravitační pole, pohyb těles a měřitelné důsledky
# ===========================================================================
# Samostatná interaktivní aplikace; nemění kód gravitační laboratoře.
# Hmotnosti, geometrie a okrajové podmínky jsou v této demonstraci vstupy.
#
# Spuštění:
#   python demo/gravitace_demo.py                 # otevře výchozí scénu
#   python demo/gravitace_demo.py --selftest      # spustí numerické autotesty
#   python demo/gravitace_demo.py --scene 2       # otevře scénu s daným indexem
#   python demo/gravitace_demo.py --scene 2 --autorun  # spustí simulaci ihned
# Požadavky: Python 3, tkinter a numpy.
#
# FYZIKÁLNÍ MODEL A JEHO PŘEDPOKLADY
# ----------------------------------
# Pole hodin L = ln Π vně bodových zdrojů je zadáno vztahem
#   B_g ln Π = κρ_s,    L(x) = −Σ_j M_j / sqrt(|x − x_j|² + ε²).
# Parametr ε numericky změkčuje singularitu. Pole se lineárně superponuje;
# nelineární překryvové korekce zdrojů nejsou zahrnuty.
# Hypotéza metriky H1: ds² = −e^(2L)dt² + e^(−2L)dx², v 1PN řádu β = γ = 1.
# Exponenciální tvar metriky mimo 1PN je volbou demonstrace, ne odvozeným výsledkem.
# Kvazistatický pohyb používá dv/dt = −e^(4L)∇L − |v|²∇L + 4v(v·∇L) a
# dτ = sqrt(e^(2L) − e^(−2L)|v|²)dt. Retardace pole se zanedbává.
# Shoda slabého pole s Newtonovou fyzikou nebo GR sama o sobě není novým důkazem.
#
# PŘEHLED SCÉN
# ------------
# 0  Precese perihelia: excentrická oběžná dráha a dlouhodobé měření precese.
# 1  Sluneční soustava (zveličená): několik planetárních drah v jednotkách G=c=1.
# 2  Ohyb světla: porovnání trajektorií fotonů s různými parametry dopadu.
# 3  Shapirovo zpoždění: porovnání zpoždění dvou paprsků v různých vzdálenostech.
# 4  Dvojhvězda a planeta: oběh dvou hvězd a vzdálenější planety.
# 5  Hodiny v poli: vlastní čas statických hodin a obíhající částice.
# 6  Sluneční soustava (skutečná): vstupní elementy J2000, vzdálenosti v AU.
# 7  Povrch Země: pád jablka, vržený míč a hodiny v různých výškách.
# 8  Země, ISS a GPS: oběžné dráhy a rozdíly rychlosti hodin.
# 9  UAP – naprogramovaný manévr: potřebné zrychlení při konvenčním pohonu.
# 10 UAP – ruční řízení: ovládání klávesnicí; „nulové pole“ je hypotéza.
# 11 Černá díra: vnější profil rovnice a kritické trajektorie, bez modelovaného horizontu.
# 12 Alfa Centauri A/B: model relativního oběhu dvojhvězdy.
# 13 Sirius A/B: dvojhvězda a gravitační rudý posuv bílého trpaslíka Siria B.
# 14 Okolí Slunce do 12 ly: prostorová mapa hvězd a extrapolace známých pohybů.
# 15 Volná scéna: pracovní plocha s pevným zdrojem, na kterou lze přidávat objekty.
#
# OVLÁDÁNÍ A ZOBRAZENÍ
# --------------------
# V ovládání lze přidávat/mazat tělesa, fotony a hodiny. Kliknutí umístí objekt;
# tažení určuje směr a rychlost. Ve scénách u Země tažení představuje vrh.
# Hmotnost se zadává jako poměr m/Mref; Mref je nejhmotnější zdroj při načtení scény a zůstává pevná.
# Tempo se zadává jednotně jako simulační časové jednotky za sekundu reálného času;
# v SI scénách se zobrazuje čas v sekundách, v geometrických scénách modelové jednotky.
# 3D náhled lze zapnout u rovinných scén. Jde pouze o projekci roviny z=0,
# nikoli o prostorovou dynamiku. Hvězdná mapa má skutečné souřadnice x,y,z.
# Ve hvězdné mapě lze objekty přidávat formulářem nebo nástrojem na pracovní ploše.
#
# SROVNÁNÍ A MĚŘENÍ
# -----------------
# Srovnávají se H1, obecná relativita (GR), Newtonova korpuskulární aproximace,
# pravidlo F3 a PPN metrika s parametry β a γ. Sledují se precese perihelia,
# ohyb světla, Shapirovo zpoždění a vlastní čas hodin a částic.
# Druhý řád ohybu v exponenciální metrice je 4π(M/b)²; GR předpovídá
# 15π(M/b)²/4. Rozdíl plyne z volby metriky mimo 1PN řád.
#
# AUTOTESTY (--selftest)
# ---------------------
# Test 1: Keplerova perioda v Newtonovské limitě.
# Test 2: precese perihelia proti vztahu 6πM/[a(1−e²)].
# Test 3: odchylka světla proti předpovědi prvního a druhého řádu.
# Test 4: Shapirovo zpoždění proti analytickému vztahu.
# Test 5: nulová podmínka fotonu |v|=e^(2L) a rychlost statických hodin.
# Test 6: doba pádu z 10 m proti sqrt(2h/g).
# Test 7: rozdíl rychlosti hodin GPS a pozemních hodin za den.
# Test 8: rozdíl rychlosti hodin ve výšce 30 m.
# Test 9: precese Merkuru pro Newton, H1, GR a PPN s β=1/2.
# Test 10: zachycení fotonu u černé díry a kritický parametr dopadu H1.
# Autotesty kontrolují implementaci a použité vzorce, nikoli platnost celé teorie.
# Karta „Stav teorie“ rozlišuje testované výsledky, spojité odvození,
# vstupní předpoklady a neověřené hypotézy. Laboratorní panel zobrazuje uložené záznamy.
# ===========================================================================

from __future__ import annotations

import math
import sys
import time
from collections import deque

import numpy as np

EPS = 1e-3  # Numerické změkčení singularity bodových zdrojů v délkových jednotkách.
C_SI = 299792458.0  # Rychlost světla [m/s]; v SI scénách je čas interně vyjádřen jako ct.
G_C2 = (
    7.42591e-28  # G/c² [m/kg]: převod kilogramů na geometrickou délkovou jednotku.
)
M_EARTH = 4.43502e-3  # G M_Earth / c^2 [m]
R_EARTH = 6.371e6  # střední poloměr Země [m]
G_SURF = (
    M_EARTH / R_EARTH**2 * C_SI**2
)  # Gravitační zrychlení na povrchu [m/s²].
AU = 1.495978707e11  # Astronomická jednotka [m].
YEAR = 365.25 * 86400.0  # Juliánský rok [s].
M_SUN = 1476.625  # GM☉/c² [m].
ARCSEC = math.pi / (180.0 * 3600.0)  # Úhlová vteřina v radiánech.
LY = 9.4607304725808e15  # Světelný rok [m].
M_SUN_KG = 1.98847e30  # Hmotnost Slunce [kg].
R_SUN = 6.957e8  # Jmenovitý poloměr Slunce [m].
TAU_NULL = (
    0.05 * C_SI
)  # Hypotéza „nulového pole“: doba odezvy řízení rychlosti 0,05 s, interně v metrech.
L_MAX = (
    0.2  # Mez slabého pole; pro |L| > L_MAX není aproximace 1PN dostačující.
)
# Chyba je řádu L²; tělesa se při vstupu sloučí a fotony jsou pohlceny.


def fmt_len(m):
    """Převede délku v metrech na vhodnou zobrazovací jednotku."""
    a = abs(m)
    if a >= 0.1 * LY:
        return f"{m / LY:.4g} ly"
    if a >= 0.01 * AU:
        return f"{m / AU:.4g} AU"
    if a >= 1000.0:
        return f"{m / 1000.0:.4g} km"
    return f"{m:.4g} m"


def fmt_time(sec):
    """Převede čas v sekundách na sekundy, minuty, hodiny, dny nebo roky."""
    a = abs(sec)
    if a < 120:
        return f"{sec:.2f} s"
    if a < 7200:
        return f"{sec / 60:.2f} min"
    if a < 2 * 86400:
        return f"{sec / 3600:.2f} h"
    if a < 2 * YEAR:
        return f"{sec / 86400:.2f} dní"
    return f"{sec / YEAR:.3f} let"


# =============================================================================
# Fyzika: pole, geodetické zrychlení a integrátor
# =============================================================================


class Body:
    """Těleso scény. Typ „mass“ je hmotný zdroj, „photon“ foton bez hmotnosti, „clock“ pevné hodiny a „craft“
    těleso s předepsaným pohybem. U režimu craft představuje vlastní zrychlení požadovaný tah pohonu."""

    def __init__(
        self,
        name,
        x,
        y,
        vx=0.0,
        vy=0.0,
        mass=0.0,
        kind="mass",
        fixed=False,
        color="#ffffff",
        radius=0.0,
    ):
        self.name, self.kind, self.color = name, kind, color
        self.radius = float(
            radius
        )  # Kladný poloměr označuje pevný povrch; ostatní tělesa na něj dosednou.
        self.z = 0.0  # Třetí souřadnice, využívá ji skutečná 3D mapa hvězd.
        self.thrust = (0.0, 0.0)  # Řízené zrychlení [m/s²] nastavené klávesnicí.
        self.info = ""  # Doplňující popis pro hvězdnou mapu a dvojhvězdy.
        self.landed_t = None  # Simulační čas dosednutí ve scénách u Země.
        self.script = None  # Předepsaný pohyb ve scéně UAP; viz CraftScript.
        self.max_g = 0.0  # 'craft': largest proper acceleration so far [m/s^2]
        self.now_g = 0.0
        self.mass = float(mass) if kind == "mass" else 0.0
        self.fixed = bool(fixed) or kind == "clock"
        self.x, self.y, self.vx, self.vy = float(x), float(y), float(vx), float(vy)
        self.tau = 0.0  # Nasčítaný vlastní čas.
        self.tau_def = 0.0  # Nasčítaný rozdíl t − τ; zachovává přesnost i pro velmi slabá pole.
        self.t_birth = None  # Čas vložení tělesa do scény, nastavuje aplikace.
        self.trail = deque(maxlen=4000)
        self.alive = True


_BH = {}


def bh_profile():
    """
    Vnější vakuové řešení rovnice laboratoře B_g ln Π = 0 pro H1, 3+1 rozměry a kulově symetrický zdroj M = 1:
        e^(2L) [(1 − L) ΔL + L |∇L|²] = 0  (tvar F1.7).
    Pro s = ln r platí dL/ds = P a dP/ds = −P − LP²/(1 − L), kde P = r dL/dr.
    Řešení se integruje metodou RK4 z r = 1e5 do r = 1e−4 s krokem ds = 1e−3. Asymptotický začátek je
    L = −1/r + 1/(6r³); druhý člen pouze předefinuje hmotnost. Výsledek (s, L, P, info) obsahuje poloměr
    fotonové sféry (r dL/dr = 1/2) a kritický parametr dopadu b = r e^(−2L). Data se ukládají do mezipaměti.
    """
    if _BH:
        return _BH["s"], _BH["L"], _BH["P"], _BH["info"]

    def rhs(L, P):
        return P, -P - L * P * P / (1.0 - L)

    s0, s1, h = math.log(1e5), math.log(1e-4), -1e-3
    r0 = 1e5
    L, P = -1.0 / r0 + 1.0 / (6 * r0**3), 1.0 / r0 - 1.0 / (2 * r0**3)
    ss, LL, PP = [s0], [L], [P]
    s_ = s0
    while s_ > s1:
        k1 = rhs(L, P)
        k2 = rhs(L + 0.5 * h * k1[0], P + 0.5 * h * k1[1])
        k3 = rhs(L + 0.5 * h * k2[0], P + 0.5 * h * k2[1])
        k4 = rhs(L + h * k3[0], P + h * k3[1])
        L += h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        P += h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        s_ += h
        ss.append(s_), LL.append(L), PP.append(P)
    s_arr, L_arr, P_arr = (
        np.array(ss[::-1]),
        np.array(LL[::-1]),
        np.array(PP[::-1]),
    )  # increasing s
    # Fotonová sféra: nejvzdálenější r, pro které P = 1/2
    k = np.where((P_arr[:-1] - 0.5) * (P_arr[1:] - 0.5) <= 0)[0][-1]
    f = (0.5 - P_arr[k]) / (P_arr[k + 1] - P_arr[k])
    s_ph = s_arr[k] + f * (s_arr[k + 1] - s_arr[k])
    L_ph = L_arr[k] + f * (L_arr[k + 1] - L_arr[k])
    r_ph = math.exp(s_ph)
    info = {
        "r_ph": r_ph,
        "R_ph": r_ph * math.exp(-L_ph),
        "b_crit": r_ph * math.exp(-2 * L_ph),
        "b_crit_gr": 3 * math.sqrt(3),
        "R_ph_gr": 3.0,
    }
    _BH.update(s=s_arr, L=L_arr, P=P_arr, info=info)
    return s_arr, L_arr, P_arr, info


def bh_L(r, M=1.0):
    """Vrátí L a dL/dr v zadaných poloměrech pro hmotnost M; pod 1e−4 M použije logaritmický chvost."""
    s_arr, L_arr, P_arr, _ = bh_profile()
    r = np.maximum(np.asarray(r, float), 1e-300)
    sx = np.log(r / M)
    L = np.interp(sx, s_arr, L_arr)
    P = np.interp(sx, s_arr, P_arr)
    lo = sx < s_arr[0]
    L = np.where(lo, L_arr[0] + P_arr[0] * (sx - s_arr[0]), L)
    hi = sx > s_arr[-1]
    L = np.where(hi, -M / r, L)
    P = np.where(hi, M / r, P)
    return L, P / r


class Law:
    """
    Pohybový zákon srovnávacího běhu. U = −L = Σ M/r je Newtonův potenciál pro G = c = 1.
    Typ „metric“ používá geodetiky statické metriky ds² = −A(U)dt² + B(U)dx² a rovnici
    dv/dt = −(A′/2B)∇U − v(v·∇U)B′/B + |v|²(B′/2B)∇U + v(v·∇U)A′/A.
    Světlo má rychlost |v| = sqrt(A/B), hodiny dτ/dt = sqrt(A − Bv²). Hodnota A−1 se počítá zvlášť
    kvůli zachování rozdílů hodin řádu 10⁻¹⁵. Typ „newton“ používá dv/dt = ∇U pro tělesa i fotony a
    nezahrnuje gravitační zpomalování hodin.
    """

    def __init__(
        self,
        key,
        name,
        color,
        kind="metric",
        A=None,
        dA=None,
        B=None,
        dB=None,
        Am1=None,
    ):
        self.key, self.name, self.color, self.kind = key, name, color, kind
        self.A, self.dA, self.B, self.dB, self.Am1 = A, dA, B, dB, Am1


def _law_h1():
    """Metrika H1 používaná v demonstraci i laboratoři: L = ln Π, A = e^(2L) = e^(−2U), B = e^(−2L) = e^(2U)."""
    return Law(
        "H1",
        "kauzální (H1, ln Π)",
        "#ffffff",
        A=lambda U: np.exp(-2 * U),
        dA=lambda U: -2 * np.exp(-2 * U),
        B=lambda U: np.exp(2 * U),
        dB=lambda U: 2 * np.exp(2 * U),
        Am1=lambda U: np.expm1(-2 * U),
    )


def _law_gr():
    """Obecná relativita: Schwarzschildova metrika v izotropních souřadnicích, A = ((1−U/2)/(1+U/2))² a B = (1+U/2)⁴.
    Pro jeden zdroj je přesná; u více zdrojů se potenciály sčítají, což je aproximace správná do 1PN řádu."""
    return Law(
        "GR",
        "obecná relativita (Schwarzschild)",
        "#ff7ad9",
        A=lambda U: ((1 - U / 2) / (1 + U / 2)) ** 2,
        dA=lambda U: -2 * (1 - U / 2) / (1 + U / 2) ** 3,
        B=lambda U: (1 + U / 2) ** 4,
        dB=lambda U: 2 * (1 + U / 2) ** 3,
        Am1=lambda U: -2 * U / (1 + U / 2) ** 2,
    )


def _law_newton():
    return Law("N", "Newton (světlo jako korpuskule)", "#7aff8a", kind="newton")


def _law_f3():
    """Transportní pravidlo F3 laboratoře odpovídá geodetikám Π²η s γ = −1: A = B = e^(−2U)."""
    return Law(
        "F3",
        "pravidlo F3 (Π²η, γ = −1)",
        "#ffb37a",
        A=lambda U: np.exp(-2 * U),
        dA=lambda U: -2 * np.exp(-2 * U),
        B=lambda U: np.exp(-2 * U),
        dB=lambda U: -2 * np.exp(-2 * U),
        Am1=lambda U: np.expm1(-2 * U),
    )


def _law_ppn(beta, gamma):
    """Parametrizovaná postnewtonovská metrika do prvního řádu: A = 1 − 2U + 2βU² a B = 1 + 2γU."""
    return Law(
        "PPN",
        f"PPN β = {beta:g}, γ = {gamma:g}",
        "#7ad0ff",
        A=lambda U: 1 - 2 * U + 2 * beta * U * U,
        dA=lambda U: -2 + 4 * beta * U,
        B=lambda U: 1 + 2 * gamma * U,
        dB=lambda U: 2 * gamma + 0 * U,
        Am1=lambda U: -2 * U + 2 * beta * U * U,
    )


class Sim:
    """Systém těles a jeho integrátor bez grafiky, sdílený aplikací a autotesty."""

    def __init__(self, law=None):
        self.law = law if law is not None else _law_h1()
        self.profile = None  # Pevný zdroj, jehož pole je přesným vnějším řešením.
        self.no_guard = (
            False  # True vypne slučování a pohlcování v silném poli černé díry.
        )
        self.bodies = []
        self.t = 0.0
        self.eta = 0.02  # Bezpečnostní faktor časového kroku.
        self.dt_max = 0.5
        self.events = []  # Záznamy o sloučení a pohlcení v silném poli.
        self._extra = None

    def remove(self, body):
        """Odstraní těleso. Odstranění zdroje profilu (černé díry) odstraní také jeho přesné pole a obnoví ochranu silného pole."""
        body.alive = False
        if body in self.bodies:
            self.bodies.remove(body)
        if body is self.profile or getattr(body, "is_profile", False):
            self.profile = None
            self.no_guard = False

    # ---- pole aktuálního stavu
    def arrays(self):
        live = [b for b in self.bodies if b.alive]
        pos = np.array([[b.x, b.y] for b in live], float).reshape(-1, 2)
        vel = np.array([[b.vx, b.vy] for b in live], float).reshape(-1, 2)
        # Pole zdroje s profilem přidává ext_field, proto se nezahrnuje do součtu bodových hmot.
        mass = np.array([0.0 if b is self.profile else b.mass for b in live], float)
        photon = np.array([b.kind == "photon" for b in live], bool)
        fixed = np.array(
            [b.fixed or b.kind == "craft" for b in live], bool
        )  # Řízené objekty se pohybují podle svého scénáře.
        return live, pos, vel, mass, photon, fixed

    def ext_field(self, points):
        """Pole zdroje s profilem: přesné vnější řešení rovnice pole pro H1. Srovnávací běh GR používá běžný bodový potenciál U."""
        if self.profile is None:
            return 0.0, 0.0
        c = self.profile
        d = points - np.array([c.x, c.y])
        r = np.sqrt(np.sum(d * d, axis=1))
        L, dL = bh_L(r, c.mass)
        with np.errstate(over="ignore", invalid="ignore"):
            g = np.nan_to_num(
                (dL / np.maximum(r, 1e-300))[:, None] * d
            )  # Ve středu zdroje je příspěvek nulový; zdroj je pevný.
        return L, g

    def L_points(self, pts):
        """Hodnota L v bodech od všech zdrojů, tedy bodových hmot i zdroje s profilem."""
        live, pos, vel, mass, photon, fixed = self.arrays()
        L = np.zeros(len(pts))
        if np.any(mass > 0):
            L, _ = self.field_at(pts, pos, mass)
        Le, _ = self.ext_field(pts)
        return L + Le

    def extra_acc(self, live):
        """Vnější zrychlení řízeného tělesa v SI převede z m/s² na interní jednotku 1/délka dělením c²."""
        return np.array(
            [
                (
                    [b.thrust[0] / C_SI**2, b.thrust[1] / C_SI**2]
                    if b.kind == "pilot"
                    else [0.0, 0.0]
                )
                for b in live
            ],
            float,
        ).reshape(-1, 2)

    @staticmethod
    def field_at(points, src_pos, src_mass, exclude=None):
        """Vypočítá L a ∇L v bodech tvaru (P, 2). Parametr exclude[p] určuje zdroj vynechaný pro bod p."""
        d = points[:, None, :] - src_pos[None, :, :]  # (P, S, 2)
        s2 = np.sum(d * d, axis=2) + EPS * EPS
        s = np.sqrt(s2)
        w = np.broadcast_to(src_mass[None, :], s.shape).copy()
        if exclude is not None:
            for p, j in enumerate(exclude):
                if j is not None and j >= 0:
                    w[p, j] = 0.0  # Vlastní pole se nezapočítává.
        L = -np.sum(w / s, axis=1)
        g = np.sum((w / (s2 * s))[:, :, None] * d, axis=1)  # grad of -M/s is +M d / s^3
        return L, g

    @staticmethod
    def geodesic_accel(vel, L, gL):
        """Zrychlení geodetiky H1: dv/dt = −e^(4L)∇L − |v|²∇L + 4v(v·∇L)."""
        v2 = np.sum(vel * vel, axis=1)[:, None]
        vg = np.sum(vel * gL, axis=1)[:, None]
        return -np.exp(4.0 * L)[:, None] * gL - v2 * gL + 4.0 * vel * vg

    def law_accel(self, vel, L, gL):
        """Vypočítá zrychlení podle zvoleného zákona; H1 používá uzavřený vztah, ostatní obecnou statickou metriku."""
        law = self.law
        if law.key == "H1":
            return self.geodesic_accel(vel, L, gL)
        U, gU = -L, -gL
        if law.kind == "newton":
            return gU.copy()
        A, dA, B, dB = law.A(U), law.dA(U), law.B(U), law.dB(U)
        v2 = np.sum(vel * vel, axis=1)
        vg = np.sum(vel * gU, axis=1)
        return (-dA / (2 * B) + v2 * dB / (2 * B))[:, None] * gU + (
            (dA / A - dB / B) * vg
        )[:, None] * vel

    def deriv(self, pos, vel, mass, fixed):
        """Časové derivace společného stavu; polohy zdrojů se mění v každé fázi integrátoru RK4."""
        n = len(mass)
        L, gL = self.field_at(pos, pos, mass, exclude=list(range(n)))
        if self.profile is not None:
            Le, ge = self.ext_field(pos)
            L, gL = L + Le, gL + ge
        acc = self.law_accel(vel, L, gL)
        if self._extra is not None:
            acc = acc + self._extra
        for i, vt in getattr(self, "_vt", []):
            acc[i] = (vt - vel[i]) / TAU_NULL
        acc[fixed] = 0.0
        dpos = vel.copy()
        dpos[fixed] = 0.0
        return dpos, acc

    def choose_dt(self, pos, vel, mass, photon):
        """
        Adaptive step: a small fraction of the shortest orbital/crossing time of the bodies that MOVE (static clocks,
        landed and fixed bodies do not limit it). With nothing moving the field is static and one long step is exact:
        returns inf (the caller takes the rest of the frame), or dt_max if a scripted craft is present. Near a solid
        surface the step is also limited by the time to impact, so a fall ends in ~one RK4 step on the surface.
        """
        live = [b for b in self.bodies if b.alive]
        moving = np.array([not (b.fixed or b.kind == "craft") for b in live], bool)
        if not np.any(moving):
            return self.dt_max if any(b.kind == "craft" for b in live) else math.inf
        speed = np.sqrt(np.sum(vel * vel, axis=1)) + 1e-12
        if self.profile is not None:  # Časové škály v okolí zdroje s profilem.
            c = self.profile
            d = np.sqrt((pos[:, 0] - c.x) ** 2 + (pos[:, 1] - c.y) ** 2) + 1e-12
            tb = np.minimum(d / speed, np.where(photon, np.inf, np.sqrt(d**3 / c.mass)))
            dt = float(np.clip(self.eta * np.min(tb[moving]), 1e-6, self.dt_max))
        else:
            src = mass > 0
            if not np.any(src) or len(mass) < 2:
                dt = self.dt_max
            else:
                d = pos[:, None, :] - pos[None, src, :]
                sdist = np.sqrt(np.sum(d * d, axis=2) + EPS * EPS)
                idx = np.where(src)[0]
                for k, j in enumerate(idx):
                    sdist[j, k] = np.inf  # Zdroj neomezuje vlastní časový krok.
                t_cross = np.min(sdist, axis=1) / speed
                t_orb = np.min(np.sqrt(sdist**3 / (mass[src][None, :] + 1e-30)), axis=1)
                t_body = np.where(photon, t_cross, np.minimum(t_cross, t_orb))
                dt = float(
                    np.clip(self.eta * np.min(t_body[moving]), 1e-6, self.dt_max)
                )
        # Odhad času do dopadu na pevný povrch; počítá se jen v blízkosti povrchu.
        for o in (b for b in live if b.radius > 0):
            for i, b in enumerate(live):
                if not moving[i] or photon[i] or b.radius > 0:
                    continue
                rx, ry = pos[i, 0] - o.x, pos[i, 1] - o.y
                r = math.hypot(rx, ry)
                h = r - o.radius
                if r <= 0.0 or h <= 0.0 or h > 0.05 * o.radius:
                    continue
                vr = (vel[i, 0] * rx + vel[i, 1] * ry) / r
                vt2 = max(float(vel[i] @ vel[i]) - vr * vr, 0.0)
                g = max(
                    o.mass / (r * r) - vt2 / r, 0.0
                )  # Gravitace po odečtení dostředivého příspěvku při oběhu.
                u = max(-vr, 0.0)  # Radiální rychlost směrem dovnitř.
                if g > 0:
                    t_imp = (-u + math.sqrt(u * u + 2.0 * g * h)) / g
                else:
                    t_imp = h / u if u > 0 else math.inf
                dt = min(
                    dt, t_imp * (1.0 + 1e-7) + 1e-9
                )  # Krok skončí těsně za povrchem, aby se zaznamenalo dosednutí.
        return dt

    def step(self, dt=None):
        """Provede jeden krok RK4 pro živá tělesa a poté upraví rychlost fotonů na nulovou podmínku e^(2L)."""
        live, pos, vel, mass, photon, fixed = self.arrays()
        if not live:
            return 0.0
        self._extra = (
            self.extra_acc(live) if any(b.kind == "pilot" for b in live) else None
        )
        # Hypotéza „nulového pole“: rychlost řízeného tělesa se během TAU_NULL přiblíží požadované hodnotě bez ohledu na gravitaci.
        self._vt = [
            (i, np.array(b.vtarget, float))
            for i, b in enumerate(live)
            if b.kind == "pilot" and getattr(b, "vtarget", None) is not None
        ]
        if dt is None:
            dt = self.choose_dt(pos, vel, mass, photon)
            if not math.isfinite(dt):
                dt = self.dt_max
        k1p, k1v = self.deriv(pos, vel, mass, fixed)
        k2p, k2v = self.deriv(pos + 0.5 * dt * k1p, vel + 0.5 * dt * k1v, mass, fixed)
        k3p, k3v = self.deriv(pos + 0.5 * dt * k2p, vel + 0.5 * dt * k2v, mass, fixed)
        k4p, k4v = self.deriv(pos + dt * k3p, vel + dt * k3v, mass, fixed)
        pos_n = pos + dt / 6.0 * (k1p + 2 * k2p + 2 * k3p + k4p)
        vel_n = vel + dt / 6.0 * (k1v + 2 * k2v + 2 * k3v + k4v)
        n = len(live)
        L, _ = self.field_at(pos_n, pos_n, mass, exclude=list(range(n)))
        if self.profile is not None:
            L = L + self.ext_field(pos_n)[0]
        for i, b in enumerate(live):
            if b.fixed:
                v2 = 0.0
            elif b.kind == "craft":
                v2 = b.vx * b.vx + b.vy * b.vy
            else:
                b.x, b.y = pos_n[i]
                if (
                    b.kind == "photon" and self.law.kind == "metric"
                ):  # Zachová nulovou podmínku |v| = sqrt(A/B).
                    sp = math.hypot(*vel_n[i])
                    U = -float(L[i])
                    vel_n[i] *= math.sqrt(
                        float(self.law.A(U)) / float(self.law.B(U))
                    ) / max(sp, 1e-300)
                b.vx, b.vy = vel_n[i]
                v2 = float(vel_n[i] @ vel_n[i])
            if b.kind != "photon" and self.law.kind == "newton":
                b.tau += dt  # Newtonův model: čas plyne stejně, bez gravitačního zpomalení hodin.
            elif b.kind != "photon" and self.law.key != "H1":
                U = -max(float(L[i]), -20.0)
                delta = (
                    float(self.law.Am1(U)) - float(self.law.B(U)) * v2
                )  # A - B v^2 - 1
                d = -math.expm1(0.5 * math.log1p(max(delta, -1.0 + 1e-16)))
                b.tau_def += d * dt
                b.tau += (1.0 - d) * dt
            elif b.kind != "photon":
                # dτ/dt = sqrt(e^(2L) − e^(−2L)|v|²); deficit 1 − dτ/dt se počítá stabilním logaritmickým zápisem.
                # Deficit se počítá pomocí expm1/log1p, aby i pole Země řádu 10⁻¹⁰ zachovalo rozdíly hodin řádu 10⁻¹⁵,
                # které by se při prostém sčítání τ ≈ t zaokrouhlily pryč.
                Li = max(
                    float(L[i]), -20.0
                )  # Ochrana proti přetečení; kontrola silného pole zasahuje dříve.
                q = min(math.exp(-4.0 * Li) * v2, 1.0 - 1e-16)
                d = -math.expm1(Li + 0.5 * math.log1p(-q))
                b.tau_def += d * dt
                b.tau += (1.0 - d) * dt
        self.t += dt
        for i, b in enumerate(live):
            if b.kind == "craft" and b.script is not None:
                self.move_craft(b, pos_n, mass, i)
        self.land_on_solids()
        if not self.no_guard:
            self.strong_field_guard()
        return dt

    def move_craft(self, b, pos, mass, i):
        """Provede předepsaný pohyb tělesa podle scénáře SI a určí vlastní zrychlení a_kin − a_grav."""
        x, y, vx, vy, ax, ay = b.script.state(self.t / C_SI)
        b.x, b.y, b.vx, b.vy = x, y, vx / C_SI, vy / C_SI
        _, g = self.field_at(np.array([[x, y]]), pos, mass, exclude=[i])
        agx, agy = (
            -g[0, 0] * C_SI**2,
            -g[0, 1] * C_SI**2,
        )  # Newtonovské gravitační zrychlení v místě řízeného tělesa [m/s²].
        b.now_g = math.hypot(ax - agx, ay - agy)
        b.max_g = max(b.max_g, b.now_g)

    def land_on_solids(self):
        """Těleso, které pronikne dovnitř pevného objektu, se zastaví na jeho povrchu; foton je pohlcen."""
        solids = [o for o in self.bodies if o.alive and o.radius > 0]
        for b in self.bodies:
            if not b.alive or b.fixed or b.radius > 0 or b.kind == "craft":
                continue
            for o in solids:
                dx, dy = b.x - o.x, b.y - o.y
                r = math.hypot(dx, dy)
                if (
                    r < o.radius and b.kind == "pilot"
                ):  # Řízené těleso se zastaví na povrchu, ale může znovu vzlétnout.
                    ux, uy = dx / max(r, 1e-12), dy / max(r, 1e-12)
                    b.x, b.y = o.x + ux * o.radius, o.y + uy * o.radius
                    vr = b.vx * ux + b.vy * uy
                    if vr < 0:
                        b.landing_speed = -vr * C_SI
                        b.vx, b.vy = b.vx - vr * ux, b.vy - vr * uy
                        b.vx *= 0.98  # some ground friction
                        b.vy *= 0.98
                    continue
                if r < o.radius:
                    if r < 1e-12 * max(
                        o.radius, 1.0
                    ):  # Přesně ve středu není definován radiální směr,
                        dx, dy, r = 0.0, 1.0, 1.0  # proto se těleso umístí na horní bod povrchu.
                    if b.kind == "photon":
                        b.alive = False
                    else:
                        b.x, b.y = o.x + dx / r * o.radius, o.y + dy / r * o.radius
                        b.vx = b.vy = 0.0
                        b.fixed = True
                        b.landed_t = self.t

    def strong_field_guard(self):
        """
        Pokud příspěvky ostatních zdrojů vytvoří u tělesa |L| > L_MAX, aproximace 1PN se nepoužije.
        Foton je pohlcen a hmotné těleso splyne se zdrojem, který dominuje jeho poli. Při sloučení se zachovává
        hmotnost a hybnost v Newtonovském přiblížení; událost se zapíše do self.events.
        """
        src = [b for b in self.bodies if b.alive and b.mass > 0]
        for b in [x for x in self.bodies if x.alive and not x.fixed]:
            if (
                not b.alive
            ):  # Těleso bylo sloučeno dříve v tomto průchodu; jeho hmotnost už je započtena.
                continue
            others = [o for o in src if o is not b and o.alive]
            if not others:
                continue
            contrib = [
                (o.mass / math.sqrt((b.x - o.x) ** 2 + (b.y - o.y) ** 2 + EPS * EPS), o)
                for o in others
            ]
            L = -sum(c for c, _ in contrib)
            if L >= -L_MAX:
                continue
            dom = max(contrib, key=lambda c: c[0])[1]
            if b.kind == "photon":
                b.alive = False
                self.events.append(
                    f"t = {self.t:.2f}: {b.name} pohlcen u {dom.name} (|L| > {L_MAX:g}, mimo slabé pole)"
                )
            elif b.kind == "mass":
                big, small = (dom, b) if dom.mass >= b.mass else (b, dom)
                if small.fixed and not big.fixed:
                    big, small = small, big
                M = big.mass + small.mass
                if not big.fixed:
                    big.x = (big.mass * big.x + small.mass * small.x) / M
                    big.y = (big.mass * big.y + small.mass * small.y) / M
                    big.vx = (big.mass * big.vx + small.mass * small.vx) / M
                    big.vy = (big.mass * big.vy + small.mass * small.vy) / M
                big.mass = M
                small.alive = False
                self.events.append(
                    f"t = {self.t:.2f}: {small.name} splynulo s {big.name} (|L| > {L_MAX:g}, mimo slabé pole)"
                )

    def advance(self, T, max_steps=5000, on_step=None):
        """Posune simulaci o čas T nebo do limitu kroků; po každém kroku zavolá on_step(sim)."""
        t_end = self.t + T
        steps = 0
        while self.t < t_end and steps < max_steps:
            live, pos, vel, mass, photon, _ = self.arrays()
            dt = min(self.choose_dt(pos, vel, mass, photon), t_end - self.t)
            self.step(max(dt, 1e-9))
            steps += 1
            if on_step:
                on_step(self)
        return steps

    def pull_split(self, body, main):
        """Velikost ∇L u tělesa zvlášť od zvoleného hlavního zdroje a od všech ostatních zdrojů."""
        gm, go = np.zeros(2), np.zeros(2)
        for o in self.bodies:
            if not o.alive or o.mass <= 0 or o is body:
                continue
            d = np.array([body.x - o.x, body.y - o.y])
            s2 = float(d @ d) + EPS * EPS
            g = o.mass * d / (s2 * math.sqrt(s2))
            if o is main:
                gm += g
            else:
                go += g
        return float(np.hypot(*gm)), float(np.hypot(*go))

    def L_at(self, x, y, skip=None):
        """Hodnota L v jednom bodě od všech zdrojů; skip označuje těleso, které se má vynechat."""
        src = [
            b
            for b in self.bodies
            if b.alive and b.mass > 0 and b is not skip and b is not self.profile
        ]
        Le = (
            float(self.ext_field(np.array([[x, y]]))[0][0])
            if self.profile is not None
            else 0.0
        )
        if not src:
            return Le
        sp = np.array([[b.x, b.y] for b in src])
        sm = np.array([b.mass for b in src])
        L, _ = self.field_at(np.array([[x, y]]), sp, sm)
        return float(L[0]) + Le


# =============================================================================
# Měření: precese perihélia, odklon světla a Shapirův časový posun
# =============================================================================


def heaviest_other(sim, body):
    others = [b for b in sim.bodies if b.alive and b.mass > 0 and b is not body]
    return max(others, key=lambda b: b.mass) if others else None


class PerihelionTracker:
    """Měří úhel perihélia tělesa vzhledem k centrálnímu zdroji a porovnává posun s předpovědí 1PN."""

    def __init__(self, body, center):
        self.body, self.center = body, center
        self.r_hist = deque(maxlen=3)
        self.peri = []  # (čas, úhel, vzdálenost)
        self.r_max = 0.0
        self.r_min = math.inf
        self.other_max = 0.0  # Největší působení ostatních těles vzhledem k hlavnímu zdroji.
        self.lrl = []  # (čas, směr vektoru Laplace–Runge–Lenz v perihéliu)

    def lrl_angle(self):
        """Směr oskulačního Newtonovského vektoru excentricity e = (v²/μ − 1/r)r − (r·v/μ)v.
        Vzorkování v každém perihéliu potlačí periodické kolísání oskulační dráhy; dlouhodobý posun udává
        precese perihélia přesněji než přímé hledání diskrétních minim vzdálenosti."""
        b, c = self.body, self.center
        rx, ry = b.x - c.x, b.y - c.y
        vx, vy = b.vx - c.vx, b.vy - c.vy
        mu = c.mass + b.mass
        r = math.hypot(rx, ry)
        v2, rv = vx * vx + vy * vy, rx * vx + ry * vy
        ex, ey = (v2 / mu - 1.0 / r) * rx - rv / mu * vx, (
            v2 / mu - 1.0 / r
        ) * ry - rv / mu * vy
        return math.atan2(ey, ex)

    def rate(self):
        """Určí sekulární rychlost precese [rad na jednotku času] lineárním fitem úhlů vektoru Laplace–Runge–Lenz."""
        if len(self.lrl) < 3:
            return None
        t = np.array([x[0] for x in self.lrl])
        a = np.unwrap([x[1] for x in self.lrl])
        return float(np.polyfit(t, a, 1)[0]) * getattr(self, "sense", 1.0)

    def update(self, sim):
        b, c = self.body, self.center
        dx, dy = b.x - c.x, b.y - c.y
        r = math.hypot(dx, dy)
        self.r_hist.append((r, math.atan2(dy, dx), sim.t))
        self.r_max, self.r_min = max(self.r_max, r), min(self.r_min, r)
        gm, go = sim.pull_split(b, c)
        if gm > 0:
            self.other_max = max(self.other_max, go / gm)
        lz = dx * (b.vy - c.vy) - dy * (b.vx - c.vx)
        self.sense = 1.0 if lz >= 0 else -1.0  # +1 counter-clockwise, -1 clockwise
        if len(self.r_hist) == 3:
            (r0, a0, t0), (r1, a1, t1), (r2, a2, t2) = self.r_hist
            near_min = (
                r1 - self.r_min <= 0.3 * (self.r_max - self.r_min) + 1e-12
            )  # Vyloučí drobné lokální kolísání mimo perihélium.
            if (
                r1 < r0 and r1 <= r2 and near_min
            ):  # Minimum vzdálenosti se zpřesní parabolickou interpolací.
                den = r0 - 2 * r1 + r2
                u = 0.5 * (r0 - r2) / den if den > 0 else 0.0
                da_m = math.atan2(math.sin(a1 - a0), math.cos(a1 - a0))
                da_p = math.atan2(math.sin(a2 - a1), math.cos(a2 - a1))
                ang = a1 + (u * (da_p if u > 0 else da_m))
                if not self.peri or t1 - self.peri[-1][0] > 1e-9:
                    self.peri.append((t1, ang, r1))
                    self.lrl.append((sim.t, self.lrl_angle()))

    def result(self):
        """Vrátí posun na oběh, předpověď GR pro β = γ = 1, velkou poloosu, excentricitu a počet oběhů.
        Bez alespoň dvou úplných oběhů nebo pro téměř kruhovou dráhu (e < 0,05) vrátí None."""
        if len(self.peri) < 3 or self.r_max <= 0:
            return None
        a = 0.5 * (self.r_min + self.r_max)
        e = (self.r_max - self.r_min) / (self.r_max + self.r_min)
        if e < 0.05:  # U téměř kruhové dráhy není směr perihélia dobře definován.
            return None
        steps = [
            math.atan2(math.sin(p2[1] - p1[1]), math.cos(p2[1] - p1[1]))
            for p1, p2 in zip(self.peri[:-1], self.peri[1:])
        ]
        adv = float(np.median(steps)) * getattr(self, "sense", 1.0)
        M = self.center.mass + self.body.mass
        pred = 6.0 * math.pi * M / (a * (1.0 - e * e))
        return adv, pred, a, e, len(steps)


class LightTracker:
    """Měří odklon fotonu a Shapirův časový posun vůči předpovědi pro přímou trajektorii."""

    def __init__(self, photon, mass_body, sim):
        self.ph, self.m = photon, mass_body
        sp = math.hypot(photon.vx, photon.vy)
        self.d0 = (photon.vx / sp, photon.vy / sp)
        rx, ry = photon.x - mass_body.x, photon.y - mass_body.y
        self.x1 = rx * self.d0[0] + ry * self.d0[1]  # poloha podél počáteční přímky
        self.b = abs(rx * self.d0[1] - ry * self.d0[0])  # parametr dopadu
        self.t0 = sim.t
        self.path = 0.0  # Eukleidovská délka skutečné dráhy.
        self.last = (photon.x, photon.y)
        self.int_main, self.int_other = (
            0.0,
            0.0,
        )  # Integrály působení hlavního zdroje a ostatních těles podél dráhy.
        self.turn = 0.0  # Nasčítané podepsané natočení směru pohybu.
        self.last_dir = math.atan2(photon.vy, photon.vx)

    def update(self, sim):
        """Sčítá délku skutečné, mírně zakřivené dráhy pro výpočet geometrické části zpoždění."""
        x, y = self.ph.x, self.ph.y
        dl = math.hypot(x - self.last[0], y - self.last[1])
        self.path += dl
        self.last = (x, y)
        gm, go = sim.pull_split(self.ph, self.m)
        self.int_main += gm * dl
        self.int_other += go * dl
        ang = math.atan2(
            self.ph.vy, self.ph.vx
        )  # Úhel se sčítá po jednotlivých krocích, takže
        self.turn += math.atan2(
            math.sin(ang - self.last_dir), math.cos(ang - self.last_dir)
        )  # foton, který
        self.last_dir = ang  # oběhne černou díru, může mít odklon větší než 180°.

    @staticmethod
    def _x_plus_r(x, b):
        """Vypočítá x + sqrt(x² + b²) bez ztráty přesnosti pro x < 0; tehdy použije ekvivalentní b²/(r − x)."""
        r = math.hypot(x, b)
        return x + r if x >= 0 else b * b / (r - x)

    def state(self, sim):
        ph, m, d0 = self.ph, self.m, self.d0
        sp = math.hypot(ph.vx, ph.vy)
        cosang = max(-1.0, min(1.0, (ph.vx * d0[0] + ph.vy * d0[1]) / max(sp, 1e-300)))
        defl = abs(self.turn) if self.turn != 0.0 else math.acos(cosang)
        rx, ry = ph.x - m.x, ph.y - m.y
        x2 = rx * d0[0] + ry * d0[1]
        b, x1, M = self.b, self.x1, m.mass
        other = self.int_other / self.int_main if self.int_main > 0 else 0.0
        # Paprsek míří přímo na zdroj nebo tak blízko, že v běžné scéně vstoupí do oblasti silného pole;
        # ochrana jej pohltí a předpovědi slabého pole zde neplatí.
        if b <= (EPS if sim.no_guard else max(EPS, M / L_MAX)):
            return {
                "defl": defl,
                "pred": None,
                "pred2": None,
                "b": b,
                "delay": (sim.t - self.t0) - self.path,
                "delay_pred": None,
                "other": other,
                "direct": True,
            }
        geom = (x2 / math.hypot(x2, b) - x1 / math.hypot(x1, b)) * M / b
        pred = {g: (1.0 + g) * geom for g in (1.0, 0.0, -1.0)}
        # Druhý řád této metriky: PPN dává π(8 + 8γ − 4β + 3ε)/4 · (M/b)²; exponenciální metrika
        # gᵢⱼ = 1 + 2U + 2U² má ε = 4/3, tedy 4π(M/b)²; GR dává 15π/4. Jde o limitu nekonečné dráhy.
        pred2 = 4.0 * math.pi * (M / b) ** 2
        # Souřadnicový čas minus doba letu po skutečné dráze v plochém prostoru; samotný ohyb přidává geometrické zpoždění
        # druhého řádu (~ α²L/2), které není součástí Shapirova vztahu, a proto se odečítá.
        delay_meas = (sim.t - self.t0) - self.path
        lr = math.log(self._x_plus_r(x2, b) / self._x_plus_r(x1, b))
        delay_pred = {g: (1.0 + g) * M * lr for g in (1.0, 0.0, -1.0)}
        return {
            "defl": defl,
            "pred": pred,
            "pred2": pred2,
            "b": b,
            "delay": delay_meas,
            "delay_pred": delay_pred,
            "other": other,
            "direct": False,
        }


def make_trackers(sim, orbits, old_p=None, old_l=None):
    """Vytvoří sledovače perihélia a fotonů. Existující sledovače zachová podle identity tělesa, takže přidání objektu
    nesmaže dosavadní měření."""
    old_p, old_l = old_p or {}, old_l or {}
    trackers, light = [], []
    for b in sim.bodies:
        if not b.alive:
            continue
        if b.kind == "mass" and not b.fixed and orbits:
            c = heaviest_other(sim, b)
            if id(b) in old_p and old_p[id(b)].center is c:
                trackers.append(old_p[id(b)])
            elif c is not None and (
                c.mass > 10 * b.mass or (orbits == "binary" and c.mass >= b.mass)
            ):
                trackers.append(PerihelionTracker(b, c))
        if b.kind == "photon":
            if id(b) in old_l:
                light.append(old_l[id(b)])
            else:
                m = heaviest_other(sim, b)
                if m is not None:
                    light.append(LightTracker(b, m, sim))
    return trackers, light


def mercury_experiment(laws, years, planets):
    """
    Generator: perihelion advance of Mercury in arcseconds per century under each law, over `years` years, with
    or without the other planets (J2000, 2D). Yields (law name, fraction done) while running and returns
    {law key: rate}. The rate is the slope of the Laplace-Runge-Lenz direction at the perihelia (PerihelionTracker).
    Expected (observations and GR): ~43''/century from relativity alone; the other planets add several
    hundred arcseconds of Newtonian perturbation (the real 3D value is ~532''/century; this 2D model differs).
    """
    out = {}
    for law in laws:
        sim = Sim(law)
        sim.bodies = solar_system(planets=planets, only=None if planets else ["Merkur"])
        sim.eta, sim.dt_max = 0.004, 2.0 * 86400 * C_SI
        sun, mer = sim.bodies[0], next(b for b in sim.bodies if b.name == "Merkur")
        tr = PerihelionTracker(mer, sun)
        T = years * YEAR * C_SI
        k = 0
        while sim.t < T:
            live, pos, vel, mass, photon, _ = sim.arrays()
            sim.step(max(min(sim.choose_dt(pos, vel, mass, photon), T - sim.t), 1e-9))
            tr.update(sim)
            k += 1
            if k % 200 == 0:
                yield law.name, sim.t / T
        r = tr.rate()
        out[law.key] = None if r is None else r * 100 * YEAR * C_SI / ARCSEC
    return out


# =============================================================================
# scenes (presets)
# =============================================================================


def kepler_peri_velocity(M, a, e):
    """Newtonovská rychlost v perihéliu dráhy s velkou poloosou a a excentricitou e."""
    return math.sqrt(M * (1 + e) / (a * (1 - e)))


# Střední elementy J2000: hmotnost [kg], a [AU], excentricita e a délky perihélia a střední délka [°].
PLANETS = [
    ("Merkur", 3.3011e23, 0.387098, 0.205630, 77.456, 252.251, "#c9a37a"),
    ("Venuše", 4.8675e24, 0.723332, 0.006772, 131.533, 181.980, "#ffe0a0"),
    ("Země+Měsíc", 6.0458e24, 1.000001, 0.016709, 102.947, 100.464, "#7ac8ff"),
    ("Mars", 6.4171e23, 1.523679, 0.093400, 336.041, 355.453, "#ff7a5a"),
    ("Jupiter", 1.89819e27, 5.20260, 0.048900, 14.753, 34.404, "#ffcf8a"),
    ("Saturn", 5.6834e26, 9.55490, 0.056500, 92.432, 49.944, "#ffe7b0"),
    ("Uran", 8.6813e25, 19.2184, 0.046381, 170.964, 313.232, "#a0f0ff"),
    ("Neptun", 1.02413e26, 30.1104, 0.009456, 44.971, 304.880, "#7a9aff"),
]


def kepler_state(mu, a, e, varpi, lam):
    """Poloha a rychlost na Keplerově elipse při G = c = 1, rychlost je vyjádřena v jednotkách c.
    Vstupy: gravitační parametr μ [m], velká poloosa a [m], excentricita e, délka perihélia varpi a střední délka lam [rad]."""
    M = (lam - varpi) % (2 * math.pi)
    E = M if e < 0.8 else math.pi
    for _ in range(50):  # Newtonovy iterace Keplerovy rovnice
        E -= (E - e * math.sin(E) - M) / (1 - e * math.cos(E))
    n = math.sqrt(mu / a**3)
    x, y = a * (math.cos(E) - e), a * math.sqrt(1 - e * e) * math.sin(E)
    k = a * n / (1 - e * math.cos(E))
    vx, vy = -k * math.sin(E), k * math.sqrt(1 - e * e) * math.cos(E)
    c, s_ = math.cos(varpi), math.sin(varpi)
    return x * c - y * s_, x * s_ + y * c, vx * c - vy * s_, vx * s_ + vy * c


def solar_system(planets=True, only=None):
    """Slunce a planety pro epochu J2000 ve 2D; Slunce se pohybuje tak, aby celková hybnost byla nulová.
    Parametr only omezuje výběr na zadaná jména."""
    sun = Body("Slunce", 0, 0, mass=M_SUN, color="#ffd27a")
    B = [sun]
    if planets or only:
        for nm, kg, a, e, vp, lam, col in PLANETS:
            if only is not None and nm not in only:
                continue
            x, y, vx, vy = kepler_state(
                M_SUN, a * AU, e, math.radians(vp), math.radians(lam)
            )
            B.append(Body(nm, x, y, vx, vy, mass=kg * G_C2, color=col))
    px = sum(b.mass * b.vx for b in B[1:])
    py = sum(b.mass * b.vy for b in B[1:])
    sun.vx, sun.vy = -px / sun.mass, -py / sun.mass
    return B


class CraftScript:
    """
    Předepsaný pohyb v SI: počáteční stav a fáze zadané délkou [s] a konstantním kinematickým zrychlením ax, ay [m/s²].
    Po poslední fázi těleso pokračuje poslední rychlostí; připravené scénáře končí v klidu.
    """

    def __init__(self, x0, y0, vx0, vy0, phases):
        self.phases = []
        t, x, y, vx, vy = 0.0, float(x0), float(y0), float(vx0), float(vy0)
        for dur, ax, ay in phases:
            self.phases.append((t, t + dur, x, y, vx, vy, ax, ay))
            x, y = (
                x + vx * dur + 0.5 * ax * dur * dur,
                y + vy * dur + 0.5 * ay * dur * dur,
            )
            vx, vy = vx + ax * dur, vy + ay * dur
            t += dur
        self.end = (t, x, y, vx, vy)

    def state(self, t):
        """Vrátí (x, y, vx, vy, ax, ay) v čase t [s]."""
        for t0, t1, x, y, vx, vy, ax, ay in self.phases:
            if t < t1:
                d = max(t - t0, 0.0)
                return (
                    x + vx * d + 0.5 * ax * d * d,
                    y + vy * d + 0.5 * ay * d * d,
                    vx + ax * d,
                    vy + ay * d,
                    ax,
                    ay,
                )
        t0, x, y, vx, vy = self.end
        d = t - t0
        return x + vx * d, y + vy * d, vx, vy, 0.0, 0.0

    def phase_name(self, t, names):
        for k, (t0, t1, *_rest) in enumerate(self.phases):
            if t < t1:
                return names[k] if k < len(names) else ""
        return names[-1] if names else ""


def scene(name):
    """
    Bodies of a preset scene, its view (center x, y, half-width) and meta:
      si        True: délky v metrech a zobrazení času v sekundách; False: geometrické jednotky G = c = 1
      dt_max    největší integrační krok
      speed     výchozí logaritmické tempo simulace na sekundu reálného času ve všech scénách
      display   pevná barevná škála pole nebo relativní odchylka od referenční hodnoty
      click     „orbit“ vytvoří kruhovou dráhu; „rest“ pustí objekt z klidu
    """
    B = []
    meta = {
        "si": False,
        "dt_max": 0.5,
        # V původním ovládání představovala hodnota 0 přibližně jednu jednotku na snímek (≈30/s).
        "speed": math.log10(30.0),
        "display": ("log", -4.0, -1.0),
        "click": "orbit",
    }
    if name == "Precese perihelia":
        M = 0.01
        B.append(Body("hvězda", 0, 0, mass=M, color="#ffd27a"))
        a, e = 1.0, 0.5
        # Newtonovská rychlost v perihéliu; účinek 1PN je zde slabší o faktor (1 − 4M/r + v²), proto je skutečná dráha
        # větší než a = 1. Hodnoty a a e se následně měří ze simulované dráhy.
        B.append(
            Body(
                "planeta",
                a * (1 - e),
                0,
                0,
                kepler_peri_velocity(M, a, e),
                mass=1e-7,
                color="#7ac8ff",
            )
        )
        view = (-0.6, 0.0, 3.0)
    elif name == "Sluneční soustava (zveličeno)":
        M = 0.005
        B.append(Body("hvězda", 0, 0, mass=M, color="#ffd27a"))
        for nm, a, e, m, col in (
            ("A", 0.5, 0.3, 2e-7, "#c9a37a"),
            ("B", 0.9, 0.05, 1e-6, "#7affb0"),
            ("C", 1.5, 0.1, 5e-6, "#7ac8ff"),
            ("D", 2.6, 0.2, 2e-5, "#ff9a7a"),
        ):
            B.append(
                Body(
                    nm,
                    a * (1 - e),
                    0,
                    0,
                    kepler_peri_velocity(M, a, e),
                    mass=m,
                    color=col,
                )
            )
        view = (0.0, 0.0, 3.2)
    elif name == "Ohyb světla":
        B.append(Body("hmota", 0, 0, mass=0.01, fixed=True, color="#ffd27a"))
        for i, bb in enumerate((0.25, 0.4, 0.7, 1.0, 1.6)):
            B.append(
                Body(
                    f"foton b={bb:g}",
                    -6.0,
                    bb,
                    1.0,
                    0.0,
                    kind="photon",
                    color="#fff27a",
                )
            )
        view = (0.0, 0.0, 6.5)
    elif name == "Shapirovo zpoždění":
        B.append(Body("hmota", 0, 0, mass=0.01, fixed=True, color="#ffd27a"))
        B.append(
            Body("foton b=0.3", -20.0, 0.3, 1.0, 0.0, kind="photon", color="#fff27a")
        )
        B.append(
            Body("foton b=3", -20.0, 3.0, 1.0, 0.0, kind="photon", color="#ffb37a")
        )
        view = (0.0, 0.0, 21.0)
    elif name == "Dvojhvězda a planeta":
        M1 = M2 = 0.005
        sep = 0.6
        v = 0.5 * math.sqrt((M1 + M2) / sep)
        B.append(Body("hvězda 1", -sep / 2, 0, 0, -v, mass=M1, color="#ffd27a"))
        B.append(Body("hvězda 2", sep / 2, 0, 0, v, mass=M2, color="#ff9a5a"))
        rp = 2.4
        B.append(
            Body(
                "planeta",
                rp,
                0,
                0,
                math.sqrt((M1 + M2) / rp),
                mass=1e-6,
                color="#7ac8ff",
            )
        )
        view = (0.0, 0.0, 3.2)
    elif name == "Hodiny v poli":
        M = 0.02
        B.append(Body("hmota", 0, 0, mass=M, fixed=True, color="#ffd27a"))
        for r in (0.3, 1.0, 4.0):
            B.append(Body(f"hodiny r={r:g}", r, 0, kind="clock", color="#b0b0ff"))
        a = 1.0
        B.append(
            Body(
                "hodiny na oběžné dráze r=1",
                0,
                a,
                -math.sqrt(M / a),
                0,
                mass=1e-9,
                color="#7affb0",
            )
        )
        view = (0.0, 0.0, 4.5)
    elif name == "Povrch Země: jablko a hodiny":
        R = R_EARTH
        B.append(
            Body("Země", 0, 0, mass=M_EARTH, fixed=True, color="#3f6b3a", radius=R)
        )
        B.append(
            Body("jablko (10 m)", -12.0, R + 10.0, mass=0.2 * G_C2, color="#ff5a4a")
        )
        v0, ang = 14.0, math.radians(
            50
        )  # vržený míč: 14 m/s pod úhlem 50° z výšky 1,5 m.
        B.append(
            Body(
                "míč (14 m/s, 50°)",
                -30.0,
                R + 1.5,
                v0 * math.cos(ang) / C_SI,
                v0 * math.sin(ang) / C_SI,
                mass=0.4 * G_C2,
                color="#ffd27a",
            )
        )
        B.append(Body("hodiny na zemi", 20.0, R + 0.01, kind="clock", color="#b0b0ff"))
        B.append(
            Body("hodiny na věži (30 m)", 20.0, R + 30.0, kind="clock", color="#b0b0ff")
        )
        view = (0.0, R + 15.0, 40.0)
        meta.update(
            si=True,
            dt_max=10.0 * C_SI,
            speed=0.0,
            display=("rel", -M_EARTH / R, M_EARTH / R**2 * 60.0),
            click="rest",
        )
    elif name == "Země: ISS a GPS":
        R = R_EARTH
        B.append(
            Body("Země", 0, 0, mass=M_EARTH, fixed=True, color="#3f6b3a", radius=R)
        )
        B.append(Body("hodiny na zemi", 0, R + 0.01, kind="clock", color="#b0b0ff"))
        for nm, r, col in (("ISS", 6.771e6, "#ffffff"), ("GPS", 2.6560e7, "#ffd27a")):
            B.append(
                Body(nm, r, 0, 0, math.sqrt(M_EARTH / r), mass=1e3 * G_C2, color=col)
            )
        view = (0.0, 0.0, 3.0e7)
        meta.update(
            si=True,
            dt_max=10.0 * C_SI,
            speed=math.log10(2000.0),
            display=("log", -10.2, -9.1),
            click="orbit",
        )
    elif name == "Sluneční soustava (skutečná)":
        B = solar_system()
        view = (0.0, 0.0, 1.8 * AU)
        meta.update(
            si=True,
            dt_max=2.0 * 86400 * C_SI,
            speed=math.log10(30 * 86400.0),
            eta=0.004,
            display=("log", -9.7, -7.0),
            click="orbit",
            orbits=True,
        )
    elif name == "UAP: manévr (hypotéza)":
        R = R_EARTH
        B.append(
            Body("Země", 0, 0, mass=M_EARTH, fixed=True, color="#3f6b3a", radius=R)
        )
        x0, y0 = -3000.0, R + 1000.0
        craft = Body("UAP", x0, y0, kind="craft", color="#9affff")
        # Ilustrační manévr: vznášení, rychlý let, prudké zastavení, vznášení, kolmý výstup a zastavení.
        craft.script = CraftScript(
            x0,
            y0,
            0.0,
            0.0,
            [
                (2.0, 0, 0),
                (1.0, 1000.0, 0),
                (2.0, 0, 0),
                (0.2, -5000.0, 0),
                (1.5, 0, 0),
                (1.0, 0, 1500.0),
                (1.0, 0, -1500.0),
                (30.0, 0, 0),
            ],
        )
        craft.phase_names = [
            "vznáší se",
            "zrychlení 1000 m/s²",
            "let 1000 m/s",
            "zastavení za 0,2 s",
            "vznáší se",
            "stoupání 1500 m/s²",
            "brzdění",
            "vznáší se",
        ]
        B.append(craft)
        B.append(
            Body("kámen (volný pád z 1000 m)", x0, y0, mass=1.0 * G_C2, color="#ff9a5a")
        )
        view = (0.0, R + 1200.0, 4000.0)
        meta.update(
            si=True,
            dt_max=0.01 * C_SI,
            speed=0.0,
            display=("rel", -M_EARTH / R, M_EARTH / R**2 * 3000.0),
            click="rest",
        )
    elif name == "Černá díra (silné pole, neověřeno)":
        # Jednotky G = c = 1 a M = 1, délky v GM/c². Hlavní scéna používá přesné vnější řešení
        # rovnice pole laboratoře (bh_profile); srovnávací běh GR používá Schwarzschildovu metriku.
        bh = Body("černá díra? (M = 1)", 0, 0, mass=1.0, fixed=True, color="#202020")
        bh.is_profile = True
        B.append(bh)
        for bb in (3.5, 4.0, 4.5, 5.0, 5.15, 5.21, 5.23, 5.26, 5.5, 6.0, 7.0, 8.0):
            # Zachovaný parametr dopadu je počáteční výška y₀ násobená indexem lomu n = sqrt(B/A) = e^(−2L), proto
            # je počáteční výška y₀ = b/n. Index GR se zde liší až ve druhém řádu, přibližně o 10⁻⁴.
            y0 = bb
            for _ in range(5):
                y0 = bb * math.exp(2.0 * float(bh_L([math.hypot(40.0, y0)])[0][0]))
            B.append(
                Body(
                    f"foton b={bb:g}",
                    -40.0,
                    y0,
                    1.0,
                    0.0,
                    kind="photon",
                    color="#fff27a",
                )
            )
        for R_ar in (
            3.0,
            4.0,
            6.0,
            10.0,
        ):  # Statické hodiny v plošných poloměrech R; izotropní r splňuje R = r e^(−L).
            lo, hi = 0.05, R_ar
            for _ in range(80):
                mid = 0.5 * (lo + hi)
                if mid * math.exp(-float(bh_L([mid])[0][0])) < R_ar:
                    lo = mid
                else:
                    hi = mid
            B.append(
                Body(
                    f"hodiny R={R_ar:g}M",
                    -0.5 * (lo + hi),
                    0.0,
                    kind="clock",
                    color="#b0b0ff",
                )
            )
        r0 = 12.0
        B.append(
            Body(
                "částice (periapsida 12 M)",
                r0,
                0,
                0,
                0.85 * math.sqrt(1.0 / r0),
                mass=1e-9,
                color="#7ac8ff",
            )
        )
        view = (0.0, 0.0, 14.0)
        meta.update(
            dt_max=1.0, eta=0.01, display=("log", -2.0, 0.7), bh=True, speed=math.log10(10**(-0.2) / 0.03)
        )
    elif name == "UAP: řízení (klávesnice)":
        R = R_EARTH
        B.append(
            Body("Země", 0, 0, mass=M_EARTH, fixed=True, color="#3f6b3a", radius=R)
        )
        pl = Body("UAP (pilot)", 0.0, R + 500.0, kind="pilot", color="#9affff")
        pl.hover = True
        B.append(pl)
        B.append(
            Body(
                "kámen (volný pád z 500 m)",
                60.0,
                R + 500.0,
                mass=1.0 * G_C2,
                color="#ff9a5a",
            )
        )
        view = (0.0, R + 500.0, 1500.0)
        meta.update(
            si=True,
            dt_max=0.005 * C_SI,
            speed=0.0,
            display=("rel", -M_EARTH / R, M_EARTH / R**2 * 3000.0),
            click="rest",
            pilot=True,
        )
    elif name in ("α Centauri A/B", "Sirius A/B (bílý trpaslík)"):
        if name.startswith("α"):  # masses [M_sun], a [AU], e, period 79.91 yr
            (n1, m1, c1), (n2, m2, c2), a_au, e = (
                ("α Cen A", 1.0788, "#fff4d0"),
                ("α Cen B", 0.9092, "#ffd29a"),
                23.33,
                0.5208,
            )
        else:  # Sirius: period 50.13 yr
            (n1, m1, c1), (n2, m2, c2), a_au, e = (
                ("Sirius A", 2.063, "#cfe0ff"),
                ("Sirius B", 1.018, "#ffffff"),
                19.8,
                0.5923,
            )
        M1, M2 = m1 * M_SUN, m2 * M_SUN
        x, y, vx, vy = kepler_state(
            M1 + M2, a_au * AU, e, 0.0, 0.0
        )  # Relativní dráha, počáteční bod je periastr.
        f1, f2 = -M2 / (M1 + M2), M1 / (M1 + M2)
        B.append(Body(n1, f1 * x, f1 * y, f1 * vx, f1 * vy, mass=M1, color=c1))
        B.append(Body(n2, f2 * x, f2 * y, f2 * vx, f2 * vy, mass=M2, color=c2))
        view = (0.0, 0.0, 1.3 * a_au * AU)
        meta.update(
            si=True,
            dt_max=20 * 86400 * C_SI,
            speed=math.log10(2 * YEAR),
            eta=0.004,
            display=("log", -9.5, -6.5),
            click="orbit",
            orbits="binary",
        )
    elif name == "Okolí Slunce do 12 ly (3D)":
        # Ekvatoriální polohy J2000 (OpenStax Astronomy, dodatek I) a přibližné hmotnosti [M☉].
        STARS = [
            ("Slunce", 0.0, 0.0, 0.0, "G2 V", 1.0),
            ("Proxima Centauri", 14.48, -62.67, 4.2, "M5.5 V", 0.12),
            ("α Centauri A", 14.65, -60.83, 4.4, "G2 V", 1.08),
            ("α Centauri B", 14.65, -60.83, 4.4, "K1 V", 0.91),
            ("Barnardova hvězda", 17.95, 4.70, 6.0, "M4 V", 0.14),
            ("Wolf 359", 10.93, 7.00, 7.8, "M6 V", 0.11),
            ("Lalande 21185", 11.05, 35.97, 8.3, "M2 V", 0.39),
            ("Sirius A", 6.75, -16.70, 8.6, "A1 V", 2.06),
            ("Sirius B", 6.75, -16.72, 8.6, "DA2 (bílý trpaslík)", 1.02),
            ("Luyten 726-8 A", 1.65, -17.95, 8.7, "M5.5 V", 0.12),
            ("Ross 154", 18.82, -23.83, 9.7, "M3.5 V", 0.17),
            ("Ross 248", 23.68, 44.17, 10.3, "M5.5 V", 0.14),
            ("ε Eridani", 3.53, -9.45, 10.5, "K2 V", 0.82),
            ("Lacaille 9352", 23.08, -35.85, 10.7, "M0.5 V", 0.48),
            ("Ross 128", 11.78, 0.80, 10.9, "M4 V", 0.17),
            ("EZ Aquarii A", 22.63, -15.28, 11.3, "M5 V", 0.11),
            ("Procyon A", 7.65, 5.22, 11.4, "F5 IV", 1.50),
            ("Procyon B", 7.65, 5.22, 11.4, "bílý trpaslík", 0.60),
            ("61 Cygni A", 21.10, 38.73, 11.4, "K5 V", 0.70),
            ("61 Cygni B", 21.10, 38.73, 11.4, "K7 V", 0.63),
            ("Struve 2398 A", 18.70, 59.62, 11.5, "M3 V", 0.33),
            ("Groombridge 34 A", 0.30, 44.02, 11.6, "M1.5 V", 0.38),
            ("ε Indi A", 22.05, -56.77, 11.8, "K5 V", 0.76),
            ("τ Ceti", 1.73, -15.93, 11.9, "G8.5 V", 0.78),
        ]
        # Vlastní pohyb μ [″/rok], poziční úhel [°, od severu přes východ] a radiální rychlost [km/s]:
        # μ a poziční úhel vycházejí z tabulky RECONS pro 100 nejbližších systémů (GJ 15 A z Kervella et al. 2019); radiální rychlosti
        # z Kervella et al. 2019, A&A 623, A72, tabulka 1; neuvedené hodnoty jsou neznámé (None znamená jen pohyb na obloze).
        # Průvodci bez vlastních údajů sdílejí pohyb hlavní hvězdy. Chybějící μ znamená, že data o pohybu nejsou k dispozici.
        MOTION = {
            "Proxima Centauri": (3.853, 281.5, None),
            "α Centauri A": (3.710, 277.5, None),
            "α Centauri B": (3.724, 284.8, None),
            "Barnardova hvězda": (10.358, 355.6, -110.40),
            "Wolf 359": (4.696, 234.6, None),
            "Lalande 21185": (4.802, 186.9, None),
            "Sirius A": (1.339, 204.1, None),
            "Sirius B": (1.339, 204.1, None),
            "Luyten 726-8 A": (None, None, None),
            "Ross 154": (0.666, 106.8, -10.41),
            "Ross 248": (1.617, 177.0, None),
            "ε Eridani": (0.977, 271.1, 16.20),
            "Lacaille 9352": (6.896, 78.9, None),
            "Ross 128": (1.361, 153.6, -30.99),
            "EZ Aquarii A": (3.254, 46.6, None),
            "Procyon A": (1.259, 214.7, None),
            "Procyon B": (1.259, 214.7, None),
            "61 Cygni A": (5.281, 51.9, -65.73),
            "61 Cygni B": (5.172, 52.6, -64.04),
            "Struve 2398 A": (None, None, None),
            "Groombridge 34 A": (2.919, 81.9, 11.81),
            "ε Indi A": (4.704, 122.7, -40.10),
            "τ Ceti": (1.922, 296.4, -16.75),
        }
        colors = {
            "O": "#9bb0ff",
            "B": "#aabfff",
            "A": "#cad7ff",
            "F": "#f8f7ff",
            "G": "#fff4ea",
            "K": "#ffd2a1",
            "M": "#ffb07a",
            "D": "#ffffff",
            "b": "#ffffff",
        }
        for nm, ra_h, dec_d, d_ly, sp, msun in STARS:
            a, dlt, d = math.radians(15 * ra_h), math.radians(dec_d), d_ly * LY
            st = Body(
                nm,
                d * math.cos(dlt) * math.cos(a),
                d * math.cos(dlt) * math.sin(a),
                mass=msun * M_SUN,
                fixed=True,
                color=colors.get(sp[0], "#ffffff"),
            )
            st.z = d * math.sin(dlt)
            st.info = f"{sp}, {msun:g} M☉, {d_ly:g} ly"
            mu, pa, rv = MOTION.get(nm, (None, None, None))
            st.v3 = (
                0.0,
                0.0,
                0.0,
            )  # Prostorová rychlost v jednotkách c vůči Slunci, které je v klidu.
            st.motion_note = (
                "Slunce (počátek)" if nm == "Slunce" else "bez údajů o pohybu"
            )
            if mu is not None:
                vt = 4.74047 * mu * d_ly / 3.26156 * 1000.0  # Tečná rychlost [m/s].
                rr = (
                    math.cos(dlt) * math.cos(a),
                    math.cos(dlt) * math.sin(a),
                    math.sin(dlt),
                )
                ee = (-math.sin(a), math.cos(a), 0.0)  # east
                nn = (
                    -math.sin(dlt) * math.cos(a),
                    -math.sin(dlt) * math.sin(a),
                    math.cos(dlt),
                )  # north
                p_ = math.radians(pa)
                vr = (rv or 0.0) * 1000.0
                st.v3 = tuple(
                    (vr * rr[k] + vt * (math.sin(p_) * ee[k] + math.cos(p_) * nn[k]))
                    / C_SI
                    for k in range(3)
                )
                st.motion_note = (
                    f"v_t {vt / 1000:.0f} km/s, v_r {rv:+.1f} km/s"
                    if rv is not None
                    else f"v_t {vt / 1000:.0f} km/s, v_r neznámá (jen pohyb po obloze)"
                )
            st.trail3 = deque(maxlen=4000)
            B.append(st)
        view = (0.0, 0.0, 13.5 * LY)
        meta.update(
            si=True,
            dt_max=YEAR * C_SI,
            speed=math.log10(1000 * YEAR),
            display=("log", -14.0, -12.0),
            map3d=True,
        )
    else:  # Volná scéna používá pevný zdroj, aby přidané oběžné dráhy zůstaly stabilní.
        B.append(Body("hmota (pevná)", 0, 0, mass=0.01, fixed=True, color="#ffd27a"))
        view = (0.0, 0.0, 3.0)
    return B, view, meta


def category_astro(m):
    """Zařadí těleso podle hmotnosti v jednotkách G = c = 1. Orientační hranice: hvězda ≥ 0,08 M☉,
    hnědý trpaslík ≥ 0,013 M☉, obří planeta ≥ 10⁻⁴ M☉, planeta ≥ 3×10⁻⁷ M☉, měsíc ≥ 10⁻¹⁰ M☉,
    jinak planetka."""
    r = m / 0.005
    for lim, name in (
        (0.08, "hvězda"),
        (0.013, "hnědý trpaslík"),
        (1e-4, "obří planeta"),
        (3e-7, "planeta"),
        (1e-10, "měsíc"),
    ):
        if r >= lim:
            return name
    return "asteroid"


def category_si(kg):
    """Pojmenuje objekt u Země podle jeho hmotnosti v kilogramech."""
    for lim, name in (
        (1.6e29, "hvězda"),
        (2.5e28, "hnědý trpaslík"),
        (1.9e26, "obří planeta"),
        (1e23, "planeta"),
        (1e19, "měsíc"),
        (1e15, "planetka"),
        (1e8, "balvan / budova"),
        (1e3, "auto / letadlo"),
        (5.0, "kámen"),
        (0.05, "jablko"),
    ):
        if kg >= lim:
            return name
    return "kamínek"


SCENES = [
    "Precese perihelia",
    "Sluneční soustava (zveličeno)",
    "Ohyb světla",
    "Shapirovo zpoždění",
    "Dvojhvězda a planeta",
    "Hodiny v poli",
    "Sluneční soustava (skutečná)",
    "Povrch Země: jablko a hodiny",
    "Země: ISS a GPS",
    "UAP: manévr (hypotéza)",
    "UAP: řízení (klávesnice)",
    "Černá díra (silné pole, neověřeno)",
    "α Centauri A/B",
    "Sirius A/B (bílý trpaslík)",
    "Okolí Slunce do 12 ly (3D)",
    "Volná scéna",
]

STATUS = [
    (
        "Newtonův zákon, Gaussův tok zdroje",
        "v laboratoři odvozeno z diskrétního operátoru a otestováno (F1.3, K.12); "
        "v demonstraci VLOŽENO jako vstup (L = −ΣM/r, hodnoty GM zadané), neodvozuje se",
    ),
    ("Vazba κ ↔ G", "vstup (kalibrace), neodvozeno"),
    (
        "Pohyb po geodetice pole",
        "otestováno ve 2+1 (F2.7: pohyb přesně na geodetice pole kauzální množiny)",
    ),
    (
        "Diskrétní operátor → spojitá rovnice",
        "otestováno: K.11, K.13a; K.12 bez stínění (r(0) = 1)",
    ),
    (
        "β = 1 (precese perihelia)",
        "spojitě odvozeno (F1.7); diskrétně F1.9a: měřidlo selhalo v kontrole kvadratury při σ = 4, "
        "popisně c → spojitá hodnota (β ≈ 1 při σ = 16), bez verdiktu",
    ),
    (
        "γ = 1 (ohyb světla, Shapiro)",
        "VSTUP: hypotéza H1, ne odvození; pravidlo F3 dává γ = −1 (otevřené, F7)",
    ),
    ("Exponenciální tvar metriky za 1PN", "volba demonstrace, neodvozeno"),
    (
        "Hmota jako bodové částice",
        "vstup demonstrace; vznik hmoty z dynamiky zatím nevyřešen",
    ),
    ("Superpozice polí", "lineární; nelineární oprava druhého řádu (F1.4) zanedbána"),
    ("Pole kvazistatické", "zpoždění šíření zanedbáno (pomalé zdroje)"),
    (
        "Vzájemná dynamika více těles",
        "každé těleso jako zkušební v poli ostatních (přiblížení)",
    ),
    ("a₀, ploché rotační křivky", "není v modelu (neodvozeno)"),
    (
        "Země: pád, hodiny, GPS",
        "slabé pole: 1PN metrika dává totéž co GR a pozorování (Pound–Rebka, GPS); nic nového",
    ),
    (
        "UAP „null field“",
        "HYPOTÉZA autora, není v teorii ani v laboratoři; scéna jen ukazuje, co by musela vysvětlit",
    ),
    (
        "Černá díra: bez horizontu, stín +0,91 %",
        "spojitá limita (přesné vakuové řešení rovnice laboratoře); diskrétně "
        "NEOVĚŘENO; předpokládá H1 i v silném poli; vnitřek (hmota) neznámý",
    ),
]

ABOUT = """Rovnice pole laboratoře: B_g ln Π = κρ_s (rozhodnutí 1),
geometrie H1: ds² = −e^{2L} dt² + e^{−2L} dx²,  L = ln Π.

Mimo zdroj (druhý řád) jen přejmenuje hmotu, takže pole
bodových hmot je L = −Σ M/r (G = c = 1).

Pohyb (geodetika, souřadnicový čas):
  dv/dt = −e^{4L}∇L − v²∇L + 4v(v·∇L)
Newtonova limita dv/dt = −∇L; světlo |v| = e^{2L},
příčné zrychlení −2∇L (dvojnásobek Newtona, γ = 1).

Vlastní čas: dτ = √(e^{2L} − e^{−2L}v²) dt.

Měření: precese ↔ (2+2γ−β)/3 · 6πM/(a(1−e²)),
ohyb ↔ (1+γ)(M/b)(sin θ₂ − sin θ₁),
Shapiro ↔ (1+γ) M ln((x₂+r₂)/(x₁+r₁)).

Jde o demonstraci spojité limity, ne o výzkumný přístroj.
V 1PN řádu (β = γ = 1) se chová jako obecná relativita;
odlišnost teorie by se ukázala až v silném poli."""


# =============================================================================
# Propojení s laboratoří: záznamy, nástěnka predikcí a spojitá reference
# =============================================================================


def default_results_dir():
    """results/ of the laboratory next to demo/ (D:\\Projekty\\causality_lab\\results when the demo lives in demo/)."""
    import os

    return os.path.normpath(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
    )


def load_records(folder, include_quick=False):
    """Načte uložené záznamy testů ve formátu JSON od nejnovějšího. Poškozené soubory přeskočí."""
    import json
    import os

    recs = []
    if not os.path.isdir(folder):
        return recs
    for tid in sorted(os.listdir(folder)):
        sub = os.path.join(folder, tid)
        if not os.path.isdir(sub):
            continue
        paths = [os.path.join(sub, f) for f in os.listdir(sub) if f.endswith(".json")]
        if include_quick and os.path.isdir(os.path.join(sub, "quick")):
            paths += [
                os.path.join(sub, "quick", f)
                for f in os.listdir(os.path.join(sub, "quick"))
                if f.endswith(".json")
            ]
        for pth in paths:
            try:
                with open(pth, encoding="utf-8") as fh:
                    rec = json.load(fh)
            except (OSError, ValueError):
                continue
            rec["_path"] = pth
            recs.append(rec)
    recs.sort(key=lambda r: r.get("time_utc", ""), reverse=True)
    return recs


# Riziková tvrzení teorie a testy, které je měří. Položka „reading“ obsahuje interpretaci uvedenou v
# docs/nalezy.md včetně data; stavy testů se načítají přímo z uložených záznamů.
PREDICTIONS = [
    (
        "Daleké pole zdroje je newtonovské (Gaussův tok)",
        ["F1.3"],
        "jiný exponent nebo síla pole",
        "drží: pole 1/r s tokem podle náboje (F1.3)",
    ),
    (
        "Diskrétní operátor → □_g − R/2 (i v H1)",
        ["K.11", "K.11c", "K.13a"],
        "nekonvergence ke spojitému tvaru",
        "drží v prvním řádu (K.11, K.13a, 0,8–2,4 % při σ = 32)",
    ),
    (
        "Operátor nestíní pole (2+1)",
        ["K.12"],
        "dno nezávislé na σ po odečtení přístroje",
        "24.–25. 10.: stínění vyloučeno analyticky (r(0) = 1); h16 potvrdil rozklad; nelineární zbytek ~0,06 otevřený",
    ),
    (
        "Pohyb po geodetice pole kauzální množiny",
        ["F2.7"],
        "dráha mimo geodetiku",
        "pohyb přesně na geodetice pole (2,027 proti 2,023); srovnání se spojitým polem čeká na kalibraci K.12",
    ),
    (
        "β = 1 v diskrétním poli",
        ["F1.9a"],
        "c mimo 0 o víc než 1/4 (β = 1/2 nebo 3/2)",
        "F1.9a MĚŘIDLO SELHALO (kvadratura σ = 4); popisně c(σ = 16) = −0,002 a +0,055 proti spojitým +0,008 a +0,058",
    ),
    (
        "γ = 1: ohyb světla ze substrátu",
        ["F7.1", "F7.4"],
        "γ ≠ 1 (Cassini: |γ − 1| < 2,3·10⁻⁵)",
        "otevřené: pravidlo F3 dává γ = −1; H1 je zatím vstup",
    ),
    (
        "Hmota vzniká z dynamiky kauzálního řádu",
        [],
        "žádná stabilní struktura",
        "neměřeno: největší mezera",
    ),
    (
        "a₀ = cH₀/2π z dynamiky (ploché rotační křivky)",
        ["F6.7"],
        "jiná škála nebo vývoj a₀(z) proti datům",
        "neodvozeno (rozbor F6.7, 23. 10.)",
    ),
    (
        "Gravitační vlny: tenzorové módy, rychlost c, kvadrupólový vzorec",
        [],
        "skalární/dipólové záření (pulsar Hulse–Taylor, LIGO)",
        "neměřeno: teorie zatím nemá dynamický tenzorový sektor",
    ),
    (
        "Silné pole: horizont, fotonová sféra, stín",
        [],
        "stín proti EHT (M87*, Sgr A*)",
        "vakuové řešení spojité rovnice (demo): žádný horizont, fotonová sféra 3,07 M, stín +0,91 % proti GR; diskrétně "
        "neměřeno; spojitá rovnice má meze existence řešení se zdrojem (F1.7 N)",
    ),
]


def continuum_a_2p1(R, sigma, q, n_grid=20001):
    """Spočítá spojité řešení rovnice (a) v rozměrech 2+1: ΔΠ = (1/3)Π²ρ_s s podmínkou Π(R) = 1 a κ = 1.
    Jde o nezávislou referenci k profilům K.12/F2.7 z laboratoře, řešenou relaxací s postupným zesilováním zdroje."""
    r = np.linspace(0.0, R, n_grid)
    rs = q * np.exp(-(r**2) / (2 * sigma**2))

    def poisson(f):  # ΔΨ = f v radiální symetrii; Ψ'(0) = 0 a Ψ(R) = 0
        dens = f * 2 * math.pi * r
        Q = np.concatenate(
            [[0.0], np.cumsum(0.5 * (dens[1:] + dens[:-1]) * np.diff(r))]
        )
        du = np.zeros_like(r)
        du[1:] = Q[1:] / (2 * math.pi * r[1:])
        U = np.concatenate([[0.0], np.cumsum(0.5 * (du[1:] + du[:-1]) * np.diff(r))])
        return U - U[-1]

    Pi = np.ones_like(r)
    for t in (0.125, 0.25, 0.5, 0.75, 1.0):
        for _ in range(3000):
            new = 1.0 + poisson(t * Pi**2 * rs / 3.0)
            step = float(np.max(np.abs(new - Pi)))
            Pi = 0.5 * Pi + 0.5 * new
            if step < 1e-11:
                break
    return r, Pi


# =============================================================================
# Autotest kódu demonstrace bez otevření okna
# =============================================================================


def selftest():
    """Ověří numerickou implementaci vzorců; nejde o důkaz fyzikální teorie."""
    ok_all = True

    def report(name, val, ref, tol):
        nonlocal ok_all
        rel = abs(val - ref) / abs(ref)
        ok = rel <= tol
        ok_all &= ok
        stav = "V POŘÁDKU" if ok else "CHYBA"
        print(
            f"  {stav:10s} {name}: {val:.6g}; očekáváno {ref:.6g} "
            f"(relativní odchylka {rel:.2e}, tolerance {tol:g})"
        )

    print("Autotesty gravitace_demo.py (jednotky G = c = 1)")
    # 1. Newtonovská limita: perioda téměř kruhové dráhy při malé hmotnosti.
    s = Sim()
    M = 1e-5
    s.bodies = [
        Body("c", 0, 0, mass=M, fixed=True),
        Body("p", 1.0, 0, 0, math.sqrt(M), mass=0.0),
    ]
    tr = PerihelionTracker(s.bodies[1], s.bodies[0])
    crossings, prev_y = [], 0.0
    p = s.bodies[1]

    def on(sm):
        nonlocal prev_y
        if prev_y < 0 <= p.y and p.x > 0:
            crossings.append(sm.t)
        prev_y = p.y

    s.advance(2.2 * 2 * math.pi / math.sqrt(M), max_steps=200000, on_step=on)
    if len(crossings) >= 2:
        report(
            "Keplerova perioda (M = 1e-5, r = 1)",
            crossings[1] - crossings[0],
            2 * math.pi / math.sqrt(M),
            1e-3,
        )
    else:
        print("  CHYBA Keplerovy periody: oběh nebyl dokončen")
        ok_all = False
    # 2. Precese perihelia proti vztahu 6πM/[a(1−e²)]; korekce vyššího řádu jsou malé.
    s = Sim()
    M, a, e = 0.002, 1.0, 0.5
    s.bodies = [
        Body("c", 0, 0, mass=M, fixed=True),
        Body("p", a * (1 - e), 0, 0, kepler_peri_velocity(M, a, e)),
    ]
    tr = PerihelionTracker(s.bodies[1], s.bodies[0])
    s.advance(
        5.3 * 2 * math.pi * math.sqrt(a**3 / M), max_steps=400000, on_step=tr.update
    )
    res = tr.result()
    if res:
        report(f"precese perihélia na oběh ({res[4]} oběhů)", res[0], res[1], 0.03)
    else:
        print("  CHYBA precese perihélia: nebylo nalezeno perihélium")
        ok_all = False
    # 3.–4. Ohyb světla a Shapirovo zpoždění proti analytickým vztahům pro γ=1.
    s = Sim()
    s.bodies = [
        Body("m", 0, 0, mass=1e-3, fixed=True),
        Body("ph", -50.0, 0.2, 1.0, 0.0, kind="photon"),
    ]
    # Předpověď: první řád na konečné přímé dráze a druhý řád exponenciální metriky.
    lt = LightTracker(s.bodies[1], s.bodies[0], s)
    s.advance(100.0, max_steps=400000, on_step=lt.update)
    st = lt.state(s)
    report(
        "odklon světla (M = 1e-3, b = 0,2)",
        st["defl"],
        st["pred"][1.0] + st["pred2"],
        2e-3,
    )
    report("Shapirovo zpoždění", st["delay"], st["delay_pred"][1.0], 0.03)
    # 5. Nulová podmínka fotonu a rychlost statických hodin.
    ph = s.bodies[1]
    Lp = s.L_at(ph.x, ph.y)
    report("rychlost fotonu |v| = e^(2L)", math.hypot(ph.vx, ph.vy), math.exp(2 * Lp), 1e-9)
    s = Sim()
    s.bodies = [
        Body("m", 0, 0, mass=0.02, fixed=True),
        Body("clk", 0.5, 0, kind="clock"),
    ]
    s.advance(10.0)
    report(
        "rychlost statických hodin dτ/dt = e^L",
        s.bodies[1].tau / s.t,
        math.exp(-0.02 / math.hypot(0.5, EPS)),
        1e-9,
    )
    # 6. Země: doba pádu z výšky 10 m proti vztahu sqrt(2h/g).
    s = Sim()
    s.dt_max = 0.002 * C_SI
    h = 10.0
    s.bodies = [
        Body("E", 0, 0, mass=M_EARTH, fixed=True, radius=R_EARTH),
        Body("a", 0, R_EARTH + h, mass=1e-28),
    ]
    while s.bodies[1].landed_t is None and s.t < 5 * C_SI:
        s.step()
    g = M_EARTH / (R_EARTH + h) ** 2 * C_SI**2
    report(
        "doba pádu z výšky 10 m [s]",
        s.bodies[1].landed_t / C_SI,
        math.sqrt(2 * h / g),
        2e-3,
    )
    # 7. GPS: rozdíl hodin na oběžné dráze a na povrchu za jeden den.
    s = Sim()
    s.dt_max = 10.0 * C_SI
    r = 2.6560e7
    s.bodies = [
        Body("E", 0, 0, mass=M_EARTH, fixed=True, radius=R_EARTH),
        Body("g", 0, R_EARTH + 0.01, kind="clock"),
        Body("GPS", r, 0, 0, math.sqrt(M_EARTH / r), mass=1e-24),
    ]
    s.advance(86400 * C_SI, max_steps=100000)
    d_meas = (s.bodies[1].tau_def - s.bodies[2].tau_def) / s.t * 86400e6
    d_pred = (M_EARTH / (R_EARTH + 0.01) - M_EARTH / r - 0.5 * M_EARTH / r) * 86400e6
    report("rozdíl hodin GPS a povrchu [μs/den]", d_meas, d_pred, 1e-2)
    # 8. Země: hodiny ve výšce 30 m a očekávaný gravitační rudý posuv.
    s = Sim()
    s.dt_max = 0.002 * C_SI
    s.bodies = [
        Body("E", 0, 0, mass=M_EARTH, fixed=True, radius=R_EARTH),
        Body("g", 0, R_EARTH + 0.01, kind="clock"),
        Body("t", 0, R_EARTH + 30.0, kind="clock"),
    ]
    s.advance(1.0 * C_SI, max_steps=10000)
    report(
        "relativní rychlost hodin ve výšce 30 m",
        (s.bodies[1].tau_def - s.bodies[2].tau_def) / s.t,
        M_EARTH / (R_EARTH + 0.01) - M_EARTH / (R_EARTH + 30.0),
        1e-3,
    )
    # 9. Merkur: relativistická precese pro jednotlivé modely bez planetárních poruch.
    gen = mercury_experiment(
        [_law_newton(), _law_h1(), _law_gr(), _law_ppn(0.5, 1.0)], 3.0, False
    )
    try:
        while True:
            next(gen)
    except StopIteration as stop:
        res = stop.value
    report("Merkur, Newton bez planet [''/století]", abs(res["N"]) + 1.0, 1.0, 0.05)
    report("Merkur, H1 [''/století]", res["H1"], 42.98, 0.01)
    report("Merkur, GR [''/století]", res["GR"], 42.98, 0.01)
    report("Merkur, PPN β = 1/2 [''/století]", res["PPN"], 42.98 * 3.5 / 3.0, 0.01)
    # 10. Černá díra: porovnání zachycení fotonu s parametrem b=5,23 M a kritického dopadu.
    # V H1 je foton zachycen; v Schwarzschildově řešení unikne.
    caught = {}
    for law in (_law_h1(), _law_gr()):
        B, _, m = scene("Černá díra (silné pole, neověřeno)")
        s = Sim(law)
        s.bodies = [B[0]] + [b for b in B if b.name == "foton b=5.23"]
        s.dt_max, s.eta, s.no_guard = m["dt_max"], m["eta"], True
        if law.key == "H1":
            s.profile = B[0]
        ph = s.bodies[1]
        while s.t < 250 and math.hypot(ph.x, ph.y) >= 1.0:
            s.step()
        caught[law.key] = math.hypot(ph.x, ph.y) < 1.0
    report(
        "černá díra: b = 5,23, zachycení v H1 a únik v GR",
        1.0 if (caught["H1"] and not caught["GR"]) else 0.0,
        1.0,
        1e-9,
    )
    report(
        "černá díra: kritický parametr dopadu (H1, M)",
        bh_profile()[3]["b_crit"],
        5.2435,
        1e-3,
    )
    print("VÝSLEDEK:", "všechny kontroly prošly" if ok_all else "některé kontroly selhaly")
    return 0 if ok_all else 1


# =============================================================================
# Aplikace a její grafické rozhraní Tkinter
# =============================================================================


def run_app(start_scene=None, autorun=False):
    import tkinter as tk
    from tkinter import ttk

    root = tk.Tk()
    root.title("Gravitace kauzální fyziky – interaktivní demonstrace")
    root.geometry("1360x860")
    root.configure(bg="#14161c")

    state = {
        "sim": Sim(),
        "view": (0.0, 0.0, 3.0),
        "running": False,
        "trackers": [],
        "light": [],
        "ghosts": [],
        "exp": None,
        "exp_result": None,
        "exp_status": "",
        "hover": None,
        "drag": None,
        "pan": None,
        "frame": 0,
        "field_img": None,
        "scene": SCENES[0],
        "meta": {"si": False, "display": ("log", -4.0, -1.0), "click": "orbit"},
        "name_seq": 0,
        "mass_reference": None,
        "view3d": False,
        "readout_hold": False,
        "readout_text": "",
        "readout_pending": "",
    }

    # ---------- layout
    left = tk.Frame(root, bg="#1c1f27", width=330)
    left.pack(side="left", fill="y")
    left.pack_propagate(False)
    right = tk.Frame(root, bg="#14161c")
    right.pack(side="left", fill="both", expand=True)
    canvas = tk.Canvas(right, bg="#0b0d12", highlightthickness=0)
    canvas.pack(fill="both", expand=True)
    out_frame = tk.Frame(right, bg="#10131a")
    out_frame.pack(fill="x")
    out_toolbar = tk.Frame(out_frame, bg="#1c1f27")
    out_toolbar.pack(side="top", fill="x")
    out_status = tk.Label(out_toolbar, text="Výpis měření", bg="#1c1f27", fg="#9aa0ad", anchor="w")
    out_status.pack(side="left", fill="x", expand=True, padx=6)
    out_text_frame = tk.Frame(out_frame, bg="#10131a")
    out_text_frame.pack(side="top", fill="both", expand=True)
    out_sb = tk.Scrollbar(out_text_frame, orient="vertical")
    out_xsb = tk.Scrollbar(out_text_frame, orient="horizontal")
    out_sb.pack(side="right", fill="y")
    out_xsb.pack(side="bottom", fill="x")
    out = tk.Text(
        out_text_frame,
        height=12,
        bg="#10131a",
        fg="#d8dce6",
        font=("Consolas", 10),
        relief="flat",
        wrap="none",
        yscrollcommand=out_sb.set,
        xscrollcommand=out_xsb.set,
        undo=False,
        exportselection=True,
    )
    out.pack(side="left", fill="both", expand=True)
    out_sb.config(command=out.yview)
    out_xsb.config(command=out.xview)
    err_lbl = tk.Label(
        right,
        text="",
        bg="#14161c",
        fg="#ff6b6b",
        anchor="w",
        justify="left",
        wraplength=1000,
    )
    err_lbl.pack(fill="x")
    seen_errors = set()

    def report_error(where, exc):
        """Zobrazí chybu pod výpisem, ale ponechá simulaci spuštěnou; trasování vypíše jen jednou."""
        import traceback

        err_lbl.config(
            text=f"Chyba ({where}): {type(exc).__name__}: {exc}  (simulace běží dál)"
        )
        key = (where, type(exc).__name__, str(exc))
        if key not in seen_errors:
            seen_errors.add(key)
            traceback.print_exc()

    def safe(where, fn):
        """Obalí obsluhu události tak, aby chyba nezastavila celé GUI."""

        def wrapped(*a, **kw):
            try:
                return fn(*a, **kw)
            except (
                Exception
            ) as exc:  # noqa: BLE001 -- shown to the user, the app goes on
                report_error(where, exc)

        return wrapped

    def _apply_readout(text):
        """Aktualizuje výpis bez zbytečného přepisování nezměněného obsahu."""
        text = str(text)
        if text == state.get("readout_text", ""):
            return
        y_top = out.yview()[0]
        byl_dole = out.yview()[1] >= 0.98
        out.configure(state="normal")
        out.delete("1.0", "end")
        out.insert("1.0", text)
        out.configure(state="disabled")
        state["readout_text"] = text
        if byl_dole:
            out.see("end")
        else:
            out.yview_moveto(y_top)

    def set_readout(text):
        """Aktualizuje data, ale během výběru textu nechá obsah stabilní pro kopírování."""
        text = str(text)
        state["readout_pending"] = text
        if state.get("readout_hold"):
            return
        _apply_readout(text)

    def copy_text_widget(widget, only_selection=False, status=None):
        """Zkopíruje vybraný text nebo celý obsah textového pole do systémové schránky."""
        try:
            bez_vyberu = False
            if only_selection:
                try:
                    content = widget.get("sel.first", "sel.last")
                except Exception:
                    content = widget.get("1.0", "end-1c")
                    bez_vyberu = True
            else:
                content = widget.get("1.0", "end-1c")
            widget.winfo_toplevel().clipboard_clear()
            widget.winfo_toplevel().clipboard_append(content)
            widget.winfo_toplevel().update_idletasks()
            if status is not None:
                if bez_vyberu:
                    status.config(text=f"Není vybrán text; zkopírován celý výpis ({len(content)} znaků).")
                elif content:
                    status.config(text=f"Zkopírováno {len(content)} znaků do schránky.")
        except Exception as exc:
            if status is not None:
                status.config(text=f"Kopírování selhalo: {exc}")
        return "break"

    def select_all_text(widget):
        """Vybere celý obsah textového pole."""
        widget.tag_add("sel", "1.0", "end-1c")
        widget.mark_set("insert", "1.0")
        widget.see("insert")
        return "break"

    def bind_text_copy(widget, status=None):
        """Přidá běžné klávesové zkratky a kontextovou nabídku pro kopírování."""
        widget.bind("<Control-a>", lambda _e, w=widget: select_all_text(w))
        widget.bind("<Control-A>", lambda _e, w=widget: select_all_text(w))
        widget.bind("<Control-c>", lambda _e, w=widget, st=status: copy_text_widget(w, True, st))
        widget.bind("<Control-C>", lambda _e, w=widget, st=status: copy_text_widget(w, True, st))
        widget.bind("<Control-Insert>", lambda _e, w=widget, st=status: copy_text_widget(w, True, st))
        menu = tk.Menu(widget, tearoff=0)
        menu.add_command(label="Kopírovat výběr", command=lambda w=widget, st=status: copy_text_widget(w, True, st))
        menu.add_command(label="Kopírovat vše", command=lambda w=widget, st=status: copy_text_widget(w, False, st))
        menu.add_command(label="Vybrat vše", command=lambda w=widget: select_all_text(w))
        def show_menu(event, m=menu):
            try:
                m.tk_popup(event.x_root, event.y_root)
            finally:
                m.grab_release()
            return "break"
        widget.bind("<Button-3>", show_menu)
        widget.bind("<Button-2>", show_menu)

    def _uvolni_vypis():
        """Po opuštění výpisu obnoví živá data, ne během označování či kopírování."""
        try:
            focus = root.focus_get()
        except tk.TclError:
            focus = None
        if focus is out:
            return
        state["readout_hold"] = False
        if state.get("readout_pending") != state.get("readout_text"):
            _apply_readout(state.get("readout_pending", ""))

    def _vypis_focus_in(_event=None):
        state["readout_hold"] = True

    def _vypis_focus_out(_event=None):
        # Krátká prodleva umožní použít tlačítko „Kopírovat výběr“ bez závodu s obnovou.
        root.after(350, _uvolni_vypis)

    out.bind("<FocusIn>", _vypis_focus_in)
    out.bind("<FocusOut>", _vypis_focus_out)
    bind_text_copy(out, out_status)
    tk.Button(out_toolbar, text="Kopírovat vše", command=lambda: copy_text_widget(out, False, out_status)).pack(side="left", padx=2, pady=2)
    tk.Button(out_toolbar, text="Kopírovat výběr", command=lambda: copy_text_widget(out, True, out_status)).pack(side="left", padx=2, pady=2)
    tk.Button(out_toolbar, text="Vybrat vše", command=lambda: select_all_text(out)).pack(side="left", padx=2, pady=2)
    tk.Button(out_toolbar, text="Obnovit výpis", command=lambda: (state.update(readout_hold=False), readouts())).pack(side="left", padx=2, pady=2)

    nb = ttk.Notebook(left)
    nb.pack(fill="both", expand=True, padx=6, pady=6)

    def vytvor_posuvnou_kartu(nazev):
        """Vytvoří kartu svisle posuvnou i v menším okně."""
        obal = tk.Frame(nb, bg="#1c1f27")
        nb.add(obal, text=nazev)
        plocha = tk.Canvas(obal, bg="#1c1f27", highlightthickness=0, bd=0)
        posuvnik = ttk.Scrollbar(obal, orient="vertical", command=plocha.yview)
        obsah = tk.Frame(plocha, bg="#1c1f27")
        okno = plocha.create_window((0, 0), window=obsah, anchor="nw")
        plocha.configure(yscrollcommand=posuvnik.set)
        plocha.pack(side="left", fill="both", expand=True)
        posuvnik.pack(side="right", fill="y")
        obsah.bind("<Configure>", lambda _e, c=plocha: c.configure(scrollregion=c.bbox("all")))
        plocha.bind("<Configure>", lambda e, c=plocha, w=okno: c.itemconfigure(w, width=e.width))
        def posun_karty(event, c=plocha, o=obal):
            widget = event.widget
            while widget is not None:
                if widget is o:
                    delta = getattr(event, "delta", 0)
                    if delta:
                        c.yview_scroll(-int(delta / 120) or (-1 if delta > 0 else 1), "units")
                    elif getattr(event, "num", None) == 4:
                        c.yview_scroll(-1, "units")
                    elif getattr(event, "num", None) == 5:
                        c.yview_scroll(1, "units")
                    return "break"
                widget = getattr(widget, "master", None)
            return None
        obal.bind_all("<MouseWheel>", posun_karty, add="+")
        obal.bind_all("<Button-4>", posun_karty, add="+")
        obal.bind_all("<Button-5>", posun_karty, add="+")
        return obal, obsah

    card_c, tab_c = vytvor_posuvnou_kartu("Ovládání")
    card_g, tab_g = vytvor_posuvnou_kartu("Srovnání")
    card_s, tab_s = vytvor_posuvnou_kartu("Stav teorie")
    card_a, tab_a = vytvor_posuvnou_kartu("O modelu")

    def label(parent, text, **kw):
        return tk.Label(
            parent,
            text=text,
            bg="#1c1f27",
            fg="#d8dce6",
            anchor="w",
            justify="left",
            **kw,
        )

    label(tab_c, "Scéna").pack(fill="x", padx=8, pady=(8, 0))
    scene_var = tk.StringVar(value=SCENES[0])
    ttk.Combobox(tab_c, textvariable=scene_var, values=SCENES, state="readonly").pack(
        fill="x", padx=8
    )

    btns = tk.Frame(tab_c, bg="#1c1f27")
    btns.pack(fill="x", padx=8, pady=8)
    run_btn = tk.Button(btns, text="Spustit", width=9)
    run_btn.pack(side="left")
    tk.Button(btns, text="Krok", width=6, command=lambda: tick(single=True)).pack(
        side="left", padx=4
    )
    tk.Button(
        btns, text="Znovu", width=7, command=lambda: load_scene(scene_var.get())
    ).pack(side="left")

    speed_lbl = label(tab_c, "Tempo simulace (log₁₀ simulačních jednotek za sekundu reálného času)")
    speed_lbl.pack(fill="x", padx=8)
    speed_var = tk.DoubleVar(value=0.0)
    speed_scale = tk.Scale(
        tab_c,
        from_=-3.0,
        to=12.5,
        resolution=0.05,
        orient="horizontal",
        variable=speed_var,
        bg="#1c1f27",
        fg="#d8dce6",
        highlightthickness=0,
    )
    speed_scale.pack(fill="x", padx=8)
    speed_val = tk.Label(tab_c, text="", bg="#1c1f27", fg="#9ad9ff", anchor="w")
    speed_val.pack(fill="x", padx=8)

    def speed_value_label(*_):
        """Zobrazuje stejné ovládací tempo ve všech scénách, nikoli tempo závislé na FPS."""
        tempo = 10.0 ** speed_var.get()
        if state["meta"].get("si"):
            txt = f"1 s reálného času = {fmt_time(tempo)} simulace"
            if state.get("achieved") is not None:
                txt += f"  (dosaženo {fmt_time(state['achieved'])}/s)"
        else:
            txt = f"{tempo:.4g} jednotek modelového času za sekundu reálného času"
            if state.get("achieved") is not None:
                txt += f"  (dosaženo {state['achieved']:.4g} jednotek/s)"
        speed_val.config(text=txt)

    speed_var.trace_add("write", speed_value_label)

    trails_var, field_var, names_var = (
        tk.BooleanVar(value=True),
        tk.BooleanVar(value=True),
        tk.BooleanVar(value=True),
    )
    for txt, var in (
        ("Stopy", trails_var),
        ("Pole Π (pevná škála scény, legenda vlevo dole)", field_var),
        ("Popisky", names_var),
    ):
        tk.Checkbutton(
            tab_c,
            text=txt,
            variable=var,
            bg="#1c1f27",
            fg="#d8dce6",
            selectcolor="#2b2f3a",
            activebackground="#1c1f27",
            anchor="w",
        ).pack(fill="x", padx=8)

    projection_var = tk.BooleanVar(value=False)
    projection_check = tk.Checkbutton(
        tab_c,
        text="3D náhled rovinné scény (pouze projekce)",
        variable=projection_var,
        command=lambda: change_projection(),
        bg="#1c1f27",
        fg="#d8dce6",
        selectcolor="#2b2f3a",
        activebackground="#1c1f27",
        anchor="w",
        wraplength=285,
        justify="left",
    )
    projection_check.pack(fill="x", padx=8, pady=(4, 0))

    label(tab_c, "Myš: kolečko = zoom, pravé tlačítko = posun").pack(
        fill="x", padx=8, pady=(10, 0)
    )
    click_lbl = label(tab_c, "")
    click_lbl.pack(fill="x", padx=8)
    tool_var = tk.StringVar(value="none")
    tool_widgets = {}
    for txt, val in (
        ("nic", "none"),
        ("přidat těleso", "mass"),
        ("přidat foton", "photon"),
        ("přidat hodiny", "clock"),
        ("smazat (klik na objekt)", "delete"),
    ):
        rb = tk.Radiobutton(
            tab_c,
            text=txt,
            value=val,
            variable=tool_var,
            bg="#1c1f27",
            fg="#d8dce6",
            selectcolor="#2b2f3a",
            activebackground="#1c1f27",
            anchor="w",
        )
        rb.pack(fill="x", padx=16)
        tool_widgets[val] = rb
    mass_lbl = label(tab_c, "Poměr hmotnosti nového tělesa m/Mref")
    mass_lbl.pack(fill="x", padx=8, pady=(6, 0))
    mass_var = tk.StringVar(value="1e-6")
    tk.Entry(tab_c, textvariable=mass_var).pack(fill="x", padx=8)
    cat_lbl = (
        label(tab_c, "", fg="#9ad9ff")
        if False
        else tk.Label(tab_c, text="", bg="#1c1f27", fg="#9ad9ff", anchor="w")
    )
    cat_lbl.pack(fill="x", padx=8)
    pilot_frame = tk.Frame(tab_c, bg="#1c1f27")
    pilot_frame.pack(fill="x", padx=8, pady=(6, 0))
    tk.Label(
        pilot_frame,
        text="Řízení UAP: šipky/WASD, Shift = ×10, mezerník zapíná/vypíná\nvznášení (1 g). Posuvník nastavuje tah [log₁₀ g] nebo\nrychlost v hypotéze „nulového pole“ [log₁₀(m/s) − 2]:",
        bg="#1c1f27",
        fg="#d8dce6",
        anchor="w",
        justify="left",
    ).pack(fill="x")
    thrust_var = tk.DoubleVar(value=0.5)
    tk.Scale(
        pilot_frame,
        from_=-1.0,
        to=3.0,
        resolution=0.05,
        orient="horizontal",
        variable=thrust_var,
        bg="#1c1f27",
        fg="#d8dce6",
        highlightthickness=0,
    ).pack(fill="x")
    nullf_var, follow_var = tk.BooleanVar(value=False), tk.BooleanVar(value=True)
    tk.Checkbutton(
        pilot_frame,
        text="Hypotéza „nulového pole“: řízení rychlosti bez setrvačnosti, posádka 0 g",
        variable=nullf_var,
        bg="#1c1f27",
        fg="#d8dce6",
        selectcolor="#2b2f3a",
        activebackground="#1c1f27",
        anchor="w",
    ).pack(fill="x")
    tk.Checkbutton(
        pilot_frame,
        text="Kamera sleduje UAP",
        variable=follow_var,
        bg="#1c1f27",
        fg="#d8dce6",
        selectcolor="#2b2f3a",
        activebackground="#1c1f27",
        anchor="w",
    ).pack(fill="x")
    pilot_frame.pack_forget()

    # Formulář pro ruční vložení objektu do skutečné prostorové mapy hvězd.
    map3d_frame = tk.LabelFrame(
        tab_c, text="Nový objekt ve 3D mapě", bg="#1c1f27", fg="#d8dce6",
        bd=1, padx=6, pady=4,
    )
    map_name_var = tk.StringVar(value="Nové těleso")
    map_xyz_vars = [tk.StringVar(value=v) for v in ("0", "0", "0")]
    map_v_vars = [tk.StringVar(value=v) for v in ("0", "0", "0")]
    tk.Label(map3d_frame, text="Název", bg="#1c1f27", fg="#d8dce6", anchor="w").grid(row=0, column=0, columnspan=3, sticky="w")
    tk.Entry(map3d_frame, textvariable=map_name_var).grid(row=1, column=0, columnspan=3, sticky="ew", pady=(0, 4))
    for j, osa in enumerate(("X", "Y", "Z")):
        tk.Label(map3d_frame, text=f"{osa} [ly]", bg="#1c1f27", fg="#d8dce6").grid(row=2, column=j, sticky="w")
        tk.Entry(map3d_frame, textvariable=map_xyz_vars[j], width=8).grid(row=3, column=j, sticky="ew", padx=(0, 3))
    for j, osa in enumerate(("Vx", "Vy", "Vz")):
        tk.Label(map3d_frame, text=f"{osa} [km/s]", bg="#1c1f27", fg="#d8dce6").grid(row=4, column=j, sticky="w")
        tk.Entry(map3d_frame, textvariable=map_v_vars[j], width=8).grid(row=5, column=j, sticky="ew", padx=(0, 3))
    for j in range(3):
        map3d_frame.columnconfigure(j, weight=1)
    tk.Button(map3d_frame, text="Přidat objekt", command=lambda: add_map3d_object()).grid(row=6, column=0, columnspan=3, sticky="ew", pady=(5, 0))
    map3d_frame.pack_forget()
    keys = set()

    def parse_mass():
        """Vrátí kladný konečný poměr m/Mref; při chybném vstupu vrátí None."""
        try:
            pomer = float(mass_var.get().strip().replace(",", "."))
        except ValueError:
            return None
        return pomer if math.isfinite(pomer) and pomer > 0 else None

    def referencni_hmotnost(sim=None):
        """Vrátí hmotnost nejhmotnějšího živého zdroje scény v interních jednotkách."""
        sim = sim or state["sim"]
        if sim is state.get("sim") and state.get("mass_reference") is not None:
            return float(state["mass_reference"])
        hodnoty = [b.mass for b in sim.bodies if b.alive and b.mass > 0 and b is not sim.profile]
        return max(hodnoty, default=(M_SUN if state.get("meta", {}).get("map3d") else 1.0))

    def update_category(*_):
        pomer = parse_mass()
        if pomer is None:
            text = "Neplatný poměr: zadej kladné konečné číslo."
        else:
            hmotnost = pomer * referencni_hmotnost()
            meta = state.get("meta", {})
            if meta.get("map3d"):
                text = f"Odhad: {hmotnost / M_SUN:.4g} M☉"
            elif meta.get("si"):
                kg = hmotnost / G_C2
                text = f"Odhad: {category_si(kg)}, {kg:.4g} kg"
            else:
                text = f"Odhad: {category_astro(hmotnost)}, {hmotnost:.4g} jednotek modelu"
        cat_lbl.config(text=text)

    mass_var.trace_add("write", update_category)
    note_lbl = label(tab_c, "")
    note_lbl.pack(fill="x", padx=8)

    for i, (k, v) in enumerate(STATUS):
        col = (
            "#ffb37a"
            if v.startswith(("VSTUP", "HYPOTÉZA")) or "neodvozeno" in v or "vstup" in v
            else "#9ad99a"
        )
        label(tab_s, k, font=("Segoe UI", 9, "bold")).pack(
            fill="x", padx=8, pady=(6 if i else 10, 0)
        )
        tk.Label(
            tab_s,
            text=v,
            bg="#1c1f27",
            fg=col,
            anchor="w",
            justify="left",
            wraplength=300,
        ).pack(fill="x", padx=16)
    tk.Label(
        tab_a,
        text=ABOUT,
        bg="#1c1f27",
        fg="#d8dce6",
        anchor="nw",
        justify="left",
        font=("Consolas", 9),
    ).pack(fill="both", padx=8, pady=8)

    # ---------- Karta srovnávacích zákonů, posuvníky PPN a experiment s Merkurem
    tk.Label(
        tab_g,
        text="Fyzikální rámce: stejná scéna souběžně v jiném\nrámci (čárkované stopy; hlavní scéna používá H1). "
        "Zaškrtnutí scénu hned spustí znovu.",
        bg="#1c1f27",
        fg="#d8dce6",
        anchor="w",
        justify="left",
    ).pack(fill="x", padx=8, pady=(8, 4))
    ghost_vars = {}
    for key, txt in (
        ("GR", "obecná relativita (Schwarzschild)"),
        ("N", "Newton (světlo jako korpuskule)"),
        ("F3", "pravidlo F3 (Π²η, γ = −1)"),
        ("PPN", "PPN s posuvníky β, γ"),
    ):
        ghost_vars[key] = tk.BooleanVar(value=False)
        tk.Checkbutton(
            tab_g,
            text=txt,
            variable=ghost_vars[key],
            bg="#1c1f27",
            fg="#d8dce6",
            selectcolor="#2b2f3a",
            activebackground="#1c1f27",
            anchor="w",
        ).pack(fill="x", padx=8)
        ghost_vars[key].trace_add(
            "write", lambda *a: safe("rámce", load_scene)(scene_var.get())
        )
    beta_var, gamma_var = tk.DoubleVar(value=1.0), tk.DoubleVar(value=1.0)
    for txt, var, lo, hi in (
        ("β (nelinearita, GR = 1)", beta_var, 0.0, 2.0),
        ("γ (zakřivení prostoru, GR = 1)", gamma_var, -1.0, 2.0),
    ):
        tk.Label(tab_g, text=txt, bg="#1c1f27", fg="#d8dce6", anchor="w").pack(
            fill="x", padx=8
        )
        tk.Scale(
            tab_g,
            from_=lo,
            to=hi,
            resolution=0.05,
            orient="horizontal",
            variable=var,
            bg="#1c1f27",
            fg="#d8dce6",
            highlightthickness=0,
        ).pack(fill="x", padx=8)
    tk.Button(
        tab_g,
        text="Použít β, γ (scéna začne znovu)",
        command=lambda: load_scene(scene_var.get()),
    ).pack(fill="x", padx=8, pady=4)
    tk.Label(
        tab_g,
        text="\nExperiment Merkur: precese perihelia v ″/století\n(Newton, H1, GR a PPN z posuvníků; počítá se\n"
        "na pozadí, okno zůstává použitelné)",
        bg="#1c1f27",
        fg="#d8dce6",
        anchor="w",
        justify="left",
    ).pack(fill="x", padx=8)
    exp_frame = tk.Frame(tab_g, bg="#1c1f27")
    exp_frame.pack(fill="x", padx=8)
    tk.Label(exp_frame, text="let:", bg="#1c1f27", fg="#d8dce6").pack(side="left")
    years_var = tk.StringVar(value="10")
    tk.Entry(exp_frame, textvariable=years_var, width=5).pack(side="left", padx=4)
    planets_var = tk.BooleanVar(value=False)
    tk.Checkbutton(
        exp_frame,
        text="s planetami",
        variable=planets_var,
        bg="#1c1f27",
        fg="#d8dce6",
        selectcolor="#2b2f3a",
        activebackground="#1c1f27",
    ).pack(side="left")
    tk.Button(
        tab_g, text="Spustit experiment Merkur", command=lambda: start_experiment()
    ).pack(fill="x", padx=8, pady=4)
    exp_lbl = tk.Label(
        tab_g,
        text="",
        bg="#1c1f27",
        fg="#9ad9ff",
        anchor="w",
        justify="left",
        wraplength=300,
    )
    exp_lbl.pack(fill="x", padx=8)
    tk.Label(
        tab_g,
        text="\nMost k laboratoři: záznamy testů, nástěnka\npredikcí, profily polí proti spojitému řešení",
        bg="#1c1f27",
        fg="#d8dce6",
        anchor="w",
        justify="left",
    ).pack(fill="x", padx=8)
    tk.Button(
        tab_g,
        text="Otevřít laboratoř a nástěnku",
        command=lambda: safe("laboratoř", open_lab)(),
    ).pack(fill="x", padx=8, pady=4)

    def open_lab():
        """Otevře okno se záznamy laboratoře, nástěnkou predikcí a prohlížečem výsledků."""
        win = tk.Toplevel(root)
        win.title("Laboratoř: nástěnka predikcí a záznamy")
        win.geometry("1250x820")
        win.configure(bg="#14161c")
        top = tk.Frame(win, bg="#14161c")
        top.pack(fill="x", padx=6, pady=6)
        tk.Label(top, text="Složka results:", bg="#14161c", fg="#d8dce6").pack(
            side="left"
        )
        path_var = tk.StringVar(value=default_results_dir())
        tk.Entry(top, textvariable=path_var, width=70).pack(side="left", padx=4)
        quick_var = tk.BooleanVar(value=False)
        tk.Checkbutton(
            top,
            text="i rychlé běhy",
            variable=quick_var,
            bg="#14161c",
            fg="#d8dce6",
            selectcolor="#2b2f3a",
            activebackground="#14161c",
        ).pack(side="left")
        lnb = ttk.Notebook(win)
        lnb.pack(fill="both", expand=True, padx=6, pady=6)
        f_board, f_rec = tk.Frame(lnb, bg="#14161c"), tk.Frame(lnb, bg="#14161c")
        lnb.add(f_board, text="Nástěnka predikcí")
        lnb.add(f_rec, text="Záznamy a profily polí")
        # board
        bsb = tk.Scrollbar(f_board)
        bsb.pack(side="right", fill="y")
        board = tk.Text(
            f_board,
            bg="#10131a",
            fg="#d8dce6",
            font=("Segoe UI", 10),
            wrap="word",
            yscrollcommand=bsb.set,
        )
        bsb.config(command=board.yview)
        board_tools = tk.Frame(f_board, bg="#14161c")
        board_tools.pack(side="top", fill="x")
        board_status = tk.Label(board_tools, text="", bg="#14161c", fg="#9aa0ad", anchor="w")
        board_status.pack(side="right", padx=6)
        tk.Button(board_tools, text="Kopírovat vše", command=lambda: copy_text_widget(board, False, board_status)).pack(side="left", padx=3, pady=3)
        tk.Button(board_tools, text="Kopírovat výběr", command=lambda: copy_text_widget(board, True, board_status)).pack(side="left", padx=3, pady=3)
        tk.Button(board_tools, text="Vybrat vše", command=lambda: select_all_text(board)).pack(side="left", padx=3, pady=3)
        board.pack(fill="both", expand=True)
        bind_text_copy(board, board_status)
        for tag, col in (
            ("h", "#ffffff"),
            ("ok", "#9ad99a"),
            ("bad", "#ff8a7a"),
            ("inst", "#ffb37a"),
            ("none", "#9aa0ad"),
            ("txt", "#d8dce6"),
        ):
            board.tag_config(tag, foreground=col)
        board.tag_config("h", font=("Segoe UI", 11, "bold"))
        # records
        left_f = tk.Frame(f_rec, bg="#14161c")
        left_f.pack(side="left", fill="y")
        lb = tk.Listbox(
            left_f,
            width=58,
            bg="#10131a",
            fg="#d8dce6",
            font=("Consolas", 9),
            selectbackground="#2b4a6a",
        )
        lsb = tk.Scrollbar(left_f, command=lb.yview)
        lb.config(yscrollcommand=lsb.set)
        lsb.pack(side="right", fill="y")
        lb.pack(side="left", fill="y")
        right_f = tk.Frame(f_rec, bg="#14161c")
        right_f.pack(side="left", fill="both", expand=True)
        plot = tk.Canvas(right_f, bg="#0b0d12", height=330, highlightthickness=0)
        plot.pack(fill="x")
        record_tools = tk.Frame(right_f, bg="#14161c")
        record_tools.pack(side="top", fill="x")
        record_status = tk.Label(record_tools, text="", bg="#14161c", fg="#9aa0ad", anchor="w")
        record_status.pack(side="right", padx=6)
        tk.Button(record_tools, text="Kopírovat vše", command=lambda: copy_text_widget(txt, False, record_status)).pack(side="left", padx=3, pady=3)
        tk.Button(record_tools, text="Kopírovat výběr", command=lambda: copy_text_widget(txt, True, record_status)).pack(side="left", padx=3, pady=3)
        tk.Button(record_tools, text="Vybrat vše", command=lambda: select_all_text(txt)).pack(side="left", padx=3, pady=3)
        tsb = tk.Scrollbar(right_f)
        tsb.pack(side="right", fill="y")
        txt = tk.Text(
            right_f,
            bg="#10131a",
            fg="#d8dce6",
            font=("Consolas", 9),
            wrap="word",
            yscrollcommand=tsb.set,
        )
        txt.pack(fill="both", expand=True)
        tsb.config(command=txt.yview)
        bind_text_copy(txt, record_status)
        lab = {"recs": []}

        def status_tag(st):
            return {"PASSED": "ok", "FAILED": "bad", "INSTRUMENT_FAILED": "inst"}.get(
                st, "none"
            )

        def fill_board():
            board.delete("1.0", "end")
            recs = lab["recs"]
            board.insert(
                "end",
                f"Nástěnka predikcí (záznamy: {len(recs)} ze složky {path_var.get()})\n"
                "U každého tvrzení: čím by se vyvrátilo, poslední záznamy příslušných testů (automaticky) "
                "a výklad z docs/nalezy.md.\n\n",
                "txt",
            )
            for claim, tests, falsifier, reading in PREDICTIONS:
                board.insert("end", claim + "\n", "h")
                board.insert("end", f"   vyvrátilo by: {falsifier}\n", "txt")
                found = False
                for tid in tests:
                    latest = {}
                    for rec in recs:  # newest record per part
                        if rec.get("id") == tid and not rec.get("quick"):
                            latest.setdefault(rec.get("part") or "", rec)
                    for part, rec in sorted(latest.items()):
                        found = True
                        st = rec.get("status", "?")
                        board.insert("end", f"   {tid} {part or ''}: ", "txt")
                        board.insert("end", f"{st}", status_tag(st))
                        board.insert("end", f"  ({rec.get('time_utc', '')})\n", "none")
                if not found:
                    board.insert("end", "   záznamy: žádné (neměřeno)\n", "none")
                board.insert("end", f"   výklad: {reading}\n\n", "txt")

        def fill_list():
            lb.delete(0, "end")
            for rec in lab["recs"]:
                q = " quick" if rec.get("quick") else ""
                lb.insert(
                    "end",
                    f"{str(rec.get('id') or '?'):7s} {str(rec.get('part') or ''):10s} "
                    f"{str(rec.get('status') or '?'):17s} {rec.get('time_utc') or ''}{q}",
                )

        def load():
            lab["recs"] = load_records(path_var.get(), quick_var.get())
            fill_board()
            fill_list()

        def show(_ev=None):
            sel = lb.curselection()
            if not sel:
                return
            rec = lab["recs"][sel[0]]
            txt.delete("1.0", "end")
            txt.insert(
                "end",
                f"{rec.get('id')} {rec.get('part') or ''}: {rec.get('status')}  ({rec.get('_path')})\n\n",
            )
            txt.insert("end", "\n".join(rec.get("lines", [])))
            plot.delete("all")
            d, prm = rec.get("data") or {}, rec.get("params") or {}
            if d.get("nodes") and d.get("Pi_cs") and int(prm.get("D", 0)) == 3:
                draw_profile(rec)
            elif rec.get("id") == "F1.9a" and d.get("table"):
                draw_f19a(rec)
            else:
                plot.create_text(
                    10,
                    10,
                    anchor="nw",
                    fill="#9aa0ad",
                    text="Pro tento záznam není graf (graf umí: pole (a) ve 2+1 z K.12/F2.7 a tabulku F1.9a).",
                )

        def axes(xmin, xmax, ymin, ymax, xl, yl, xticks=None):
            w, h = max(plot.winfo_width(), 400), 330
            X = lambda x: 60 + (x - xmin) / (xmax - xmin) * (w - 80)
            Y = lambda y: h - 40 - (y - ymin) / (ymax - ymin) * (h - 70)
            plot.create_line(60, h - 40, w - 20, h - 40, fill="#9aa0ad")
            plot.create_line(60, 30, 60, h - 40, fill="#9aa0ad")
            for k in range(6):
                xv = xmin + (xmax - xmin) * k / 5
                yv = ymin + (ymax - ymin) * k / 5
                if xticks is None:
                    plot.create_text(
                        X(xv),
                        h - 28,
                        text=f"{xv:.3g}",
                        fill="#9aa0ad",
                        font=("Segoe UI", 8),
                    )
                plot.create_text(
                    52,
                    Y(yv),
                    text=f"{yv:.3g}",
                    fill="#9aa0ad",
                    anchor="e",
                    font=("Segoe UI", 8),
                )
            for xv, lab_ in xticks or []:
                plot.create_text(
                    X(xv), h - 28, text=lab_, fill="#9aa0ad", font=("Segoe UI", 8)
                )
            plot.create_text(
                w / 2, h - 12, text=xl, fill="#c8ccd6", font=("Segoe UI", 9)
            )
            plot.create_text(
                10, 14, text=yl, fill="#c8ccd6", anchor="w", font=("Segoe UI", 9)
            )
            return X, Y

        def draw_profile(rec):
            d, prm = rec["data"], rec["params"]
            nodes, pic = np.array(d["nodes"], float), np.array(d["Pi_cs"], float)
            R, sig, q = float(prm["R"]), float(prm["sigma"]), float(prm["q"])
            r, pa = continuum_a_2p1(R, sig, q)
            X, Y = axes(0, R, min(pic.min(), pa.min()) - 0.01, 1.005, "r [ℓ_k]", "Π(r)")
            pts = []
            for rr, pp in zip(r[::50], pa[::50]):
                pts += [X(rr), Y(pp)]
            plot.create_line(*pts, fill="#7ad0ff", width=2)
            pts = []
            for rr, pp in zip(nodes, pic):
                pts += [X(rr), Y(pp)]
                plot.create_oval(
                    X(rr) - 2,
                    Y(pp) - 2,
                    X(rr) + 2,
                    Y(pp) + 2,
                    fill="#ffd27a",
                    outline="",
                )
            plot.create_line(*pts, fill="#ffd27a")
            w = max(plot.winfo_width(), 400)
            plot.create_text(
                w - 20,
                250,
                anchor="ne",
                fill="#ffd27a",
                text="● pole kauzální množiny (uzly záznamu)",
            )
            plot.create_text(
                w - 20,
                266,
                anchor="ne",
                fill="#7ad0ff",
                text=f"— spojité řešení (a) se stejnou stěnou, σ = {sig:g}, q = {q:g}, R = {R:g}",
            )

        def draw_f19a(rec):
            tab = rec["data"]["table"]
            rows = [t for t in tab if t.get("gamma") in ("h1", "conf")]
            vals = [v for t in rows for v in t["c_disc"] + t["c_cont"]]
            ns = len(rows[0]["c_disc"]) if rows else 1
            sig = (rec.get("params") or {}).get("sigmas", list(range(ns)))
            X, Y = axes(
                -0.15,
                max(ns - 1, 1) + 0.15,
                min(vals + [-0.5]) - 0.1,
                max(vals + [2.0]) + 0.1,
                "σ / ℓ_k (vždy oba řádky r/σ)",
                "c = dB_L[L] / |∇L|²",
                xticks=[(i, f"σ = {sig[i]:g}") for i in range(min(ns, len(sig)))],
            )
            cols = {"h1": "#ffffff", "conf": "#ffb37a"}
            for t in rows:
                pts = []
                for i, v in enumerate(t["c_disc"]):
                    pts += [X(i), Y(v)]
                    plot.create_oval(
                        X(i) - 3,
                        Y(v) - 3,
                        X(i) + 3,
                        Y(v) + 3,
                        fill=cols[t["gamma"]],
                        outline="",
                    )
                if len(pts) >= 4:
                    plot.create_line(*pts, fill=cols[t["gamma"]])
                plot.create_line(
                    X(0),
                    Y(t["c_cont"][0]),
                    X(max(ns - 1, 1)),
                    Y(t["c_cont"][0]),
                    fill=cols[t["gamma"]],
                    dash=(3, 3),
                )
            plot.create_text(
                70,
                40,
                anchor="nw",
                fill="#ffffff",
                text="● H1 (předpověď c = 0, čárkovaně spojitá hodnota směsi)",
            )
            plot.create_text(
                70, 56, anchor="nw", fill="#ffb37a", text="● kontrola γ = −1 (c = 2)"
            )

        tk.Button(top, text="Načíst", command=safe("laboratoř", load)).pack(
            side="left", padx=6
        )
        lb.bind("<<ListboxSelect>>", safe("laboratoř", show))
        safe("laboratoř", load)()

    def selected_laws():
        laws = {"GR": _law_gr, "N": _law_newton, "F3": _law_f3}
        out = [laws[k]() for k in ("GR", "N", "F3") if ghost_vars[k].get()]
        if ghost_vars["PPN"].get():
            out.append(_law_ppn(beta_var.get(), gamma_var.get()))
        return out

    def start_experiment():
        try:
            years = max(1.0, min(200.0, float(years_var.get().replace(",", "."))))
        except ValueError:
            years = 10.0
        laws = [
            _law_newton(),
            _law_h1(),
            _law_gr(),
            _law_ppn(beta_var.get(), gamma_var.get()),
        ]
        state["exp"] = mercury_experiment(laws, years, planets_var.get())
        state["exp_meta"] = (
            years,
            planets_var.get(),
            beta_var.get(),
            gamma_var.get(),
            laws,
        )
        state["exp_result"] = None
        exp_lbl.config(text="běží …")

    # ---------- coordinates
    def W2S(x, y):
        """Převod souřadnic scény na obrazovku, případně přes nakloněnou 3D projekci."""
        cx, cy, hw = state["view"]
        w, h = max(canvas.winfo_width(), 10), max(canvas.winfo_height(), 10)
        sc = 0.5 * w / hw
        if state.get("view3d") and not state["meta"].get("map3d"):
            yaw, pitch = state.get("rot", [0.25, 0.9])
            dx, dy = x - cx, y - cy
            xr = math.cos(yaw) * dx - math.sin(yaw) * dy
            yr = math.sin(yaw) * dx + math.cos(yaw) * dy
            return w / 2 + xr * sc, h / 2 - math.sin(pitch) * yr * sc
        return w / 2 + (x - cx) * sc, h / 2 - (y - cy) * sc

    def S2W(sx, sy):
        """Inverzní převod kliknutí na souřadnice scény včetně 3D projekce roviny z=0."""
        cx, cy, hw = state["view"]
        w, h = max(canvas.winfo_width(), 10), max(canvas.winfo_height(), 10)
        sc = 0.5 * w / hw
        xr = (sx - w / 2) / sc
        if state.get("view3d") and not state["meta"].get("map3d"):
            yaw, pitch = state.get("rot", [0.25, 0.9])
            yr = (h / 2 - sy) / (sc * max(abs(math.sin(pitch)), 0.15))
            dx = math.cos(yaw) * xr + math.sin(yaw) * yr
            dy = -math.sin(yaw) * xr + math.cos(yaw) * yr
            return cx + dx, cy + dy
        return cx + xr, cy - (sy - h / 2) / sc

    # ---------- správa scén
    def orbits_tracked():
        return (not state["meta"]["si"]) or state["meta"].get("orbits", False)

    def setup_trackers():
        """Připraví měřicí sledovače hlavní scény i srovnávacích rámců; stávající průběh zachová."""
        sim = state["sim"]
        old_p = {
            id(t.body): t for t in state["trackers"] if t.body.alive and t.center.alive
        }
        old_l = {id(t.ph): t for t in state["light"] if t.ph.alive and t.m.alive}
        state["trackers"], state["light"] = make_trackers(
            sim, orbits_tracked(), old_p, old_l
        )
        for g in state["ghosts"]:
            op = {
                id(t.body): t for t in g["trackers"] if t.body.alive and t.center.alive
            }
            ol = {id(t.ph): t for t in g["light"] if t.ph.alive and t.m.alive}
            g["trackers"], g["light"] = make_trackers(
                g["sim"], orbits_tracked(), op, ol
            )

    def build_ghosts(name, meta):
        """Vytvoří samostatné kopie scény pro každý vybraný srovnávací zákon."""
        state["ghosts"] = []
        if meta.get("map3d"):  # Hvězdná mapa nemá gravitační dynamiku, kterou by bylo možné srovnávat.
            return
        for law in selected_laws():
            gs = Sim(law)
            gs.bodies = [
                b for b in scene(name)[0] if b.kind != "pilot"
            ]  # Klávesnicí se řídí jen hlavní scéna.
            gs.dt_max, gs.eta = meta["dt_max"], meta.get("eta", 0.02)
            gs.no_guard = bool(meta.get("bh"))
            if meta.get("bh") and law.key == "H1":
                gs.profile = next(
                    b for b in gs.bodies if getattr(b, "is_profile", False)
                )
            for b in gs.bodies:
                b.t_birth, b.start = 0.0, (b.x, b.y, b.vx, b.vy)
            state["ghosts"].append({"sim": gs, "trackers": [], "light": [], "law": law})

    def load_scene(name):
        sim = Sim()
        sim.bodies, state["view"], meta = scene(name)
        sim.dt_max, sim.eta = meta["dt_max"], meta.get("eta", 0.02)
        if meta.get("bh"):  # Přesné vnější řešení rovnice pole bez ochrany silného pole.
            sim.profile = next(b for b in sim.bodies if getattr(b, "is_profile", False))
            sim.no_guard = True
        state["sim"], state["scene"], state["meta"] = sim, name, meta
        state["mass_reference"] = max(
            (b.mass for b in sim.bodies if b.alive and b.mass > 0 and b is not sim.profile),
            default=(M_SUN if meta.get("map3d") else 1.0),
        )
        tool_var.set("none")
        state["rot"] = [0.6, 0.5]  # otočení 3D mapy: azimut a sklon [rad]
        for b in sim.bodies:  # Počáteční stav pro měření ve scénách u Země: výšky a rychlosti.
            b.t_birth = 0.0
            b.start = (b.x, b.y, b.vx, b.vy)
        build_ghosts(name, meta)
        # Všechny scénáře používají stejný rozsah tempa a stejný poměr hmotnosti m/Mref.
        speed_scale.config(from_=-3.0, to=12.5)
        speed_lbl.config(text="Tempo simulace (log₁₀ jednotek simulačního času za sekundu reálného času)")
        mass_lbl.config(text="Hmotnost nového tělesa: poměr m/Mref (Mref je pevná pro načtenou scénu)")
        ref_mass = referencni_hmotnost(sim)
        if meta.get("map3d"):
            mass_var.set("0.05")
            note_lbl.config(text="\nSkutečná 3D mapa hvězd. Poloha se zadává v ly a rychlost v km/s.\nKliknutím s nástrojem „přidat těleso“ vložíš objekt do roviny pohledu;\nformulář dovoluje zadat i hloubku Z a prostorovou rychlost.")
            click_lbl.config(text="Levým tlačítkem otáčíš mapou. Vyber „přidat těleso“ nebo „smazat“ pro úpravu objektů.")
            map3d_frame.pack(fill="x", padx=8, pady=(6, 0))
            projection_check.pack_forget()
            for key in ("photon", "clock"):
                tool_widgets[key].pack_forget()
            if not tool_widgets["mass"].winfo_manager():
                tool_widgets["mass"].pack(fill="x", padx=16, before=tool_widgets["delete"])
            tool_var.set("none")
            state["view3d"] = False
            projection_var.set(False)
        else:
            map3d_frame.pack_forget()
            if not projection_check.winfo_manager():
                projection_check.pack(fill="x", padx=8, pady=(4, 0), before=click_lbl)
            for key in ("photon", "clock"):
                if not tool_widgets[key].winfo_manager():
                    tool_widgets[key].pack(fill="x", padx=16, before=tool_widgets["delete"])
            if meta.get("si"):
                if state["scene"].startswith("Sluneční soustava (skutečná)"):
                    default_ratio = 6.0458e24 * G_C2 / max(ref_mass, 1e-300)
                elif state["scene"].startswith("Povrch Země") or state["scene"].startswith("UAP:"):
                    default_ratio = 0.2 * G_C2 / max(ref_mass, 1e-300)
                else:
                    default_ratio = 1e3 * G_C2 / max(ref_mass, 1e-300)
                note_lbl.config(text="\nScéna používá vstupy SI; výpočty mají interní geometrické jednotky.\nHmotnost je poměr k referenční hmotnosti Mref stanovené při načtení scény.")
                click_lbl.config(text=("Kliknutí: pustit v klidu; tažení: vrh (80 px = 20 m/s)." if meta["click"] == "rest" else "Kliknutí: kruhová dráha; tažení určuje směr a násobek kruhové rychlosti (80 px = 1×)."))
            else:
                default_ratio = 1e-6
                note_lbl.config(text="\nScéna používá geometrické jednotky G = c = 1. Hmotnost je poměr\nk referenční hmotnosti Mref stanovené při načtení scény.")
                click_lbl.config(text="Kliknutí: kruhová dráha; tažení určuje směr a rychlost (80 px = 1× kruhové rychlosti).")
            mass_var.set(f"{default_ratio:.6g}")
            state["view3d"] = bool(projection_var.get())
        if meta.get("pilot"):
            pilot_frame.pack(fill="x", padx=8, pady=(6, 0))
        else:
            pilot_frame.pack_forget()
        speed_var.set(meta["speed"])
        speed_value_label()
        update_category()
        state["trackers"], state["light"], state["n_events"] = [], [], 0
        state["achieved"] = None
        err_lbl.config(text="")
        setup_trackers()
        canvas.focus_set()  # Klávesy obsluhuje plátno, nikoli seznam scén.
        draw(force_field=True)

    def change_projection():
        """Zapne/vypne pouze prostorové zobrazení rovinných dat."""
        state["view3d"] = bool(projection_var.get()) and not state["meta"].get("map3d")
        state["rot"] = [0.25, 0.95] if state["view3d"] else [0.0, 0.0]
        draw(force_field=True)

    scene_var.trace_add("write", lambda *a: load_scene(scene_var.get()))

    def toggle_run():
        state["running"] = not state["running"]
        run_btn.config(text="Pauza" if state["running"] else "Spustit")

    run_btn.config(command=toggle_run)

    # ---------- obraz pole s pevnou barevnou škálou
    def field_image():
        sim = state["sim"]
        src = [b for b in sim.bodies if b.alive and b.mass > 0]
        w, h = canvas.winfo_width(), canvas.winfo_height()
        k = 6
        gw, gh = max(w // k, 2), max(h // k, 2)
        if not src:
            return None
        if state.get("view3d") and not state["meta"].get("map3d"):
            sx, sy = np.meshgrid(
                np.arange(gw) * k + k / 2,
                np.arange(gh) * k + k / 2,
            )
            cx, cy, hw = state["view"]
            scale = 0.5 * max(w, 10) / hw
            yaw, pitch = state.get("rot", [0.25, 0.95])
            xr = (sx - w / 2) / scale
            yr = (h / 2 - sy) / (scale * max(abs(math.sin(pitch)), 0.15))
            X = cx + math.cos(yaw) * xr + math.sin(yaw) * yr
            Y = cy - math.sin(yaw) * xr + math.cos(yaw) * yr
        else:
            xs = np.array([S2W(i * k + k / 2, h / 2)[0] for i in range(gw)])
            ys = np.array([S2W(w / 2, j * k + k / 2)[1] for j in range(gh)])
            X, Y = np.meshgrid(xs, ys)
        pts = np.stack([X.ravel(), Y.ravel()], axis=1)
        L = sim.L_points(pts)  # Bodové hmoty a případný přesný vnější profil.
        disp = state["meta"]["display"]
        if disp[0] == "log":  # Pevná logaritmická škála |L|.
            lo, hi = disp[1], disp[2]
            t = np.clip(
                (np.log10(np.maximum(-L, 1e-300)) - lo) / (hi - lo), 0.0, 1.0
            ).reshape(gh, gw)
        else:  # Režim „rel“ zobrazuje odchylku L od referenční hodnoty.
            t = np.clip(1.0 - (L - disp[1]) / disp[2], 0.0, 1.0).reshape(gh, gw)
            t = np.where(t >= 1.0, 1.0, t)
        # Barevná škála: tmavě modrá → fialová → oranžová → světle žlutá; rozsah je pevný pro danou scénu.
        r = np.clip(40 + 230 * t**0.9, 0, 255)
        g = np.clip(18 + 200 * t**2.2, 0, 255)
        bl = np.clip(
            60 + 120 * np.sin(np.pi * np.minimum(t * 1.2, 1.0)) - 40 * t, 0, 255
        )
        rgb = np.stack([r, g, bl], axis=2).astype(np.uint8)
        dim = ((t > 0) | (disp[0] == "rel")).astype(np.uint8)[
            ..., None
        ]  # 'log': L ~ 0 far away is dark
        rgb = (rgb * dim + np.array([11, 13, 18], np.uint8) * (1 - dim)).astype(
            np.uint8
        )
        ppm = b"P6 %d %d 255\n" % (gw, gh) + rgb.tobytes()
        img = tk.PhotoImage(data=ppm, format="PPM")
        return img.zoom(k)

    # ---------- drawing
    def draw_map3d():
        """Vykreslí skutečnou 3D mapu hvězd ortografickou projekcí po otočení (azimut, sklon)."""
        sim = state["sim"]
        w, h = max(canvas.winfo_width(), 10), max(canvas.winfo_height(), 10)
        cx, cy, hw = state["view"]
        sc = 0.5 * w / hw
        yaw, pitch = state["rot"]
        cyw, syw, cp, sp = (
            math.cos(yaw),
            math.sin(yaw),
            math.cos(pitch),
            math.sin(pitch),
        )

        def proj(x, y, z):
            x1, y1 = (
                cyw * x - syw * y,
                syw * x + cyw * y,
            )  # Otočení kolem nebeského pólu.
            y2, z2 = cp * y1 - sp * z, sp * y1 + cp * z  # pitch
            return w / 2 + x1 * sc, h / 2 - z2 * sc, y2  # Hloubka y2 určuje pořadí vykreslení.

        for rr in (4, 8, 12):  # Kružnice v rovníkové rovině.
            pts = []
            for k in range(73):
                a = 2 * math.pi * k / 72
                X, Y, _ = proj(rr * LY * math.cos(a), rr * LY * math.sin(a), 0.0)
                pts += [X, Y]
            canvas.create_line(*pts, fill="#33415a", dash=(2, 3))
            X, Y, _ = proj(rr * LY, 0, 0)
            canvas.create_text(
                X + 4,
                Y,
                text=f"{rr} ly",
                fill="#5a6a88",
                anchor="w",
                font=("Segoe UI", 8),
            )
        for x, y, z, lab_ in (
            (13 * LY, 0, 0, "RA 0h"),
            (0, 13 * LY, 0, "RA 6h"),
            (0, 0, 13 * LY, "sever"),
        ):
            X0, Y0, _ = proj(0, 0, 0)
            X1, Y1, _ = proj(x, y, z)
            canvas.create_line(X0, Y0, X1, Y1, fill="#46506a")
            canvas.create_text(X1, Y1, text=lab_, fill="#7a86a0", font=("Segoe UI", 8))
        if trails_var.get():  # Trajektorie hvězd při přímočarém prostorovém pohybu.
            for b in sim.bodies:
                tr = getattr(b, "trail3", None)
                if tr and len(tr) > 1:
                    pts = []
                    for x, y, z in list(tr)[::2] + [(b.x, b.y, b.z)]:
                        X, Y, _ = proj(x, y, z)
                        pts += [X, Y]
                    canvas.create_line(*pts, fill=b.color, width=1)
        items = []
        projected_points = {}
        for b in sim.bodies:
            X, Y, dep = proj(b.x, b.y, b.z)
            items.append((dep, X, Y, b))
            projected_points[id(b)] = (X, Y, b)
        state["map3d_points"] = projected_points
        for dep, X, Y, b in sorted(items, key=lambda t: -t[0]):  # nejvzdálenější hvězdy se kreslí první
            rad = max(
                2.0, min(9.0, 3.0 + 3.0 * math.log10(max(b.mass / M_SUN, 0.05) / 0.05))
            )
            X0, Y0, _ = proj(b.x, b.y, 0.0)
            canvas.create_line(
                X, Y, X0, Y0, fill="#2a3346"
            )  # Spojnice k rovníkové rovině.
            canvas.create_oval(
                X - rad, Y - rad, X + rad, Y + rad, fill=b.color, outline=""
            )
            if names_var.get():
                canvas.create_text(
                    X + rad + 3,
                    Y - 3,
                    text=f"{b.name} ({b.info.split(', ')[-1]})",
                    fill="#c8ccd6",
                    anchor="w",
                    font=("Segoe UI", 8),
                )
        canvas.create_text(
            10,
            10,
            anchor="nw",
            fill="#9aa0ad",
            font=("Segoe UI", 9),
            text=(f"{state['scene']}   t = {fmt_time(sim.t / C_SI)} od J2000   "
                  f"(nástroj: {tool_var.get()}, levé tlačítko = {('otáčení' if tool_var.get() == 'none' else 'úprava objektu')}, "
                  "kolečko = zoom; pohyb hvězd je přímočarý)"),
        )
        scale_bar()
        readouts()

    def draw(force_field=False, refresh_field=True):
        canvas.delete("all")
        sim = state["sim"]
        if state["meta"].get("map3d"):
            draw_map3d()
            return
        if field_var.get() and (
            force_field
            or (refresh_field and state["frame"] % 4 == 0)
            or state["field_img"] is None
        ):
            try:
                state["field_img"] = field_image()
            except tk.TclError:
                state["field_img"] = None
        if field_var.get() and state["field_img"] is not None:
            canvas.create_image(0, 0, image=state["field_img"], anchor="nw")
        for b in sim.bodies:  # Pevná tělesa, například Země, se kreslí jako první.
            if b.alive and b.radius > 0:
                draw_solid(b)
        if (
            state["meta"].get("bh") and sim.profile is not None
        ):  # Fotonové sféry a kritické parametry dopadu.
            info = bh_profile()[3]
            c = sim.profile
            for rr, col, dash in (
                (info["r_ph"], "#ffffff", (4, 3)),
                (1.866, "#ff7ad9", (2, 3)),
            ):
                x0, y0 = W2S(c.x - rr, c.y + rr)
                x1, y1 = W2S(c.x + rr, c.y - rr)
                canvas.create_oval(x0, y0, x1, y1, outline=col, dash=dash)
            x0, y0 = W2S(c.x - 0.5, c.y + 0.5)
            x1, y1 = W2S(c.x + 0.5, c.y - 0.5)
            canvas.create_oval(
                x0, y0, x1, y1, outline="#ff7ad9", width=2
            )  # Schwarzschild horizon (isotropic r = M/2)
            for bb, col in (
                (info["b_crit"], "#ffffff"),
                (info["b_crit_gr"], "#ff7ad9"),
            ):
                xa, ya = W2S(-40, bb)
                xb, yb = W2S(-6, bb)
                canvas.create_line(xa, ya, xb, yb, fill=col, dash=(1, 4))
        for g in state["ghosts"]:  # ghosts: thin dashed trails, hollow markers
            col = g["law"].color
            for b in g["sim"].bodies:
                if not b.alive or b.radius > 0 or b.fixed:
                    continue
                if trails_var.get() and len(b.trail) > 1:
                    pts = []
                    for x, y in list(b.trail)[::2]:
                        pts.extend(W2S(x, y))
                    if len(pts) >= 4:
                        canvas.create_line(*pts, fill=col, width=1, dash=(4, 3))
                sx, sy = W2S(b.x, b.y)
                canvas.create_oval(sx - 4, sy - 4, sx + 4, sy + 4, outline=col, width=1)
        for b in sim.bodies:
            if not b.alive or b.radius > 0:
                continue
            if trails_var.get() and len(b.trail) > 1:
                pts = []
                for x, y in list(b.trail)[::2]:
                    pts.extend(W2S(x, y))
                if len(pts) >= 4:
                    canvas.create_line(*pts, fill=b.color, width=1)
            sx, sy = W2S(b.x, b.y)
            if b.kind == "mass":
                rad = max(
                    2.5,
                    min(12.0, 3.0 + 2.2 * math.log10(max(b.mass, 1e-9) / 1e-9) / 2.0),
                )
                canvas.create_oval(
                    sx - rad, sy - rad, sx + rad, sy + rad, fill=b.color, outline=""
                )
            elif b.kind == "photon":
                canvas.create_oval(
                    sx - 2, sy - 2, sx + 2, sy + 2, fill=b.color, outline=""
                )
            elif b.kind in ("craft", "pilot"):  # a saucer-shaped marker
                canvas.create_oval(
                    sx - 12, sy - 4, sx + 12, sy + 4, fill=b.color, outline="#ffffff"
                )
                canvas.create_oval(
                    sx - 5, sy - 8, sx + 5, sy - 1, fill="#d0ffff", outline=""
                )
            else:
                canvas.create_rectangle(
                    sx - 4, sy - 4, sx + 4, sy + 4, outline=b.color, width=2
                )
            if names_var.get():
                canvas.create_text(
                    sx + 8,
                    sy - 8,
                    text=b.name,
                    fill="#c8ccd6",
                    anchor="w",
                    font=("Segoe UI", 8),
                )
        if state["drag"] is not None:
            (x0, y0), (x1, y1) = state["drag"]
            canvas.create_line(
                x0, y0, x1, y1, fill="#ffffff", arrow="last", dash=(3, 2)
            )
            if tool_var.get() == "mass" and math.hypot(x1 - x0, y1 - y0) >= 6:
                k = math.hypot(x1 - x0, y1 - y0) / 80.0
                txt = (
                    f"v = {20 * k:.1f} m/s"
                    if state["meta"]["si"] and state["meta"]["click"] == "rest"
                    else f"v = {k:.2f} × kruhová rychlost"
                )
                canvas.create_text(
                    x1 + 10,
                    y1,
                    anchor="w",
                    fill="#ffffff",
                    font=("Segoe UI", 9),
                    text=txt,
                )
        if state["meta"]["si"]:
            ttxt = "t = " + fmt_time(sim.t / C_SI)
            canvas.create_text(
                10,
                10,
                anchor="nw",
                fill="#9aa0ad",
                font=("Segoe UI", 9),
                text=f"{state['scene']}   {ttxt}   (SI; velké vzdálenosti v AU)",
            )
        else:
            canvas.create_text(
                10,
                10,
                anchor="nw",
                fill="#9aa0ad",
                font=("Segoe UI", 9),
                text=f"{state['scene']}   t = {sim.t:.2f}   (G = c = 1)",
            )
        scale_bar()
        legend()
        if state["ghosts"]:
            w = max(canvas.winfo_width(), 10)
            y = 10
            canvas.create_text(
                w - 10,
                y,
                anchor="ne",
                fill="#ffffff",
                font=("Segoe UI", 9),
                text="— kauzální H1 (plné stopy); fyzikální rámce:",
            )
            for g in state["ghosts"]:
                y += 15
                canvas.create_text(
                    w - 10,
                    y,
                    anchor="ne",
                    fill=g["law"].color,
                    font=("Segoe UI", 9),
                    text=f"- - {g['law'].name}",
                )
        if state["hover"] is not None:
            hover_text()
        readouts()

    def hover_text():
        """Zobrazí hodnoty v místě kurzoru: L, Π a v SI také rychlost hodin za den."""
        sim = state["sim"]
        hx, hy = state["hover"]
        wx, wy = S2W(hx, hy)
        L = sim.L_at(wx, wy)
        h = max(canvas.winfo_height(), 10)
        if state["meta"]["si"]:
            src = [b for b in sim.bodies if b.alive and b.mass > 0]
            c = max(src, key=lambda b: b.mass) if src else None
            dist = f"{fmt_len(math.hypot(wx - c.x, wy - c.y))} od {c.name}" if c else ""
            txt = (
                f"{dist} | L = ln Π = {L:.4e}, Π = 1{math.expm1(L):+.3e} | hodiny v klidu: "
                f"{math.expm1(L) * 86400e6:+.3f} μs/den proti nekonečnu"
            )
        else:
            txt = f"x = {wx:.3f}, y = {wy:.3f} | L = ln Π = {L:.4e}, Π = {math.exp(L):.6f} (rychlost hodin v klidu)"
        canvas.create_text(
            10, h - 28, anchor="w", text=txt, fill="#ffffff", font=("Segoe UI", 9)
        )

    def draw_solid(b):
        """Pevné těleso vykreslí jako kruh, pokud se vejde do pohledu, jinak jako viditelnou část povrchu."""
        cx, cy, hw = state["view"]
        w = max(canvas.winfo_width(), 10)
        sc = 0.5 * w / hw
        if b.radius * sc < 5000:
            x0, y0 = W2S(b.x - b.radius, b.y + b.radius)
            x1, y1 = W2S(b.x + b.radius, b.y - b.radius)
            canvas.create_oval(x0, y0, x1, y1, fill=b.color, outline="#7fb07a")
        else:
            phi0 = math.atan2(
                cy - b.y, cx - b.x
            )  # Směr od středu scény k pozorovateli.
            dphi = 4.0 * hw / b.radius
            pts = []
            for k in range(81):
                ph = phi0 - dphi + 2 * dphi * k / 80
                pts.extend(
                    W2S(b.x + b.radius * math.cos(ph), b.y + b.radius * math.sin(ph))
                )
            deep = b.radius - 10.0 * hw
            for ph in (phi0 + dphi, phi0 - dphi):
                pts.extend(W2S(b.x + deep * math.cos(ph), b.y + deep * math.sin(ph)))
            canvas.create_polygon(*pts, fill=b.color, outline="#7fb07a")
        if names_var.get():
            sx, sy = W2S(b.x, b.y)
            if b.radius * sc < 5000:
                canvas.create_text(
                    sx, sy, text=b.name, fill="#e0f0e0", font=("Segoe UI", 9)
                )

    def scale_bar():
        """Vykreslí měřítko o zaokrouhlené délce 1, 2 nebo 5 × 10^k, přibližně pětinu šířky pohledu."""
        cx, cy, hw = state["view"]
        w, h = max(canvas.winfo_width(), 10), max(canvas.winfo_height(), 10)
        sc = 0.5 * w / hw
        target = 0.4 * hw
        unit = 1.0
        if state["meta"]["si"]:  # Při velkém zobrazení se používají světelné roky nebo astronomické jednotky.
            unit = LY if target >= 0.1 * LY else (AU if target >= 0.01 * AU else 1.0)
        p10 = 10 ** math.floor(math.log10(target / unit))
        length = max(m * p10 for m in (1, 2, 5) if m * p10 <= target / unit) * unit
        px = length * sc
        x0, y0 = w - 30 - px, h - 22
        canvas.create_line(x0, y0, x0 + px, y0, fill="#ffffff", width=2)
        for xx in (x0, x0 + px):
            canvas.create_line(xx, y0 - 5, xx, y0 + 5, fill="#ffffff", width=2)
        if state["meta"]["si"]:
            txt = (
                f"{length / LY:g} ly"
                if unit == LY
                else (
                    f"{length / AU:g} AU"
                    if unit == AU
                    else (f"{length/1000:g} km" if length >= 1000 else f"{length:g} m")
                )
            )
        else:
            txt = f"{length:g} jednotek délky"
        canvas.create_text(
            x0 + px / 2, y0 - 10, text=txt, fill="#ffffff", font=("Segoe UI", 9)
        )

    def legend():
        """Vykreslí legendu pevné barevné škály gravitačního pole."""
        if not field_var.get():
            return
        h = max(canvas.winfo_height(), 10)
        disp = state["meta"]["display"]
        if disp[0] == "log":
            txt = f"barva: log10|L| od {disp[1]:g} (tmavá) do {disp[2]:g} (světlá), L = ln Π"
        else:
            txt = (
                f"barva: zpomalení hodin vůči výšce; světlá = povrch, tmavá = o {disp[2]:.2e} rychlejší "
                f"(ΔL = gh/c² na rozsahu {disp[2]*C_SI**2/G_SURF:.0f} m)"
            )
        canvas.create_text(
            10, h - 12, anchor="w", text=txt, fill="#c8ccd6", font=("Segoe UI", 8)
        )

    def readouts_si():
        """Vypíše výsledky scén u Země: pády, vrhy, rozdíly hodin a naprogramovaný manévr UAP."""
        sim = state["sim"]
        earth = next((b for b in sim.bodies if b.alive and b.radius > 0), None)
        lines = [
            "Země: M = GM/c² = 4,435 mm, R = 6371 km, g = %.3f m/s². Pole slabé (L ≈ −7·10⁻¹⁰): 1PN = GR, Newton "
            "na povrchu stačí; demonstrace ukazuje, že teorie tu dává totéž co známá fyzika."
            % G_SURF
        ]
        if earth is None:
            return lines
        for b in sim.bodies:
            if (
                not b.alive
                or b.kind not in ("mass",)
                or b.radius > 0
                or not hasattr(b, "start")
            ):
                continue
            x0, y0, vx0, vy0 = b.start
            r0 = math.hypot(x0 - earth.x, y0 - earth.y)
            if r0 <= 0.0:
                continue
            h0 = max(r0 - earth.radius, 0.0)
            g = earth.mass / r0**2 * C_SI**2
            ux, uy = (x0 - earth.x) / r0, (y0 - earth.y) / r0  # local vertical
            vv, vh = (vx0 * ux + vy0 * uy) * C_SI, (vx0 * uy - vy0 * ux) * C_SI
            if h0 > 2e5:  # Obíhající těleso: údaj o volném pádu se nezobrazuje.
                continue
            t_pred = (
                vv + math.sqrt(vv * vv + 2 * g * h0)
            ) / g  # Výpočet vrhu v homogenním gravitačním poli z výšky h0.
            if b.landed_t is not None:
                t_meas = (b.landed_t - (b.t_birth or 0.0)) / C_SI
                dist = math.hypot(b.x - x0, b.y - y0)
                lines.append(
                    f"DOPAD {b.name}: z výšky {h0:.1f} m za {t_meas:.3f} s (Newton {t_pred:.3f} s), "
                    f"vzdálenost {dist:.2f} m (Newton {abs(vh) * t_pred:.2f} m)"
                )
            else:
                lines.append(
                    f"LET {b.name}: výška {math.hypot(b.x - earth.x, b.y - earth.y) - earth.radius:.2f} m, "
                    f"rychlost {math.hypot(b.vx, b.vy) * C_SI:.2f} m/s"
                )
        ground = next(
            (
                b
                for b in sim.bodies
                if b.alive and b.kind == "clock" and "na zemi" in b.name
            ),
            None,
        )
        if ground is None:  # Pokud chybí pozemní hodiny, porovnávají se všechny s hodinami v nekonečnu.
            lines += clock_lines()
        if ground is not None:
            Lg = sim.L_at(ground.x, ground.y)
            age_g = sim.t - (ground.t_birth or 0.0)
            def_g = (
                ground.tau_def / age_g if age_g > 0 else float("nan")
            )  # 1 minus rychlost pozemních hodin.
            for b in sim.bodies:
                if (
                    not b.alive
                    or b is ground
                    or b.kind not in ("clock", "mass")
                    or b.radius > 0
                ):
                    continue
                if (
                    b.kind == "mass"
                    and math.hypot(b.x - earth.x, b.y - earth.y) - earth.radius < 2e5
                ):
                    continue
                age = sim.t - (b.t_birth or 0.0)
                if age <= 0:
                    continue
                d_meas = def_g - b.tau_def / age  # rate(b) - rate(ground)
                L = sim.L_at(b.x, b.y, skip=b)
                d_pred = (L - Lg) - 0.5 * (
                    b.vx * b.vx + b.vy * b.vy
                )  # 1PN: Delta L - v^2/2 (c = 1)
                f, unit = (
                    (86400e6, "μs/den")
                    if abs(d_pred) * 86400e6 >= 0.1
                    else (86400e9, "ns/den")
                )
                lines.append(
                    f"HODINY {b.name} proti hodinám na zemi: {d_meas * f:+.4f} {unit} "
                    f"(předpověď ΔL − v²/2: {d_pred * f:+.4f} {unit}; relativně {d_meas:+.3e})"
                )
        for b in sim.bodies:
            if b.alive and b.kind == "craft" and b.script is not None:
                names = getattr(b, "phase_names", [])
                ph = b.script.phase_name(sim.t / C_SI, names)
                lines.append(
                    f"UAP: {ph}; rychlost {math.hypot(b.vx, b.vy) * C_SI:.0f} m/s; vlastní zrychlení teď "
                    f"{b.now_g / 9.81:.1f} g, nejvíc {b.max_g / 9.81:.0f} g (tolik by musel dodat tah; "
                    f"vznášení = 1 g, volný pád = 0 g)"
                )
                lines.append(
                    "   Hypotéza autora „null field“ (NENÍ součást teorie ani laboratoře): loď by měnila pole Π "
                    "kolem sebe a padala by v něm volně, posádka by cítila 0 g. Potřebný gradient "
                    f"|∇L| = a/c² ≈ {b.max_g / C_SI ** 2:.1e} 1/m. Rovnice pole zná jen kladný zdroj κρ_s, který přitahuje; "
                    "zdroj řízený lodí ani zachování hybnosti bez reakční hmoty v ní nejsou vyřešené."
                )
        return lines

    def clock_lines():
        """Pro každé hodiny vypíše měřenou rychlost vůči hodinám v nekonečnu a statickou předpověď e^L.
        Ve scénách SI se rozdíly uvádějí v mikrosekundách nebo nanosekundách za den."""
        sim = state["sim"]
        si = state["meta"]["si"]
        lines = []
        for b in sim.bodies:
            if not (b.alive and b.kind == "clock"):
                continue
            age = sim.t - (b.t_birth or 0.0)
            if age <= 0:
                lines.append(f"HODINY {b.name}: čekám, až uběhne čas …")
                continue
            d = b.tau_def / age  # Naměřené 1 − dτ/dt.
            pred = -math.expm1(max(sim.L_at(b.x, b.y), -50.0))  # Statická předpověď 1 − e^L.
            if si:
                f, unit = (
                    (86400e6, "μs/den")
                    if abs(pred) * 86400e6 >= 0.1
                    else (86400e9, "ns/den")
                )
                lines.append(
                    f"HODINY {b.name}: proti hodinám v nekonečnu {-d * f:+.4f} {unit} (předpověď e^L teď: "
                    f"{-pred * f:+.4f} {unit}; relativně {-d:+.3e})"
                )
            else:
                lines.append(
                    f"HODINY {b.name}: dτ/dt = {1 - d:.6f} (předpověď e^L teď = {1 - pred:.6f})"
                )
        return lines

    def experiment_lines():
        """Vypíše výsledky experimentu s Merkurem, pokud jsou dostupné, a porovná je s pozorováním a vztahem PPN."""
        r = state.get("exp_result")
        if not r:
            return []
        years, planets, beta, gamma, laws = state["exp_meta"]
        names = {law.key: law.name for law in laws}
        lines = [
            f"EXPERIMENT MERKUR ({years:g} let, {'s planetami' if planets else 'jen Slunce + Merkur'}), precese ″/století:"
        ]
        base = r.get("N")
        for k in ("N", "H1", "GR", "PPN"):
            if r.get(k) is None:
                lines.append(f"   {names.get(k, k)}: málo perihelií")
                continue
            extra = (
                ""
                if k == "N" or base is None
                else f"   (rozdíl proti Newtonovi {r[k] - base:+.2f}″)"
            )
            lines.append(f"   {names[k]}: {r[k]:+.2f}″{extra}")
        lines.append(
            f"   pozorovaný relativistický příspěvek 42,98 ± 0,04″/století (GR 42,98″); PPN (2+2γ−β)/3 · 42,98 = "
            f"{(2 + 2 * gamma - beta) / 3 * 42.98:.2f}″ pro β = {beta:g}, γ = {gamma:g}"
        )
        if planets:
            lines.append(
                "   newtonovská část od planet je v této 2D soustavě jiná než skutečných ~532″ (sklony drah chybí); "
                "srovnává se proto jen rozdíl proti Newtonovi"
            )
        return lines

    def comparison_lines():
        """Porovná stejnou měřenou veličinu v hlavní scéně H1 a ve všech srovnávacích bězích."""
        if not state["ghosts"]:
            return []
        si = state["meta"]["si"]
        lines = []
        for tr in [t for t in state["trackers"] if t.body.alive]:
            vals = [("H1", tr)] + [
                (
                    g["law"].key,
                    next(
                        (t for t in g["trackers"] if t.body.name == tr.body.name), None
                    ),
                )
                for g in state["ghosts"]
            ]
            parts = []
            for key, t in vals:
                if t is None:
                    continue
                if si:
                    rr = t.rate()
                    parts.append(
                        f"{key} {rr * 100 * YEAR * C_SI / ARCSEC:+.2f}″/stol."
                        if rr is not None
                        else f"{key} …"
                    )
                else:
                    res = t.result()
                    parts.append(
                        f"{key} {math.degrees(res[0]):+.3f}°/oběh"
                        if res
                        else f"{key} …"
                    )
            lines.append(f"RÁMCE precese {tr.body.name}: " + " | ".join(parts))
        for lt in [t for t in state["light"] if t.ph.alive]:
            parts = [f"H1 {math.degrees(lt.state(state['sim'])['defl']):.4f}°"]
            for g in state["ghosts"]:
                t = next((x for x in g["light"] if x.ph.name == lt.ph.name), None)
                if t is not None:
                    parts.append(
                        f"{g['law'].key} {math.degrees(t.state(g['sim'])['defl']):.4f}°"
                    )
            lines.append(f"RÁMCE ohyb {lt.ph.name}: " + " | ".join(parts))
        sim = state["sim"]

        def twin(g, b):
            return next((x for x in g["sim"].bodies if x.name == b.name), None)

        for b in [
            x for x in sim.bodies if x.alive and x.kind == "clock"
        ]:  # Rychlosti hodin vůči hodinám v nekonečnu.
            parts = []
            for key, sm, bb in [("H1", sim, b)] + [
                (g["law"].key, g["sim"], twin(g, b)) for g in state["ghosts"]
            ]:
                if bb is None:
                    continue
                age = sm.t - (bb.t_birth or 0.0)
                if age <= 0:
                    continue
                d = bb.tau_def / age
                parts.append(
                    f"{key} {-d * 86400e6:+.4f} μs/den" if si else f"{key} {1 - d:.6f}"
                )
            if parts:
                lines.append(
                    f"RÁMCE hodiny {b.name} ({'proti nekonečnu' if si else 'dτ/dt'}): "
                    + " | ".join(parts)
                    + (
                        "  (Newton: hodiny pole necítí)"
                        if any(g["law"].key == "N" for g in state["ghosts"])
                        else ""
                    )
                )
        for b in [
            x
            for x in sim.bodies
            if x.alive and x.kind == "mass" and x.landed_t is not None
        ]:  # fall times
            parts = [
                f"H1 {(b.landed_t - (b.t_birth or 0.0)) / (C_SI if si else 1.0):.6f}"
            ]
            for g in state["ghosts"]:
                bb = twin(g, b)
                if bb is not None and bb.landed_t is not None:
                    parts.append(
                        f"{g['law'].key} {(bb.landed_t - (bb.t_birth or 0.0)) / (C_SI if si else 1.0):.6f}"
                    )
            lines.append(
                f"RÁMCE dopad {b.name} [{'s' if si else 'čas'}]: " + " | ".join(parts)
            )
        return lines

    def readouts_solar():
        """Vypíše pro skutečnou Sluneční soustavu posun perihélia planet v úhlových vteřinách za století."""
        lines = [
            "Sluneční soustava J2000 (2D, sklony drah zanedbány), skutečné hmotnosti; pole slabé (L ≈ −10⁻⁸): teorie "
            "= GR, od Newtona se liší o 1PN efekty. Precese z vektoru Laplace–Runge–Lenz v periheliích."
        ]
        for tr in [t for t in state["trackers"] if t.body.alive]:
            rr = tr.rate()
            M = tr.center.mass + tr.body.mass
            a = 0.5 * (tr.r_min + tr.r_max) if tr.r_max > 0 else 0.0
            e = (tr.r_max - tr.r_min) / (tr.r_max + tr.r_min) if tr.r_max > 0 else 0.0
            if rr is None or a <= 0:
                lines.append(
                    f"PRECESE {tr.body.name}: čekám na tři průchody perihéliem …"
                )
                continue
            period = 2 * math.pi * math.sqrt(a**3 / M)
            gr = (
                6
                * math.pi
                * M
                / (a * (1 - e * e))
                * (100 * YEAR * C_SI / period)
                / ARCSEC
            )
            lines.append(
                f"PRECESE {tr.body.name}: {rr * 100 * YEAR * C_SI / ARCSEC:+.2f}″/století celkem "
                f"({len(tr.lrl)} perihelií; samotná relativita by dala {gr:.2f}″, zbytek jsou poruchy planet)"
            )
        return lines + clock_lines()

    def readouts_bh():
        """Porovná vnější řešení modelu černé díry se Schwarzschildovým řešením."""
        sim = state["sim"]
        info = bh_profile()[3]
        if sim.profile is None:
            return [
                "ČERNÁ DÍRA SMAZÁNA: její pole zmizelo, fotony a částice se dál pohybují jen v poli zbylých těles "
                "(tlačítko Znovu scénu obnoví)."
            ]
        lines = [
            "ČERNÁ DÍRA? Spojitá limita rovnice laboratoře B_g ln Π = 0 vně zdroje, geometrie H1 (předpoklad i v silném "
            "poli); diskrétně NEOVĚŘENO; co je uvnitř, teorie neříká. Délky v GM/c².",
            f"   horizont: ŽÁDNÝ (Π → 0 až v r → 0) | GR: horizont R = 2M (růžový kruh)",
            f"   fotonová sféra: R = {info['R_ph']:.3f} M (bílý kruh) | GR 3,000 M (růžový)",
            f"   kritický parametr (velikost stínu): b = {info['b_crit']:.4f} M | GR {info['b_crit_gr']:.4f} M "
            f"→ stín o {100 * (info['b_crit'] / info['b_crit_gr'] - 1):+.2f} % větší (tečkované čáry)",
        ]
        gr = next((g for g in state["ghosts"] if g["law"].key == "GR"), None)
        for b in sim.bodies:
            if b.kind != "photon":
                continue
            st = (
                "zachycen (zamrzá u středu)"
                if getattr(b, "captured", False)
                else "letí / unikl"
            )
            lt = next((t for t in state["light"] if t.ph is b), None)
            dfl = (
                f", ohyb {math.degrees(lt.state(sim)['defl']):.1f}°"
                if lt and not getattr(b, "captured", False)
                else ""
            )
            line = f"FOTON {b.name}: H1 {st}{dfl}"
            if gr is not None:
                gb = next((x for x in gr["sim"].bodies if x.name == b.name), None)
                if gb is not None:
                    line += f" | GR {'zachycen (horizont)' if getattr(gb, 'captured', False) else 'letí / unikl'}"
            lines.append(line)
        for b in sim.bodies:
            if b.kind == "clock":
                L = max(sim.L_at(b.x, b.y), -50.0)
                R_ar = math.hypot(b.x - sim.profile.x, b.y - sim.profile.y) * math.exp(
                    -L
                )
                if R_ar <= 0.0:
                    continue
                z_h1 = math.exp(-L) - 1
                z_gr = (1 - 2 / R_ar) ** -0.5 - 1 if R_ar > 2 else float("inf")
                lines.append(
                    f"HODINY {b.name}: rychlost Π = {math.exp(L):.4f} (GR √(1−2M/R) = "
                    f"{math.sqrt(max(1 - 2 / R_ar, 0)):.4f}); rudý posuv světla odtud z = {z_h1:.3f} (GR {z_gr:.3f})"
                )
        if gr is None:
            lines.append(
                "Tip: v záložce Srovnání zapni GR: fotony b = 5,21 a 5,23 tady zachycené, v GR uniknou."
            )
        return lines

    def readouts_pilot():
        sim = state["sim"]
        earth = next((b for b in sim.bodies if b.radius > 0), None)
        lines = [
            "UAP ŘÍZENÍ: šipky/WASD tah, Shift ×10, mezerník vznášení. Gravitace Země skutečná (SI); tah = vlastní "
            "zrychlení, které by loď musela dodat."
        ]
        for b in sim.bodies:
            if b.kind != "pilot":
                continue
            h_ = (
                math.hypot(b.x - earth.x, b.y - earth.y) - earth.radius
                if earth
                else 0.0
            )
            v = math.hypot(b.vx, b.vy) * C_SI
            g_now, g_max = b.now_g / 9.81, b.max_g / 9.81
            limit = (
                "bezpečné"
                if g_max < 9
                else (
                    "jen pro piloty v kombinéze (≤ 9 g)"
                    if g_max < 15
                    else (
                        "krátce přežitelné (≲ 45 g)"
                        if g_max < 45
                        else (
                            "SMRTELNÉ pro posádku (> 45 g)"
                            if g_max < 100
                            else "SMRTELNÉ, i stroj by měl problém (> 100 g)"
                        )
                    )
                )
            )
            lines.append(
                f"UAP: výška {fmt_len(h_)}, rychlost {v:.1f} m/s ({v * 3.6:.0f} km/h), vznášení "
                f"{'ZAP' if getattr(b, 'hover', False) else 'vyp'}"
            )
            if nullf_var.get():
                lines.append(
                    f"   režim null field (HYPOTÉZA): klávesy zadávají rychlost {getattr(b, 'v_cmd', 0.0):.0f} m/s, "
                    f"dosažená za {TAU_NULL / C_SI:.2f} s bez ohledu na gravitaci (okamžité zastavení a zatáčky); "
                    "posádka cítí 0 g."
                )
                lines.append(
                    f"   běžná loď by na to potřebovala tah {g_now:.0f} g (teď), nejvíc {g_max:.0f} g ({limit}); "
                    f"vlastní pole lodi by muselo mít gradient |∇L| = a/c² ≈ {b.now_g / C_SI ** 2:.2e} 1/m, "
                    f"max {b.max_g / C_SI ** 2:.2e} 1/m. V teorii ani v laboratoři pro to není zdroj ani mechanismus."
                )
            else:
                lines.append(
                    f"   přetížení posádky teď {g_now:.1f} g, nejvíc {g_max:.1f} g → {limit}"
                )
            if getattr(b, "landing_speed", None):
                lines.append(
                    f"   poslední dosednutí rychlostí {b.landing_speed:.1f} m/s"
                    + (" (havárie!)" if b.landing_speed > 10 else "")
                )
        return lines

    def readouts_binary():
        sim = state["sim"]
        lines = [
            f"{state['scene']}: skutečné hmotnosti a dráha (2D). Pole slabé (L ~ 10⁻⁸): teorie = GR; stáčení periastra "
            "je malé, ukazuje se proti předpovědi 6πM/(a(1−e²))."
        ]
        for tr in state["trackers"]:
            rr = tr.rate()
            M = tr.center.mass + tr.body.mass
            a = 0.5 * (tr.r_min + tr.r_max) if tr.r_max > 0 else 0
            e = (tr.r_max - tr.r_min) / (tr.r_max + tr.r_min) if tr.r_max > 0 else 0
            if a > 0:
                period = 2 * math.pi * math.sqrt(a**3 / M)
                gr = (
                    6
                    * math.pi
                    * M
                    / (a * (1 - e * e))
                    * (100 * YEAR * C_SI / period)
                    / ARCSEC
                )
                lines.append(
                    f"PERIASTRUM {tr.body.name} kolem {tr.center.name}: oběžná doba {fmt_time(period / C_SI)}, "
                    f"a = {fmt_len(a)}, e = {e:.3f}; stáčení "
                    + (
                        f"{rr * 100 * YEAR * C_SI / ARCSEC:.4f}″/století"
                        if rr is not None
                        else "… (čekám na 3 periastra)"
                    )
                    + f" | 1PN předpověď {gr:.4f}″/století"
                )
        if state["scene"].startswith("Sirius"):
            M_B, R_B = 1.018 * M_SUN, 0.0081 * R_SUN
            z = math.expm1(
                M_B / R_B
            )  # H1: z = e^(−L) − 1, kde L = −M/R; korekce (M/R)³/6 je zanedbatelná.
            lines.append(
                f"SIRIUS B (bílý trpaslík): gravitační rudý posuv z povrchu c·z = {z * C_SI / 1000:.1f} km/s pro "
                "M = 1,018 M☉, R = 0,0081 R☉ (teorie = GR v tomto slabém poli); pozorováno 80,0 ± 1,6 km/s "
                "(Joyce et al. 2018, HST). Výsledek závisí hlavně na poloměru hvězdy."
            )
        return lines + clock_lines()

    def readouts_map():
        sim = state["sim"]
        sun = sim.bodies[0]
        yrs = sim.t / C_SI / YEAR
        lines = [
            f"OKOLÍ SLUNCE, {yrs:,.0f} let od J2000 (počáteční polohy J2000). Hvězdy se pohybují přímočaře svou skutečnou "
            "prostorovou rychlostí vůči Slunci (vlastní pohyb RECONS, radiální rychlost Kervella et al. 2019, kde je "
            "známá); vzájemná gravitace ani pole Galaxie se nezapočítávají, na tisíce let je to dobrá aproximace. "
            "Pole hvězd je tu nepatrné (L ~ 10⁻¹³), o pohybu rozhodují jejich rychlosti, ne gravitace."
        ]
        lines.append(
            'Úprava mapy: objekt lze přidat formulářem nebo kliknutím po výběru nástroje „přidat těleso“; '
            'nástroj „smazat“ odstraní nejbližší označenou hvězdu. Nové objekty se pohybují podle zadané prostorové rychlosti.'
        )
        rows = []
        for b in sim.bodies[1:]:
            r = (b.x - sun.x, b.y - sun.y, b.z - sun.z)
            d = math.sqrt(sum(c * c for c in r))
            v = getattr(b, "v3", (0.0, 0.0, 0.0))
            v2 = sum(c * c for c in v)
            extra = ""
            if (
                v2 > 0 and "neznámá" not in b.motion_note
            ):  # Nejmenší vzdálenost na přímočaré dráze, pokud je známa 3D rychlost.
                t_min = -sum(r[k] * v[k] for k in range(3)) / v2
                d_min = math.sqrt(sum((r[k] + v[k] * t_min) ** 2 for k in range(3)))
                when = (sim.t + t_min) / C_SI / YEAR
                extra = f"; nejblíž {d_min / LY:.2f} ly v roce {2000 + when:,.0f}"
            rows.append((d, b, extra))
        for d, b, extra in sorted(rows, key=lambda t: t[0]):
            lines.append(
                f"   {b.name:20s} {d / LY:6.2f} ly   {b.info.split(', ')[0]:22s} {b.motion_note}{extra}"
            )
        return lines

    def readouts():
        sim = state["sim"]
        meta = state["meta"]
        if meta.get("bh"):
            set_readout("\n".join(readouts_bh() + comparison_lines()))
            return
        if meta.get("map3d"):
            set_readout("\n".join(readouts_map()))
            return
        if meta.get("pilot"):
            note = (
                [
                    "Fyzikální rámce: řízený UAP se v nich nepočítá (klávesnice řídí jen hlavní scénu), srovnává se zbytek scény."
                ]
                if state["ghosts"]
                else []
            )
            set_readout(
                "\n".join(
                    readouts_pilot() + readouts_si()[1:] + note + comparison_lines()
                )
            )
            return
        if meta.get("orbits") == "binary":
            set_readout(
                "\n".join(readouts_binary() + comparison_lines() + experiment_lines())
            )
            return
        if state["meta"]["si"]:
            body = readouts_solar() if state["meta"].get("orbits") else readouts_si()
            set_readout(
                "\n".join(
                    body + comparison_lines() + experiment_lines() + sim.events[-3:]
                )
            )
            return
        lines = []
        for tr in [t for t in state["trackers"] if t.body.alive and t.center.alive]:
            r = tr.result()
            if r:
                adv, pred, a, e, n = r
                line = (
                    f"PRECESE {tr.body.name}: změřeno {math.degrees(adv):+.3f}°/oběh ({n} oběhů, a = {a:.3f}, "
                    f"e = {e:.3f}) | 1PN β=γ=1: {math.degrees(pred):.3f}° | β=1/2: {math.degrees(pred*7/6):.3f}° | "
                    f"β=3/2: {math.degrees(pred*5/6):.3f}°"
                )
                # Efekt 1PN tvoří jen přibližně 3M/p Newtonovského působení, proto jej může ovlivnit i slabé působení dalších těles.
                # Podobně jako planetární poruchy u Merkuru; upozornění se zobrazí při překročení 10 % hodnoty 3M/p.
                rel_1pn = 3.0 * (tr.center.mass + tr.body.mass) / (a * (1.0 - e * e))
                if tr.other_max > 1.0:  # Jiné těleso působí silněji než hlavní zdroj.
                    line += (
                        f"  !! RUŠENO: jiné těleso táhne silněji než {tr.center.name} (až {tr.other_max:.0f}×); "
                        f"nejde o oběh kolem jedné hmoty, předpověď neplatí"
                    )
                elif tr.other_max > 0.1 * rel_1pn:
                    line += (
                        f"  !! RUŠENO: ostatní tělesa táhnou až {100*tr.other_max:.2f} % síly {tr.center.name} "
                        f"(1PN efekt je jen {100*rel_1pn:.2f} %); newtonovské poruchy mohou převážit"
                    )
                lines.append(line)
            else:
                lines.append(
                    f"PRECESE {tr.body.name}: čekám na dva oběhy (u téměř kruhové dráhy, e < 0,05, se precese nehlásí) …"
                )
        for lt in [t for t in state["light"] if t.ph.alive and t.m.alive]:
            st = lt.state(sim)
            if st.get("direct"):
                lines.append(
                    f"SVĚTLO {lt.ph.name}: míří přímo na {lt.m.name} (b = {st['b']:.3g}, blíž než {lt.m.mass / L_MAX:.3g}, "
                    f"kde |L| > {L_MAX:g}: bude pohlcen); předpověď ohybu ani Shapirova zpoždění se nedává"
                )
                continue
            lines.append(
                f"SVĚTLO b={st['b']:.2f}: ohyb {math.degrees(st['defl']):.3f}° | H1 γ=1: "
                f"{math.degrees(st['pred'][1.0] + st['pred2']):.3f}° (1. řád {math.degrees(st['pred'][1.0]):.3f}°) | "
                f"γ=0: {math.degrees(st['pred'][0.0]):.3f}° | F3 γ=−1: 0° | Shapiro {st['delay']:+.4f} vs {st['delay_pred'][1.0]:+.4f}"
                + (
                    f"  !! RUŠENO: ostatní tělesa {100*st['other']:.0f} % tahu {lt.m.name}; předpověď počítá jen s ním"
                    if st["other"] > 0.02
                    else ""
                )
            )
        for b in sim.bodies:
            if b.alive and (
                b.kind == "clock"
                or (b.kind == "mass" and b.mass < 1e-3 and not b.fixed)
            ):
                L = sim.L_at(b.x, b.y, skip=b)
                age = sim.t - (
                    b.t_birth or 0.0
                )  # Vlastní čas od vložení tělesa do scény.
                rate = b.tau / age if age > 0 else float("nan")
                lines.append(
                    f"HODINY {b.name}: τ = {b.tau:.3f}, τ/t od vložení = {rate:.6f} | statické hodiny na tomto místě "
                    f"teď: e^L = {math.exp(L):.6f}"
                )
        if sim.events:
            lines.extend(sim.events[-4:])
        head = (
            "Předpovědi jsou 1PN pro zkušební těleso u JEDNÉ hmoty (u světla s 2. řádem metriky); rozdíly řádu M/r jsou "
            "vyšší řády, ne chyba. Ruší-li ostatní tělesa, je to vyznačeno.\n"
        )
        lines += comparison_lines() + experiment_lines()
        set_readout(
            head
            + (
                "\n".join(lines)
                if lines
                else "Žádná měření v této scéně (přidej těleso nebo foton)."
            )
        )

    # ---------- časový krok simulace
    def tick_experiment_only():
        t2 = time.perf_counter()
        try:
            while time.perf_counter() - t2 < 0.025:
                lawname, frac = next(state["exp"])
                exp_lbl.config(text=f"běží: {lawname}, {100 * frac:.0f} %")
        except StopIteration as stop:
            state["exp_result"], state["exp"] = stop.value, None
            exp_lbl.config(text="hotovo, výsledky ve výpisu dole")
            readouts()

    def pilot_control(sim):
        """
        The pilot from the keys. Normal mode: thrust = direction x 10^slider g (Shift x10) plus hovering (1 g up); the
        craft has inertia. 'Null field' mode (the author's HYPOTHESIS): the keys command a VELOCITY (10^(slider+2) m/s,
        Shift x10; no key = stop), reached within TAU_NULL = 0.05 s whatever the gravity -- instant stops and turns, crew
        at 0 g; shown is what a conventional thrust (or a field gradient of the craft) would have to supply for it.
        """
        for b in sim.bodies:
            if b.kind != "pilot" or not b.alive:
                continue
            dx = (1 if ({"Right", "d", "D"} & keys) else 0) - (
                1 if ({"Left", "a", "A"} & keys) else 0
            )
            dy = (1 if ({"Up", "w", "W"} & keys) else 0) - (
                1 if ({"Down", "s", "S"} & keys) else 0
            )
            if nullf_var.get():
                V = 10 ** (thrust_var.get() + 2.0) * (
                    10.0 if ({"Shift_L", "Shift_R"} & keys) else 1.0
                )  # m/s
                n = math.hypot(dx, dy)
                vt = (V * dx / n / C_SI, V * dy / n / C_SI) if n else (0.0, 0.0)
                b.vtarget, b.thrust = vt, (0.0, 0.0)
                gx = gy = 0.0  # Newtonovské gravitační zrychlení v místě řízeného tělesa [m/s²].
                for o in sim.bodies:
                    if o.alive and o.mass > 0 and o is not b:
                        ddx, ddy = b.x - o.x, b.y - o.y
                        r3 = (ddx * ddx + ddy * ddy) ** 1.5
                        gx -= o.mass * ddx / r3 * C_SI**2
                        gy -= o.mass * ddy / r3 * C_SI**2
                ax = (
                    (vt[0] - b.vx) * C_SI**2 / TAU_NULL
                )  # kinematic acceleration [m/s^2]
                ay = (vt[1] - b.vy) * C_SI**2 / TAU_NULL
                b.v_cmd = V if n else 0.0
                b.now_g = math.hypot(
                    ax - gx, ay - gy
                )  # what a thrust would have to supply
                b.max_g = max(b.max_g, b.now_g)
                continue
            b.vtarget = None
            mag = (
                9.81
                * 10 ** thrust_var.get()
                * (10.0 if ({"Shift_L", "Shift_R"} & keys) else 1.0)
            )
            n = math.hypot(dx, dy)
            ax, ay = (mag * dx / n, mag * dy / n) if n else (0.0, 0.0)
            if getattr(
                b, "hover", False
            ):  # Kompenzace místní gravitace v Newtonovské aproximaci [m/s²].
                g = sim.pull_split(b, None)[1] * C_SI**2
                ex = next((o for o in sim.bodies if o.radius > 0), None)
                if ex is not None:
                    r = math.hypot(b.x - ex.x, b.y - ex.y)
                    ax += g * (b.x - ex.x) / r
                    ay += g * (b.y - ex.y) / r
            b.thrust = (ax, ay)
            b.now_g = math.hypot(ax, ay)
            b.max_g = max(b.max_g, b.now_g)

    def tick(single=False):
        sim = state["sim"]
        now = time.perf_counter()
        real_dt = min(max(now - state.get("last_tick", now - 0.03), 0.001), 0.1)
        state["last_tick"] = now
        # Stejné ovládání ve všech scénách: čas simulace za sekundu reálného času.
        # SI scény používají interní délku ct, proto se sekundy násobí rychlostí světla.
        jednotky_ct = C_SI if state["meta"].get("si") else 1.0
        T = 10.0 ** speed_var.get() * real_dt * jednotky_ct
        if state["meta"].get(
            "map3d"
        ):  # Hvězdná mapa: každá hvězda se pohybuje svou prostorovou rychlostí.
            hw = state["view"][2]
            for b in sim.bodies:
                vx, vy, vz = getattr(b, "v3", (0.0, 0.0, 0.0))
                if vx or vy or vz:
                    b.x, b.y, b.z = b.x + vx * T, b.y + vy * T, b.z + vz * T
                    tr = b.trail3
                    if not tr or math.dist(tr[-1], (b.x, b.y, b.z)) > 0.004 * hw:
                        tr.append((b.x, b.y, b.z))
            sim.t += T
            state["achieved"] = T / C_SI / real_dt
            speed_value_label()
            draw(refresh_field=False)
            return
        if state["meta"].get("pilot"):
            pilot_control(sim)
        t_start = sim.t

        def on(sm):
            for tr in state["trackers"]:
                tr.update(sm)
            for lt in state["light"]:
                lt.update(sm)

        t0 = time.perf_counter()
        if single:
            sim.step()
            on(sim)
        else:
            steps = 0
            t_end = sim.t + T
            while (
                sim.t < t_end and time.perf_counter() - t0 < 0.025
            ):  # Zachová odezvu okna.
                live, pos, vel, mass, photon, _ = sim.arrays()
                dt = min(sim.choose_dt(pos, vel, mass, photon), t_end - sim.t)
                sim.step(max(dt, 1e-9))
                on(sim)
                steps += 1
        if not single:  # Skutečně dosažené tempo může být nižší kvůli výpočetní náročnosti.
            state["achieved"] = (sim.t - t_start) / (C_SI if state["meta"].get("si") else 1.0) / real_dt
            speed_value_label()
        if state["meta"].get("pilot") and follow_var.get():
            pl = next((b for b in sim.bodies if b.kind == "pilot" and b.alive), None)
            if pl is not None:
                state["view"] = (pl.x, pl.y, state["view"][2])
        if state["meta"].get(
            "bh"
        ):  # Foton hluboko v oblasti nekonečného rudého posuvu se zastaví.
            for sm in [sim] + [g["sim"] for g in state["ghosts"]]:
                c = next(
                    (b for b in sm.bodies if getattr(b, "is_profile", False)), None
                )
                for b in sm.bodies:
                    if b.kind == "photon" and b.alive and not b.fixed and c is not None:
                        if math.hypot(b.x - c.x, b.y - c.y) < 1.0:
                            b.fixed, b.captured = True, True
        # Srovnávací běhy: stejná scéna podle jiných zákonů, posunutá ke stejnému času vlastními kroky.
        t1 = time.perf_counter()
        for g in state["ghosts"]:
            gs = g["sim"]
            n = 0
            while gs.t < sim.t and n < 20000 and time.perf_counter() - t1 < 0.05:
                live, pos, vel, mass, photon, _ = gs.arrays()
                dt = min(gs.choose_dt(pos, vel, mass, photon), sim.t - gs.t)
                gs.step(max(dt, 1e-9))
                for tr in g["trackers"]:
                    tr.update(gs)
                for lt in g["light"]:
                    lt.update(gs)
                n += 1
        # Experiment s Merkurem: na jeden snímek připadá nejvýše přibližně 15 ms výpočtu.
        if state["exp"] is not None:
            t2 = time.perf_counter()
            try:
                while time.perf_counter() - t2 < 0.015:
                    lawname, frac = next(state["exp"])
                    exp_lbl.config(text=f"běží: {lawname}, {100 * frac:.0f} %")
            except StopIteration as stop:
                state["exp_result"], state["exp"] = stop.value, None
                exp_lbl.config(text="hotovo, výsledky ve výpisu dole")
        if sim.events and state.get("n_events", 0) != len(sim.events):
            state["n_events"] = len(sim.events)
            setup_trackers()  # Došlo ke sloučení tělesa nebo pohlcení fotonu.
        for sm in [sim] + [g["sim"] for g in state["ghosts"]]:
            for b in sm.bodies:
                if b.alive and not b.fixed:
                    if (
                        not b.trail
                        or math.hypot(b.x - b.trail[-1][0], b.y - b.trail[-1][1])
                        > state["view"][2] * 0.003
                    ):
                        b.trail.append((b.x, b.y))
        state["frame"] += 1
        draw()

    def loop():
        try:
            if state["running"]:
                tick()
            elif (
                state["exp"] is not None
            ):  # Experiment běží i při pozastavené scéně.
                tick_experiment_only()
        except (
            Exception
        ) as exc:  # noqa: BLE001 -- chyba jednoho snímku nesmí zastavit animaci.
            report_error("snímek", exc)
        finally:
            root.after(30, loop)  # always rescheduled

    # ---------- pomocné nástroje 3D mapy a ovládání myší
    def map3d_souradnice_z_kliknuti(sx, sy):
        """Převede kliknutí na bod v rovině pohledu 3D mapy; hloubku určí projekční rovina."""
        w, h = max(canvas.winfo_width(), 10), max(canvas.winfo_height(), 10)
        hw = state["view"][2]
        scale = 0.5 * w / hw
        yaw, pitch = state["rot"]
        x1 = (sx - w / 2) / scale
        z2 = (h / 2 - sy) / scale
        # Zvolíme bod s hloubkovou souřadnicí y2 = 0 v aktuální projekci.
        y1 = math.sin(pitch) * z2
        z = math.cos(pitch) * z2
        x = math.cos(yaw) * x1 + math.sin(yaw) * y1
        y = -math.sin(yaw) * x1 + math.cos(yaw) * y1
        return x, y, z

    def smaz_neblizsi_3d_objekt(sx, sy):
        """Odstraní nejbližší bod mapy do 16 obrazových bodů od kliknutí."""
        body = state.get("map3d_points", {})
        nejlepsi, vzdalenost = None, 16.0
        for px, py, b in body.values():
            d = math.hypot(px - sx, py - sy)
            if d < vzdalenost:
                nejlepsi, vzdalenost = b, d
        if nejlepsi is not None:
            if nejlepsi.name == "Slunce":
                report_error("úprava 3D mapy", ValueError("Slunce je počátek souřadnicové soustavy a nelze jej odstranit."))
                return
            state["sim"].remove(nejlepsi)
            update_category()
            draw(refresh_field=False)
            readouts()

    # ---------- ovládání myší
    def on_press(ev):
        if state["meta"].get("map3d"):
            nastroj = tool_var.get()
            if nastroj == "delete":
                smaz_neblizsi_3d_objekt(ev.x, ev.y)
            elif nastroj == "mass":
                state["map3d_add_click"] = (ev.x, ev.y)
            else:
                state["rot_drag"] = (ev.x, ev.y, list(state["rot"]))
            return
        if state.get("view3d") and tool_var.get() == "none":
            state["rot_drag"] = (ev.x, ev.y, list(state.get("rot", [0.25, 0.95])))
            return
        if tool_var.get() == "delete":  # Odstraní nejbližší objekt do vzdálenosti 15 px.
            sim = state["sim"]
            best, bd = None, 15.0
            for b in sim.bodies:
                if b.alive:
                    sx, sy = W2S(b.x, b.y)
                    d = math.hypot(sx - ev.x, sy - ev.y)
                    if d < bd:
                        best, bd = b, d
            if best is not None:
                sim.remove(best)  # Při odstranění černé díry se smaže i její přesné pole.
                for g in state[
                    "ghosts"
                ]:  # Odstraní těleso se stejným unikátním názvem ve všech srovnávacích bězích.
                    for gb in [b for b in g["sim"].bodies if b.name == best.name]:
                        g["sim"].remove(gb)
                setup_trackers()
                draw(force_field=True)
            return
        if tool_var.get() != "none":
            state["drag"] = ((ev.x, ev.y), (ev.x, ev.y))

    def on_motion(ev):
        if (state["meta"].get("map3d") or state.get("view3d")) and state.get("rot_drag"):
            x0, y0, (yw, pt) = state["rot_drag"]
            state["rot"] = [
                yw + (ev.x - x0) * 0.01,
                max(-1.45, min(1.45, pt + (ev.y - y0) * 0.01)),
            ]
            draw(refresh_field=False)
            return
        if state["drag"] is not None:
            state["drag"] = (state["drag"][0], (ev.x, ev.y))
            if not state["running"]:
                draw()

    def on_release(ev):
        if state["meta"].get("map3d"):
            zacatek = state.get("map3d_add_click")
            state["map3d_add_click"] = None
            if zacatek is not None:
                x3, y3, z3 = map3d_souradnice_z_kliknuti(ev.x, ev.y)
                map_xyz_vars[0].set(f"{x3 / LY:.8g}")
                map_xyz_vars[1].set(f"{y3 / LY:.8g}")
                map_xyz_vars[2].set(f"{z3 / LY:.8g}")
                add_map3d_object(poloha_m=(x3, y3, z3))
            state["rot_drag"] = None
            return
        if state.get("view3d") and state.get("rot_drag") is not None:
            state["rot_drag"] = None
            return
        if state["drag"] is None:
            return
        (x0, y0), (x1, y1) = state["drag"]
        state["drag"] = None
        sim = state["sim"]
        wx, wy = S2W(x0, y0)
        ex, ey = S2W(x1, y1)
        dvx, dvy = ex - wx, ey - wy
        tool = tool_var.get()
        n = len(sim.bodies)
        state["name_seq"] += 1
        seq = state["name_seq"]  # Pořadové číslo se po smazání nerecykluje, názvy zůstávají jedinečné.
        if tool in ("mass", "clock"):
            for o in [
                b for b in sim.bodies if b.alive and b.radius > 0
            ]:  # Při kliknutí dovnitř Země se objekt přesune na povrch.
                rx, ry = wx - o.x, wy - o.y
                r = math.hypot(rx, ry)
                if r < o.radius:
                    ux, uy = (rx / r, ry / r) if r > 0 else (0.0, 1.0)
                    wx, wy = o.x + ux * o.radius * (1 + 1e-9), o.y + uy * o.radius * (
                        1 + 1e-9
                    )
            if (
                sim.profile is not None
            ):  # Objekt nesmí být přesně ve středu černé díry, kde L → −∞.
                c = sim.profile
                rx, ry = wx - c.x, wy - c.y
                r = math.hypot(rx, ry)
                r_min = 0.05 * c.mass
                if r < r_min:
                    ux, uy = (rx / r, ry / r) if r > 0 else (1.0, 0.0)
                    wx, wy = c.x + ux * r_min, c.y + uy * r_min
        if tool == "mass":
            pomer = parse_mass()
            if pomer is None:
                report_error("vstup hmotnosti", ValueError("Poměr m/Mref musí být kladné konečné číslo."))
                state["name_seq"] -= 1
                return
            meta = state["meta"]
            m_len = pomer * referencni_hmotnost(sim)
            if meta["si"]:
                kg = m_len / G_C2
                name = f"{category_si(kg)} {seq}"
                if meta["click"] == "rest":  # Kliknutí pustí těleso; tažení mu udělí počáteční rychlost.
                    dpx, dpy = x1 - x0, y1 - y0
                    dl = max(math.hypot(dpx, dpy), 1e-9)
                    k = dl / 80.0 * 20.0 / C_SI
                    vx, vy = (k * dpx / dl, -k * dpy / dl) if dl >= 6 else (0.0, 0.0)
                else:
                    vx, vy, _ = new_body_velocity(wx, wy, m_len, x0, y0, x1, y1)
            else:
                name = f"{category_astro(m_len)} {seq}"
                vx, vy, _ = new_body_velocity(wx, wy, m_len, x0, y0, x1, y1)
            sim.bodies.append(Body(name, wx, wy, vx, vy, mass=m_len, color="#9ad9ff"))
            sim.bodies[-1].start = (wx, wy, vx, vy)
        elif tool == "photon":
            sp = math.hypot(dvx, dvy)
            if sp == 0:
                dvx, dvy, sp = 1.0, 0.0, 1.0
            sim.bodies.append(
                Body(
                    f"foton {seq}",
                    wx,
                    wy,
                    dvx / sp,
                    dvy / sp,
                    kind="photon",
                    color="#fff27a",
                )
            )
        elif tool == "clock":
            sim.bodies.append(
                Body(f"hodiny {seq}", wx, wy, kind="clock", color="#b0b0ff")
            )
            sim.bodies[-1].start = (wx, wy, 0.0, 0.0)
        sim.bodies[-1].t_birth = sim.t
        if len(sim.bodies) > n:  # Zkopíruje nový objekt do všech srovnávacích běhů.
            nb_ = sim.bodies[-1]
            for g in state["ghosts"]:
                cp = Body(
                    nb_.name,
                    nb_.x,
                    nb_.y,
                    nb_.vx,
                    nb_.vy,
                    mass=nb_.mass,
                    kind=nb_.kind,
                    color=nb_.color,
                )
                cp.t_birth, cp.start = g["sim"].t, getattr(
                    nb_, "start", (nb_.x, nb_.y, nb_.vx, nb_.vy)
                )
                g["sim"].bodies.append(cp)
        update_category()
        setup_trackers()
        draw(force_field=True)

    def add_map3d_object(poloha_m=None):
        """Přidá objekt do 3D mapy; souřadnice se zadávají v ly a rychlosti v km/s."""
        if not state.get("meta", {}).get("map3d"):
            return
        try:
            if poloha_m is None:
                xyz = tuple(float(v.get().strip().replace(",", ".")) * LY for v in map_xyz_vars)
            else:
                xyz = tuple(float(v) for v in poloha_m)
            rychlosti = tuple(float(v.get().strip().replace(",", ".")) * 1000.0 / C_SI for v in map_v_vars)
            if not all(math.isfinite(v) for v in (*xyz, *rychlosti)):
                raise ValueError("Souřadnice i rychlosti musí být konečná čísla.")
        except ValueError as exc:
            report_error("vstup objektu 3D mapy", exc)
            return
        pomer = parse_mass()
        if pomer is None:
            report_error("vstup hmotnosti", ValueError("Poměr m/Mref musí být kladné konečné číslo."))
            return
        m_len = pomer * referencni_hmotnost(state["sim"])
        nazev = map_name_var.get().strip() or f"Nové těleso {state['name_seq'] + 1}"
        if any(b.alive and b.name == nazev for b in state["sim"].bodies):
            report_error("vstup objektu 3D mapy", ValueError(f"Objekt s názvem ‚{nazev}‘ už existuje."))
            return
        b = Body(nazev, xyz[0], xyz[1], mass=m_len, fixed=True, color="#9ad9ff")
        b.z = xyz[2]
        b.v3 = rychlosti
        b.trail3 = deque(maxlen=4000)
        vzdalenost = math.sqrt(sum(v * v for v in xyz)) / LY
        b.info = f"uživatelský objekt, {m_len / M_SUN:.4g} M☉, {vzdalenost:.4g} ly"
        v_kms = tuple(v * C_SI / 1000.0 for v in rychlosti)
        b.motion_note = f"v = ({v_kms[0]:.4g}, {v_kms[1]:.4g}, {v_kms[2]:.4g}) km/s"
        b.t_birth = state["sim"].t
        b.start = (b.x, b.y, 0.0, 0.0)
        state["sim"].bodies.append(b)
        state["name_seq"] += 1
        map_name_var.set(f"Nové těleso {state['name_seq'] + 1}")
        update_category()
        draw(refresh_field=False)
        readouts()

    def new_body_velocity(wx, wy, m, x0, y0, x1, y1):
        """
        Velocity of a new body: a click without dragging gives a circular orbit (counter-clockwise) around the
        nejhmotnějšího zdroje. Tažení určí směr a násobek kruhové rychlosti (80 px = 1×).
        Bez zdroje se použije desetina délky tažení v souřadnicích scény.
        Vrací (vx, vy, násobek kruhové rychlosti nebo None).
        """
        sim = state["sim"]
        src = [b for b in sim.bodies if b.alive and b.mass > 0]
        dpx, dpy = x1 - x0, y1 - y0
        dlen = math.hypot(dpx, dpy)
        if not src:
            ex, ey = S2W(x1, y1)
            return 0.1 * (ex - wx), 0.1 * (ey - wy), None
        c = max(src, key=lambda b: b.mass)
        rx, ry = wx - c.x, wy - c.y
        r = max(math.hypot(rx, ry), 1e-6)
        vc = math.sqrt((c.mass + m) / r)
        if dlen < 6:  # Samotné kliknutí vytvoří kruhový oběh.
            return c.vx - vc * ry / r, c.vy + vc * rx / r, 1.0
        k = dlen / 80.0
        return (
            c.vx + vc * k * dpx / dlen,
            c.vy - vc * k * dpy / dlen,
            k,
        )  # Osa y obrazovky míří dolů.

    def on_wheel(ev):
        delta = (
            ev.delta
            if hasattr(ev, "delta") and ev.delta
            else (120 if ev.num == 4 else -120)
        )
        cx, cy, hw = state["view"]
        wx, wy = S2W(ev.x, ev.y)
        f = 0.85 if delta > 0 else 1 / 0.85
        state["view"] = (wx + (cx - wx) * f, wy + (cy - wy) * f, hw * f)
        draw(force_field=True)

    def on_pan_press(ev):
        state["pan"] = (ev.x, ev.y, state["view"])

    def on_pan_motion(ev):
        if state["pan"] is None:
            return
        x0, y0, (cx, cy, hw) = state["pan"]
        w = max(canvas.winfo_width(), 10)
        sc = 0.5 * w / hw
        state["view"] = (cx - (ev.x - x0) / sc, cy + (ev.y - y0) / sc, hw)
        draw(force_field=True)

    canvas.bind(
        "<ButtonPress-1>", safe("myš", lambda ev: (canvas.focus_set(), on_press(ev)))
    )
    canvas.bind("<B1-Motion>", safe("myš", on_motion))
    canvas.bind("<ButtonRelease-1>", safe("myš", on_release))
    canvas.bind("<MouseWheel>", safe("myš", on_wheel))
    canvas.bind("<Button-4>", safe("myš", on_wheel))
    canvas.bind("<Button-5>", safe("myš", on_wheel))
    canvas.bind("<ButtonPress-3>", safe("myš", on_pan_press))
    canvas.bind("<B3-Motion>", safe("myš", on_pan_motion))

    CONTROL_KEYS = {
        "Up",
        "Down",
        "Left",
        "Right",
        "w",
        "a",
        "s",
        "d",
        "W",
        "A",
        "S",
        "D",
        "space",
    }

    def on_key(ev):
        k = ev.keysym
        if not state["meta"].get("pilot"):
            return None
        if (
            k in CONTROL_KEYS and not state["running"]
        ):  # Scéna řízení se spustí prvním stiskem ovládací klávesy.
            toggle_run()
        if k == "space":
            if (
                "space" not in keys
            ):  # U podržené klávesy opakování nesmí přepnout režim znovu.
                for b in state["sim"].bodies:
                    if b.kind == "pilot":
                        b.hover = not getattr(b, "hover", False)
            keys.add("space")
            return "break"  # Stisk nesmí zároveň aktivovat právě zaostřené tlačítko.
        keys.add(k)
        return (
            "break" if k in CONTROL_KEYS else None
        )  # Šipky nesmějí měnit právě zaostřený posuvník ani seznam.

    def on_key_up(ev):
        keys.discard(ev.keysym)
        return (
            "break"
            if state["meta"].get("pilot") and ev.keysym in CONTROL_KEYS
            else None
        )

    # Fokus má plátno (po načtení scény nebo přejetí myší); vazby okna zachytí ostatní vstupy.
    canvas.bind("<KeyPress>", safe("klávesnice", on_key))
    canvas.bind("<KeyRelease>", safe("klávesnice", on_key_up))
    root.bind(
        "<KeyPress>",
        safe("klávesnice", lambda ev: on_key(ev) if ev.widget is not canvas else None),
    )
    root.bind(
        "<KeyRelease>",
        safe(
            "klávesnice", lambda ev: on_key_up(ev) if ev.widget is not canvas else None
        ),
    )
    root.bind(
        "<FocusOut>", lambda ev: keys.clear()
    )  # Při ztrátě fokusu okna se žádná klávesa nesmí tvářit jako stále stisknutá.
    canvas.bind("<Enter>", lambda ev: canvas.focus_set())

    def on_hover(ev):
        state["hover"] = (ev.x, ev.y)
        if not state["running"]:
            draw(refresh_field=False)

    canvas.bind("<Motion>", safe("myš", on_hover))
    canvas.bind("<Leave>", lambda e: state.update(hover=None))
    canvas.bind("<Configure>", safe("okno", lambda e: draw(force_field=True)))

    if start_scene is not None:
        scene_var.set(SCENES[start_scene])  # Změna proměnné zavolá load_scene přes sledování změny.
    else:
        load_scene(SCENES[0])
    if autorun:
        toggle_run()
    loop()
    root.mainloop()


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        raise SystemExit(selftest())
    # Volby: --scene N (0–15) otevře scénu podle indexu, --autorun ihned spustí simulaci.
    sc = None
    if "--scene" in sys.argv:
        try:
            sc = int(sys.argv[sys.argv.index("--scene") + 1]) % len(SCENES)
        except (IndexError, ValueError):
            sc = None
    run_app(sc, "--autorun" in sys.argv)
