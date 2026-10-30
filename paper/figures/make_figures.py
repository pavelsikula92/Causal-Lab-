"""
Figures of the paper, drawn from the saved records of the laboratory
=====================================================================

Every point in every figure is read from a record in results/<ID>/ (the newest full run of the test, quick runs are
ignored), never typed in by hand. The only numbers written here are those of analytic curves that the paper states
as formulas (the window of the non-locality scale, Sec. "The window of the non-locality scale"), and they are marked
as such. Error bars are +-1 standard error as stored in the record.

Run from the repository root (Windows cmd):
    cd /d D:\\Projekty\\causality_lab & python paper\\figures\\make_figures.py
Writes paper/figures/fig_*.pdf (vector, for LaTeX) and fig_*.png (preview), and prints the numbers each figure uses,
so that they can be compared with the text of the paper.

Figures (the section of the paper in brackets):
  fig_blind      the blind instrument: gamma read against gamma true, 2+1 and 3+1; potentials by radius (sec:reading)
  fig_response   gamma of the three responses of the order against the predictions (sec:response)
  fig_counter    counter-streaming flows: clock potential by radius, amplitude ratio against the rules (sec:rule)
  fig_conformal  the conformal sector: K_loc(k) with its fit; the age term collapses onto c_T k T (sec:bdaction)
  fig_sectors    the three sectors approach the Einstein-Hilbert values 6, 8, 2 as k -> 0 (sec:optical)
  fig_optical    the age term of the optical sector grows faster than linearly; the conformal one does not (sec:optical)
  fig_window     the window of the non-locality scale: L* against l_k for four readings of T (sec:window)
  fig_fall       free fall from the apex (F5.1d): displacement against time, and R = a/a_pred with the criterion (sec:fall)
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import matplotlib                                   # noqa: E402  (after the path set-up)

matplotlib.use("Agg")                               # files only, no window
import matplotlib.pyplot as plt                     # noqa: E402
import matplotlib.transforms                        # noqa: E402

from core.provenance import load_all                # noqa: E402

OUT = Path(__file__).resolve().parent

# ----------------------------------------------------------------------------------------------- style
# One text width of the paper (A4, margins 2.6 cm) is 15.8 cm = 6.22 in.
WIDTH = 6.22
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
# Categorical slots in fixed order (blue, orange, aqua); identity is never carried by colour alone: every series
# also has its own marker and a legend entry.
C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"
plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "cm", "font.size": 8.5, "axes.titlesize": 9,
    "axes.labelsize": 8.5, "legend.fontsize": 7.5, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.4, "lines.markersize": 5,
    "legend.frameon": False, "savefig.bbox": "tight", "savefig.pad_inches": 0.03,
})
ERR = dict(capsize=2, elinewidth=0.9, capthick=0.9)


def record(test_id: str, part: str | None = None, seed: int | None = None) -> dict:
    """The newest full-run record of a test (optionally of one part / seed); a clear error if there is none."""
    recs = [r for r in load_all(test_id) if not r.get("quick")]
    if part is not None:
        recs = [r for r in recs if r.get("part") == part]
    if seed is not None:
        recs = [r for r in recs if r.get("seed") == seed]
    if not recs:
        raise FileNotFoundError(f"no full-run record of {test_id} (part {part}, seed {seed}) in results/{test_id}/")
    rec = max(recs, key=lambda r: r.get("time_utc", ""))
    print(f"  {test_id}: {rec['_file']} ({rec['status']})")
    return rec


def save(fig, name: str):
    """Vector PDF for the paper and a PNG preview."""
    fig.savefig(OUT / f"{name}.pdf")
    fig.savefig(OUT / f"{name}.png", dpi=200)
    plt.close(fig)
    print(f"  -> paper/figures/{name}.pdf")


def mids(bins):
    """Centres of radial bins."""
    b = np.asarray(bins, float)
    return 0.5 * (b[1:] + b[:-1])


def kfit(k, K0, c1, c0):
    """The fit form of the paper: K = K0 - k^2 (c1 ln(1/k) + c0)."""
    k = np.asarray(k, float)
    return K0 - k * k * (c1 * np.log(1.0 / k) + c0)


# ===================================================================================== fig_blind (K.16b, K.16c)

def fig_blind():
    print("fig_blind: the blind instrument")
    r2, r3 = record("K.16b"), record("K.16c")
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(WIDTH, 2.7), gridspec_kw={"width_ratios": [1, 1.25], "wspace": 0.38})
    # (a) gamma read against gamma true; small horizontal offsets keep the two paths apart
    ax.plot([-1.3, 1.3], [-1.3, 1.3], color=INK2, lw=0.8, ls=(0, (4, 3)), label="read = true")
    for rec, colour, lab, cfgs in ((r2, C1, "2+1", ("g-1", "g0", "g+1")), (r3, C2, "3+1", ("g-1", "g1/3", "g+1"))):
        D = rec["data"]
        gam = {c[0]: float(c[2]) for c in rec["params"]["configs"]}     # name, amplitude a, gamma of the metric
        for path, mk, dx in (("gA", "o", -0.03), ("gB", "s", 0.03)):
            g = [gam[c] for c in cfgs]
            v = [D[c][path][0] for c in cfgs]
            e = [D[c][path][1] for c in cfgs]
            ax.errorbar(np.array(g) + dx, v, yerr=e, fmt=mk, color=colour, mfc=colour if path == "gA" else "white",
                        label=f"{lab}, path {'A' if path == 'gA' else 'B'}", **ERR)
            print(f"    {lab} {path}: " + ", ".join(f"{c} {x:+.3f} +- {s:.3f}" for c, x, s in zip(cfgs, v, e)))
    wk = r2["data"]["g+1_weak"]                       # the weak field, reported as unreliable in the paper
    for path, mk in (("gA", "o"), ("gB", "s")):
        ax.errorbar(1.0 + (0.09 if path == "gA" else 0.15), wk[path][0], yerr=wk[path][1], fmt=mk, color=INK2,
                    mfc="white", **ERR)
    ax.annotate(r"weak field $|\phi|=0.03$" + "\n(unreliable)", (1.12, wk["gB"][0]), xytext=(0.30, -0.75),
                textcoords="data", fontsize=7, color=INK2, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6))
    ax.set_xlabel(r"$\gamma$ of the metric")
    ax.set_ylabel(r"$\gamma$ read from the order")
    ax.set_xlim(-1.3, 1.35)
    ax.set_ylim(-1.3, 1.35)
    ax.set_title("(a) calibration", loc="left")
    ax.legend(loc="upper left", handlelength=1.6)
    # (b) potentials by radius for gamma = +1 in 2+1: true profile and the two reading paths
    D, P = r2["data"]["g+1"], r2["params"]
    r = mids(P["bins"])
    rr = np.linspace(0.0, P["bins"][-1], 200)
    cfg = {c[0]: (float(c[1]), float(c[2])) for c in P["configs"]}["g+1"]
    a, w, g = cfg[0], P["width"], cfg[1]
    phi = -a * np.exp(-rr ** 2 / (2 * w * w))
    bx.plot(rr, phi, color=INK, lw=1.1, label=r"$\phi=\psi/\gamma$ of the metric")
    for key, colour, nm in (("phi", C1, r"\phi"), ("psi", C2, r"\psi")):
        for path, mk, dx in (("A", "o", -0.15), ("B", "s", 0.15)):
            m, s = D[f"{key}{path}"]
            bx.errorbar(r + dx, m, yerr=s, fmt=mk, color=colour, mfc=colour if path == "A" else "white", ms=4,
                        label=f"${nm}$ read, path {path}", **ERR)
    bx.set_xlabel(r"$r/\ell$")
    bx.set_ylabel("potential")
    bx.set_title(r"(b) 2+1, $\gamma=+1$, $a=0.1$", loc="left")
    bx.legend(loc="upper left", ncol=1, handlelength=1.4, fontsize=6.8)
    save(fig, "fig_blind")


# ===================================================================================== fig_response (F7.1a/b/c)

def fig_response():
    print("fig_response: gamma of the three responses")
    a = record("F7.1a")
    b = [record("F7.1b", seed=s) for s in (0, 1)]
    c = [record("F7.1c", seed=s) for s in (0, 1)]
    groups = [  # label, prediction, [(record, config)], control (record, config) or None
        ("cones only\n2+1 (F7.1a)", 0.5, [(a, "capacity")], None),
        ("throughput\n2+1 (F7.1b)", 2.0, [(r_, "throughput") for r_ in b], [(r_, "ctrl_g2") for r_ in b]),
        ("cones only\n3+1 (F7.1c)", 1.0 / 3.0, [(r_, "cones_only") for r_ in c], [(r_, "ctrl_g13") for r_ in c]),
        ("throughput\n3+1 (F7.1c)", 1.0, [(r_, "throughput") for r_ in c], [(r_, "ctrl_g1") for r_ in c]),
    ]
    fig, ax = plt.subplots(figsize=(WIDTH, 2.7))
    seen = set()
    for i, (lab, pred, runs, ctrls) in enumerate(groups):
        ax.plot([i - 0.38, i + 0.38], [pred, pred], color=INK, lw=1.6, label="prediction" if "p" not in seen else None)
        seen.add("p")
        pts = [(r_, cfg, path) for r_, cfg in runs for path in ("gA", "gB")]
        xs = np.linspace(i - 0.27, i + (0.02 if ctrls else 0.27), len(pts)) if len(pts) > 1 else [i - 0.1]
        for x, (r_, cfg, path) in zip(xs, pts):
            m, s = r_["data"][cfg][path]
            key = f"order, path {'A' if path == 'gA' else 'B'}"
            ax.errorbar(x, m, yerr=s, fmt="o" if path == "gA" else "s", color=C1, mfc=C1 if path == "gA" else "white",
                        label=key if key not in seen else None, **ERR)
            seen.add(key)
            print(f"    {lab.split(chr(10))[0]} {lab.split(chr(10))[1]} seed {r_['seed']} {path}: {m:+.3f} +- {s:.3f}")
        if ctrls:
            cp = [(r_, cfg, path) for r_, cfg in ctrls for path in ("gA", "gB")]
            for x, (r_, cfg, path) in zip(np.linspace(i + 0.1, i + 0.32, len(cp)), cp):
                m, s = r_["data"][cfg][path]
                key = f"exact metric (control), path {'A' if path == 'gA' else 'B'}"
                ax.errorbar(x, m, yerr=s, fmt="o" if path == "gA" else "s", color=C2,
                            mfc=C2 if path == "gA" else "white", label=key if key not in seen else None, **ERR)
                seen.add(key)
                print(f"      control {cfg} seed {r_['seed']} {path}: {m:+.3f} +- {s:.3f}")
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([g_[0] for g_ in groups])
    ax.set_ylabel(r"$\gamma$ read from the order")
    ax.set_xlim(-0.5, len(groups) - 0.5)
    ax.grid(axis="x", visible=False)
    ax.legend(loc="upper right", ncol=1, fontsize=7)
    save(fig, "fig_response")


# ===================================================================================== fig_counter (F7.4a)

def fig_counter():
    print("fig_counter: counter-streaming flows")
    rec = record("F7.4a")
    D, P = rec["data"], rec["params"]
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(WIDTH, 2.6), gridspec_kw={"width_ratios": [1.1, 1], "wspace": 0.75})
    # (a) clock potential across the flows: read (points) against the tensor rule (line)
    r = mids(P["bins"])
    for cfg, colour, mk, lab in (("two_v0", C1, "o", r"two at rest"), ("two_v06", C2, "s", r"two at $\pm0.6c$"),
                                 ("single_v06", C3, "D", r"one at $0.6c$")):
        m, s = D["configs"][cfg]["phic"]
        ax.errorbar(r, m, yerr=s, fmt=mk, color=colour, ms=4, label=lab, **ERR)
        ax.plot(r, D["predicted"][cfg]["phic"], color=colour, lw=1.1)
        print(f"    {cfg}: A = {D['configs'][cfg]['A'][0]:.3f} +- {D['configs'][cfg]['A'][1]:.3f}")
    ax.plot([], [], color=INK2, lw=1.1, label="tensor rule")
    ax.set_xlabel(r"distance from the axis of the flows, $\rho/\ell$")
    ax.set_ylabel(r"clock potential $\phi$")
    ax.set_title("(a) clock field across the flows", loc="left")
    ax.legend(loc="lower right", fontsize=6.8, handlelength=1.4)
    # (b) the amplitude ratio against the candidate rules
    cands = D["candidates"]
    names = list(cands.keys())
    vals = [cands[n_] for n_ in names]
    short = {"tensor (exact, this rule)": "tensor, exact", "tensor (linear)": "tensor, linear"}
    y = np.arange(len(names))
    R, se = D["R"]
    bx.axvspan(R - se, R + se, color=C2, alpha=0.22, lw=0)
    bx.axvline(R, color=C2, lw=1.4, label=f"read: $R={R:.3f}\\pm{se:.3f}$")
    bx.plot(vals, y, "o", color=INK, ms=5, label="prediction of the rule")
    bx.set_yticks(y)
    bx.set_yticklabels([short.get(n_, n_) for n_ in names])
    bx.invert_yaxis()
    bx.set_xlabel(r"$R=A(0.6c)/A(0)$")
    bx.set_title("(b) which rule", loc="left")
    bx.grid(axis="y", visible=False)
    bx.legend(loc="lower right", fontsize=7)
    print(f"    R = {R:.3f} +- {se:.3f}; candidates " + ", ".join(f"{n_} {v:.4f}" for n_, v in zip(names, vals)))
    save(fig, "fig_counter")


# ===================================================================================== fig_conformal (F6.8a)

def fig_conformal():
    print("fig_conformal: the conformal sector")
    rec = record("F6.8a")
    D, P = rec["data"], rec["params"]
    ks = sorted(float(k) for k in D["per_k"])
    Kloc = np.array([_row(D, k)["K_loc"] for k in ks])
    F = D["fit"]
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(WIDTH, 2.6))
    # (a) the local part K_loc(k) with the fit over k <= fit_kmax and its continuation
    kk = np.geomspace(0.004, 0.6, 300)
    inside = kk <= float(P["fit_kmax"])
    ax.plot(kk[inside], kfit(kk[inside], F["K0"], F["c1"], F["c0"]), color=C1, lw=1.3,
            label=fr"fit, $k\ell_k\leq{P['fit_kmax']}$")
    ax.plot(kk[~inside], kfit(kk[~inside], F["K0"], F["c1"], F["c0"]), color=C1, lw=1.0, ls=(0, (3, 2)),
            label="fit continued (not fitted)")
    ax.plot(ks, Kloc, "o", color=INK, ms=4.5, label=r"$K_{\rm loc}$ from the action")
    # the sign change, by linear interpolation between the two computed points around it
    j = int(np.flatnonzero(np.diff(np.sign(Kloc)) != 0)[0])
    kz = ks[j] + (ks[j + 1] - ks[j]) * Kloc[j] / (Kloc[j] - Kloc[j + 1])
    ax.axhline(0.0, color=INK2, lw=0.7)
    ax.axhline(1.0, color=INK2, lw=0.7, ls=(0, (1, 2)))
    ax.annotate(fr"sign change $k\ell_k\simeq{kz:.2f}$", (kz, 0.0), xytext=(0.035, 0.42), fontsize=7, color=INK2,
                arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6))
    ax.set_xscale("log")
    ax.set_xlabel(r"$k\ell_k$")
    ax.set_ylabel(r"$K_{\rm loc}$ (Einstein--Hilbert $=1$)")
    ax.set_ylim(-0.12, 1.08)
    ax.set_title(fr"(a) local part, $K_0={F['K0']:.5f}$", loc="left")
    ax.legend(loc="lower left", fontsize=6.8)
    print(f"    K0 = {F['K0']:.5f} (6 K0 = {6 * F['K0']:.4f}), c1 = {F['c1']:.3f}, c0 = {F['c0']:.3f}, rms {F['rms']:.1e}; "
          f"sign change at k = {kz:.3f}")
    # (b) the age term: (K_loc - K_hat) / k^2 against k T collapses onto c_T k T for all k
    cT = []
    for i, k in enumerate(ks):
        row = _row(D, k)
        T, Kh = np.array(row["T"]), np.array(row["Khat"])
        x, yv = k * T, (row["K_loc"] - Kh) / k ** 2
        bx.plot(x, yv, "o", color=plt.cm.Blues(0.35 + 0.6 * i / (len(ks) - 1)), ms=4,
                label=(fr"$k\ell_k={ks[0]:g}$ (light) to ${ks[-1]:g}$ (dark)" if i == 0 else None))
        cT.append(row["c_T"])
    cTm = cT_conformal()
    xx = np.linspace(0.0, 42.0, 10)
    bx.plot(xx, cTm * xx, color=C2, lw=1.2, label=fr"$c_T\,kT$, $c_T={cTm:.4f}$")
    bx.set_xlabel(r"$kT$")
    bx.set_ylabel(r"$(K_{\rm loc}-\hat K)/k^2$")
    bx.set_title("(b) the age term is linear in $T$", loc="left")
    bx.legend(loc="upper left", fontsize=6.8)
    print(f"    c_T (mean of the per-k slopes, 0.05 <= k <= 0.2) = {cTm:.4f}; per k: " + ", ".join(f"{c:.4f}" for c in cT))
    save(fig, "fig_conformal")


def _row(D, k):
    """per_k entries are stored with the wavenumber as a string key."""
    for kk in (str(k), repr(k), f"{k:g}"):
        if kk in D["per_k"]:
            return D["per_k"][kk]
    raise KeyError(k)


def _per_k(D, key, k):
    return _row(D, k)[key]


def cT_conformal():
    """c_T of the conformal sector: the mean of the per-k slopes of F6.8a over 0.05 <= k l_k <= 0.2 (the two smallest
    k carry a finite-age effect, the largest the discreteness; the paper quotes c_T = 0.1035 +- 0.0003)."""
    D = record("F6.8a")["data"]
    ks = sorted(float(k) for k in D["per_k"])
    return float(np.mean([_row(D, k)["c_T"] for k in ks if 0.05 <= k <= 0.2]))


# ===================================================================================== fig_sectors (F6.8a, b2, c2)

def fig_sectors():
    print("fig_sectors: the three sectors")
    A, B, C = record("F6.8a"), record("F6.8b2"), record("F6.8c2")
    fig, ax = plt.subplots(figsize=(WIDTH, 2.9))
    kk = np.geomspace(0.003, 0.32, 300)
    # conformal: local part and its fit
    Da, Fa = A["data"], A["data"]["fit"]
    ka = sorted(float(k) for k in Da["per_k"] if float(k) <= 0.3)
    ya = [_per_k(Da, "K_loc", k) for k in ka]
    # mixed: local part and its fit
    Db, Fb = B["data"], B["data"]["fit"]
    kb = sorted(float(k) for k in Db["per_k"] if float(k) <= 0.3)
    yb = [_per_k(Db, "Kloc_ab", k) for k in kb]
    # optical: K_hat at fixed n = kT = judge_n, and the fit at that n
    Dc, n0 = C["data"], float(C["params"]["judge_n"])
    rows = sorted((r_ for r_ in Dc["rows"] if abs(r_["n"] - n0) < 1e-9), key=lambda r_: r_["k"])
    kc, yc = [r_["k"] for r_ in rows], [r_["Kbb"] for r_ in rows]
    Fc = Dc["fits"][str(n0)]["bb"] if str(n0) in Dc["fits"] else Dc["fits"][f"{n0:g}"]["bb"]
    series = [
        (r"conformal $\mathcal{A}$ (EH 6)", ka, ya, (Fa["K0"], Fa["c1"], Fa["c0"]), float(A["params"]["fit_kmax"]), C1, "o", 6),
        (r"mixed $\mathcal{B}$ (EH 8)", kb, yb, (Fb["K0_ab"], Fb["c1"], Fb["c0"]), float(B["params"]["fit_kmax"]), C2, "s", 8),
        (fr"optical $\mathcal{{C}}$ (EH 2), at $kT={n0:g}$", kc, yc, tuple(Fc[:3]), max(kc), C3, "D", 2),
    ]
    for lab, k_, y_, (K0, c1, c0), kmax, colour, mk, eh in series:
        inside = kk <= kmax
        ax.plot(kk[inside], kfit(kk[inside], K0, c1, c0), color=colour, lw=1.2)
        ax.plot(kk[~inside], kfit(kk[~inside], K0, c1, c0), color=colour, lw=0.9, ls=(0, (3, 2)))
        ax.plot(k_, y_, mk, color=colour, ms=5, mec="white", mew=0.6,
                label=fr"{lab}: $k\to0$ gives {eh}$\times${K0:.5f} $=$ {eh * K0:.4f}")
        print(f"    {lab}: K0 = {K0:.5f} -> {eh * K0:.4f} (EH {eh}); points " +
              ", ".join(f"k {a_:g}: {b_:.4f}" for a_, b_ in zip(k_, y_)))
    ax.axhline(1.0, color=INK2, lw=0.8, ls=(0, (1, 2)))
    ax.set_xscale("log")
    ax.set_xlabel(r"$k\ell_k$")
    ax.set_ylabel(r"stiffness / Einstein--Hilbert value")
    ax.set_ylim(0.55, 1.03)
    ax.legend(loc="lower left", fontsize=7.2)
    save(fig, "fig_sectors")


# ===================================================================================== fig_optical (F6.8c2)

def fig_optical():
    print("fig_optical: the optical age term")
    rec = record("F6.8c2")
    D = rec["data"]
    rows = list(D["rows"]) + [dict(r_, n=r_["k"] * r_["T"]) for r_ in D["extra"]]
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(WIDTH, 2.6))
    # (a) one wavenumber over the whole range of ages: 1 - K_hat for both sectors
    kx = float(D["extra"][0]["k"])
    sel = sorted((r_ for r_ in rows if abs(r_["k"] - kx) < 1e-12), key=lambda r_: r_["n"])
    n = np.array([r_["n"] for r_ in sel])
    for key, colour, mk, lab in (("Kbb", C3, "D", r"optical $\mathcal{C}$"), ("Kc", C1, "o", r"conformal $\mathcal{A}$")):
        y = 1.0 - np.array([r_[key] for r_ in sel])
        ax.plot(n, y, mk + "-", color=colour, lw=0.9, ms=5, label=lab)
        print(f"    k = {kx:g} {key}: " + ", ".join(f"kT {a_:g}: 1-K {b_:.4f}" for a_, b_ in zip(n, y)))
    ax.set_xlabel(r"$kT$")
    ax.set_ylabel(r"$1-\hat K$")
    ax.set_title(fr"(a) $k\ell_k={kx:g}$", loc="left")
    ax.legend(loc="upper left")
    # (b) the age part per k^2, measured from the shortest age: conformal collapses onto c_T, optical does not
    ks = sorted({r_["k"] for r_ in rows})
    for key, colour, mk, lab in (("Kbb", C3, "D", r"optical $\mathcal{C}$"), ("Kc", C1, "o", r"conformal $\mathcal{A}$")):
        first = True
        for k in ks:
            s_ = sorted((r_ for r_ in rows if r_["k"] == k), key=lambda r_: r_["n"])
            n_ = np.array([r_["n"] for r_ in s_])
            K_ = np.array([r_[key] for r_ in s_])
            bx.plot(n_ - n_[0], (K_[0] - K_) / k ** 2, mk, color=colour, ms=4, mec="white", mew=0.5,
                    label=lab if first else None)
            first = False
    cT = cT_conformal()
    xx = np.linspace(0, 18.5, 10)
    bx.plot(xx, cT * xx, color=INK2, lw=1.0, ls=(0, (4, 2)), label=fr"$c_T\,\Delta(kT)$, $c_T={cT:.4f}$ (F6.8a)")
    bx.set_xlabel(r"$\Delta(kT)$ from the shortest age")
    bx.set_ylabel(r"$[\hat K(kT_{\min})-\hat K]/k^2$")
    bx.set_title("(b) all wavenumbers", loc="left")
    bx.legend(loc="upper left", fontsize=6.8)
    print(f"    reference slope c_T (F6.8a) = {cT:.4f}")
    save(fig, "fig_optical")


# ===================================================================================== fig_window (analytic)

def fig_window():
    print("fig_window: the window of the non-locality scale (analytic, Eq. Lstar)")
    # Numbers stated in the paper (Sec. window): c_T of the conformal sector and the four readings of T in metres.
    c_T = 0.1035
    readings = [("age of the universe", 1.3e26, C1, "-"), ("age of Earth's field", 4.3e25, C2, (0, (5, 2))),
                ("one year", 9.5e15, C3, (0, (2, 1.5))), ("one day", 2.6e13, INK2, (0, (5, 1.5, 1, 1.5)))]
    L_test, delta = 52e-6, 1e-2                         # shortest distance at which gravity is measured; 1 % deviation
    l_P, l_LHC, l_opto = 1.616e-35, 1e-19, 2e-15
    lk = np.geomspace(1e-36, 1e-13, 400)
    fig, ax = plt.subplots(figsize=(WIDTH, 3.0))
    L_bound = L_test * delta ** (1.0 / 3.0)            # delta = (L*/L)^3 < 1 % at L = 52 um  <=>  L* < L_bound
    ax.axhspan(L_bound, 1e-2, color=C2, alpha=0.10, lw=0)
    ax.axhline(L_bound, color=C2, lw=0.9)
    ax.text(2e-36, L_bound * 2.0, fr"excluded: deviation $>1\%$ at {L_test * 1e6:.0f} $\mu$m ($L_*>{L_bound * 1e6:.0f}\,\mu$m)",
            fontsize=7, color=INK2, va="bottom")
    top = matplotlib.transforms.blended_transform_factory(ax.transData, ax.transAxes)
    for x, lab in ((l_P, r"$\ell_P$"), (l_LHC, r"LHC"), (l_opto, r"opto-mech.")):
        ax.axvline(x, color=INK2, lw=0.7, ls=(0, (1, 2)))
        ax.text(x, 1.01, lab, transform=top, fontsize=7, color=INK2, ha="center", va="bottom")
    for lab, T, colour, ls in readings:
        Ls = (c_T * lk ** 2 * T) ** (1.0 / 3.0)
        ax.plot(lk, Ls, color=colour, ls=ls, lw=1.3, label=fr"$T$ = {lab} ({T:.1e} m)".replace("e+", r"$\times10^{").replace(" m)", "}$ m)"))
        print(f"    {lab}: L*(l_P) = {(c_T * l_P ** 2 * T) ** (1 / 3):.2e} m, L*(1e-19 m) = {(c_T * l_LHC ** 2 * T) ** (1 / 3):.2e} m, "
              f"bound l_k < {math.sqrt(delta * L_test ** 3 / (c_T * T)):.1e} m")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e-36, 1e-13)
    ax.set_ylim(1e-21, 1e-2)
    ax.set_xlabel(r"non-locality scale $\ell_k$ (m)")
    ax.set_ylabel(r"$L_*=(c_T\,\ell_k^2\,T)^{1/3}$ (m)")
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.83), fontsize=7)
    save(fig, "fig_window")


def fig_fall():
    """F5.1d: (a) radial displacement from the apex against time, loaded and flat order (mean of the per-sprinkling
    series, band +-1 SE across sprinklings), with the Newtonian fall from rest at the apex radius (the prediction of
    the test) and the same fall shifted by the fitted start velocity v0; (b) R = a/a_pred and the control R_max with
    the criterion band 1/2 < R < 2."""
    from phases.F5_gravity.F5_01_free_fall_order import fall_curve
    rec = record("F5.1d")
    D, p = rec["data"], rec["params"]
    tau = np.linspace(0.0, float(p["span"]), 23)               # the time grid of the test (F5_01d, seed_task)
    keys = sorted(D["res"], key=float)                          # "0.0" flat, "0.45" loaded
    flat_k, load_k = keys[0], keys[-1]
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(WIDTH, 2.6), gridspec_kw={"width_ratios": [1.6, 1]})
    for key, col, mk, lab in [(load_k, C1, "o", "loaded order"), (flat_k, C2, "s", "flat order")]:
        rows = [r for r in D["rows"] if abs(float(r["sigma0"]) - float(key)) < 1e-9]
        S = np.array([r["series"] for r in rows], float)
        m = np.nanmean(S, axis=0)
        se = np.nanstd(S, axis=0, ddof=1) / np.sqrt(np.sum(np.isfinite(S), axis=0))
        ax.fill_between(tau, m - se, m + se, color=col, alpha=0.18, lw=0)
        ax.plot(tau, m, color=col, marker=mk, ms=3, lw=1.1, label=lab)
        print(f"    {lab}: displacement at tau = {tau[-1]:g}: {m[-1]:+.3f} +- {se[-1]:.3f} ({len(rows)} sprinklings)")
    rows = [r for r in D["rows"] if abs(float(r["sigma0"]) - float(load_k)) < 1e-9]
    pred = np.mean([fall_curve(dict(p, r0=float(r["r_apex"]), t0=0.0), tau) for r in rows], axis=0)
    v0 = float(D["res"][load_k]["v0"][0])
    ax.plot(tau, pred, color=INK2, ls=":", lw=1.2, label="Newtonian fall from rest (prediction)")
    ax.plot(tau, pred + v0 * tau, color=INK, ls="--", lw=1.0, label=f"same, with fitted $v_0$ = {v0:+.3f}")
    print(f"    prediction at tau = {tau[-1]:g}: {pred[-1]:+.3f}; with v0: {pred[-1] + v0 * tau[-1]:+.3f}")
    ax.axhline(0, color=INK2, lw=0.6)
    ax.set_xlabel(r"time from the apex $\tau$ [$\ell$]")
    ax.set_ylabel(r"radial displacement [$\ell$]")
    ax.set_title("(a) fall from the apex", loc="left")
    ax.legend(loc="lower left")
    # (b) R and the control
    bx.axhspan(0.5, 2.0, color=C3, alpha=0.12, lw=0)
    bx.axhline(1.0, color=INK2, lw=0.8, ls=":")
    items = [("$R$ loaded", D["res"][load_k]["R"], C1, "o"), ("$R$ flat", D["res"][flat_k]["R"], C2, "s"),
             (r"$R_{\max}$ loaded", D["res"][load_k]["R_max"], C1, "D"),
             (r"$R_{\max}$ flat", D["res"][flat_k]["R_max"], C2, "D")]
    for i, (lab, (v, e), col, mk) in enumerate(items):
        bx.errorbar([i], [v], yerr=[e], color=col, marker=mk, ls="none", **ERR)
        print(f"    {lab}: {v:+.3f} +- {e:.3f}")
    bx.set_xticks(range(len(items)))
    bx.set_xticklabels([it[0] for it in items], rotation=30, ha="right")
    bx.set_ylabel(r"$a/a_{\mathrm{pred}}$")
    bx.set_ylim(-0.4, 2.2)
    bx.text(3.35, 1.95, "criterion band", ha="right", va="top", fontsize=7, color=INK2)
    bx.set_title("(b) acceleration ratio", loc="left")
    save(fig, "fig_fall")


# ===================================================================================== main

FIGURES = [fig_blind, fig_response, fig_counter, fig_conformal, fig_sectors, fig_optical, fig_window, fig_fall]

if __name__ == "__main__":
    only = set(sys.argv[1:])                           # e.g. "python make_figures.py fig_window" draws one figure
    failed = []
    for f in FIGURES:
        if only and f.__name__ not in only:
            continue
        try:
            f()
        except (FileNotFoundError, KeyError) as e:     # a missing record skips that figure, never the others
            failed.append((f.__name__, e))
            print(f"  SKIPPED {f.__name__}: {e}")
    if failed:
        print("figures skipped: " + ", ".join(n_ for n_, _ in failed))
