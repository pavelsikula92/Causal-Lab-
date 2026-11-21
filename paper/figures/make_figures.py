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
  fig_fall       free fall from the apex, F5.1d (2+1) and F5.1e (3+1): displacement against time, R = a/a_pred (sec:fall)
  fig_matter     toward matter: the width of a bundle in its own field (F3.12b) and the cost of a locked defect (F3.13b2)
                 (sec:matter)
  fig_defect     matter as a defect of a field of the order: inertia (F3.14b, F3.15), swerves (F3.14c), fall towards a
                 source in proportion to its energy (F5.2b) (sec:fielddefect)
  fig_chain      the chain of the field equation of the order: what is computed, from which tests, with which status
                 (a diagram; the numbers in it are read from the records) (sec:variation)
  fig_fieldeq    the linear field equation of the order: stiffness at kT = 4 (F6.8g), source form factors (F6.10b), and
                 the response computed from them (diag/ov3_lin.py) (sec:variation)
"""

from __future__ import annotations

import math
import re
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


def _fall_panel(ax, rec, fall, title, legend):
    """One displacement panel of fig_fall: loaded and flat order (mean of the per-sprinkling series, band +-1 SE across
    sprinklings), the Newtonian fall from rest at the apex radius (the prediction of the test) and the same fall
    shifted by the fitted start velocity v0."""
    D, p = rec["data"], rec["params"]
    keys = sorted(D["res"], key=float)                          # "0.0" flat, the other loaded
    flat_k, load_k = keys[0], keys[-1]
    rows_l = [r for r in D["rows"] if abs(float(r["sigma0"]) - float(load_k)) < 1e-9]
    tau = np.linspace(0.0, float(p["span"]), len(rows_l[0]["series"]))   # the time grid of the test
    for key, col, mk, lab in [(load_k, C1, "o", "loaded order"), (flat_k, C2, "s", "flat order")]:
        rows = [r for r in D["rows"] if abs(float(r["sigma0"]) - float(key)) < 1e-9]
        S = np.array([r["series"] for r in rows], float)
        m = np.nanmean(S, axis=0)
        se = np.nanstd(S, axis=0, ddof=1) / np.sqrt(np.sum(np.isfinite(S), axis=0))
        ax.fill_between(tau, m - se, m + se, color=col, alpha=0.18, lw=0)
        ax.plot(tau, m, color=col, marker=mk, ms=2.5, lw=1.1, label=lab)
        print(f"    {title} {lab}: displacement at tau = {tau[-1]:g}: {m[-1]:+.3f} +- {se[-1]:.3f} ({len(rows)} sprinklings)")
    pred = np.mean([fall(dict(p, r0=float(r["r_apex"]), t0=0.0), tau) for r in rows_l], axis=0)
    v0 = float(D["res"][load_k]["v0"][0])
    ax.plot(tau, pred, color=INK2, ls=":", lw=1.2, label="prediction: fall from rest")
    ax.plot(tau, pred + v0 * tau, color=INK, ls="--", lw=1.0, label=r"prediction + fitted $v_0\tau$")
    print(f"    {title} prediction at tau = {tau[-1]:g}: {pred[-1]:+.3f}; with v0 = {v0:+.3f}: {pred[-1] + v0 * tau[-1]:+.3f}")
    ax.axhline(0, color=INK2, lw=0.6)
    ax.set_xlabel(r"time from the apex $\tau$ [$\ell$]")
    ax.set_title(title, loc="left")
    return D["res"][load_k], D["res"][flat_k]


def fig_fall():
    """Free fall from the apex: (a) F5.1d in 2+1, (b) F5.1e in 3+1, displacement from the apex against time (see
    _fall_panel); (c) R = a/a_pred and the control R_max of both tests with the criterion band 1/2 < R < 2."""
    from phases.F5_gravity.F5_01_free_fall_order import fall_curve
    from phases.F5_gravity.F5_01e_free_fall_apex_3d import fall_curve3
    r2, r3 = record("F5.1d", seed=0), record("F5.1e")   # panel (a): the first run; the replication is quoted in the text
    fig, (ax, ay, bx) = plt.subplots(1, 3, figsize=(WIDTH, 2.9), gridspec_kw={"width_ratios": [1.25, 1.25, 1]})
    fig.subplots_adjust(wspace=0.42, bottom=0.30)
    l2, f2 = _fall_panel(ax, r2, fall_curve, "(a) 2+1 (F5.1d)", True)
    l3, f3 = _fall_panel(ay, r3, fall_curve3, "(b) 3+1 (F5.1e)", False)
    ax.set_ylabel(r"radial displacement [$\ell$]")
    h, lab = ax.get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=4, fontsize=7, bbox_to_anchor=(0.42, -0.02))
    bx.axhspan(0.5, 2.0, color=C3, alpha=0.12, lw=0)
    bx.axhline(1.0, color=INK2, lw=0.8, ls=":")
    bx.axhline(0.0, color=INK2, lw=0.6)
    items = [("$R$ 2+1", l2["R"], C1, "o"), (r"$R_{\max}$ 2+1", l2["R_max"], C1, "D"),
             ("$R$ 3+1", l3["R"], C1, "^"), (r"$R_{\max}$ 3+1", l3["R_max"], C1, "D"),
             ("$R$ flat 2+1", f2["R"], C2, "o"), ("$R$ flat 3+1", f3["R"], C2, "^")]
    for i, (lab, (v, e), col, mk) in enumerate(items):
        bx.errorbar([i], [v], yerr=[e], color=col, marker=mk, ls="none", **ERR)
        print(f"    {lab}: {v:+.3f} +- {e:.3f}")
    bx.set_xticks(range(len(items)))
    bx.set_xticklabels([it[0] for it in items], rotation=40, ha="right", fontsize=6.5)
    bx.set_ylabel(r"$a/a_{\mathrm{pred}}$")
    bx.set_ylim(-0.5, 2.2)
    bx.text(5.4, 1.95, "criterion band", ha="right", va="top", fontsize=6.5, color=INK2)
    bx.set_title("(c) ratio", loc="left")
    save(fig, "fig_fall")


def fig_matter():
    """(a) F3.12b: the width R_rms of a bundle against time (means over sprinklings, as stored) for the reference, the
    cold bundle in its field, warm streams without merging and with capacity; (b) F3.13b2 at kappa = 0.3: the core cost
    of the hedgehog, of the hedgehog whose core was randomised, and of the twist (band +-1 SE across sprinklings);
    (c) the mean cost per carrying event of the hedgehog and of the twist (the twist relaxes, the hedgehog does not)."""
    rb, rc = record("F3.12b"), record("F3.13b2")
    st = rb["data"]["stats"]
    tb = np.arange(5.0, 5.0 + 2.0 * len(st["free_0"]["R"]), 2.0)
    fig, (ax, bx, cx) = plt.subplots(1, 3, figsize=(WIDTH, 2.5))
    fig.subplots_adjust(wspace=0.42)
    for key, col, ls, lab in [("free_0", INK2, ":", "no field"), ("field_0", C2, "-", "cold, own field"),
                              ("cl_0.06", C1, "-", r"warm streams $\eta=0.06$"),
                              ("cap4_0.06", C3, "--", r"same with capacity $C=4$")]:
        ax.plot(tb, st[key]["R"], color=col, ls=ls, lw=1.2, label=lab)
        print(f"    F3.12b {key}: R {st[key]['R'][0]:.2f} -> {st[key]['R'][-1]:.2f}; slopes W1 {st[key]['s1'][0]:+.4f} "
              f"W2 {st[key]['s2'][0]:+.4f} ({st[key]['c1']}, {st[key]['c2']})")
    ax.set_xlabel("time")
    ax.set_ylabel(r"width $R_{\rm rms}$ [$\ell$]")
    ax.set_title("(a) bundle in its own field", loc="left")
    ax.set_ylim(top=5.0)
    ax.legend(loc="upper left", fontsize=5.8)
    rows = rc["data"]["rows"]
    t = np.asarray(rows[0]["arms"]["defect_0.3"]["t"], float)
    for key, col, mk, lab in [("defect_0.3", C1, "o", "hedgehog"), ("erase_0.3", C3, "D", "hedgehog, core erased"),
                              ("twist_0.3", C2, "s", "twist (no charge)")]:
        for q, axx in (("core_cost", bx), ("mean_cost", cx)):
            if q == "mean_cost" and key == "erase_0.3":
                continue
            A = np.array([r_["arms"][key][q] for r_ in rows], float)
            m = np.nanmean(A, axis=0)
            se = np.nanstd(A, axis=0, ddof=1) / np.sqrt(np.sum(np.isfinite(A), axis=0))
            axx.fill_between(t, m - se, m + se, color=col, alpha=0.18, lw=0)
            axx.plot(t, m, color=col, marker=mk, ms=3, lw=1.1, label=lab)
            print(f"    F3.13b2 {key} {q}: {m[0]:.4f} -> {m[-1]:.4f} (+- {se[-1]:.4f})")
    bx.axvspan(6.0, 7.0, color=GRID, alpha=0.8, lw=0)
    bx.set_xlabel("time")
    bx.set_ylabel("core cost")
    bx.set_title(r"(b) core, $\kappa=0.3$", loc="left")
    bx.set_ylim(0.0, 0.75)
    bx.legend(loc="upper right", fontsize=5.8)
    cx.set_xlabel("time")
    cx.set_ylabel("mean cost per event")
    cx.set_title("(c) all events", loc="left")
    cx.set_ylim(bottom=0)
    save(fig, "fig_matter")


# ===================================================================================== fig_defect (F3.14b, F3.15, F3.14c, F5.2b)

def _line_value(rec, pattern):
    """A number from the printed summary of a record (rec["lines"]), by a regular expression with groups (value, SE)."""
    for ln in rec["lines"]:
        m = re.search(pattern, ln)
        if m:
            return tuple(float(g) for g in m.groups())
    raise KeyError(f"{rec['_file']}: no line matches {pattern!r}")


def fig_defect():
    """(a) inertia r_M = M_in / M_rest: F3.14b per density (1+1, both seeds) and F3.15 (2+1 wall, both seeds);
    (b) swerves: R = D_p rho / ((32/15) mu^(5/2)) per cell of F3.14c (the M 16 cell is reported, not judged);
    (c) fall under the growth rule (F5.2b): r_B = displacement of the test kink / continuum, against the source mass, for
    a test kink of mass 8 and, at source mass 8, of mass 4."""
    print("fig_defect: matter as a defect of a field of the order")
    fig, (ax, bx, cx) = plt.subplots(1, 3, figsize=(WIDTH, 2.45), gridspec_kw={"width_ratios": [1.15, 1.0, 1.0]})
    fig.subplots_adjust(wspace=0.45)
    # (a) inertia
    for seed, mk, off in ((0, "o", 0.97), (1, "s", 1.03)):
        pd = record("F3.14b", seed=seed)["data"]["per_density"]
        rho = sorted(float(r_) for r_ in pd)
        y = [pd[f"{r_:.1f}" if f"{r_:.1f}" in pd else str(r_)]["rM"] for r_ in rho]
        ax.errorbar([r_ * off for r_ in rho], [v[0] for v in y], [v[1] for v in y], fmt=mk, color=C1, ms=4,
                    mec="white", mew=0.5, label=f"1+1 kink, seed {seed}", **ERR)
        print(f"    F3.14b seed {seed}: " + ", ".join(f"rho {r_:g}: {v[0]:.4f} +- {v[1]:.4f}" for r_, v in zip(rho, y)))
    for seed, mk, off in ((0, "D", 0.95), (1, "^", 1.05)):
        rec = record("F3.15", seed=seed)
        rM, se = _line_value(rec, r"r_M = M_in / M_rest ([0-9.]+) \+- ([0-9.]+)")
        rho = float(rec["params"].get("rho", 30))
        ax.errorbar([rho * off], [rM], [se], fmt=mk, color=C2, ms=4, mec="white", mew=0.5,
                    label=f"2+1 wall, seed {seed}", **ERR)
        print(f"    F3.15 seed {seed}: rho {rho:g}: r_M {rM:.4f} +- {se:.4f}")
    ax.axhline(1.0, color=INK2, lw=0.8, ls=(0, (1, 2)))
    ax.set_xscale("log")
    ax.set_xticks([30, 100, 200, 400, 800])
    ax.set_xticklabels(["30", "100", "200", "400", "800"])
    ax.minorticks_off()
    ax.set_xlabel(r"density $\rho$")
    ax.set_ylabel(r"$r_M=M_{\rm in}/M_{\rm rest}$")
    ax.set_ylim(0.80, 1.08)
    ax.set_title("(a) inertia", loc="left")
    ax.legend(loc="lower right", fontsize=5.4, handletextpad=0.2)
    # (b) swerves
    cells = record("F3.14c")["data"]["cells"]
    keys = sorted(cells, key=lambda c: (float(c.split("_")[1]), float(c.split("_")[0])))
    for i, c in enumerate(keys):
        rho, mu = (float(v) for v in c.split("_"))
        judged = not (mu > 1.0)
        bx.errorbar([i], [cells[c]["R"]], [cells[c]["R_se"]], fmt="o", color=C1 if judged else INK2, ms=4,
                    mfc=C1 if judged else "white", mec=C1 if judged else INK2, mew=0.8, **ERR)
        print(f"    F3.14c {c}: R {cells[c]['R']:.3f} +- {cells[c]['R_se']:.3f}{'' if judged else ' (reported)'}")
    bx.set_xticks(range(len(keys)))
    bx.set_xticklabels([f"$M\\,{8 * float(c.split('_')[1]) ** 0.5:g}$\n$\\rho\\,{float(c.split('_')[0]):g}$" for c in keys],
                       fontsize=6.2)
    bx.axhline(1.0, color=INK2, lw=0.8, ls=(0, (1, 2)))
    bx.set_ylim(0.7, 1.2)
    bx.set_xlim(-0.6, len(keys) - 0.4)
    bx.set_ylabel(r"$D_p\,\rho\,/\,(\frac{32}{15}\mu^{5/2})$")
    bx.set_title("(b) swerves", loc="left")
    bx.text(len(keys) - 1, 0.74, "reported", ha="center", fontsize=5.8, color=INK2)
    # (c) fall under the growth rule
    for seed, mk, off in ((0, "o", -0.25), (1, "s", 0.25)):
        rec = record("F5.2b", seed=seed)
        pts = [(ma, *_line_value(rec, rf"^\s*M_A {ma}, M_B 8: .*?r_B ([0-9.]+) \+- ([0-9.]+)")) for ma in (4, 8, 16)]
        cx.errorbar([p_[0] * (1 + off / 8) for p_ in pts], [p_[1] for p_ in pts], [p_[2] for p_ in pts], fmt=mk + "-",
                    color=C1, ms=4, lw=0.8, mec="white", mew=0.5, label=f"test mass 8, seed {seed}", **ERR)
        r4, s4 = _line_value(rec, r"^\s*M_A 8, M_B 4: .*?r_B ([0-9.]+) \+- ([0-9.]+)")
        cx.errorbar([8 * (1 + (off + 0.5) / 8)], [r4], [s4], fmt="D" if seed == 0 else "^", color=C2, ms=4, mec="white",
                    mew=0.5, label=f"test mass 4, seed {seed}", **ERR)
        print(f"    F5.2b seed {seed}: " + ", ".join(f"M_A {a_}: r_B {b_:.4f} +- {c_:.4f}" for a_, b_, c_ in pts) +
              f"; test mass 4 at M_A 8: {r4:.4f} +- {s4:.4f}")
    cx.axhline(1.0, color=INK2, lw=0.8, ls=(0, (1, 2)))
    cx.set_xscale("log")
    cx.set_xticks([4, 8, 16])
    cx.set_xticklabels(["4", "8", "16"])
    cx.minorticks_off()
    cx.set_xlabel(r"source mass $M_A$")
    cx.set_ylabel(r"fall / continuum, $r_B$")
    cx.set_ylim(0.94, 1.08)
    cx.set_title("(c) fall", loc="left")
    cx.legend(loc="lower left", fontsize=5.4, handletextpad=0.2)
    save(fig, "fig_defect")


# ===================================================================================== fig_fieldeq (F6.8g, F6.10b, OV3-lin-a)

def _ov3_inputs():
    """Stiffness at n = kT = 4 (F6.8g) and source form factors (F6.10b) at the common wavenumbers, from the records."""
    g, b = record("F6.8g"), record("F6.10b")
    ov3 = {float(k): v for k, v in g["data"]["ov3"].items()}
    rat = b["data"]["ratios"]
    F = {s: {float(k): v["ratio"] for k, v in rat[s].items()} for s in ("omega", "nu")}
    dF = {s: {float(k): v["delta"] for k, v in rat[s].items()} for s in ("omega", "nu")}
    return g, ov3, F, dF


def fig_fieldeq():
    """(a) the three static stiffnesses at the fixed age kT = 4 (F6.8g rows and fits); (b) the source form factors of the
    pairwise matter action (F6.10b; bars: the stated Delta of the record); (c) the static linear response computed from
    (a) and (b) by stationarity of the two expected actions (diag/ov3_lin.py): the field of a dust slab against
    linearised GR, gamma(k) = Psi/Phi of dust, and the lapse of a static BPS wall against that of dust of the same
    energy (sign reversed: GR gives -1). At k -> 0 the solver is fed the Einstein--Hilbert values and form factors 1."""
    from diag.ov3_lin import response
    print("fig_fieldeq: the linear field equation of the order")
    g, ov3, F, dF = _ov3_inputs()
    fig, (ax, bx, cx) = plt.subplots(1, 3, figsize=(WIDTH, 2.45))
    fig.subplots_adjust(wspace=0.45)
    D = g["data"]
    kk = np.linspace(0.0, 0.13, 200)
    kk[0] = 1e-9
    for key, rows, col, mk, lab in (("c@4", "optical", C1, "o", r"conformal $\mathcal{A}/6$"),
                                    ("ab@4", "mixed", C2, "s", r"mixed $\mathcal{B}/8$"),
                                    ("bb@4", "optical", C3, "D", r"optical $\mathcal{C}/2$")):
        K0, c1, c0 = D["fits"][key][:3]
        field = {"c@4": "Kc", "ab@4": "Kab", "bb@4": "Kbb"}[key]
        pts = sorted((r_["k"], r_[field]) for r_ in D[rows] if abs(r_.get("n", 4.0) - 4.0) < 1e-9)
        ax.plot(kk, kfit(kk, K0, c1, c0), color=col, lw=1.0)
        ax.plot([p_[0] for p_ in pts], [p_[1] for p_ in pts], mk, color=col, ms=4 if key != "ab@4" else 6,
                mfc=col if key != "ab@4" else "none", mew=0.9, label=f"{lab}, $k\\to0$: {K0:.5f}")
        print(f"    F6.8g {key}: K0 {K0:.5f}; " + ", ".join(f"k {a_:g}: {b_:.5f}" for a_, b_ in pts))
    ax.axhline(1.0, color=INK2, lw=0.8, ls=(0, (1, 2)))
    ax.set_xlim(0, 0.13)
    ax.set_ylim(0.76, 1.02)
    ax.set_xlabel(r"$k\ell_k$")
    ax.set_ylabel("stiffness / Einstein--Hilbert")
    ax.set_title(r"(a) order, $kT=4$", loc="left")
    ax.legend(loc="lower left", fontsize=5.6, handletextpad=0.3)
    for s, col, mk, lab in (("omega", C1, "o", r"$F_\omega$ (trace)"), ("nu", C3, "D", r"$F_\nu$ ($\sum_iT_{ii}$, cones)")):
        ks = sorted(F[s])
        bx.errorbar(ks, [F[s][k] for k in ks], [dF[s][k] for k in ks], fmt=mk + "-", color=col, ms=4, lw=0.9,
                    mec="white", mew=0.5, label=lab, **ERR)
        print(f"    F6.10b {s}: " + ", ".join(f"k {k:g}: {F[s][k]:.4f}" for k in ks))
    bx.axhline(1.0, color=INK2, lw=0.8, ls=(0, (1, 2)))
    bx.set_xlim(-0.005, 0.13)
    bx.set_ylim(0.76, 1.02)
    bx.set_xlabel(r"$k\ell_k$")
    bx.set_ylabel(r"source / $T_{ab}$")
    bx.set_title("(b) matter, wall", loc="left")
    bx.legend(loc="lower left", fontsize=5.6, handletextpad=0.3)
    ks = [0.0] + sorted(ov3)
    out = {"Phi": [], "gamma": [], "wall": []}
    for k in ks:
        if k == 0.0:
            r = response(1.0, 1.0, 1.0, 1.0, 1.0)
        else:
            r = response(ov3[k]["K_c"], ov3[k]["K_ab"], ov3[k]["K_bb"], F["omega"][k], F["nu"][k])
        out["Phi"].append(r["dust"]["Phi/GR"])
        out["gamma"].append(r["dust"]["Psi/Phi"])
        out["wall"].append(-r["wall"]["Phi"] / r["dust"]["Phi"])
    for key, col, mk, lab in (("Phi", C1, "o", r"dust $\Phi/\Phi_{\rm GR}$"), ("gamma", C3, "D", r"dust $\gamma(k)=\Psi/\Phi$"),
                              ("wall", C2, "s", r"$-$wall/dust lapse (repels)")):
        cx.plot(ks, out[key], mk + "-", color=col, ms=4, lw=0.9, mec="white", mew=0.5, label=lab)
        print(f"    response {key}: " + ", ".join(f"k {k:g}: {v:.4f}" for k, v in zip(ks, out[key])))
    cx.axhline(1.0, color=INK2, lw=0.8, ls=(0, (1, 2)))
    cx.set_xlim(-0.005, 0.13)
    cx.set_ylim(0.76, 1.16)
    cx.set_xlabel(r"$k\ell_k$")
    cx.set_ylabel("response / linearised GR")
    cx.set_title("(c) response", loc="left")
    cx.legend(loc="lower left", fontsize=5.6, handletextpad=0.3)
    save(fig, "fig_fieldeq")


# ===================================================================================== fig_chain (diagram)

def fig_chain():
    """The chain of the field equation of the order as a diagram: two actions of one order, their stationarity, the linear
    field equation and its consequences, with the tests that carry each step and their status. The status is written in
    each box (passed / development / derived / open) and also shown by the edge (solid: frozen test or derived from
    frozen records; dashed: development; dotted: open), so it is not carried by colour alone. Numbers are read from the
    records."""
    from matplotlib.patches import FancyBboxPatch
    from diag.ov3_lin import response
    print("fig_chain: the chain of the field equation")
    g, ov3, F, _ = _ov3_inputs()
    K0 = {k: g["data"]["fits"][k][0] for k in ("c@4", "ab@4", "bb@4")}
    r12 = response(ov3[0.12]["K_c"], ov3[0.12]["K_ab"], ov3[0.12]["K_bb"], F["omega"][0.12], F["nu"][0.12])
    v9 = record("F6.9a", seed=0)["data"]["eps"]
    worst = max(abs(v["c4eps"]) + 2 * v["se"] for v in v9.values())
    a10 = record("F6.10a")["data"]["rows"]
    rkin = max(abs(v["R"]) for v in a10.values())
    fig = plt.figure(figsize=(WIDTH, 3.9))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 64)
    ax.axis("off")
    ax.grid(False)
    style = {"frozen": dict(ls="-", lw=1.0), "dev": dict(ls=(0, (4, 2)), lw=0.9), "open": dict(ls=(0, (1, 1.5)), lw=0.9)}
    fill = {"frozen": "#eef4fc", "dev": "#fdf1ea", "open": "#f4f4f2", "plain": "white"}

    def box(x, y, w, h, title, body, kind, status=None):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3,rounding_size=1.2",
                                    fc=fill[kind], ec=INK2, **style.get(kind, style["frozen"])))
        ax.text(x + 1.2, y + h - 1.4, title, ha="left", va="top", fontsize=6.9, weight="bold", color=INK)
        ax.text(x + 1.2, y + h - 4.6, body, ha="left", va="top", fontsize=6.0, color=INK, linespacing=1.25)
        if status:
            ax.text(x + w - 1.0, y + 0.9, status, ha="right", va="bottom", fontsize=5.6, style="italic", color=INK2)

    def arrow(x0, y0, x1, y1):
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="-|>", color=INK2, lw=0.9, shrinkA=0, shrinkB=0, mutation_scale=8))

    box(30, 56, 40, 7, "a causal order of events", "Poisson sprinkling: order $+$ number $=$ geometry", "plain")
    box(1, 38, 45, 14, r"the action of the order $\langle S_{\rm BD}\rangle$",
        "static stiffness, three sectors, at $kT=4$:\n"
        f"$\\mathcal{{A}}/6\\to{K0['c@4']:.5f}$, $\\mathcal{{B}}/8\\to{K0['ab@4']:.5f}$, $\\mathcal{{C}}/2\\to{K0['bb@4']:.5f}$\n"
        f"flat order: no volume term ($<{worst * 1e4:.0f}\\times10^{{-4}}$ per element, 2D)",
        "frozen", "passed: F6.8a, b2, c2, g; F6.9a")
    box(54, 38, 45, 14, r"matter: a field on all events, $S_m^\chi$",
        "a kink carried by the Green function of the order;\n"
        "pairwise action with the BD kernel $\\to T_{ab}\\times F_\\omega, F_\\nu\\to1$\n"
        f"2D trace: kinetic share $\\leq{rkin:.3f}$, falls with density",
        "frozen", "passed: F3.14b, F3.15; F6.10a, b")
    box(1, 24, 35, 11, "selection of the order",
        "the action alone does not hold\nthe flat order; saddle of $\\langle S\\rangle$ (linear),\naction + entropy beyond",
        "dev", "development (OV1)")
    box(64, 24, 35, 11, "conservation",
        "energy of $S_m^\\chi$ as a Noether charge,\nkept by the order's dynamics\nto $0.3$--$0.5\\,\\%$ (2D)",
        "dev", "development (T3)")
    box(12, 10, 76, 11, r"stationarity $\delta(\langle S_{\rm BD}\rangle+\langle S_m\rangle)/\delta g_{ab}=0$, first order",
        "$G_{ab}+2(\\eta_{ab}\\partial^2-\\partial_a\\partial_b)F_0(\\partial^2)R=8\\pi G\\,T^{\\rm eff}_{ab}$\n"
        f"$k\\ell_k=0.12$: dust $\\Phi/\\Phi_{{\\rm GR}}={r12['dust']['Phi/GR']:.3f}$, $\\gamma={r12['dust']['Psi/Phi']:.3f}$, "
        f"wall/dust $={r12['wall']['Phi'] / r12['dust']['Phi']:+.3f}$ (repels)",
        "frozen", "derived from frozen records (OV3-lin-a)")
    box(1, 0.4, 98, 8.2, "open",
        r"second order ($\beta$); old fields ($\varepsilon=b(kT)^2$, F6.8f); entropic correction to the stiffness;" "\n"
        r"time-dependent and tensor sectors; the value of $G$; a defect in 3+1", "open")
    for x0, x1 in ((42, 23.5), (58, 76.5)):
        arrow(x0, 55.6, x1, 52.6)
    arrow(41, 37.6, 41, 21.6)
    arrow(59, 37.6, 59, 21.6)
    print(f"    K0 {K0}; F6.9a worst |c/4eps|+2SE {worst:.2e}; F6.10a max |R| {rkin:.4f}; "
          f"k 0.12: Phi/GR {r12['dust']['Phi/GR']:.4f}, gamma {r12['dust']['Psi/Phi']:.4f}, "
          f"wall/dust {r12['wall']['Phi'] / r12['dust']['Phi']:+.4f}")
    save(fig, "fig_chain")


# ===================================================================================== main

FIGURES = [fig_blind, fig_response, fig_counter, fig_conformal, fig_sectors, fig_optical, fig_window, fig_fall, fig_matter,
           fig_defect, fig_chain, fig_fieldeq]

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
