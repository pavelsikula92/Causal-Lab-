"""
C0, development: the nonlinear modes of diag/c0_ct.py against the mean field with the source renormalized by the noise
(14 Nov; decides nothing)
=====================================================================================================================
For J = mu sin phi and a field phi = f + eta with event-wise noise eta of sd s (Gaussian), the mean source is
E[mu sin(f + eta)] = mu exp(-s^2/2) sin f: the noise softens the potential (the factor of the self-consistent harmonic
approximation of sine-Gordon). The mean field then obeys
    f - delta mu q sin f = A - int F_0(D) mu q(t - D) sin f(t - D) dD,   q(t) = exp(-s(t)^2/2),
with F_0 the kernel's expectation ((A3) of diag/c0_apex.py) and s(t) = the MEASURED sd of the field of the same row
(input, not predicted; s = 0 at t = 0, linear interpolation between the slabs). No other parameter. Compared with the
measured mean deviation from the pendulum.

    python -m diag.c0_ct_renorm diag/out/c0_ct_<stamp>.json
"""
from __future__ import annotations

import json
import math
import sys

import numpy as np

from diag.c0_apex import F_k, pendulum
from diag.c0_ct import MODES, rows


def mean_field(w, ct, rho, A, tg, s, dt, F=None):
    F = F_k(w, rho, 0, tg) if F is None else F
    q = np.exp(-0.5 * s ** 2)
    f = np.zeros_like(tg)
    J = np.zeros_like(tg)
    for i in range(tg.size):
        acc = 0.5 * F[i] * J[0] if i > 0 else 0.0
        if i > 1:
            acc += np.dot(F[1:i], J[i - 1:0:-1])
        rhs = A - dt * acc
        p = rhs
        for _ in range(50):
            dp = (p - ct * q[i] * math.sin(p) - rhs) / (1 - ct * q[i] * math.cos(p))
            p -= dp
            if abs(dp) < 1e-14:
                break
        f[i] = p
        J[i] = q[i] * math.sin(p)
    return f


def main(fn):
    res = json.load(open(fn))
    for rho in sorted({r["rho"] for r in res}):
        R = [r for r in res if r["rho"] == rho]
        tb = np.array(R[0]["t"])
        X = R[0]["X"]
        dt = 0.005
        tg = np.arange(0.0, tb[-1] + 0.3, dt)
        print(f"######## rho {rho:g} (t = " + " ".join(f"{v:.2f}" for v in tb[1::2]) + ")")
        for name, w, ct in rows(rho):
            F = F_k(w, rho, 0, tg)
            print(f"== {name}")
            for m in MODES:
                if m[1] or m[2] != 0:
                    continue
                A = m[3]
                D = np.array([r["rows"][name]["modes"][m[0]] for r in R])
                amp, sd = D[:, :, 0].mean(0), D[:, :, 1].mean(0)
                s = np.interp(tg, np.r_[0.0, tb], np.r_[0.0, sd * A])
                ex = pendulum(A, tg)
                lines = {}
                for lab, ss in (("no noise", 0 * s), ("renorm", s)):
                    dev = (mean_field(w, ct, rho, A, tg, ss, dt, F) - ex) / A
                    v = []
                    for b in tb:
                        mm = np.abs(tg - b) < 0.25
                        wt = np.clip(X - tg[mm] - 0.3, 0, None) ** 3
                        v.append(np.sum(dev[mm] * wt) / np.sum(wt))
                    lines[lab] = np.array(v)
                print(f"   {m[0]:9s} meas     " + " ".join(f"{v:+.2f}" for v in amp[1::2]))
                for lab, v in lines.items():
                    print(f"   {'':9s} {lab:8s} " + " ".join(f"{x:+.2f}" for x in v[1::2])
                          + f"   max |meas - this| {np.max(np.abs(amp - v)):.3f}")


if __name__ == "__main__":
    main(sys.argv[1])
