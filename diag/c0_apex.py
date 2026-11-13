"""
C0, development: the continuum expectation of a layer kernel at finite density, and the apex term (13 Nov; decides
nothing; card C0 NOT frozen)
====================================================================================================================
A kernel w_n (n = the number of events in the interval of a pair u < v) has, on a Poisson sprinkling of density rho, the
exact expectation
    E[c sum_{u < v} w(n) J(u)] = c rho int d^4u K(sigma) J(u),   K(sigma) = sum_n w_n Pois(n; rho pi sigma^2 / 24),
sigma = tau^2 the squared proper time of the pair. With the spatial distance r of the pair, t_v - t_u = sqrt(r^2 + sigma)
and d^4u = 2 pi r^2 dr dsigma / sqrt(r^2 + sigma) exactly; the near-cone form used to derive the moments (and Johnston's
c) is 2 pi r dr dsigma. Hence, per unit measure r dr of the cone,
    1 + eps(r) = [int dsigma K(sigma) r / sqrt(r^2 + sigma)] / [int dsigma K(sigma)]                      (A1)
with no free parameter. The deficit at the apex (r of a few rho^-1/4) is this factor, not discreteness noise.

The integrated apex term. For a source smooth on the scale of the apex, the kernel acts as G*J + delta J(v) with
    delta = int_0^inf eps(r) r dr = <sigma ln sigma>/4 = (a / (8 sqrt pi)) sum_n w_n H_n,   a = sqrt(24/(pi rho)),   (A2)
provided sum_n w_n = 0 (else the integral diverges logarithmically: the slow 1/r^2 tail of the links); H_n is the
harmonic number (the moment <s^1/2 ln s> of layer n is psi(n + 1)). delta scales as rho^-1/2. In the linear rule
(J = mu phi) it shifts the mass, mu_eff = mu/(1 + delta mu), and the initial amplitude by 1/(1 + delta mu).

Two ways to remove it, both from the continuum, no measured number:
  MVA (constraint)    the minimum-variance kernel with the third constraint sum_n w_n H_n = 0;
  ct (counterterm)    the rule phi(v) = (phi0(v) - c sum w J(u)) / (1 - delta mu) with delta from (A2) (linear J;
                      for J = mu sin phi the same local term solved per event).

Prediction of the dynamics: the field of the linear rule, spatially homogeneous (k = 0) or ~ cos(k x), obeys the
Volterra equation f(t) = f0(t) - mu int_0^t F_k(D) f(t - D) dD with
    F_k(D) = c rho 2 pi int_0^{min(D^2, sigma_c)} dsigma K(sigma) sqrt(D^2 - sigma) sinc(k sqrt(D^2 - sigma)),     (A3)
(F_k(D) -> D for the continuum Green function at k = 0), solved here on a grid. sigma_c is the enumeration cut
(rho V <= lam). (A1) to (A3) are the expectation over sprinklings; the mean of a linear functional of the field differs
from it only by the feedback of the noise (second order).
"""
from __future__ import annotations

import math

import numpy as np
from scipy import integrate, stats


def harmonic(n):
    return sum(1.0 / k for k in range(1, n + 1))


def delta(w, rho):
    """(A2); nan if sum w != 0."""
    w = np.asarray(w, float)
    if abs(w.sum()) > 1e-9:
        return float("nan")
    a = math.sqrt(24.0 / (math.pi * rho))
    return a / (8.0 * math.sqrt(math.pi)) * sum(wn * harmonic(n) for n, wn in enumerate(w))


def _K(sig, w, rho):
    s = rho * math.pi * np.asarray(sig) ** 2 / 24.0
    return sum(wn * stats.poisson.pmf(n, s) for n, wn in enumerate(w) if wn != 0.0)


def eps_pred(w, rho, rbin=0.25, nrb=10, lam=22.0):
    """(A1) averaged over the shells [i rbin, (i + 1) rbin) with the measure r dr, as the laboratory bins it."""
    a = math.sqrt(24.0 / (math.pi * rho))
    sigc = a * math.sqrt(lam)
    c = 1.0 / (2.0 * math.pi * math.sqrt(6.0 * rho))
    out = []
    for i in range(nrb):
        r0, r1 = i * rbin, (i + 1) * rbin
        f = lambda sig, r: c * rho * 2 * math.pi * float(_K(sig, w, rho)) * r * r / math.sqrt(r * r + sig)
        val, _ = integrate.dblquad(f, r0, r1, 0.0, sigc, epsabs=1e-8)
        out.append(val / (0.5 * (r1 * r1 - r0 * r0)) - 1.0)
    return out


def F_k(w, rho, k, tg, lam=22.0, nsig=1501):
    """(A3) on the grid tg."""
    a = math.sqrt(24.0 / (math.pi * rho))
    sg = np.linspace(0.0, a * math.sqrt(lam), nsig)
    Kw = _K(sg, w, rho)
    c = 1.0 / (2.0 * math.pi * math.sqrt(6.0 * rho))
    F = np.zeros_like(tg)
    for i, D in enumerate(tg):
        m = sg < D * D
        if m.sum() < 2:
            continue
        r = np.sqrt(D * D - sg[m])
        fac = np.sinc(k * r / math.pi) if k else 1.0
        F[i] = c * rho * 2 * math.pi * np.trapezoid(Kw[m] * r * fac, sg[m])
    return F


def volterra(F, f0, dt, scale=1.0):
    """f = scale (f0 - int F f), trapezoid (F(0) = 0)."""
    f = np.zeros_like(f0)
    for i in range(f0.size):
        s = 0.5 * F[i] * f[0] if i > 0 else 0.0
        if i > 1:
            s += np.dot(F[1:i], f[i - 1:0:-1])
        f[i] = scale * (f0[i] - dt * s)
    return f


def dyn_pred(w, rho, k, tb, X=7.0, ct=0.0, mu=1.0, lam=22.0, dt=0.005):
    """Predicted deviation (f - exact) of the amplitude of mode k (0 or 1, data cos(k x) cos(k t)) averaged over the
    slabs |t - b| < 0.25 with the inner volume (X - t - 0.3)^3 as weight; ct = delta of the counterterm (0 = none)."""
    T = max(tb) + 0.3
    tg = np.arange(0.0, T + dt, dt)
    F = F_k(w, rho, k, tg, lam)
    f = volterra(mu * F, np.cos(k * tg), dt, 1.0 / (1.0 - ct * mu))
    exact = np.cos(math.sqrt(k * k + mu) * tg)
    dev = f - exact
    out = []
    for b in tb:
        m = np.abs(tg - b) < 0.25
        wt = np.clip(X - tg[m] - 0.3, 0, None) ** 3
        out.append(float(np.sum(dev[m] * wt) / np.sum(wt)))
    return out


def norm_pred(w, rho, tb, X=7.0, lam=22.0, dt=0.005):
    """Predicted c sum w / (t^2/2) for the constant source, slab-averaged as in c0_minvar."""
    T = max(tb) + 0.3
    tg = np.arange(0.0, T + dt, dt)
    F = F_k(w, rho, 0, tg, lam)
    S = np.concatenate([[0.0], np.cumsum(0.5 * (F[1:] + F[:-1]) * dt)])
    out = []
    for b in tb:
        m = np.abs(tg - b) < 0.25
        wt = np.clip(X - tg[m] - 0.3, 0, None) ** 3
        out.append(float(np.sum(S[m] * wt) / np.sum(0.5 * tg[m] ** 2 * wt)))
    return out


def dyn_pred_nl0(w, rho, A, tb, X=7.0, ct=0.0, mu=1.0, lam=22.0, dt=0.005):
    """Nonlinear (J = mu sin phi), spatially homogeneous: f - ct mu sin f = A - int F_0(D) mu sin f(t - D) dD with F_0
    of (A3), against the pendulum f'' = -mu sin f, f(0) = A, f'(0) = 0. Returns (f - exact)/A slab-averaged as dyn_pred.
    This is the expectation of the kernel without the feedback of the noise; for a nonlinear J the noise also
    renormalizes the source (E sin(phi + eta) = exp(-var/2) sin phi for Gaussian eta), which is measured, not predicted."""
    T = max(tb) + 0.3
    tg = np.arange(0.0, T + dt, dt)
    F = F_k(w, rho, 0, tg, lam)
    d = ct * mu
    f = np.zeros_like(tg)
    J = np.zeros_like(tg)
    for i in range(tg.size):
        s = 0.5 * F[i] * J[0] if i > 0 else 0.0
        if i > 1:
            s += np.dot(F[1:i], J[i - 1:0:-1])
        rhs = A - dt * s
        p = rhs
        for _ in range(50):
            dp = (p - d * math.sin(p) - rhs) / (1.0 - d * math.cos(p))
            p -= dp
            if abs(dp) < 1e-14:
                break
        f[i] = p
        J[i] = mu * math.sin(p)
    ex = pendulum(A, tg, mu)
    dev = (f - ex) / A
    out = []
    for b in tb:
        m = np.abs(tg - b) < 0.25
        wt = np.clip(X - tg[m] - 0.3, 0, None) ** 3
        out.append(float(np.sum(dev[m] * wt) / np.sum(wt)))
    return out


def pendulum(A, tg, mu=1.0):
    from scipy.integrate import solve_ivp
    sol = solve_ivp(lambda t, y: [y[1], -mu * math.sin(y[0])], (0.0, float(tg[-1])), [A, 0.0], t_eval=tg,
                    rtol=1e-11, atol=1e-12, method="DOP853")
    return sol.y[0]
