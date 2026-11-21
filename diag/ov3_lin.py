"""
OV3-lin-a (development; decisions 41 and 44): the static linear response of the order to a source, per k, in the
young-field reading (n = kT = 4) and for the bare ensemble action.
==================================================================================================================
Static metric of family M: g = e^{2 omega}(-dt^2 + e^{2 nu} dx^2), omega = Phi (the lapse), Psi = -(omega + nu).
Stationarity of <S_BD> + <S_m> per k (OV3 sec. 1; signs verified with sympy there):
    alpha k^2 [ 6 K  omega + 4 K_ab nu ] = F_omega tau,      alpha k^2 [ 4 K_ab omega + 2 K_bb nu ] = F_nu sigma,
tau = T^a_a and sigma = Sum_i T_ii of the source mode; K, K_ab, K_bb the stiffness factors (F6.8a/b2/c-type rows at
n = 4: F6.8g), F_omega, F_nu the source form factors (F6.10b). alpha = -1/(8 pi G) is a common factor: every output is
a ratio to linearised GR (all factors 1), so alpha drops out.
Sources (per unit amplitude of the mode): dust tau = -rho, sigma = 0; static BPS wall tau = -3 rho_E,
sigma = -2 rho_E (OV3 sec. 2).
Outputs per k: Phi/Phi_GR and Psi/Psi_GR for dust; for the wall the sign of Phi relative to dust of the same energy,
Phi_wall/Phi_GR,wall, Psi/Phi (GR: +1 dust, -1 wall); the lapse source of the wall in units of rho_E.

    python -m diag.ov3_lin                   (built-in inputs: F6.10b and the F6.8g quick rows; development)
"""
import numpy as np


def solve(K, Kab, Kbb, Fw, Fn, tau, sig):
    A = np.array([[6.0 * K, 4.0 * Kab], [4.0 * Kab, 2.0 * Kbb]])
    b = np.array([Fw * tau, Fn * sig])
    om, nu = np.linalg.solve(A, b)           # in units of 1/(alpha k^2)
    return om, nu, -(om + nu)


def response(K, Kab, Kbb, Fw, Fn):
    out = {}
    for name, tau, sig in (("dust", -1.0, 0.0), ("wall", -3.0, -2.0)):
        om, nu, psi = solve(K, Kab, Kbb, Fw, Fn, tau, sig)
        om0, nu0, psi0 = solve(1, 1, 1, 1, 1, tau, sig)
        out[name] = {"Phi": om, "Psi": psi, "Phi/GR": om / om0, "Psi/GR": psi / psi0,
                     "Psi/Phi": psi / om, "Psi/Phi GR": psi0 / om0}
    return out


# development inputs (21 Nov): F6.10b full run (sha 14801b85), F6.8g quick (rev 1, not frozen)
F610B = {0.03: (0.9847, 0.9898), 0.06: (0.9545, 0.9495), 0.12: (0.8547, 0.8179)}
F68G_QUICK = {0.06: (0.95245, 0.95241, 0.93697), 0.12: (0.84918, 0.84902, 0.80017)}


def main():
    r = response(1, 1, 1, 1, 1)
    print("local limit (all factors 1): dust Phi/GR %.3f Psi/Phi %+.3f; wall Phi/GR %.3f Psi/Phi %+.3f; "
          "wall Phi / dust Phi = %+.3f" % (r["dust"]["Phi/GR"], r["dust"]["Psi/Phi"], r["wall"]["Phi/GR"],
                                            r["wall"]["Psi/Phi"], r["wall"]["Phi"] / r["dust"]["Phi"]))
    print("k     | dust Phi/GR  Psi/GR  Psi/Phi | wall Phi/GR  Psi/GR  Psi/Phi | wall/dust lapse (GR -1)")
    for k, (K, Kab, Kbb) in sorted(F68G_QUICK.items()):
        Fw, Fn = F610B[k]
        r = response(K, Kab, Kbb, Fw, Fn)
        d, w = r["dust"], r["wall"]
        print(f"{k:<5g} | {d['Phi/GR']:.4f}      {d['Psi/GR']:.4f}  {d['Psi/Phi']:+.4f} | {w['Phi/GR']:.4f}      "
              f"{w['Psi/GR']:.4f}  {w['Psi/Phi']:+.4f} | {w['Phi'] / d['Phi']:+.4f}")


if __name__ == "__main__":
    main()
