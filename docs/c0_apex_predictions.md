# C0: predictions written before the run (13 Nov 2026, 20:05; development, decide nothing; card C0 NOT frozen)

Source: `diag/c0_apex.py` (md5 2f21957fcc9ea781fee486acaa27893f), formulas (A1)–(A3): the exact continuum expectation of a
layer kernel at density ρ (Poisson classes in the interval, the exact 4-volume element 2π r² dr dσ / √(r² + σ), the
enumeration cut ρV ≤ 22). There is no free parameter and no number taken from a measurement. Printed by
`python -m diag.c0_minvar --predict --rhos 30,60` (`diag/c0_minvar.py` md5 5b31241ec688e2eceec4bd326a11f63f).

## What v1 (ρ 30) already showed against these formulas (computed after that run)

| | measured (seeds 0–3) | predicted |
|---|---|---|
| ε(r), all 7 kernels, r 0.125–2.375 | e.g. MV1 −0.503 −0.206 −0.082 …; MV16 −0.649 −0.361 −0.194 … | MV1 −0.503 −0.201 −0.079 …; MV16 −0.648 −0.358 −0.188 … (all within ±0.006) |
| wave0 mean deviation, all 7 kernels | e.g. MV16 +0.17 −0.07 −0.47 −0.20 +0.63 | MV16 +0.18 −0.08 −0.48 −0.19 +0.63 (all within ±0.02) |
| normalisation t ≥ 3 | 0.953 (K1), 0.990–0.997 (MV), 0.904 (S16) | 0.947, 0.983–0.991, 0.899: the measurement is 0.5 % higher for every kernel (open) |

wave1 in v1 was the mean over x, diluted by ⟨cos x⟩ of the slab; with that factor it agrees too. v2 reads the projection.

The integrated apex term δ = (a/8√π) Σ w_n H_n (a = √(24/πρ), requires Σ w_n = 0) acts in the linear rule as the mass
μ_eff = μ/(1 + δμ); it scales as ρ^−1/2.

## Predictions for the v2 run (ρ 30 and ρ 60), amp-dev = deviation of the mode amplitude from the exact solution, in A

```
== K1 rho 30: sum g w^2 1.77; delta (A2) n/a; normalisation (t >= 3) 0.9471
   eps(r) -0.637 -0.363 -0.212 -0.134 -0.091 -0.065 -0.048 -0.037 -0.029 -0.024   (r = 0.125 .. 2.375)
   wave0 amp-dev +0.156 +0.052 -0.381 -0.479 +0.155  (t = 0.75 1.75 2.75 3.75 4.75)
   wave1 amp-dev +0.136 -0.085 -0.302 +0.178 +0.511  (t = 0.75 1.75 2.75 3.75 4.75)
== MV1 rho 30: sum g w^2 10.63; delta (A2) -0.0712; normalisation (t >= 3) 0.9912
   eps(r) -0.503 -0.201 -0.079 -0.035 -0.017 -0.009 -0.005 -0.003 -0.002 -0.001   (r = 0.125 .. 2.375)
   wave0 amp-dev +0.078 -0.066 -0.138 -0.011 +0.194  (t = 0.75 1.75 2.75 3.75 4.75)
   wave1 amp-dev +0.065 -0.092 -0.036 +0.136 +0.052  (t = 0.75 1.75 2.75 3.75 4.75)
== MV2 rho 30: sum g w^2 6.18; delta (A2) -0.0840; normalisation (t >= 3) 0.9897
   eps(r) -0.532 -0.228 -0.095 -0.043 -0.021 -0.011 -0.007 -0.004 -0.003 -0.002   (r = 0.125 .. 2.375)
   wave0 amp-dev +0.094 -0.075 -0.168 -0.018 +0.235  (t = 0.75 1.75 2.75 3.75 4.75)
   wave1 amp-dev +0.079 -0.108 -0.049 +0.164 +0.068  (t = 0.75 1.75 2.75 3.75 4.75)
== MV4 rho 30: sum g w^2 3.81; delta (A2) -0.1044; normalisation (t >= 3) 0.9871
   eps(r) -0.568 -0.265 -0.119 -0.057 -0.029 -0.016 -0.009 -0.006 -0.004 -0.003   (r = 0.125 .. 2.375)
   wave0 amp-dev +0.118 -0.085 -0.221 -0.035 +0.302  (t = 0.75 1.75 2.75 3.75 4.75)
   wave1 amp-dev +0.100 -0.132 -0.074 +0.210 +0.099  (t = 0.75 1.75 2.75 3.75 4.75)
== MV8 rho 30: sum g w^2 2.48; delta (A2) -0.1355; normalisation (t >= 3) 0.9833
   eps(r) -0.608 -0.310 -0.151 -0.077 -0.041 -0.024 -0.014 -0.009 -0.006 -0.004   (r = 0.125 .. 2.375)
   wave0 amp-dev +0.146 -0.089 -0.310 -0.075 +0.411  (t = 0.75 1.75 2.75 3.75 4.75)
   wave1 amp-dev +0.126 -0.160 -0.123 +0.280 +0.163  (t = 0.75 1.75 2.75 3.75 4.75)
== MV16 rho 30: sum g w^2 1.67; delta (A2) -0.1815; normalisation (t >= 3) 0.9869
   eps(r) -0.648 -0.358 -0.188 -0.099 -0.052 -0.027 -0.012 -0.004 +0.002 +0.005   (r = 0.125 .. 2.375)
   wave0 amp-dev +0.173 -0.073 -0.469 -0.196 +0.606  (t = 0.75 1.75 2.75 3.75 4.75)
   wave1 amp-dev +0.151 -0.183 -0.229 +0.389 +0.325  (t = 0.75 1.75 2.75 3.75 4.75)
== MVA8 rho 30: sum g w^2 20.86; delta (A2) +0.0000; normalisation (t >= 3) 1.0000
   eps(r) -0.421 -0.104 -0.007 +0.012 +0.012 +0.009 +0.006 +0.004 +0.003 +0.002   (r = 0.125 .. 2.375)
   wave0 amp-dev +0.019 -0.039 +0.019 +0.058 +0.027  (t = 0.75 1.75 2.75 3.75 4.75)
   wave1 amp-dev +0.010 -0.017 +0.043 +0.011 -0.052  (t = 0.75 1.75 2.75 3.75 4.75)
== MVA16 rho 30: sum g w^2 10.41; delta (A2) -0.0000; normalisation (t >= 3) 0.9760
   eps(r) -0.470 -0.150 -0.036 -0.009 -0.008 -0.012 -0.017 -0.021 -0.024 -0.026   (r = 0.125 .. 2.375)
   wave0 amp-dev +0.052 -0.080 +0.072 +0.121 +0.001  (t = 0.75 1.75 2.75 3.75 4.75)
   wave1 amp-dev +0.040 -0.049 +0.088 -0.010 -0.118  (t = 0.75 1.75 2.75 3.75 4.75)
== S16 rho 30: sum g w^2 2.03; delta (A2) n/a; normalisation (t >= 3) 0.8987
   eps(r) -0.645 -0.368 -0.218 -0.150 -0.120 -0.107 -0.102 -0.101 -0.101 -0.101   (r = 0.125 .. 2.375)
   wave0 amp-dev +0.162 -0.061 -0.202 +0.012 +0.218  (t = 0.75 1.75 2.75 3.75 4.75)
   wave1 amp-dev +0.141 -0.154 -0.100 +0.191 +0.048  (t = 0.75 1.75 2.75 3.75 4.75)
== MV2+ct rho 30: sum g w^2 6.18; delta (A2) -0.0840; counterterm -0.0840; normalisation (t >= 3) 0.9897
   eps(r) -0.532 -0.228 -0.095 -0.043 -0.021 -0.011 -0.007 -0.004 -0.003 -0.002   (r = 0.125 .. 2.375)
   wave0 amp-dev +0.043 +0.020 -0.028 -0.046 -0.006  (t = 0.75 1.75 2.75 3.75 4.75)
   wave1 amp-dev +0.047 -0.000 -0.037 +0.005 +0.045  (t = 0.75 1.75 2.75 3.75 4.75)
== MV8+ct rho 30: sum g w^2 2.48; delta (A2) -0.1355; counterterm -0.1355; normalisation (t >= 3) 0.9833
   eps(r) -0.608 -0.310 -0.151 -0.077 -0.041 -0.024 -0.014 -0.009 -0.006 -0.004   (r = 0.125 .. 2.375)
   wave0 amp-dev +0.054 +0.062 -0.051 -0.111 -0.035  (t = 0.75 1.75 2.75 3.75 4.75)
   wave1 amp-dev +0.065 +0.020 -0.085 -0.006 +0.105  (t = 0.75 1.75 2.75 3.75 4.75)
== MV16+ct rho 30: sum g w^2 1.67; delta (A2) -0.1815; counterterm -0.1815; normalisation (t >= 3) 0.9869
   eps(r) -0.648 -0.358 -0.188 -0.099 -0.052 -0.027 -0.012 -0.004 +0.002 +0.005   (r = 0.125 .. 2.375)
   wave0 amp-dev +0.047 +0.116 -0.074 -0.208 -0.080  (t = 0.75 1.75 2.75 3.75 4.75)
   wave1 amp-dev +0.064 +0.056 -0.141 -0.027 +0.190  (t = 0.75 1.75 2.75 3.75 4.75)
== K1 rho 60: sum g w^2 1.77; delta (A2) n/a; normalisation (t >= 3) 0.9605
   eps(r) -0.589 -0.308 -0.169 -0.103 -0.068 -0.048 -0.035 -0.027 -0.021 -0.017   (r = 0.125 .. 2.375)
   wave0 amp-dev +0.129 +0.015 -0.302 -0.318 +0.175  (t = 0.75 1.75 2.75 3.75 4.75)
   wave1 amp-dev +0.111 -0.086 -0.215 +0.167 +0.354  (t = 0.75 1.75 2.75 3.75 4.75)
== MV1 rho 60: sum g w^2 10.63; delta (A2) -0.0503; normalisation (t >= 3) 0.9938
   eps(r) -0.446 -0.153 -0.053 -0.021 -0.010 -0.005 -0.003 -0.002 -0.001 -0.001   (r = 0.125 .. 2.375)
   wave0 amp-dev +0.051 -0.049 -0.091 -0.002 +0.132  (t = 0.75 1.75 2.75 3.75 4.75)
   wave1 amp-dev +0.040 -0.064 -0.019 +0.093 +0.030  (t = 0.75 1.75 2.75 3.75 4.75)
== MV2 rho 60: sum g w^2 6.18; delta (A2) -0.0594; normalisation (t >= 3) 0.9927
   eps(r) -0.476 -0.176 -0.065 -0.027 -0.013 -0.006 -0.004 -0.002 -0.001 -0.001   (r = 0.125 .. 2.375)
   wave0 amp-dev +0.063 -0.057 -0.110 -0.005 +0.159  (t = 0.75 1.75 2.75 3.75 4.75)
   wave1 amp-dev +0.051 -0.076 -0.026 +0.112 +0.038  (t = 0.75 1.75 2.75 3.75 4.75)
== MV4 rho 60: sum g w^2 3.81; delta (A2) -0.0738; normalisation (t >= 3) 0.9909
   eps(r) -0.513 -0.209 -0.083 -0.036 -0.017 -0.009 -0.005 -0.003 -0.002 -0.001   (r = 0.125 .. 2.375)
   wave0 amp-dev +0.082 -0.069 -0.143 -0.011 +0.203  (t = 0.75 1.75 2.75 3.75 4.75)
   wave1 amp-dev +0.068 -0.096 -0.038 +0.142 +0.054  (t = 0.75 1.75 2.75 3.75 4.75)
== MV8 rho 60: sum g w^2 2.48; delta (A2) -0.0958; normalisation (t >= 3) 0.9882
   eps(r) -0.556 -0.251 -0.109 -0.051 -0.026 -0.014 -0.008 -0.005 -0.003 -0.002   (r = 0.125 .. 2.375)
   wave0 amp-dev +0.109 -0.082 -0.199 -0.026 +0.274  (t = 0.75 1.75 2.75 3.75 4.75)
   wave1 amp-dev +0.093 -0.123 -0.062 +0.191 +0.084  (t = 0.75 1.75 2.75 3.75 4.75)
== MV16 rho 60: sum g w^2 1.67; delta (A2) -0.1283; normalisation (t >= 3) 0.9946
   eps(r) -0.598 -0.297 -0.138 -0.065 -0.030 -0.012 -0.002 +0.003 +0.006 +0.008   (r = 0.125 .. 2.375)
   wave0 amp-dev +0.140 -0.094 -0.311 -0.073 +0.425  (t = 0.75 1.75 2.75 3.75 4.75)
   wave1 amp-dev +0.121 -0.158 -0.117 +0.284 +0.163  (t = 0.75 1.75 2.75 3.75 4.75)
== MVA8 rho 60: sum g w^2 20.86; delta (A2) +0.0000; normalisation (t >= 3) 0.9999
   eps(r) -0.358 -0.062 +0.006 +0.013 +0.009 +0.006 +0.004 +0.003 +0.002 +0.001   (r = 0.125 .. 2.375)
   wave0 amp-dev -0.009 -0.016 +0.014 +0.032 +0.011  (t = 0.75 1.75 2.75 3.75 4.75)
   wave1 amp-dev -0.014 -0.004 +0.025 +0.003 -0.030  (t = 0.75 1.75 2.75 3.75 4.75)
== MVA16 rho 60: sum g w^2 10.41; delta (A2) -0.0000; normalisation (t >= 3) 0.9733
   eps(r) -0.409 -0.102 -0.018 -0.007 -0.012 -0.017 -0.022 -0.025 -0.027 -0.029   (r = 0.125 .. 2.375)
   wave0 amp-dev +0.010 -0.020 +0.063 +0.070 -0.040  (t = 0.75 1.75 2.75 3.75 4.75)
   wave1 amp-dev +0.002 -0.002 +0.049 -0.031 -0.075  (t = 0.75 1.75 2.75 3.75 4.75)
== S16 rho 60: sum g w^2 2.03; delta (A2) n/a; normalisation (t >= 3) 0.8973
   eps(r) -0.596 -0.312 -0.179 -0.128 -0.109 -0.102 -0.101 -0.101 -0.102 -0.103   (r = 0.125 .. 2.375)
   wave0 amp-dev +0.127 -0.040 -0.067 +0.013 +0.032  (t = 0.75 1.75 2.75 3.75 4.75)
   wave1 amp-dev +0.109 -0.101 -0.043 +0.063 +0.008  (t = 0.75 1.75 2.75 3.75 4.75)
== MV2+ct rho 60: sum g w^2 6.18; delta (A2) -0.0594; counterterm -0.0594; normalisation (t >= 3) 0.9927
   eps(r) -0.476 -0.176 -0.065 -0.027 -0.013 -0.006 -0.004 -0.002 -0.001 -0.001   (r = 0.125 .. 2.375)
   wave0 amp-dev +0.029 +0.008 -0.017 -0.024 -0.000  (t = 0.75 1.75 2.75 3.75 4.75)
   wave1 amp-dev +0.030 -0.003 -0.020 +0.005 +0.024  (t = 0.75 1.75 2.75 3.75 4.75)
== MV8+ct rho 60: sum g w^2 2.48; delta (A2) -0.0958; counterterm -0.0958; normalisation (t >= 3) 0.9882
   eps(r) -0.556 -0.251 -0.109 -0.051 -0.026 -0.014 -0.008 -0.005 -0.003 -0.002   (r = 0.125 .. 2.375)
   wave0 amp-dev +0.048 +0.026 -0.034 -0.058 -0.010  (t = 0.75 1.75 2.75 3.75 4.75)
   wave1 amp-dev +0.053 +0.002 -0.046 +0.004 +0.057  (t = 0.75 1.75 2.75 3.75 4.75)
== MV16+ct rho 60: sum g w^2 1.67; delta (A2) -0.1283; counterterm -0.1283; normalisation (t >= 3) 0.9946
   eps(r) -0.598 -0.297 -0.138 -0.065 -0.030 -0.012 -0.002 +0.003 +0.006 +0.008   (r = 0.125 .. 2.375)
   wave0 amp-dev +0.054 +0.051 -0.066 -0.113 -0.005  (t = 0.75 1.75 2.75 3.75 4.75)
   wave1 amp-dev +0.064 +0.013 -0.083 +0.011 +0.111  (t = 0.75 1.75 2.75 3.75 4.75)

[exited with code 0]
```

How the run is read (development): ε(r), the normalisation and the mean dynamics of every row against the lines above
at ρ 30 and ρ 60, with no adjustment; the noise per kernel and its scaling with ρ (the diagonal model predicts sd ∝ ρ^−1/4
at fixed kernel); and the choice between MVA (constraint) and MV+ct (counterterm) is the author's.
