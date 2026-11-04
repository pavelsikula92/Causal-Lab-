# Causal Gravity: the simulator in plain words

`causal_gravity.py` is an interactive simulator of the gravity law of Causal Physics. You can watch planets,
light, clocks and probes move. The same scene is also computed in other theories of gravity, so you can see
where each of them succeeds and where it fails.

The simulator runs on its own and does not change anything in the laboratory. The formulas behind it are in the
appendix at the end of this file. They are derived in the paper *Physics of Causality*
([doi:10.5281/zenodo.23853689](https://doi.org/10.5281/zenodo.23853689)).

## What the simulator shows, honestly

In Causal Physics, gravity is not a force that sits in space. It is the way the causal order of events responds
to matter:

- near a mass, clocks run slower;
- light cones (the directions in which influence can travel) narrow;
- events become denser.

The laboratory shows three things:

- this response has exactly the structure of Einstein's gravity at the first order (tests F7.1c and F7.4a);
- it follows from the discrete action of the order (F6.8a, F6.8b2, F6.8c2);
- in the order of our own dimension, it gives the same bending of light as Einstein (γ = 1).

Beyond the first order the theory **adopts** general relativity as its completion. Deriving that step from the
causal order is still open.

So in scenes 0–21 causal gravity and Einstein agree. This is expected, not a coincidence, and the simulator does
not hide it. What these scenes add is a comparison with the descriptions that **fail**:

- Newton;
- the scalar first stage of the theory itself;
- a fixed causal order.

Scene 22 is different. It shows the one place where the theory predicts a **different number** from Einstein:
an old static field pulls more weakly at very short distances.

Some things are inputs of the simulator, not results of the theory:

- the masses and their positions;
- the strength of gravity G;
- the frame dragging of spinning bodies and the gravitational waves, both of which come from the adopted
  completion.

The *Theory status* tab lists, item by item, what is derived, what is tested and what is only assumed.

## Running it

    python Demo/causal_gravity.py                       the default scene
    python Demo/causal_gravity.py --scene 22            a scene by its number
    python Demo/causal_gravity.py --scene 22 --autorun  and start at once
    python Demo/causal_gravity.py --selftest            numerical checks of the code (no window)

It needs Python 3.9 or newer, with tkinter and numpy.

## The panels

- **Controls.** Choose the scene, the pace (simulated time per second of real time), and what is drawn. You can
  add or delete bodies, light and clocks:
  - a click places an object;
  - a drag sets its direction and speed.
- **Comparison.** Switch on other theories. Each one computes the same scene in parallel and is drawn with
  dashed trails. This tab also has two buttons:
  - *Run classical tests* measures the classical tests of gravity in every theory and prints a table next to
    the observations;
  - *Mercury experiment* measures Mercury's perihelion.
- **Theory status** and **About.** What is derived, what is tested, what is an input, and the formulas in short.
- **Measurement readout** (below the picture). What the scene measures, in numbers, as it runs.

### The theories that can be compared

| Theory | In plain words | Where it fails |
|---|---|---|
| Causal gravity (solid trails) | the gravity of the causal order; Einstein's beyond the first order | nowhere in scenes 0–21; differs only in scene 22 |
| Einstein (GR) | general relativity | the reference |
| Newton | force at a distance; light as a fast corpuscle | half the bending of light, no perihelion advance, no clock effects, no waves |
| Scalar stage (H1) | the first version of the theory: one clock field in a fixed frame | no frame dragging; radiates 1/6 of the observed power; moving matter does not gravitate correctly; no horizon |
| Fixed order (F3) | matter moves, but the causal order does not respond | light is not bent at all (γ = −1) |
| PPN | a general weak-field metric with two sliders, β and γ | depends on the sliders |

## The scenes

**0 Perihelion precession.** A planet on an elongated orbit around a star, exaggerated so that the effect is
visible. In Einstein's gravity the orbit slowly turns; this is the famous advance of Mercury's perihelion.
Newton gives no turning. The readout measures the angle per orbit in every theory.

**1 Solar system (exaggerated).** Four planets around a star, with strong gravity so that the relativistic
effects are visible within a few orbits.

**2 Light bending.** Light passes a mass at different distances. Its path bends towards the mass, and the closer
it passes, the more it bends. Einstein gives twice Newton's bending, which was confirmed in 1919. A fixed causal
order (F3) does not bend light at all. This is the clearest sign that the order itself must respond to matter.

**3 Shapiro delay.** Light passing close to a mass arrives later than light far from it, because clocks run
slower and space is stretched near the mass. The delay is measured with radar echoes from planets and spacecraft
(Cassini).

**4 Binary star and planet.** Two stars orbit each other and a planet orbits both. This is an example of
motion in the field of several bodies.

**5 Clocks in a field.** Clocks at different distances from a mass, and one clock on an orbit. The closer a clock
is to the mass, the slower it runs. The orbiting clock is also slowed by its motion. In Causal Physics a clock
rate is literally how fast the order runs at that place.

**6 Solar system (real, J2000).** The real planets with their real orbits, starting from the year 2000.

**7 Earth surface: apple and clocks.** An apple falls, a ball is thrown, and two clocks sit 30 m apart in
height. Everything is as everyday experience and precise measurement say: in a weak field every theory that
works agrees with Newton for falling bodies. The clocks differ by about 3 parts in 10¹⁵, which is the
Pound–Rebka effect.

**8 Earth: ISS and GPS.** The space station and a GPS satellite. A GPS clock gains about 38 microseconds a day
against a clock on the ground, because it sits higher (it runs faster) and moves (it runs slower). Without this
correction GPS would drift by kilometres a day.

**9 and 10 UAP.** A scripted manoeuvre and keyboard control of a craft that accelerates abruptly. These scenes
only show what such a motion would require: the accelerations a crew would feel. The "null field" option is a
**hypothesis of the author**, not part of the theory or of the laboratory.

What such a field would have to be: in the scripted manoeuvre the craft stops from 1000 m/s within 0.2 s, an
acceleration of a = 5000 m/s², about 510 g. A crew driven by an engine would feel all of it. A body that falls
freely in a gravitational field feels nothing, so a craft moved by gravity would carry its crew at 0 g. A field
that gives this acceleration needs a steepness of the potential (U in units of c²) of

  |∇U| = a/c² ≈ 5.6 · 10⁻¹⁴ m⁻¹,

about 510 times the steepness of the Earth's field at its surface (g/c² ≈ 1.1 · 10⁻¹⁶ m⁻¹). Nothing in the theory or
in the laboratory produces such a field around a craft; the number only states what the hypothesis would require.

**11 Black hole (strong field).** Light passing a black hole at different distances, clocks near the horizon, and
a particle on a tight orbit. Light that comes too close is captured, which produces the "shadow" seen by the
Event Horizon Telescope. Here causal gravity is Einstein's by adoption. The scalar stage has no horizon (a throat
instead) and is shown only for comparison; it is **not** a prediction of the theory.

**12 and 13 α Centauri A/B and Sirius A/B.** Real binary stars with their real orbits. Sirius B is a white dwarf.

**14 Solar neighbourhood to 12 ly (3D).** The real nearest stars, with their measured motions where they are
known. You can rotate the map.

**15 Free scene.** One fixed mass. Add your own bodies, light and clocks.

**16 Fast probes and light.** Probes pass a mass at 0.3–0.9 of the speed of light, together with a beam of
light. The faster a probe moves, the more it is bent: the bending grows from Newton's value at low speed to
twice Newton's for light. The fixed order bends fast probes *less*, which is wrong.

**17 Frame dragging (exaggerated).** A fast-spinning star drags space around with it. A planet orbiting with the
spin moves differently from one orbiting against it, and a gyroscope slowly turns. In the laboratory the local
mechanism is tested: a moving flow tilts the causal cones (F7.4a). The dragging by a rotating body is **not**
computed in the order; the scene uses Einstein's formula.

**18 Gravity Probe B (gyroscope).** The real experiment. Gyroscopes on a polar orbit around the Earth turn because of the curvature of space
and because the Earth's spin drags space. Einstein predicts 6606 and 39 milliarcseconds per year, and the experiment
measured 6602 and 37. The scalar stage cannot give the second number.

**19 Binary pulsar PSR B1913+16.** Two neutron stars lose energy by gravitational waves and their orbit shrinks.
The observed shrinking agrees with Einstein to a fraction of a percent; this was the first evidence of
gravitational waves. A scalar field would radiate only one sixth of the energy. In Causal Physics the waves come
from the adopted completion; how they propagate in the order has not been computed.

**20 Inspiral (strong radiation).** Two black holes spiral into each other as they radiate. This is the kind of
event LIGO observes.

**21 What gravitates: counter-streaming beams (F7.4a).** Pairs of beams stream in opposite directions at
different speeds. Their momentum cancels, but their energy and their pressure remain. A probe above each pair
measures how strongly it attracts. In Einstein's gravity, and in the laboratory test F7.4a, moving matter pulls
more strongly than the same matter at rest, by exactly Γ²(1 + v²). Newton and the scalar stage count only the
rest mass. This is the scene where the tensor nature of gravity is seen most directly, and it is a result of the
laboratory, not an input.

**22 Short range: the age of a field (prediction P2).** A tiny tungsten sphere of radius 20 µm and five
probes released at rest at distances from 60 to 250 µm. The solid trails are causal gravity and the dashed
trails are Einstein.

In Causal Physics the order is discrete, and its response is slightly non-local. A static field that has
existed for a long time T therefore pulls a little more weakly at very short distances. The weakening grows fast
as the distance shrinks: tenfold closer, a thousand times stronger. The readout shows three things:

- the predicted weakening at each probe;
- the weakening that the simulation actually measures;
- whether your choice of parameters is already **excluded** by the best laboratory test of gravity, which holds
  down to 52 µm (torsion balance, Lee et al. 2020).

You choose two things with the controls on the left:

- **The non-locality scale ℓ_k**, a free parameter of the theory.
- **What the age T of the field is.** This is open in the theory: it could be the age of the universe, of the
  Earth's field, or something much shorter.

With ℓ_k = 10⁻¹⁹ m and a cosmological T, the weakening at 52 µm would be almost total, so that choice is excluded.
Move the slider to about 10⁻²⁰ m and the theory passes the test, with the weakening setting in just below the
tested range.

This is the one place where the theory says something new and testable about gravity: **where** the deviation
from Einstein should appear. The scene marks its limits:

- the form of the term has been computed only while the weakening is small (up to about 30 %);
- the dark disc in the field map is where the simulator simply caps the weakening at 100 %;
- the meaning of T is not settled.

## Appendix: the formulas and where they come from

Units: G = c = 1 unless a scene says SI. The section numbers are those of the paper (November 2026).

- **Metric of causal gravity and of Einstein** (static weak field, superposed point masses, U = Σ M/r):
  ds² = −A dt² + B dx², with A = ((1 − U/2)/(1 + U/2))² and B = (1 + U/2)⁴. Paper: Sec. 38 "The static section to first
  post-Newtonian order" and Sec. 39 "Thermodynamic completion"; the completion is adopted.
- **Motion** along the geodesics of this metric:
  dv/dt = −(A′/2B)∇U − v(v·∇U)B′/B + v²(B′/2B)∇U + v(v·∇U)A′/A + v × (∇ × g₀ᵢ). Paper: Sec. 21 "The action of a
  flow" (proper time as the action).
- **γ = 1:**
  - the throughput principle gives γ = 2/(d − 1);
  - pure tensor coupling gives γ = 1/(d − 2);
  - the two agree for d = 3 spatial dimensions;
  - from the action of the order, γ = −c_ΦΨ/(2c_ΨΨ) = 1.

  Paper: Sec. 22 and 33.2, Sec. 35 "The rule from the action of the order", Sec. 36 "The Benincasa–Dowker action
  to second order".
- **Moving matter as a source (scene 21):** the active mass of counter-streaming beams relative to rest is
  Γ²(1 + v²), from T₀₀ + Tᵢᵢ. Paper: Sec. 34 "The substrate rule and its test".
- **Frame dragging (scenes 17, 18):** g₀ᵢ = −2(J × r)ᵢ/r³; gyroscope precession (γ + ½) v × ∇U (geodetic) and
  ((1 + γ)/2)[3(J·n)n − J]/r³ (dragging). These come from general relativity. Paper: Sec. 33.2 (status) and
  Sec. 34.2 (the tilt of the cones by a moving flow, tested).
- **Radiation (scenes 19, 20):** the 2.5PN radiation-reaction force of a binary, and the Peters formula
  dP/dt = −(192π/5)(2πM_c/P)^(5/3) f(e). These come from general relativity. Paper: Sec. 33.3 (why a scalar field
  radiates 1/6) and Sec. 40.
- **Classical tests at first post-Newtonian order:**
  - perihelion: (2 + 2γ − β)/3 · 6πM/(a(1 − e²));
  - light bending: (1 + γ)(M/b)(sin θ₂ − sin θ₁) + π(8 + 8γ − 4β + 3ε)/4 (M/b)²;
  - fast probes: bending (1 + γv²) times Newton;
  - Shapiro delay: (1 + γ) M ln((x₂ + r₂)/(x₁ + r₁)).

  Paper: Sec. 38 "The static section to first post-Newtonian order", Sec. 45 "Consistency with general
  relativity".
- **The age term (scene 22):**
  - weakening at wavenumber k: δ(k) ≃ c_T k³ ℓ_k² T, with c_T = 0.1035;
  - the simulator uses k = 1/r, so δ(r) = c_T ℓ_k² cT/r³, capped at 1;
  - the weakening reaches 1 at L* = (c_T ℓ_k² cT)^(1/3).

  Paper: Sec. 46 "The window of the non-locality scale" and Sec. 47.2 "P2: the age term". The term has been
  computed in the conformal sector (F6.8a). It is 1.78 times larger in the mixed sector, and its form at large kT
  in the optical sector is open. The scene uses the conformal coefficient as an estimate of the size.
