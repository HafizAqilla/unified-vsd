# Theory Notes — Unified VSD Lumped-Parameter Model

## 1. Model overview

This model represents the cardiovascular system as a network of electrical
circuit analogues — the *lumped-parameter* or *0-D* approach.  Each anatomical
compartment is modelled by a resistor (R), capacitor (C), and/or inductor (L):

| Analogue | Cardiovascular meaning | Unit |
|---|---|---|
| Voltage (V) | Pressure P | mmHg |
| Current (I) | Volumetric flow Q | mL/s |
| Capacitor C | Vascular compliance | mL/mmHg |
| Resistor R | Vascular resistance | mmHg·s/mL |
| Inductor L | Blood inertance | mmHg·s²/mL |

---

## 2. State vector (14 states)

The ODE is implemented in `models/system_rhs.m`.

| Index | State | Unit | Compartment |
|---|---|---|---|
| 1 | V_RA | mL | Right atrium volume |
| 2 | V_RV | mL | Right ventricle volume |
| 3 | V_LA | mL | Left atrium volume |
| 4 | V_LV | mL | Left ventricle volume |
| 5 | V_SAR | mL | Systemic arterial volume |
| 6 | Q_SAR | mL/s | Systemic arterial flow (inductor state) |
| 7 | V_SC | mL | Systemic capillary volume |
| 8 | V_SVEN | mL | Systemic venous volume |
| 9 | Q_SVEN | mL/s | Systemic venous flow (inductor state) |
| 10 | V_PAR | mL | Pulmonary arterial volume |
| 11 | Q_PAR | mL/s | Pulmonary arterial flow (inductor state) |
| 12 | P_PC | mmHg | Combined pulmonary capillary pressure (Valenti Eq. 2.7) |
| 13 | V_PVEN | mL | Pulmonary venous volume |
| 14 | Q_PVEN | mL/s | Pulmonary venous flow (inductor state) |

---

## 3. Governing equations

### 3.1 Chamber pressure  (Valenti Eq. 2.1–2.2)

$$P_{\text{ch}}(t) = E_{\text{ch}}(t) \cdot \left(V_{\text{ch}} - V_0\right)$$

where $E_{\text{ch}}(t)$ is the time-varying elastance (Section 3.2) and $V_0$
is the unstressed (zero-pressure) volume.

### 3.2 Elastance model  (Valenti Eqs. 2.2–2.4)

$$E_{\text{ch}}(t) = E_A \, e(t) + E_B$$

The normalised activation $e(t) \in [0, 1]$ follows a piecewise cosine:

**Ventricular** ($e_v$):

$$e_v(\phi) = \begin{cases}
\tfrac{1}{2}\left(1 - \cos\tfrac{\pi \phi}{T_c}\right) & 0 \le \phi \le T_c \\
\tfrac{1}{2}\left(1 + \cos\tfrac{\pi (\phi - T_c)}{T_r}\right) & T_c < \phi \le T_c + T_r \\
0 & \text{otherwise}
\end{cases}$$

**Atrial** ($e_a$): identical piecewise form but time-shifted by $t_{ac}$.

### 3.3 Valve flow  (Valenti Eq. 2.6)

$$Q = \frac{P_{\text{up}} - P_{\text{dn}}}{R_{\text{valve}}}$$

$$R_{\text{valve}} = \begin{cases} R_{\min} & P_{\text{up}} > P_{\text{dn}} \\ R_{\max} & \text{otherwise} \end{cases}$$

### 3.4 RLC vascular compartments  (Valenti Eq. 2.5)

$$L \frac{dQ}{dt} = P_{\text{in}} - P_{\text{out}} - R \, Q$$

$$\frac{dV}{dt} = Q_{\text{in}} - Q_{\text{out}}, \quad P = \frac{V - V_0}{C}$$

### 3.5 Pulmonary capillary  (Valenti Eq. 2.7)

$$\frac{dP_{PC}}{dt} = \frac{Q_{PAR} - Q_{COX} - Q_{CNO}}{C_{PC,\text{total}}}$$

$$Q_{COX} = \frac{P_{PC} - P_{PVEN}}{R_{PCOX}}, \quad Q_{CNO} = \frac{P_{PC} - P_{PVEN}}{R_{PCNO}}$$

### 3.6 VSD shunt

$$Q_{VSD} = \frac{P_{LV} - P_{RV}}{R_{VSD}}$$

- $Q_{VSD} > 0$: left-to-right (physiologically typical, pressure-driven)
- Post-surgery: $R_{VSD} \to 10^6$ mmHg·s/mL, so $Q_{VSD} \approx 0$

---

## 4. Allometric scaling

Scale factor: $s = \mathrm{BSA}_{\text{patient}} / \mathrm{BSA}_{\text{ref}}$,
where $\mathrm{BSA}_{\text{ref}} = 1.73\,\mathrm{m}^2$ (adult male, Mosteller).

| Parameter class | Exponent | Source | Physical derivation |
|---|---|---|---|
| Resistance R | $s^{-1}$ | Lundquist (2025) Table 1 | Poiseuille: $R \propto L/r^4$; $r \propto s^{1/2}$ |
| Compliance C | $s^{+1}$ | Lundquist (2025) Table 2 | $C \propto r^2 \cdot L \propto s$ |
| Inertance L | $s^{-1}$ | Geometry analogous to R | $L \propto \rho L / r^2$ |
| Unstressed V₀ | $s^{+1}$ | Lundquist (2025) Table 2 | Compartment volume ∝ body size |
| Heart rate HR | $s^{-0.33}$ | Lundquist (2025) Table 2 | Metabolic (basal) scaling |
| E_LV/LA | $s^{-1}$ | Lundquist (2025) Table 2 | Systemic pressure load |
| E_RV/RA | $s^{-1.5}$ | Lundquist (2025) Table 2 | Pulmonary (lower pressure load) |
| Valve open area | $s^{+1}$ → $R_{\text{valve}} \propto s^{-1}$ | Lundquist (2025) Table 3 | Orifice flow: $Q \propto A$ at const $\Delta P$ |
| VSD shunt area | $s^{+1}$ → $R_{\text{VSD}} \propto s^{-1}$ | Lundquist (2025) Table 3 | Patient-specific R_VSD from `params_from_clinical.m` |

> **BSA reference note.** The code uses $\mathrm{BSA}_{\text{ref}} = 1.73\,\mathrm{m}^2$, which corresponds to the standard 70 kg / 175 cm adult implicit in the Valenti (2023) parameter set. The Lundquist (2025) paper's own reference subject is a 30-year-old, 85 kg, 180 cm male ($\mathrm{BSA} \approx 2.06\,\mathrm{m}^2$) — this difference is intentional: only the **scaling exponents** are borrowed from Lundquist; the **absolute parameter values** derive from Valenti.

> **Parameters not scaled.** (a) `epsilon_valve` (tanh smoothing pressure, 0.5 mmHg): Lundquist Table 3 lists valve opening/closing constants at $s^{-0.5}$, but those are physiological pressure thresholds in the Lundquist model; `epsilon_valve` is a numerical continuity parameter in our tanh valve model and is not equivalent. (b) `R_VSD`: pathological parameter set per patient by `params_from_clinical.m` — not scaled allometrically.

---

## 5. Pre-surgery R_VSD computation

### Option A — Ohm's law (preferred)

$$R_{VSD} = \frac{\Delta P_{\text{VSD}}}{Q_{\text{shunt}}}$$

where $\Delta P_{\text{VSD}}$ is the peak gradient (mmHg) from echo Doppler
and $Q_{\text{shunt}}$ is the net shunt flow (mL/s) from catheterisation.

### Option B — Gorlin orifice equation

$$Q = C_c \cdot A \cdot \sqrt{\frac{2 \Delta P}{\rho}}, \quad R_{VSD} = \frac{\Delta P_{\text{ref}}}{Q}$$

with:
- $C_c = 0.7$ (contraction / discharge coefficient) — Gorlin & Gorlin (1951); Baumgartner et al. (2009)
- $\rho = 1060\,\mathrm{kg/m^3}$ (blood density) — Levick (2010), p. 13
- $A = \pi (D/2)^2$ determined from echo-measured defect diameter $D$
- $\Delta P_{\text{ref}} = 20\,\mathrm{mmHg}$ reference gradient (typical peak LV–RV gradient)

> **Note:** Option A is preferred when direct catheterisation shunt-flow data are available;
> Option B introduces additional uncertainty from $C_c$ and $\Delta P_{\text{ref}}$ choices.

---

## 6. Calibration

### 6.1 Two-phase architecture

**Phase 0 (deterministic, zero optimiser calls):** `utils/params_from_clinical.m`
maps the 28 clinical inputs to model parameters via Ohm's law and conservation
principles. This ensures the starting point is physically self-consistent.

**Phase 1 (Nelder-Mead polish):** `calibration/run_calibration.m` refines
elastances and `R.vsd` using `fminsearch` (Nelder-Mead simplex, Torczon 1991).
Nelder-Mead was chosen because:
1. The problem is smooth and low-dimensional ($d \le 5$).
2. `fmincon` (gradient-based) fails here — finite-difference step sizes
   (~$10^{-8}$) are swamped by ODE solver noise (~$10^{-6}$), producing
   random gradient estimates and zero-length steps.
3. `particleswarm` converges to the same solution in 4× more evaluations.

### 6.2 Objective function

$$J(x) = \sum_k w_k \left(\frac{y_k(x) - y_k^{\text{clin}}}{y_k^{\text{clin}}}\right)^2$$

No Tikhonov regularisation (`regLambda = 0`): the analytic pre-conditioning in
Phase 0 already provides a good starting point; regularisation was found to
prevent parameters from moving far enough to reach clinical targets.

Bound enforcement: quadratic penalty added to $J$ for parameter violations,
allowing `fminsearch` (which has no native bound support) to operate within
physiological ranges:
$$J_{\text{pen}}(x) = J(x) + \lambda_{\text{pen}} \sum_i \left[\max(0,\, l_i - x_i)^2 + \max(0,\, x_i - u_i)^2\right]$$

where $\lambda_{\text{pen}} = 100$.

### 6.3 Free parameters and bounds

**Pre-surgery** (5 free parameters):

| Parameter | Physical role | Bounds |
|---|---|---|
| `R.vsd` | Shunt resistance →控制 Qp/Qs | 0.05× – 20× $x_0$ |
| `E.LV.EA` | LV active elastance → LVEF, LVESV | 0.3× – 3.0× $x_0$ |
| `E.LV.EB` | LV diastolic stiffness → LVEDV | 0.3× – 3.0× $x_0$ |
| `E.RV.EA` | RV active elastance → RVEF, RVESV | 0.3× – 3.0× $x_0$ |
| `E.RV.EB` | RV diastolic stiffness → RVEDV | 0.3× – 3.0× $x_0$ |

**Post-surgery** (11 free parameters): R.SAR, R.SVEN, R.PAR, R.PCOX, R.PVEN,
C.SAR, C.PAR, + 4 elastances. See `calibration_param_sets.m` for full details.

### 6.4 Speed settings inside `objective_calibration.m`

Calibration overrides steady-state settings for speed per evaluation:

| Setting | Calibration evals | Baseline / Validation |
|---|---|---|
| `nCyclesSteady` | 40 | 80 |
| `ss_tol_P` | 0.5 mmHg | 0.1 mmHg |
| `ss_tol_V` | 0.5 mL | 0.1 mL |

---

## 7. Sobol sensitivity analysis

First-order ($S_1^i$) and total ($S_T^i$) Sobol indices via Jansen (1999) estimator:

$$S_T^i = \frac{\mathbb{E}\left[(f_A - f_{A_B^i})^2\right]}{2\,\mathrm{Var}(Y)}$$

$$S_1^i = 1 - \frac{\mathbb{E}\left[(f_B - f_{A_B^i})^2\right]}{2\,\mathrm{Var}(Y)}$$

Saltelli (2010) quasi-random sampling with Sobol sequences.

**Sample size:** $N = 256$ per direction, giving $N(d+2) \geq 4608$ model evaluations
for $d = 18$ parameters. This satisfies the Saltelli (2010) recommendation for
precision of $S_T$ within $\pm 0.05$ (95 % confidence).

---

## Assumptions

1. Incompressible, Newtonian flow throughout.
2. Valve inertia neglected (two-state resistance model).
3. Pulmonary capillary represented by a single lumped pressure state P_PC
   (Valenti Eq. 2.7), not separate oxygenated/non-oxygenated volume states.
4. No coronary circulation, no pericardium, no ventricular interaction.
5. Allometric scaling uses BSA-based geometric power laws; metabolic scaling
   for HR only.
6. Diastolic chamber stiffness is modelled as a **linear** passive elastance $E_B$
   (Valenti Eq. 2.2), not as the exponential pressure-volume law used in some
   models (Lundquist Table 2 lists a `Lambda` parameter for exponential fitting;
   this is not applicable to the Valenti formulation).
7. `epsilon_valve` (valve smoothing constant, 0.5 mmHg) is **not** scaled with
   body size. Lundquist (2025) Table 3 scales valve opening/closure constants at
   $s^{-0.5}$, but those refer to physiological pressure thresholds, not the
   numerical tanh continuity parameter used here.

---

## References

1. Valenti (2023). Thesis: *Full-order 0D cardiovascular model*. Eqs. 2.1–2.7, Table 3.3.
2. Lundquist et al. (2025). *Patient-Specific Pediatric Cardiovascular Lumped Parameter Modeling — Scaling Across Ages and Sizes*. ASAIO J. DOI: 10.1097/MAT.0000000000002528. **Table 2** (cardiac: HR, E, V₀, blood volume) and **Table 3** (valve open/closed area, shunt area, ventilation).
3. Saltelli et al. (2010). *Variance based sensitivity analysis of model output*. CPC 181:259–270.
4. Jansen (1999). *Analysis of variance designs for model output*. CPC 117:35–43.
5. Gorlin R & Gorlin SG (1951). *Hydraulic formula for calculation of the area of the stenotic mitral valve*. Am Heart J 41(1):1–29.  (orifice Cd = 0.7)
6. Baumgartner H et al. (2009). *Echocardiographic assessment of valve stenosis: EAE/ASE recommendations*. J Am Soc Echocardiogr 22(1):1–23.  (Cd confirmation)
7. Levick JR (2010). *An Introduction to Cardiovascular Physiology*, 5th ed. Hodder Arnold, London. p. 13.  (rho_blood = 1060 kg/m³)
8. Nichols WW & O’Rourke MF (2005). *McDonald’s Blood Flow in Arteries*, 5th ed. Hodder Arnold, London.  (physiological reference ranges)
