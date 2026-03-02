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

| Parameter class | Exponent | Physical derivation |
|---|---|---|
| Resistance R | $s^{-1}$ | Poiseuille: $R \propto L/r^4$; $r \propto s^{1/2}$ |
| Compliance C | $s^{+1}$ | $C \propto r^2 \cdot L \propto s$ |
| Inertance L | $s^{-1}$ | $L \propto \rho L / r^2$ |
| Unstressed V₀ | $s^{+1}$ | Compartment volume ∝ body size |
| Heart rate HR | $s^{-0.33}$ | Metabolic (basal) scaling |
| E_LV/LA | $s^{-1}$ | Systemic pressure load |
| E_RV/RA | $s^{-1.5}$ | Pulmonary (lower pressure load) |

Reference: Lundquist et al. (2025), Table 3.1.

---

## 5. Pre-surgery R_VSD computation

### Option A — Ohm's law (preferred)

$$R_{VSD} = \frac{\Delta P_{\text{VSD}}}{Q_{\text{shunt}}}$$

where $\Delta P_{\text{VSD}}$ is the peak gradient (mmHg) from echo Doppler
and $Q_{\text{shunt}}$ is the net shunt flow (mL/s) from catheterisation.

### Option B — Gorlin orifice equation

$$Q = C_c \cdot A \cdot \sqrt{\frac{2 \Delta P}{\rho}}, \quad R_{VSD} = \frac{\Delta P_{\text{ref}}}{Q}$$

with $C_c = 0.7$ (contraction coefficient), $\rho = 1060\,\mathrm{kg/m^3}$
(blood density), and $A = \pi (D/2)^2$ determined from echo-measured defect
diameter $D$.

---

## 6. Calibration objective

$$J(x) = \sum_k w_k \left(\frac{y_k(x) - y_k^{\text{clin}}}{y_k^{\text{clin}}}\right)^2
+ \lambda \sum_i \left(\frac{x_i - x_0^i}{x_0^i}\right)^2$$

Free parameters per scenario: see `calibration/calibration_param_sets.m`.

---

## 7. Sobol sensitivity analysis

First-order ($S_1^i$) and total ($S_T^i$) Sobol indices via Jansen (1999) estimator:

$$S_T^i = \frac{\mathbb{E}\left[(f_A - f_{A_B^i})^2\right]}{2\,\mathrm{Var}(Y)}$$

$$S_1^i = 1 - \frac{\mathbb{E}\left[(f_B - f_{A_B^i})^2\right]}{2\,\mathrm{Var}(Y)}$$

Saltelli (2010) quasi-random sampling with Sobol sequences.

---

## Assumptions

1. Incompressible, Newtonian flow throughout.
2. Valve inertia neglected (two-state resistance model).
3. Pulmonary capillary represented by a single lumped pressure state P_PC
   (Valenti Eq. 2.7), not separate oxygenated/non-oxygenated volume states.
4. No coronary circulation, no pericardium, no ventricular interaction.
5. Allometric scaling uses BSA-based geometric power laws; metabolic scaling
   for HR only.

---

## References

1. Valenti (2023). Thesis: *Full-order 0D cardiovascular model*. Eqs. 2.1–2.7, Table 3.3.
2. Lundquist et al. (2025). *Allometric Scaling in Pediatric Cardiology*. Eqs. 3.2, 4.1, 4.3.
3. Saltelli et al. (2010). *Variance based sensitivity analysis of model output*. CPC 181:259–270.
4. Jansen (1999). *Analysis of variance designs for model output*. CPC 117:35–43.
