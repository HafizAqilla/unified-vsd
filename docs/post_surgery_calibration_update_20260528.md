# Post-Surgery Calibration Update - 2026-05-28

This note documents the post-surgery calibration changes for the Reyna VSD
model. The goal is to keep the post-operative model physiologically traceable
when echo-derived ventricular volumes, ejection fraction, cardiac output, and
catheter-style pressure targets do not all come from one perfectly identical
measurement state.

## Why The Calibration Changed: Resolving Protocol Inconsistency

The post-surgery calibration was revised because the handwritten clinical protocol ("protokol") contains severe, mathematically impossible inconsistencies. Attempting to fit the lumped-parameter model to these conflicting targets forced the optimizer into a physiological deadlock, resulting in poor model performance.

### 1. Mathematical Inconsistencies in the Handwritten Protocol

The handwritten protocol reports the following values for post-surgery LV function:
*   $\text{LVEDV} = 32.0\text{ mL}$ (Row 23)
*   $\text{LVESV} = 23.6\text{ mL}$ (Row 24)
*   $\text{Stroke Volume (SV)} = 26.4\text{ mL}$ (Row 27)
*   $\text{Ejection Fraction (EF)} = 61.8\%$ (Row 29)
*   $\text{Cardiac Output (CO)} = 3.0\text{ L/min}$ (Row 28)

These values contain two direct mathematical contradictions:
1.  **Stroke Volume Discrepancy**: Standard physiology requires $\text{SV} = \text{EDV} - \text{ESV}$. Using the protocol's raw volumes yields an implied stroke volume of $32.0 - 23.6 = \mathbf{8.4\text{ mL}}$. This is a **$214\%$ discrepancy** compared to the reported SV of $26.4\text{ mL}$.
2.  **Ejection Fraction Discrepancy**: Standard physiology requires $\text{EF} = \text{SV} / \text{EDV}$. Using the protocol's raw volumes yields an implied ejection fraction of $8.4 / 32.0 = 26.25\%$. This is a **$135\%$ discrepancy** compared to the reported EF of $61.8\%$.

### 2. Mathematical Elaboration: Derivation of EDV, SV, and CO

To ensure absolute traceability and scientific integrity under the **Vibecoding Guardrails**, the exact mathematical derivation of the post-operative target parameters from the two raw, directly measured Biplane Simpson (BP) values is detailed below.

```mermaid
graph TD
    %% Direct measurements styling
    subgraph "Direct Measurements (Echo Screen)"
        ESV["LVESV (BP) = 13.8 mL"]
        EF["LVEF (BP) = 61.8%"]
        HR["HR = 114 bpm"]
    end

    %% Derived values styling
    subgraph "Derived Targets"
        EDV["LVEDV = 36.1 mL"]
        SV["LVSV = 22.3 mL"]
        CO["LVCO = 2.54 L/min"]
    end

    %% Mathematical flows
    ESV -->| "LVEDV = ESV / (1 - EF)" | EDV
    EF -->| "LVEDV = ESV / (1 - EF)" | EDV
    
    EDV -->| "LVSV = EDV - ESV" | SV
    ESV -->| "LVSV = EDV - ESV" | SV
    
    SV -->| "CO = SV * HR / 1000" | CO
    HR -->| "CO = SV * HR / 1000" | CO
    
    %% Styles
    classDef measured fill:#e1f5fe,stroke:#039be5,stroke-width:2px,color:#01579b;
    classDef derived fill:#efebe9,stroke:#5d4037,stroke-width:2px,color:#3e2723;
    
    class ESV,EF,HR measured;
    class EDV,SV,CO derived;
```

#### A. Parameter Summary Table
This table maps our direct clinical source variables to our derived mathematical calibration targets:

| Parameter | Source Type | Clinical Value | Mathematical Role in Model |
|---|---|---|---|
| `post.LVESV_mL` | **Direct Measurement** | **13.8 mL** (Biplane Simpson ESV) | Measured Source / Validation Holdout |
| `post.EF` | **Direct Measurement** | **61.8%** (Biplane Simpson EF) | Measured Source / Calibration Target |
| `post.HR` | **Hemodynamic Cath** | **108 bpm** (Override from catheter) | Sim Timing / Cycle Period ($0.556\text{ s}$) |
| `post.LVEDV_mL` | *Derived* | **36.1 mL** | Calibration Target |
| `post.SV_lv_mL` | *Derived* | **22.3 mL** | Calibration Target |
| `post.CO_Lmin` | *Derived* | **2.54 L/min** (at $114\text{ bpm}$) | Validation Target / Echo Benchmark |

---

#### B. Step-by-Step Mathematical Derivations

> [!IMPORTANT]
> **1. End-Diastolic Volume (LVEDV)**
> Ejection Fraction ($\text{EF}$) represents the fraction of blood ejected from the left ventricle during systole:
> $$\text{EF} = \frac{\text{LVEDV} - \text{LVESV}}{\text{LVEDV}} = 1 - \frac{\text{LVESV}}{\text{LVEDV}}$$
> Isolating and solving for the unknown $\text{LVEDV}$:
> $$\text{LVEDV} = \frac{\text{LVESV}}{1 - \text{EF}}$$
> Substituting the clinical baseline values:
> $$\text{LVEDV}_{\text{BP}} = \frac{13.8\text{ mL}}{1 - 0.618} = \frac{13.8\text{ mL}}{0.382} \approx \mathbf{36.13\text{ mL}}$$

---

> [!IMPORTANT]
> **2. Stroke Volume (LVSV)**
> Stroke Volume ($\text{SV}$) is the volume of blood pumped out of the left ventricle per contraction, computed as the difference between maximum ($\text{EDV}$) and minimum ($\text{ESV}$) ventricular volumes:
> $$\text{LVSV} = \text{LVEDV} - \text{LVESV}$$
> Substituting the derived $\text{LVEDV}$ and measured $\text{LVESV}$:
> $$\text{LVSV}_{\text{BP}} = 36.13\text{ mL} - 13.8\text{ mL} = \mathbf{22.33\text{ mL}}$$

---

> [!IMPORTANT]
> **3. Cardiac Output (LVCO)**
> Cardiac Output ($\text{CO}$) is the total volume of blood pumped by the heart per minute. To scale the flow rate from milliliters per minute ($\text{mL/min}$) to clinical reporting units ($\text{L/min}$), the product of Stroke Volume ($\text{SV}$) and Heart Rate ($\text{HR}$) is divided by $1000\text{ mL/L}$:
> $$\text{CO}\text{ [L/min]} = \frac{\text{LVSV}\text{ [mL]} \times \text{HR}\text{ [bpm]}}{1000}$$
> 
> *   **Echo-Derived Cardiac Output (at Measurement HR = 114 bpm)**:
>     $$\text{CO}_{\text{echo}} = \frac{22.33\text{ mL} \times 114\text{ bpm}}{1000} \approx \mathbf{2.54\text{ L/min}}$$
>     *This represents the static clinical comparator target.*
> 
> *   **Simulation-Aligned Cardiac Output (at Model HR = 108 bpm)**:
>     $$\text{CO}_{\text{simulation}} = \frac{22.33\text{ mL} \times 108\text{ bpm}}{1000} \approx \mathbf{2.41\text{ L/min}}$$
>     *This represents the steady-state flow rate that the cardiovascular simulation naturally converges to.*

---

### 3. The Solution: Calibration to Pure Echo Screen Data

To resolve these contradictions, we bypass the handwritten protocol and calibrate directly to a single, self-consistent dataset from the raw echo screen. The **Biplane Simpson (BP)** dataset is selected because it is the clinical gold standard for reconstructing the 3D LV cavity:

| Echo Screen Field (Biplane) | Model Variable | Value | Role |
|---|---|---|---|
| `ESV (BP)` | `post.LVESV_mL` | 13.8 mL | Measured Source |
| `EF (BP)` | `post.EF` | 0.618 | Measured Source |
| Derived EDV | `post.LVEDV_mL` | 36.1 mL | Derived Calibration Target |
| Derived SV | `post.SV_lv_mL` | 22.3 mL | Derived Calibration Target |
| Derived CO (at 114 bpm) | `post.CO_Lmin` | 2.54 L/min | Derived Calibration Target |

**Note**: To maintain absolute consistency, the Cardiac Output target is also derived **purely from the Biplane Simpson (BP)** dataset ($22.3\text{ mL} \times 114\text{ bpm} = 2.54\text{ L/min}$).

Biplane Simpson values are preferable to a single-plane estimate because they are based on two apical views and better represent LV volume. Therefore, the post-operative LV block now treats BP ESV and BP EF as the measured source values, then derives EDV and stroke volume from them.


## How The LV Target Block Works Now

The active post-surgery target block is in `run_post_surgery.m` inside
`apply_post_surgery_targets`.

```matlab
post.LVESV_mL      = 13.8;     % [mL] post-op echo, biplane Simpson ESV
post.EF            = 0.618;    % [-] post-op echo, biplane Simpson EF = 61.8%
post.LVEDV_mL      = post.LVESV_mL / (1 - post.EF); % [mL] derived from ESV/EF
post.SV_lv_mL      = post.LVEDV_mL - post.LVESV_mL; % [mL] derived LV stroke volume
post.LVESV_reported_mL = post.LVESV_mL; % [mL] direct BP echo entry
```

The calculation is:

```text
LVEDV = LVESV / (1 - EF)
      = 13.8 / (1 - 0.618)
      = 36.1 mL

LVSV  = LVEDV - LVESV
      = 36.1 - 13.8
      = 22.3 mL
```

This ordering matters. `LVESV_mL` must be assigned before `LVEDV_mL` and `SV_lv_mL` are derived. If the derived values are computed before assigning `LVESV_mL`, MATLAB uses the `NaN` initialized in `patient_reyna.m`, causing the post-operative LV volume targets to become unavailable.

## Heart Rate And Cardiac Output Handling

The model uses one heart rate for the post-surgery steady-state simulation:

```matlab
post.HR = 108; % [bpm]
```

`params_from_clinical.m` now applies a post-surgery HR override. The sequence is:

1. Load the default/scaled parameter set.
2. Apply `clinical.common.HR` if available.
3. For `post_surgery`, override with `clinical.post_surgery.HR` when finite.
4. Recompute chamber timing from the selected HR.

This means the post-operative steady state runs at:

```text
T_HB = 60 / 108 = 0.556 s/cycle
```

The current `post.CO_Lmin` is set to:

```matlab
post.CO_Lmin = 2.54; % [L/min] echo-derived LVCO = BP SV * 114 / 1000
```

This value corresponds to the echo screen heart rate of 114 bpm:

```text
CO_echo = 22.3 mL/beat * 114 bpm / 1000
        = 2.54 L/min
```

By deriving the Cardiac Output target directly from the Biplane Stroke Volume ($\text{SV}_{\text{BP}} = 22.3\text{ mL}$), we ensure that all ventricular target parameters (volumes, EF, and flows) are defined using a single, self-consistent Biplane Simpson dataset.

If cardiac output should be fully aligned to the hemodynamic model HR of 108 bpm, then the HR-aligned value would be:

```text
CO_HR108 = 22.3 mL/beat * 108 bpm / 1000
         = 2.41 L/min
```

The 2.54 L/min target is therefore interpreted as a pure Biplane echo-derived clinical comparator at the measurement heart rate of 114 bpm. If a direct catheter/Fick cardiac output is available from the same SAP/PAP/RAP measurement state, that value should supersede the echo-derived CO target.

## Target-Tier Governance

The calibration target system in `build_target_tiers.m` exposes post-operative ventricular metrics while using robust clinical data governance to avoid double-counting. 

Because `LVEDV` and `CO_Lmin` are algebraically derived from the direct source measurements of `LVESV` and `EF`, including all of them in active calibration would artificially overweight the echo block. To enforce data integrity, the system applies the following roles:

| Metric | Clinical field | Role | Inclusion in Calibration | Inclusion in Primary RMSE |
|---|---|---|---|---|
| `LVESV` | `post.LVESV_mL` | **Direct Measured Source** (Hard Target) | **YES** | **YES** |
| `LVEF` | `post.EF` | **Direct Measured Source** (Hard Target) | **YES** | **YES** |
| `LVEDV` | `post.LVEDV_mL` | *Algebraic Derivative* (Derived Validation) | NO | NO |
| `CO_Lmin` | `post.CO_Lmin` | *Algebraic Derivative* (Derived Validation) | NO | NO |
| `RVEDV` | `post.RVEDV_mL` | Soft Post-Surgery Target | **YES** | **YES** |
| `RVESV` | `post.RVESV_mL` | Soft Post-Surgery Target | **YES** | **YES** |
| `RVEF` | `post.RVEF` | Validation Target | NO | **YES** (Reporting Only) |

### Active Calibration Policy:
*   **Direct Measured Pair Fit**: The optimizer directly fits **`LVESV`** and **`LVEF`** (the independent raw echo measurements). 
*   **Algebraic Derivations Excluded**: The derived algebraic consequences, **`LVEDV`** and **`CO_Lmin`**, are set to `derived_validation` (retained as validation benchmarks for reporting, but excluded from active calibration and primary RMSE computation to avoid double-counting).
*   **Right-Sided Soft Targets**: **`RVEDV`** and **`RVESV`** are actively calibrated as softer targets to guide right ventricular sizing.

## Interpretation

The post-surgery model should now be read as a pressure-flow-volume calibration
with explicit source awareness:

- SAP/PAP/RAP targets describe the hemodynamic pressure state.
- `post.HR = 108 bpm` defines the simulation cycle timing.
- BP Simpson `LVESV = 13.8 mL` and `EF = 0.618` define the LV echo volume block.
- CO is currently retained as an echo-derived comparator at 114 bpm.
- Validation reports should distinguish active fitting targets from holdout
  consistency checks.

This prevents the model from being judged against silently inconsistent targets.
When the model misses a target, the validation table should make clear whether
the mismatch reflects physiology, calibration error, or source inconsistency
between echo-derived and hemodynamic measurements.
