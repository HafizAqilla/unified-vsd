# Vibecoding Guardrails
## World-Class MATLAB Standards for Cardiovascular Simulation Teams

---

> **The Prime Directive**
>
> *"If someone reads this code 3 years from now — without access to your brain, your thesis, or your supervisor — can they (a) understand the physiology, (b) reproduce your results, and (c) trust the numbers?"*
>
> If the answer to any of these is **no**, the code is not ready to commit.

---

## Table of Contents

1. [Philosophy](#1-philosophy)
2. [Project Structure](#2-project-structure)
3. [Naming Conventions](#3-naming-conventions)
4. [Function Standards](#4-function-standards)
5. [Parameter Handling](#5-parameter-handling)
6. [Units Policy](#6-units-policy)
7. [State Vector Management](#7-state-vector-management)
8. [Numerical Practices](#8-numerical-practices)
9. [Clinical Data Handling](#9-clinical-data-handling)
10. [Validation Requirements](#10-validation-requirements)
11. [Plotting Standards](#11-plotting-standards)
12. [Version Control](#12-version-control)
13. [Forbidden Patterns](#13-forbidden-patterns)
14. [Pre-Commit Checklist](#14-pre-commit-checklist)

---

## 1. Philosophy

### 1.1 What This Codebase Is

This is a **scientific simulation codebase** whose outputs may inform clinical understanding of cardiovascular physiology. The code is not a prototype. It is not a script. It is a **computational instrument** — and like any instrument used in medicine or engineering, it must be calibrated, documented, and trusted.

Code quality here is not about style preference. It is about **scientific integrity**.

### 1.2 The Four Pillars

Every file, function, and variable must serve all four of these properties simultaneously:

**Physically traceable** — every quantity must map to a real anatomical structure, a measurable clinical quantity, or a documented physical parameter. If it cannot be found in a cardiology textbook or a journal paper methods section, it does not belong without explanation.

**Numerically defensible** — every computational choice (solver type, tolerances, time step, initialization) must be made deliberately and documented. If someone challenges your numerical method, you must be able to justify it from first principles.

**Reproducible** — given the same inputs, the code must always produce the same outputs. No hidden state. No undocumented dependencies. No magic initialization values.

**Peer-reviewable** — the code should be embeddable, without embarrassment, as supplementary material in a journal submission. Write it that way from the beginning.

### 1.3 The Vibecoding Test

Before committing any code, ask:

> *Does this read like a paper in executable form — or does it read like a draft written at 2am?*

If the answer is the latter, refactor before committing.

---

## 2. Project Structure

Every project must follow a consistent layout that separates **physics**, **numerics**, **data**, and **presentation** into distinct layers. These layers must never bleed into each other.

```
project_root/
│
├── main.m                          ← Entry point only. No physics. No plotting.
├── README.md                       ← What the model does, how to run it, who wrote it
│
├── config/
│   ├── default_parameters.m        ← Reference (healthy / baseline) parameter set
│   └── patient_XXX_parameters.m   ← Per-case overrides, one file per patient/scenario
│
├── models/
│   ├── system_rhs.m                ← ODE right-hand side — the core physics
│   ├── elastance_model.m           ← Time-varying chamber mechanics
│   ├── valve_model.m               ← Valve and shunt flow logic
│   └── [other_subsystem].m        ← One file per physical subsystem
│
├── solvers/
│   └── integrate_system.m          ← Wraps numerical integration with standard options
│
├── utils/
│   ├── unit_conversion.m           ← All unit conversion factors in one place
│   ├── compute_clinical_indices.m  ← Derived quantities: EF, CO, SV, Qp/Qs, etc.
│   ├── plotting_tools.m
│   └── validation_tools.m
│
├── tests/
│   ├── test_baseline.m             ← Must pass before any patient or disease model
│   ├── test_valve_logic.m
│   └── test_[subsystem].m
│
├── results/
│   ├── figures/
│   └── tables/
│
└── docs/
    ├── theory_notes.md             ← Governing equations, assumptions, references
    ├── clinical_data_dictionary.md ← Maps clinical records to MATLAB variables
    └── VIBECODING_GUARDRAILS.md   ← This file
```

### Structure Rules

- `main.m` contains **no physics** — only loading parameters, calling solvers, and calling plotting functions
- Each `.m` file has **one responsibility** — one physical subsystem or one utility task
- `models/` contains physics only — no file I/O, no plotting, no hardcoded numerical values
- `results/` is excluded from version control for raw data, but figure-generation scripts are tracked
- Patient-specific data files are **never committed** without documented authorization

---

## 3. Naming Conventions

### 3.1 The Core Rule

Every variable name must answer three questions simultaneously:

1. **What** physical quantity is it?
2. **Where** in the anatomy does it belong?
3. **What phase or modifier** applies?

**General format:**
```
<quantity>_<anatomical_location>_<phase_or_modifier>
```

The test for any variable name: *Could a cardiologist read this in a Methods section without a legend?* If yes, it passes.

---

### 3.2 Pressure Variables

Format: `P_<location>_<phase>`

```matlab
% Instantaneous (no phase suffix needed inside the ODE)
P_lv        % Left ventricular pressure          [mmHg]
P_rv        % Right ventricular pressure         [mmHg]
P_la        % Left atrial pressure               [mmHg]
P_ra        % Right atrial pressure              [mmHg]
P_ao        % Aortic pressure                    [mmHg]
P_pa        % Pulmonary artery pressure          [mmHg]
P_pv        % Pulmonary venous pressure          [mmHg]
P_svc       % Superior vena cava pressure        [mmHg]
P_ivc       % Inferior vena cava pressure        [mmHg]

% With phase suffixes (for clinical reporting and validation)
P_lv_ed     % LV end-diastolic pressure          [mmHg]
P_lv_es     % LV end-systolic pressure           [mmHg]
P_ao_sys    % Aortic systolic pressure           [mmHg]
P_ao_dia    % Aortic diastolic pressure          [mmHg]
P_pa_sys    % Pulmonary artery systolic          [mmHg]
P_pa_dia    % Pulmonary artery diastolic         [mmHg]
P_pa_mean   % Mean pulmonary artery pressure     [mmHg]

% Pressure gradients (always: upstream minus downstream)
DeltaP_lv_ao    % LV-to-aorta gradient           [mmHg]
DeltaP_la_lv    % LA-to-LV gradient (mitral)     [mmHg]
```

**Phase suffix reference:**

| Phase | Suffix |
|---|---|
| End-diastolic | `_ed` |
| End-systolic | `_es` |
| Systolic peak | `_sys` |
| Diastolic minimum | `_dia` |
| Time-averaged mean | `_mean` |

---

### 3.3 Volume Variables

Format: `V_<chamber>_<phase>`

```matlab
% Instantaneous
V_lv        % LV volume                          [mL]
V_rv        % RV volume                          [mL]
V_la        % LA volume                          [mL]
V_ra        % RA volume                          [mL]

% Phase-specific (for clinical reporting)
V_lv_ed     % LVEDV                              [mL]
V_lv_es     % LVESV                              [mL]
V_rv_ed     % RVEDV                              [mL]
V_rv_es     % RVESV                              [mL]

% Model parameters
V0_lv       % LV unstressed (zero-pressure) volume  [mL]
V0_rv       % RV unstressed volume                  [mL]
V0_la       % LA unstressed volume                  [mL]
V0_ra       % RA unstressed volume                  [mL]
```

---

### 3.4 Flow Variables

Format: `Q_<from>_<to>`

```matlab
Q_lv_ao     % LV outflow through aortic valve        [mL/s]
Q_rv_pa     % RV outflow to pulmonary artery         [mL/s]
Q_la_lv     % Mitral valve flow (LA to LV)           [mL/s]
Q_ra_rv     % Tricuspid valve flow (RA to RV)        [mL/s]
Q_ao_sys    % Flow from aorta into systemic bed      [mL/s]
Q_pa_pul    % Flow from PA into pulmonary bed        [mL/s]

% Clinical aggregate flows (for reporting only — use L/min)
Q_systemic      % Total systemic flow (Qs)           [L/min]
Q_pulmonary     % Total pulmonary flow (Qp)          [L/min]
Q_ratio         % Qp/Qs ratio                        [dimensionless]

% Shunt flows — name by anatomical connection
% Sign convention must be documented in each function that computes shunt flows
Q_shunt_vsd     % Example: VSD shunt flow            [mL/s]
Q_shunt_asd     % Example: ASD shunt flow            [mL/s]
```

**Sign convention:** positive flow = physiologically forward direction. Any deviation must be documented in the function header.

---

### 3.5 Derived Clinical Indices

```matlab
SV_lv       % Left ventricular stroke volume         [mL]
SV_rv       % Right ventricular stroke volume        [mL]
EF_lv       % LV ejection fraction                   [fraction, NOT percent]
EF_rv       % RV ejection fraction                   [fraction]
CO          % Cardiac output                         [L/min]
CI          % Cardiac index (CO normalized by BSA)   [L/min/m²]
HR          % Heart rate                             [bpm]
SW_lv       % LV stroke work (area of PV loop)       [mmHg·mL]
```

---

### 3.6 Resistances

Format: `R_<location_or_element>`

```matlab
R_systemic      % Systemic vascular resistance (SVR)    [mmHg·s/mL]
R_pulmonary     % Pulmonary vascular resistance (PVR)   [mmHg·s/mL]
R_ao            % Aortic valve resistance               [mmHg·s/mL]
R_pv_valve      % Pulmonary valve resistance            [mmHg·s/mL]
R_mv            % Mitral valve resistance               [mmHg·s/mL]
R_tv            % Tricuspid valve resistance            [mmHg·s/mL]
R_shunt_vsd     % VSD shunt resistance (example)        [mmHg·s/mL]
```

Internal unit: **mmHg·s/mL** at all times. Convert to Wood units only for clinical reporting outputs.

---

### 3.7 Compliances

Format: `C_<vessel_or_compartment>`

```matlab
C_ao        % Aortic compliance                  [mL/mmHg]
C_pa        % Pulmonary artery compliance        [mL/mmHg]
C_sys       % Systemic venous compliance         [mL/mmHg]
C_pul       % Pulmonary venous compliance        [mL/mmHg]
```

---

### 3.8 Elastance

Format: `E_<chamber>` for instantaneous values; named parameters for model bounds

```matlab
E_lv        % LV time-varying elastance          [mmHg/mL]
E_rv        % RV time-varying elastance          [mmHg/mL]
E_la        % LA time-varying elastance          [mmHg/mL]
E_ra        % RA time-varying elastance          [mmHg/mL]

% Elastance model parameters
Emax_lv     % Peak (end-systolic) elastance      [mmHg/mL]
Emin_lv     % Minimum (diastolic) elastance      [mmHg/mL]
Emax_rv
Emin_rv
Emax_la
Emin_la
```

---

### 3.9 Defect and Shunt Geometry

When a model includes a structural defect or shunt, name its geometry explicitly:

```matlab
D_shunt_<n>     % Defect diameter                [mm]
A_shunt_<n>     % Defect cross-sectional area    [mm²]
r_shunt_<n>     % Defect radius                  [mm]
```

---

### 3.10 Demographic and Scaling Variables

| Clinical Term | MATLAB Variable | Unit |
|---|---|---|
| Age | `age_years` | years |
| Sex | `sex` | 0 = female, 1 = male |
| Height | `height_cm` | cm |
| Weight | `weight_kg` | kg |
| Body Surface Area | `BSA` | m² |

`BSA` must be computed from a documented formula with the formula name cited in a comment.

**Never use:** `AGE`, `BB`, `TB`, `BW`, `ht`, `wt`, or bare `age`, `height`, `weight` without a unit suffix.

---

### 3.11 Clinical vs. Model Namespace — Strictly Separate

Raw clinical measurements and model outputs live in separate structs. They must **never be mixed**.

```matlab
% ✅ Correct — raw measurements from clinical records
clinical.P_pa_mean      % [mmHg] — from right heart catheterization
clinical.V_lv_ed        % [mL]   — from echocardiography
clinical.EF_lv          % [fraction] — converted from reported percentage at load time

% ✅ Correct — simulated quantities from the ODE solution
model.P_pa_mean
model.V_lv_ed
model.EF_lv

% ❌ Forbidden — mixing source and simulation silently
P_pa_mean = 35;   % Is this clinical data or model output? No one knows.
```

---

## 4. Function Standards

### 4.1 Mandatory Header Template

Every `.m` function file must begin with this structure. No exceptions.

```matlab
function output = function_name(input1, input2, params)
% FUNCTION_NAME
% -----------------------------------------------------------------------
% One sentence describing what this function computes.
%
% Longer description: what physical system it represents, what
% assumptions are made, what the governing equations are, and
% what simplifications have been applied.
%
% INPUTS:
%   input1  - description                              [unit]
%   input2  - description                              [unit]
%   params  - struct of model parameters (see config/)
%
% OUTPUTS:
%   output  - description                              [unit]
%
% ASSUMPTIONS:
%   - List all non-obvious physical assumptions here
%   - e.g., "Incompressible, Newtonian flow assumed at valve orifice"
%
% SIGN CONVENTIONS (if flows or gradients are computed):
%   - Positive Q_X = forward physiological direction
%   - Document any exceptions explicitly
%
% REFERENCES:
%   [1] Author (Year). Title. Journal. Section / Equation number.
%
% AUTHOR:   [Name]
% DATE:     [YYYY-MM-DD]
% VERSION:  [X.Y]
% -----------------------------------------------------------------------
```

### 4.2 All Local Functions Need Abbreviated Headers

Any helper function longer than 10 lines inside a file needs at minimum:

```matlab
% <FUNCTION_NAME> — what it computes, inputs [units], outputs [units]
```

### 4.3 Function Length Limits

| File | Max Lines |
|---|---|
| ODE right-hand side (`system_rhs.m`) | 200 |
| Chamber mechanics / elastance | 80 |
| Valve or shunt model | 80 |
| Any utility function | 60 |
| `main.m` | 60 |

If a function exceeds its limit, split it. Add a comment explaining why the split was made and how the pieces connect.

### 4.4 One Function, One Physical Role

A function computes either:
- a time derivative (ODE RHS), or
- a single derived quantity (valve flow, elastance, clinical index), or
- a transformation or unit conversion

It does not do more than one of these.

---

## 5. Parameter Handling

### 5.1 No Global Variables — Ever

```matlab
% ❌ Absolutely forbidden
global R_systemic C_ao Emax_lv

% ✅ Always pass parameters as a struct
params.R_systemic
params.C_ao
params.Emax_lv
```

Global variables hide dependencies, make isolated testing impossible, and cause silent corruption when multiple scenarios run in the same MATLAB session.

### 5.2 Parameter Struct Must Be Organized by Physiological Subsystem

```matlab
% Ventricular mechanics
params.Emax_lv
params.Emin_lv
params.V0_lv
params.T_sys_lv     % Systolic duration                   [s]
params.Emax_rv
params.Emin_rv
params.V0_rv
params.T_sys_rv

% Atrial mechanics
params.Emax_la
params.Emin_la
params.V0_la
params.Emax_ra
params.Emin_ra
params.V0_ra

% Vascular resistances
params.R_systemic
params.R_pulmonary
params.R_ao
params.R_pv_valve
params.R_mv
params.R_tv

% Vascular compliances
params.C_ao
params.C_pa
params.C_sys
params.C_pul

% Simulation control
params.T_cardiac    % Cardiac cycle period                 [s]
params.idx          % State vector index struct (Section 7)

% Scenario / disease specific (document each field added)
% params.R_shunt_vsd = Inf;   % Inf = shunt absent; document default
```

### 5.3 Every Default Parameter Must Have a Source Comment

```matlab
% In config/default_parameters.m

params.R_systemic = 1.0;    % [mmHg·s/mL] — Source: Author (Year), Table X
params.C_ao       = 1.2;    % [mL/mmHg]   — Source: Author (Year), Eq. Y
params.Emax_lv    = 2.5;    % [mmHg/mL]   — Source: Author (Year), healthy adult
```

If you cannot cite a source, write: `% Source: assumed — needs literature validation`

---

## 6. Units Policy

### 6.1 The Project Unit System

Adopt this system for all internal computations. Never deviate silently.

| Quantity | Internal Unit | Clinical Reporting Unit |
|---|---|---|
| Pressure | mmHg | mmHg |
| Volume | mL | mL |
| Flow (ODE internal) | mL/s | L/min |
| Time | s | s |
| Resistance | mmHg·s/mL | Wood units (divide by 80 to convert) |
| Compliance | mL/mmHg | mL/mmHg |
| Elastance | mmHg/mL | mmHg/mL |
| Diameter | mm | mm |
| Area | mm² | mm² |
| BSA | m² | m² |

Conversions between internal and clinical reporting units must only happen in `utils/unit_conversion.m` or in clearly marked reporting blocks — never silently inside physics functions.

### 6.2 Unit Comments Are Mandatory on Every Variable Declaration

```matlab
P_lv    = X(idx.P_lv);       % [mmHg]
V_lv    = X(idx.V_lv);       % [mL]
Q_ao    = X(idx.Q_ao);       % [mL/s]
R_sys   = params.R_systemic; % [mmHg·s/mL]
```

### 6.3 Unit Conversions Must Be Named and Explicit

```matlab
% ❌ Forbidden — what is 60? What is 1000? What unit was Q before?
CO = Q * 60 / 1000;

% ✅ Required — explicit, traceable, self-documenting
mLs_to_Lmin  = 60 / 1000;                        % Conversion factor [L/min per mL/s]
CO_Lmin      = Q_systemic_mLs * mLs_to_Lmin;    % Cardiac output     [L/min]
```

### 6.4 Ejection Fraction Is Always Stored as a Fraction

```matlab
% Internal computation and storage — always fraction
EF_lv = (V_lv_ed - V_lv_es) / V_lv_ed;     % [fraction], range 0–1

% Convert to % only at display time
fprintf('LVEF = %.1f%%\n', EF_lv * 100);
```

### 6.5 Never Mix SI and Physiological Units Silently

If a formula requires SI internally (e.g., fluid mechanics using Pascals and m²), convert explicitly at entry and exit, with comments:

```matlab
DeltaP_Pa   = DeltaP_mmHg * 133.322;    % [mmHg] → [Pa]
A_m2        = A_mm2 * 1e-6;             % [mm²]  → [m²]
% ... compute Q in SI units ...
Q_mLs       = Q_m3s * 1e6;             % [m³/s] → [mL/s]
```

---

## 7. State Vector Management

### 7.1 Never Hardcode State Vector Indices

```matlab
% ❌ Forbidden — fragile, breaks silently when the state vector changes
P_lv = X(1);
V_lv = X(2);

% ✅ Required — always access via the index struct
P_lv = X(idx.P_lv);
V_lv = X(idx.V_lv);
```

### 7.2 Canonical Index Struct Definition

Define `idx` exactly once in `config/default_parameters.m` and pass it inside `params`. Every state variable must have an entry with a unit comment.

```matlab
% State vector layout — define once, never modify across files
% Units: pressures [mmHg], volumes [mL], flows [mL/s]

idx.P_lv    = 1;    % LV pressure              [mmHg]
idx.V_lv    = 2;    % LV volume                [mL]
idx.P_rv    = 3;    % RV pressure              [mmHg]
idx.V_rv    = 4;    % RV volume                [mL]
idx.P_la    = 5;    % LA pressure              [mmHg]
idx.V_la    = 6;    % LA volume                [mL]
idx.P_ra    = 7;    % RA pressure              [mmHg]
idx.V_ra    = 8;    % RA volume                [mL]
idx.P_ao    = 9;    % Aortic pressure          [mmHg]
idx.P_pa    = 10;   % Pulmonary artery press.  [mmHg]
idx.P_sys   = 11;   % Systemic venous press.   [mmHg]
idx.P_pul   = 12;   % Pulmonary venous press.  [mmHg]
% Scenario-specific states appended below with full comments
```

### 7.3 Scenario Extensions Must Be Documented

When extending the state vector for a disease model or scenario, add a block comment in `system_rhs.m`:

```matlab
% SCENARIO EXTENSION: [Disease / Scenario Name]
% State idx.<variable> = N represents [physical meaning] [unit]
% This state is only active when params.scenario == '[name]'
% Reference: docs/theory_notes.md, Section X
```

---

## 8. Numerical Practices

### 8.1 Solver Choice Must Be Deliberate and Documented

The ODE system may be stiff or non-stiff depending on the scenario. The solver choice must be documented in `solvers/integrate_system.m` with a written justification:

```matlab
% SOLVER: [solver name]
% JUSTIFICATION: [e.g., "System exhibits stiffness from near-discontinuous
% valve switching. Verified by monitoring step size — explicit solvers
% require excessively small steps. Implicit solver selected."]
%
% TOLERANCES:
%   RelTol = 1e-6  — sufficient for hemodynamic quantities; validated by
%                    halving and confirming < 0.01 mmHg change in P_ao_sys
%   AbsTol = 1e-8  — prevents drift accumulation in volume states
```

### 8.2 Tolerance Choices Must Be Validated

Before fixing tolerances for production use, run a convergence study:

```matlab
% Halve both RelTol and AbsTol.
% If the peak P_ao changes by less than 0.01 mmHg, original tolerances sufficient.
% Document this result as a comment in integrate_system.m.
```

### 8.3 Always Warm-Start to Periodic Steady State

Never report results from the first cardiac cycle. Always run to steady state first.

```matlab
% Acceptance criterion for steady state:
%   max change in any state peak value between consecutive cycles
%   < 0.1 mmHg (pressures) and < 0.1 mL (volumes)
%
% Document the number of warm-up cycles required for each scenario.
```

### 8.4 Valve and Shunt Logic Must Be Numerically Robust

Hard on/off switching creates near-discontinuities that can cause integration failures. Use regularized logic:

```matlab
% ❌ Fragile — hard switching collapses solver step size
if P_upstream > P_downstream
    Q = (P_upstream - P_downstream) / R;
else
    Q = 0;
end

% ✅ Robust — smooth max prevents discontinuities
Q = max(0, P_upstream - P_downstream) / R;

% For highly stiff scenarios, use soft switching (document epsilon):
% Q = (P_upstream - P_downstream) / R ...
%     * (0.5 + 0.5 * tanh((P_upstream - P_downstream) / params.epsilon_valve));
% where params.epsilon_valve [mmHg] is documented in config/
```

### 8.5 Initial Conditions Must Be Physiologically Motivated

Never initialize with `X0 = zeros(n, 1)`. The system will start from a non-physiological state.

```matlab
% In config/default_parameters.m — document every initial condition
X0(idx.P_lv)  = 8;      % [mmHg] — approximate LV end-diastolic pressure
X0(idx.V_lv)  = 120;    % [mL]   — approximate LVEDV, healthy adult male
X0(idx.P_ao)  = 80;     % [mmHg] — approximate diastolic aortic pressure
% Cite source for each initial value
```

### 8.6 Mass Conservation Must Be Checked at Steady State

At periodic steady state, net volume accumulation across all compartments must be zero within numerical tolerance. Assert this in `tests/test_baseline.m`.

---

## 9. Clinical Data Handling

### 9.1 A Data Dictionary Is Mandatory for Every Patient or Dataset

For every patient or clinical dataset used, create or update `docs/clinical_data_dictionary.md`:

```markdown
| Clinical Record Field | MATLAB Variable     | Unit     | Method           | Reliability | Date    |
|-----------------------|---------------------|----------|------------------|-------------|---------|
| LVEDP (cath)         | clinical.P_lv_ed   | mmHg     | Left heart cath  | High        | YYYY-MM |
| LVEDV (echo)         | clinical.V_lv_ed   | mL       | 2D echo biplane  | Moderate    | YYYY-MM |
| LVEF (echo)          | clinical.EF_lv     | fraction | Derived          | Moderate    | YYYY-MM |
| mPAP (RHC)           | clinical.P_pa_mean | mmHg     | Right heart cath | High        | YYYY-MM |
```

Reliability must be one of: `High`, `Moderate`, `Low`, or `Derived` (not directly measured).

### 9.2 Patient Data Is Never Hardcoded in `.m` Files

```matlab
% ❌ Forbidden — patient value hardcoded in source code
clinical.P_pa_mean = 42;

% ✅ Required — loaded from external, access-controlled file
data     = load('data/patient_01_clinical.mat');
clinical = data.clinical;
```

### 9.3 Annotate Measurement Uncertainty Where Known

```matlab
clinical.V_lv_ed             = 145;        % [mL]
clinical.V_lv_ed_uncertainty = 0.15;       % ±15% typical for 2D echo biplane
clinical.V_lv_ed_method      = '2D echo';  % Note: 3D echo or MRI is more accurate
```

### 9.4 Convert Clinical Record Units Explicitly at Load Time

Clinical records often report values in conventions that differ from the model's internal system. Convert once, at load time, and document it:

```matlab
clinical.EF_lv = EF_lv_percent / 100;     % [%] → [fraction] at load time
clinical.CO    = CO_Lmin;                  % [L/min] — no conversion needed
```

---

## 10. Validation Requirements

### 10.1 The Baseline Must Pass Before Any Disease or Patient Model

Run `tests/test_baseline.m` before developing or running any patient-specific model. If the baseline fails, nothing downstream can be trusted.

Baseline outputs must fall within documented physiological reference ranges, with citations:

```matlab
% Reference ranges for a healthy adult at rest (cite your source):
assert(P_ao_sys   >= 100  && P_ao_sys  <= 140,  'Aortic systolic OOR');
assert(P_ao_dia   >= 60   && P_ao_dia  <= 90,   'Aortic diastolic OOR');
assert(P_pa_mean  >= 10   && P_pa_mean <= 20,   'Mean PAP OOR');
assert(EF_lv      >= 0.55 && EF_lv    <= 0.75,  'LVEF OOR');
assert(CO         >= 4.0  && CO       <= 8.0,   'Cardiac output OOR');
assert(abs(Q_ratio - 1.0) < 0.05,               'Qp/Qs not unity at baseline');
% OOR = Out Of Range
```

### 10.2 Pressure-Volume Loop Validation

After any simulation producing ventricular pressure and volume traces, generate and inspect the P-V loop. Confirm:

- Shape is approximately rectangular
- Valve opening and closing occur at physiologically expected pressures
- End-systolic pressure-volume relationship slope is consistent with Emax
- Stroke work (loop area) is consistent with cardiac output and heart rate

### 10.3 Mass Conservation Check

At periodic steady state, net fluid accumulation across all compartments must be zero. Assert this programmatically in `test_baseline.m`:

```matlab
% Integrate dV/dt for all volume states over one complete cardiac cycle
% Result must be < numerical tolerance (e.g., 1e-3 mL)
```

### 10.4 Sensitivity to Initial Conditions

Perturb X0 by ±10% and confirm the model converges to the same limit cycle within a bounded number of cycles. Document the number of cycles required for each scenario.

---

## 11. Plotting Standards

### 11.1 All Figures Must Be Publication-Ready From Day One

Do not maintain a "quick" version and a "paper" version. Write publication-ready plotting code from the start:

```matlab
figure('Units', 'centimeters', 'Position', [0 0 16 12]);

plot(t_sol, P_lv, 'b-', 'LineWidth', 1.5, 'DisplayName', 'P_{LV}');
hold on;
plot(t_sol, P_ao, 'r-', 'LineWidth', 1.5, 'DisplayName', 'P_{Ao}');

xlabel('Time [s]',           'FontSize', 11, 'FontName', 'Arial');
ylabel('Pressure [mmHg]',    'FontSize', 11, 'FontName', 'Arial');
title('[Descriptive title]', 'FontSize', 12, 'FontName', 'Arial');
legend('FontSize', 10, 'Location', 'best');
grid on;
set(gca, 'FontSize', 10, 'FontName', 'Arial', 'Box', 'on');
```

### 11.2 Clinical Overlays Must Be Visually Distinct

```matlab
% Clinical measurement — marker, dashed line
plot(t_ref, clinical.P_pa_mean * ones(size(t_ref)), 'ks--', ...
    'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Clinical');

% Model output — solid line
plot(t_sol, model.P_pa_mean_trace, 'b-', ...
    'LineWidth', 1.5, 'DisplayName', 'Model');
```

### 11.3 Always Export as Vector Format

```matlab
exportgraphics(gcf, 'results/figures/figure_name.pdf', ...
    'ContentType', 'vector', 'Resolution', 300);
```

Never commit `.png` or `.fig` files as final outputs.

### 11.4 Figure File Names Must Be Descriptive

```
% ❌ Forbidden
figure1.pdf
plot_final.pdf
test_fig.pdf

% ✅ Required — scenario, content, patient context in the name
PV_loop_LV_baseline_healthy.pdf
pressure_traces_patient01.pdf
Qp_Qs_sensitivity_defect_radius.pdf
```

---

## 12. Version Control

### 12.1 Commit Messages Must Reference Physiology or Numerics

Every commit message must tell a reader **what changed physically or numerically**:

```bash
# ❌ Bad — no physical meaning conveyed
git commit -m "fix bug"
git commit -m "updated model"
git commit -m "changes"

# ✅ Good — physically traceable
git commit -m "Add time-varying atrial elastance model for LA and RA"
git commit -m "Fix sign error in pulmonary venous return flow"
git commit -m "Tighten AbsTol to 1e-8 to resolve drift in volume states"
git commit -m "Separate healthy baseline from patient config files"
git commit -m "Add test_baseline.m with physiological range assertions"
```

### 12.2 Tag Stable Milestones

Before sharing with a collaborator, supervisor, or submitting to a journal, tag the exact commit:

```bash
git tag -a v1.0-baseline-validated   -m "Baseline passes all physiological range checks"
git tag -a v1.1-patient01-calibrated -m "Patient 01 matches clinical data within stated RMSE"
```

### 12.3 What Must Never Be Committed

- Patient clinical data files — even anonymized — without documented authorization
- `.m` files with hardcoded patient-specific numerical values
- Default parameters without source comments
- Large binary result files (`.mat` > 10 MB) — store externally, document path in `README.md`
- Broken code — never commit a repository state that does not run

---

## 13. Forbidden Patterns

These patterns are banned. Any code review finding them will require revision before merge.

### 13.1 Ambiguous Variable Names

```matlab
% ❌ All forbidden
x1, x2, x3
data, Data
temp, tmp
val, val2
res, comp
p1, p2, p3
q1, q_in, q_out
flow1, flow_left
test, test2
```

### 13.2 Global and Persistent State

```matlab
% ❌ Forbidden
global params
global R_systemic
persistent prev_state   % Hides bugs, destroys reproducibility
```

### 13.3 Magic Numbers

```matlab
% ❌ Forbidden — what is 133.322? What is 1060? What is 0.7?
Q = 0.7 * A * sqrt(2 * 1060 * abs(dP * 133.322));

% ✅ Required — named, commented, sourced
Cd          = 0.7;        % Discharge coefficient [dimensionless] — Source: [ref]
rho_blood   = 1060;       % Blood density [kg/m^3]               — Source: [ref]
mmHg_to_Pa  = 133.322;    % Unit conversion factor [Pa/mmHg]
Q           = Cd * A * sqrt(2 * rho_blood * abs(dP * mmHg_to_Pa));
```

### 13.4 Physics Inside main.m

```matlab
% ❌ Forbidden — main.m must never compute physiological quantities
P_lv = E_lv * (V_lv - V0_lv);    % This belongs in models/, not main.m
```

### 13.5 Silent Unit Mixing

```matlab
% ❌ Forbidden — Pa or mmHg? m² or mm²?
Q = Cd * A * sqrt(2 * rho * abs(dP));   % Units are ambiguous
```

### 13.6 Undocumented Solver Choices

```matlab
% ❌ Forbidden — no justification, no tolerances documented
[t, X] = ode45(@rhs, tspan, X0);
```

### 13.7 Mixing clinical and model Structs

```matlab
% ❌ Forbidden — overwrites model output with clinical data silently
model.EF_lv = clinical.EF_lv;
```

---

## 14. Pre-Commit Checklist

Run through every item before committing. If any box cannot be checked, do not commit.

### Physics and Naming
- [ ] Every new variable follows the `<quantity>_<location>_<phase>` convention
- [ ] Every variable has a unit comment on the same line
- [ ] All pressure gradients are defined as upstream minus downstream
- [ ] Sign conventions for all flows and shunts are documented in function headers
- [ ] No magic numbers — all constants are named, commented, and sourced

### Code Structure
- [ ] No physics in `main.m`
- [ ] No global or persistent variables
- [ ] Every function has a complete header (Section 4.1)
- [ ] All functions are within length limits (Section 4.3)
- [ ] One physical role per function

### Parameter Handling
- [ ] All parameters passed as a struct — no globals
- [ ] Every default parameter has a source comment in `config/`
- [ ] Patient-specific values are loaded from external files — not hardcoded

### Units
- [ ] Unit comments present on every variable declaration
- [ ] All unit conversions are named and explicit
- [ ] EF stored as fraction (0–1), not percentage
- [ ] No silent mixing of SI and physiological units

### Numerics
- [ ] Solver choice is documented and justified in `solvers/`
- [ ] Tolerances are documented and validated by convergence study
- [ ] Warm-start to steady state confirmed before results are extracted
- [ ] Initial conditions are physiologically motivated and documented
- [ ] Valve and shunt logic uses smooth or regularized switching

### Validation
- [ ] `test_baseline.m` passes before any scenario-specific code is run
- [ ] P-V loop shape has been visually inspected
- [ ] Mass conservation check passes at steady state

### Clinical Data (if applicable)
- [ ] Data dictionary updated in `docs/clinical_data_dictionary.md`
- [ ] No patient values hardcoded in `.m` files
- [ ] Measurement uncertainty annotated where known
- [ ] Unit conversions from clinical records performed at load time

### Version Control
- [ ] Commit message describes what changed physically or numerically
- [ ] No patient data committed without authorization
- [ ] Repository runs without errors in its current state

---

*This document is a living standard. When a new failure mode is discovered, a new forbidden pattern identified, or a new best practice established — add it here with a date and author name.*

*Last updated: [YYYY-MM-DD] — [Author]*
