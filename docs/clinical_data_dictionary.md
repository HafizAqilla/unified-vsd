# Clinical Data Dictionary

Maps clinical record terminology to MATLAB variable names in the unified VSD model.

All variables live inside the `clinical` struct returned by `patient_template.m`.
**Raw clinical data and model outputs are strictly separate namespaces** —
never assign clinical values directly to model output variables.

---

## patient_template.m sub-sections

```matlab
clinical.common       % demographics; applies to both scenarios
clinical.pre_surgery  % measurements before VSD closure
clinical.post_surgery % measurements after VSD closure
```

---

## COMMON  —  `clinical.common`

| MATLAB field | Clinical term | Source | Unit | Reliability | Date |
|---|---|---|---|---|---|
| `age_years` | Patient age | Medical record | years | High | YYYY-MM |
| `weight_kg` | Body weight | Anthropometry | kg | High | YYYY-MM |
| `height_cm` | Standing height | Anthropometry | cm | High | YYYY-MM |
| `sex` | Biological sex | Medical record | 'M' / 'F' | High | YYYY-MM |
| `BSA` | Body surface area | Mosteller formula | m² | Derived | YYYY-MM |
| `HR` | Resting heart rate | ECG / pulse ox | bpm | High | YYYY-MM |

---

## PRE-SURGERY  —  `clinical.pre_surgery`

### Shunt geometry

| MATLAB field | Clinical term | Source | Unit | Reliability | Date |
|---|---|---|---|---|---|
| `VSD_diameter_mm` | VSD defect diameter | 2D echo | mm | Moderate | YYYY-MM |
| `VSD_gradient_mmHg` | Peak instantaneous VSD pressure gradient | CW Doppler | mmHg | High | YYYY-MM |
| `Q_shunt_Lmin` | Net shunt flow (Qp − Qs) | Fick / thermodilution | L/min | Moderate | YYYY-MM |
| `QpQs` | Pulmonary-to-systemic flow ratio | Oximetry / Fick | dimensionless | Moderate | YYYY-MM |

### Pulmonary haemodynamics

| MATLAB field | Clinical term | Source | Unit | Reliability | Date |
|---|---|---|---|---|---|
| `PAP_sys_mmHg` | Systolic pulmonary artery pressure | RHC | mmHg | High | YYYY-MM |
| `PAP_dia_mmHg` | Diastolic pulmonary artery pressure | RHC | mmHg | High | YYYY-MM |
| `PAP_mean_mmHg` | Mean pulmonary artery pressure (mPAP) | RHC | mmHg | High | YYYY-MM |
| `PVR_WU` | Pulmonary vascular resistance | RHC / Fick | Wood units | Moderate | YYYY-MM |

### Systemic haemodynamics

| MATLAB field | Clinical term | Source | Unit | Reliability | Date |
|---|---|---|---|---|---|
| `SAP_sys_mmHg` | Systolic arterial pressure | Arterial line / cuff | mmHg | High | YYYY-MM |
| `SAP_dia_mmHg` | Diastolic arterial pressure | Arterial line / cuff | mmHg | High | YYYY-MM |
| `SAP_mean_mmHg` | Mean arterial pressure | Arterial line | mmHg | High | YYYY-MM |
| `SVR_WU` | Systemic vascular resistance | Fick calculation | Wood units | Moderate | YYYY-MM |

### Atrial pressures

| MATLAB field | Clinical term | Source | Unit | Reliability | Date |
|---|---|---|---|---|---|
| `RAP_mean_mmHg` | Right atrial pressure (mean) | RHC | mmHg | High | YYYY-MM |
| `LAP_mean_mmHg` | Left atrial pressure (mean) OR PCWP | RHC | mmHg | Moderate | YYYY-MM |

### Ventricular volumes and function

| MATLAB field | Clinical term | Source | Unit | Reliability | Date |
|---|---|---|---|---|---|
| `LVEDV_mL` | Left ventricular end-diastolic volume | Echo / MRI | mL | Moderate | YYYY-MM |
| `LVESV_mL` | Left ventricular end-systolic volume | Echo / MRI | mL | Moderate | YYYY-MM |
| `RVEDV_mL` | Right ventricular end-diastolic volume | Echo / MRI | mL | Moderate | YYYY-MM |
| `RVESV_mL` | Right ventricular end-systolic volume | Echo / MRI | mL | Moderate | YYYY-MM |
| `LVEF` | LV ejection fraction | Echo / MRI | fraction (not %) | Moderate | YYYY-MM |
| `CO_Lmin` | Systemic cardiac output | Fick / thermodilution | L/min | Moderate | YYYY-MM |

---

## POST-SURGERY  —  `clinical.post_surgery`

| MATLAB field | Clinical term | Source | Unit | Reliability | Date |
|---|---|---|---|---|---|
| `QpQs` | Residual shunt ratio (expected ~1.0) | Echo / oximetry | dimensionless | Moderate | YYYY-MM |
| `PAP_sys_mmHg` | Systolic PA pressure (post-op) | Echo / RHC | mmHg | High | YYYY-MM |
| `PAP_dia_mmHg` | Diastolic PA pressure (post-op) | Echo / RHC | mmHg | High | YYYY-MM |
| `PAP_mean_mmHg` | Mean PA pressure (post-op) | Echo / RHC | mmHg | High | YYYY-MM |
| `PVR_WU` | PVR (post-op, should be lower) | Fick | Wood units | Moderate | YYYY-MM |
| `SAP_sys_mmHg` | Systolic arterial pressure | Cuff / arterial line | mmHg | High | YYYY-MM |
| `SAP_dia_mmHg` | Diastolic arterial pressure | Cuff / arterial line | mmHg | High | YYYY-MM |
| `MAP_mmHg` | Mean arterial pressure | Arterial line | mmHg | High | YYYY-MM |
| `SVR_WU` | SVR (post-op) | Fick | Wood units | Moderate | YYYY-MM |
| `RAP_mean_mmHg` | Right atrial pressure | RHC | mmHg | High | YYYY-MM |
| `LAP_mean_mmHg` | Left atrial pressure | RHC / PCWP | mmHg | Moderate | YYYY-MM |
| `LVEDV_mL` | LVEDV (expected smaller post-op) | Echo / MRI | mL | Moderate | YYYY-MM |
| `LVESV_mL` | LVESV | Echo / MRI | mL | Moderate | YYYY-MM |
| `RVEDV_mL` | RVEDV (volume unloaded post-op) | Echo / MRI | mL | Moderate | YYYY-MM |
| `RVESV_mL` | RVESV | Echo / MRI | mL | Moderate | YYYY-MM |
| `LVEF` | LVEF | Echo / MRI | fraction | Moderate | YYYY-MM |
| `RVEF` | RVEF | Echo / MRI | fraction | Moderate | YYYY-MM |
| `CO_Lmin` | Systemic cardiac output | Fick / thermodilution | L/min | Moderate | YYYY-MM |

> **How to fill in the Date column:** Enter the date of the clinical measurement or
> record in YYYY-MM format. Write `N/A` for derived quantities with no direct measurement date.
>
> **Reliability levels (AGENTS.md §9.1):**
> - `High` — direct invasive measurement (RHC, arterial line)
> - `Moderate` — non-invasive imaging or indirect calculation (echo, Fick)
> - `Low` — indirect estimate or single measurement without repeat
> - `Derived` — computed from other measured fields (BSA, EF from volumes)

---

## Units conversion reference

| From | To | Factor |
|---|---|---|
| Wood units (WU) | mmHg·s/mL | × 0.06 |
| mmHg·s/mL | Wood units | ÷ 0.06 |
| L/min | mL/s | × 1000/60 |
| mL/s | L/min | × 60/1000 |
| EF (%) | EF (fraction) | ÷ 100 |

---

## Source abbreviations

| Abbreviation | Meaning |
|---|---|
| RHC | Right heart catheterisation |
| CW | Continuous-wave |
| Echo | Two-dimensional transthoracic or transoesophageal echocardiography |
| MRI | Cardiac magnetic resonance imaging |
| Fick | Fick oxygen method for cardiac output |
| TD | Thermodilution |
| PCWP | Pulmonary capillary wedge pressure (surrogate for LAP) |
