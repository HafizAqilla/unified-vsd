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

| MATLAB field | Clinical term | Source | Unit |
|---|---|---|---|
| `age_years` | Patient age | Medical record | years |
| `weight_kg` | Body weight | Anthropometry | kg |
| `height_cm` | Standing height | Anthropometry | cm |
| `sex` | Biological sex | Medical record | 'M' / 'F' |
| `BSA` | Body surface area | Mosteller formula | m² |
| `HR` | Resting heart rate | ECG / pulse ox | bpm |

---

## PRE-SURGERY  —  `clinical.pre_surgery`

### Shunt geometry

| MATLAB field | Clinical term | Source | Unit |
|---|---|---|---|
| `VSD_diameter_mm` | VSD defect diameter | 2D echo | mm |
| `VSD_gradient_mmHg` | Peak instantaneous VSD pressure gradient | CW Doppler | mmHg |
| `Q_shunt_Lmin` | Net shunt flow (Qp − Qs) | Fick / thermodilution | L/min |
| `QpQs` | Pulmonary-to-systemic flow ratio | Oximetry / Fick | dimensionless |

### Pulmonary haemodynamics

| MATLAB field | Clinical term | Source | Unit |
|---|---|---|---|
| `PAP_sys_mmHg` | Systolic pulmonary artery pressure | RHC | mmHg |
| `PAP_dia_mmHg` | Diastolic pulmonary artery pressure | RHC | mmHg |
| `PAP_mean_mmHg` | Mean pulmonary artery pressure (mPAP) | RHC | mmHg |
| `PVR_WU` | Pulmonary vascular resistance | RHC / Fick | Wood units |

### Systemic haemodynamics

| MATLAB field | Clinical term | Source | Unit |
|---|---|---|---|
| `SAP_sys_mmHg` | Systolic arterial pressure | Arterial line / cuff | mmHg |
| `SAP_dia_mmHg` | Diastolic arterial pressure | Arterial line / cuff | mmHg |
| `SAP_mean_mmHg` | Mean arterial pressure | Arterial line | mmHg |
| `SVR_WU` | Systemic vascular resistance | Fick calculation | Wood units |

### Atrial pressures

| MATLAB field | Clinical term | Source | Unit |
|---|---|---|---|
| `RAP_mean_mmHg` | Right atrial pressure (mean) | RHC | mmHg |
| `LAP_mean_mmHg` | Left atrial pressure (mean) OR PCWP | RHC | mmHg |

### Ventricular volumes and function

| MATLAB field | Clinical term | Source | Unit |
|---|---|---|---|
| `LVEDV_mL` | Left ventricular end-diastolic volume | Echo / MRI | mL |
| `LVESV_mL` | Left ventricular end-systolic volume | Echo / MRI | mL |
| `RVEDV_mL` | Right ventricular end-diastolic volume | Echo / MRI | mL |
| `RVESV_mL` | Right ventricular end-systolic volume | Echo / MRI | mL |
| `LVEF` | LV ejection fraction | Echo / MRI | fraction (not %) |
| `CO_Lmin` | Systemic cardiac output | Fick / thermodilution | L/min |

---

## POST-SURGERY  —  `clinical.post_surgery`

| MATLAB field | Clinical term | Source | Unit |
|---|---|---|---|
| `QpQs` | Residual shunt ratio (expected ~1.0) | Echo / oximetry | dimensionless |
| `PAP_sys_mmHg` | Systolic PA pressure (post-op) | Echo / RHC | mmHg |
| `PAP_dia_mmHg` | Diastolic PA pressure (post-op) | Echo / RHC | mmHg |
| `PAP_mean_mmHg` | Mean PA pressure (post-op) | Echo / RHC | mmHg |
| `PVR_WU` | PVR (post-op, should be lower) | Fick | Wood units |
| `SAP_sys_mmHg` | Systolic arterial pressure | Cuff / arterial line | mmHg |
| `SAP_dia_mmHg` | Diastolic arterial pressure | Cuff / arterial line | mmHg |
| `MAP_mmHg` | Mean arterial pressure | Arterial line | mmHg |
| `SVR_WU` | SVR (post-op) | Fick | Wood units |
| `RAP_mean_mmHg` | Right atrial pressure | RHC | mmHg |
| `LAP_mean_mmHg` | Left atrial pressure | RHC / PCWP | mmHg |
| `LVEDV_mL` | LVEDV (expected smaller post-op) | Echo / MRI | mL |
| `LVESV_mL` | LVESV | Echo / MRI | mL |
| `RVEDV_mL` | RVEDV (volume unloaded post-op) | Echo / MRI | mL |
| `RVESV_mL` | RVESV | Echo / MRI | mL |
| `LVEF` | LVEF | Echo / MRI | fraction |
| `RVEF` | RVEF | Echo / MRI | fraction |
| `CO_Lmin` | Systemic cardiac output | Fick / thermodilution | L/min |

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
