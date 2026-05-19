# Reyna Targeted Grid Polish - 2026-05-15

## Purpose

This note documents the targeted Reyna pre-surgery polish run that started
from the Lundquist-scaled best-candidate package:

`results/runs/20260514_155253_reyna_pre_surgery/mat/params_best_candidate_pre_surgery.mat`

The goal was to reduce the governed primary RMSE below 0.10 while preserving
physiological validity and avoiding hidden changes to clinical target
governance.

## Method

Two local search passes were added:

- `scripts/scan_reyna_targeted_sensitivity.m` performs one-at-a-time
  perturbations around the saved best candidate.
- `scripts/run_reyna_targeted_grid_polish.m` performs a deterministic local
  grid over the strongest residual directions: `C.SAR`, `E.LV.EB`, `V0.LV`,
  and `E.LV.EA`.

The grid selector now includes gate awareness, so it favors candidates with
fewer primary targets above 10 percent absolute error instead of choosing only
the lowest RMSE.

## Accepted Registry-Bounded Candidate

Run folder:

`results/runs/20260515_173742_reyna_targeted_grid_polish`

Summary:

| Quantity | Start | Best |
|---|---:|---:|
| Governed primary RMSE | 0.106544 | 0.088354 |
| Full transparent RMSE | 0.123790 | 0.109307 |
| Hard-only RMSE | 0.093174 | 0.082640 |
| Soft-only RMSE | 0.093367 | 0.089138 |
| Primary targets above 10 percent | 6 | 3 |

Best parameter movement:

| Parameter | Start | Best | Ratio |
|---|---:|---:|---:|
| `C.SAR` | 0.631064 | 0.501710 | 0.795022 |
| `E.LV.EB` | 0.159684 | 0.178846 | 1.120000 |
| `V0.LV` | 6.510851 | 7.137572 | 1.096258 |
| `E.LV.EA` | 7.422022 | 7.125141 | 0.960000 |

The direct pressure-flow block is within the 10 percent gate:

| Metric | Error |
|---|---:|
| `RAP_mean` | -0.441 percent |
| `PAP_mean` | +0.356 percent |
| `SAP_min` | -7.435 percent |
| `SAP_max` | -9.440 percent |
| `SAP_mean` | -1.344 percent |
| `QpQs` | -5.342 percent |
| `SVR` | +9.473 percent |
| `CO_Lmin` | -9.948 percent |

Remaining misses:

| Metric | Error | Interpretation |
|---|---:|---|
| `LAP_mean` | -12.274 percent | Estimated target; not directly fitted |
| `LVEDV` | +14.653 percent | Conflicts with flow-derived stroke volume |
| `LVEF` | +13.909 percent | Coupled to LVEDV/LVESV tradeoff |
| `RVEDV` | +25.649 percent | Consistency-check only, excluded from primary RMSE |

## Exploratory Relaxed-Bound Probe

Run folder:

`results/runs/20260515_173544_reyna_targeted_grid_polish`

This probe relaxed selected registry bounds for `C.SAR`, `V0.LV`, and
`E.LV.EA`. It improved governed primary RMSE to 0.086275 and made `LAP_mean`
pass the 10 percent gate, but `CO_Lmin` moved outside the gate at -10.439
percent. This candidate should remain exploratory unless the relaxed bounds are
clinically justified.

## Clinical Consistency Finding

The consistency audit exported:

`results/tables/clinical_co_consistency_reyna_pre_surgery_20260515_172856.csv`

Key values:

| Quantity | Value |
|---|---:|
| Fick `Qs` target | 3.423 L/min |
| `Qp` from `Qs * QpQs` | 4.087 L/min |
| LV output from echo LVSV and HR | 2.582 L/min |
| LV output / Fick-implied `Qp` | 0.632 |

This is why all-target closure is difficult: the catheter-derived
pressure-flow block and echo-derived ventricular volume block are not mutually
consistent as exact point targets.

## Literature Framing

- Zhang et al. support allometric age/size scaling for infant, child, and
  adolescent cardiovascular models using parameter-type scaling exponents:
  https://pubmed.ncbi.nlm.nih.gov/31005012/
- Lundquist et al. support patient-specific 0D scaling using age, weight,
  height, and sex, and report realistic physiology across body sizes and ages:
  https://researchinformation.umcutrecht.nl/en/publications/patient-specific-size-and-age-scaling-in-a-zero-dimensional-cardi/
- Pediatric ventricular volume measurements are modality-sensitive; 3D echo
  can be feasible/reproducible, but volume agreement should still be treated
  with explicit uncertainty when it conflicts with catheter flow:
  https://pubmed.ncbi.nlm.nih.gov/22150697/

## Recommendation

Use the registry-bounded candidate as the current best trustworthy simulation
state. The relaxed-bound result is useful evidence that the remaining error is
not purely numerical, but it should not replace the bounded candidate without a
documented physiological reason to relax the parameter registry.
