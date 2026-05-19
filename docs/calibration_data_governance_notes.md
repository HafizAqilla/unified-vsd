# Calibration Data Governance Notes

Date: 2026-05-14

## Purpose

This note documents the calibration-governance change introduced for the Reyna pre-surgery VSD model. The goal is not to hide difficult targets, but to prevent the optimizer from treating internally inconsistent mixed-modality measurements as equally hard physiological constraints.

## Clinical Consistency Rule

Before calibration, the workflow now computes:

| Quantity | Formula | Interpretation |
|---|---:|---|
| `SV_Qs` | `CO_Lmin * 1000 / HR` | systemic stroke volume implied by catheter/Fick systemic flow |
| `SV_Qp` | `CO_Lmin * QpQs * 1000 / HR` | pulmonary stroke volume implied by catheter/Fick flow ratio |
| `SV_LV` | `LVEDV - LVESV` | echo-derived LV stroke volume |
| `SV_RV` | `RVEDV - RVESV` | echo-derived RV stroke volume |

The thresholds are explicit in `audit_clinical_consistency.m`:

| Severity | Relative difference |
|---|---:|
| mild | > 20% |
| strong | > 35% |
| critical | > 50% |

For Reyna, the audit reports:

| Quantity | Value |
|---|---:|
| `SV_Qs` | 28.765 mL/beat |
| `SV_Qp` | 34.345 mL/beat |
| `SV_LV` | 21.700 mL/beat |
| `SV_RV` | 18.500 mL/beat |
| `SV_Qp` vs `SV_RV` relative difference | 59.968% |
| Audit severity | critical |

## RVEDV Decision

The PDF record shows `RVEDV = 30.5 mL` and `RVESV = 12 mL` as echocardiography-derived absolute volumes, not BSA-indexed values. However, those values imply `SV_RV = 18.5 mL`, while the catheter-derived pulmonary flow implies `SV_Qp = 34.345 mL/beat`.

Therefore, `RVEDV` is now:

| Role | Setting |
|---|---|
| Shown in validation table | yes |
| Included in full RMSE | yes |
| Included in calibration loss | no |
| Included in primary governed RMSE | no |
| Tier | `consistency_check_only` |
| Flag | `inconsistent_or_unverified` |

This supports the thesis/paper statement:

> The model identified inconsistency between catheter-derived pulmonary/systemic flow estimates and echo-derived RV volume targets. Therefore, RVEDV was retained for transparent validation reporting but was not used as a hard fitting target. Primary calibration metrics exclude the flagged inconsistent target, while full RMSE is reported for completeness.

## Current Target Tiers

Hard targets:

`CO_Lmin`, `QpQs`, `PAP_mean`, `SAP_mean`, `RAP_mean`, `LVEDV`, `LVESV`, `LVEF`

Soft targets:

`SAP_max`, `SAP_min`, `SVR`, `RVESV`

Consistency-check only:

`RVEDV`

Validation-only:

`LAP_mean` remains reported but is not used for calibration because it is estimated rather than directly measured.

## Literature Support

The governance logic is based on conservation identities and measurement uncertainty, not on claiming that one modality is automatically wrong.

Relevant support:

- Zhang et al. introduced allometric scaling laws for cardiovascular parameters across infants, children, and adolescents in a multiscale cardiovascular model, supporting the general need for age/size-aware pediatric parameter scaling: [PubMed 31005012](https://pubmed.ncbi.nlm.nih.gov/31005012/).
- Lundquist et al. describe patient-specific size and age scaling for a 0D cardiovascular model using anthropometrics, supporting the active comparison of size-scaling methods rather than tuning from an adult baseline alone: [University of Twente record](https://research.utwente.nl/en/publications/patient-specific-size-and-age-scaling-in-a-zero-dimensional-cardi/).
- Pediatric RV volume measurement by 3D echo can correlate well with MRI, but RV geometry makes 2D estimates less robust; this supports treating RV volume targets as modality-sensitive rather than blindly hard constraints: [PubMed 17628408](https://pubmed.ncbi.nlm.nih.gov/17628408/) and [PubMed 21460148](https://pubmed.ncbi.nlm.nih.gov/21460148/).

## Implementation Files

- `src/validation/audit_clinical_consistency.m`
- `src/calibration/build_target_tiers.m`
- `src/calibration/build_case_calibration_profile.m`
- `src/calibration/objective_calibration.m`
- `src/utils/validation_report.m`
- `src/utils/params_from_clinical.m`
- `src/utils/build_initial_conditions.m`
- `tests/test_clinical_consistency_target_tiers.m`
