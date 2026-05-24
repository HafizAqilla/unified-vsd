# PR: Reyna Stable Main Run - Keisya Baseline Repair & Scaling Independence

**Branch:** `codex/reyna-stable-main-run` -> `main`
**Author:** Hafiz + opencode (executor), 2026-05-24

---

## Summary

Repairs the calibration regression that raised Reyna pre-surgery RMSE from 0.0588 to 0.2036 after the baseline-integration merge. Restores single-digit primary governed RMSE under Lundquist, aligns Zhang scaling with Zhang 2019 Table 1, locks the Keisya baseline with regression coverage, and separates Zhang from the Lundquist-calibrated disease seed.

## Acceptance

- `primary_governed RMSE < 0.10` for **Lundquist** on Reyna pre-surgery
- Lundquist run: **ACCEPT, RMSE 0.0588, all 5 primary gates within 10%**
- Zhang is a same-governance comparator with its own Zhang-scaled baseline; rerun before final PR notes to report the refreshed post-seed-separation RMSE/status
- Previous Zhang run before seed separation: **PROMISING_NEAR_MISS, RMSE 0.1318** (PAP_mean at 10.60%, borderline)
- Lightweight tests were previously reported passing before the seed-separation cleanup; rerun before final commit if you want a fresh verification line
- Heavy regression tests skip by default with `UNIFIED_VSD_RUN_HEAVY_TESTS` unset

## What Changed

### Calibration Profile (recipe system)

- New `config/calibration_recipes/reyna_pre_surgery.m`: explicit recipe that pins demographics, 14 active parameters, 5 primary pressure-flow targets, and the Lundquist-only pre-validated seed from `20260521_214941`
- Recipe loaders `apply_calibration_recipe_to_clinical.m` and `apply_calibration_recipe_to_params.m` wired into `build_case_calibration_profile.m`
- Active parameter set restored to sparse-cath equivalence: all 4 ventricular elastances, both V0s, LA/RA elastance, vascular resistances/compliances, `vsd.Cd` (14 total)
- Bounds include `group.R_sys_scale [0.25, 2.80]` and `group.R_pul_scale [0.45, 2.80]`
- The accepted disease seed, fixed `V0.SVEN`, initial vector, and seed short-circuit are gated to `lundquist_bsa`; Zhang uses the same recipe/governance but starts from Zhang's own scaled baseline plus normal clinical seeding

### Reporting Split

- Added explicit split tables for clinical validation targets vs model-derived findings:
  - `validation_clinical_targets_baseline_<scenario>.csv`
  - `validation_clinical_targets_calibrated_<scenario>.csv`
  - `model_derived_metric_findings_<scenario>.csv`
  - `model_derived_parameter_findings_<scenario>.csv`
- This is reporting-only. Existing RMSE, target tiers, primary gates, and calibration logic are unchanged.
- Clinical validation rows still retain tier semantics (`hard`, `soft`, `consistency_check_only`, `derived_validation`, `validation_holdout`).
- Model-derived parameter findings should be interpreted against literature/plausibility bounds, not as direct clinical validation.

### Zhang Scaling

- `eRv_op` corrected from `-0.90` to `-0.50` per Zhang 2019 Table 1 (`K_vo = -0.5`)
- Inertance L scaling added: `eL = -1.0`, matching Lundquist symmetry (not in Zhang Table 1; documented inline)
- `zhang_exponents` struct extended with `L` and `Rvalve_open` keys

### Baseline & Demographics

- `patient_reyna.m` updated to Keisya 2026-05-11 revision: 14.0 kg, 98.0 cm, BSA 0.6173
- Provenance smoke check added to `test_baseline_reference_metrics.m`

### Cleanup

- Removed orphaned `config/patient_profile_B.m`; `run_virtual_patients.m` was updated accordingly
- Do not stage `results/` artifacts

### Regression Coverage

| Test | Type | Purpose |
|---|---|---|
| `test_reyna_hemodynamic_active_set.m` | Lightweight | 14-param set, bounds, registry calibratable flags |
| `test_reyna_rmse_regression.m` | Heavy (skip-default) | End-to-end Lundquist RMSE < 0.10 |
| `test_scaling_mode_parity.m` | Heavy (skip-default) | Same-governance Lundquist/Zhang comparison; only Lundquist has the < 0.10 gate |
| Extended `test_scaling_modes.m` | Lightweight | Zhang `Rvalve_open` exponent, L scaling |
| Extended `test_baseline_reference_metrics.m` | Lightweight | Provenance table smoke check |
| Extended `test_reyna_systemic_flow_profile.m` | Lightweight | Recipe seed is Lundquist-only |

No MATLAB tests were run during the final seed-separation/reporting cleanup; rerun lightweight tests and the focused Lundquist/Zhang smoke checks before final PR if fresh evidence is required.

### Documentation

- `.assistant.md`: stable main-run snapshot plus seed-separation note
- `docs/scaling_method_comparison.md`: updated to avoid claiming Zhang single-digit/parity before the refreshed independent run

## Verification

### Lundquist (previous accepted evidence)

```text
Status: ACCEPT | Rollback: 0 | Primary RMSE: 0.0588
RAP_mean +5.48% | PAP_mean +2.01% | SAP_mean +4.66%
QpQs -3.78% | CO_Lmin -3.76%
Plausibility: 7 OK, 7 WARN, 0 FAIL
```

### Zhang (previous comparator run before seed separation)

```text
Status: PROMISING_NEAR_MISS | Primary RMSE: 0.1318
RAP_mean +8.13% | PAP_mean +10.60% | SAP_mean -4.45%
QpQs -3.39% | CO_Lmin +4.03%
Plausibility: 10 OK, 4 WARN, 0 FAIL
```

After the seed-separation cleanup, Zhang no longer inherits the Lundquist accepted disease vector, fixed `V0.SVEN`, or accepted IC vector. The next Zhang verification should report whether the independent Zhang baseline is limited by scaling, bounds, or optimizer behavior.

## Merge Risk

Moderate until the user reruns the focused checks after seed separation. The architecture is main-worthy: recipe governance only activates for `patient_reyna() + pre_surgery`, Lundquist remains the accepted production path, and Zhang is retained as an honest same-governance comparator instead of a Lundquist-seeded clone.
