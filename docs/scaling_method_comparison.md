# Scaling Method Comparison: Zhang vs Lundquist

Date: 2026-05-14

## Purpose

This document compares the two active pediatric scaling methods in the repository:

- `scaling.method = "zhang"`
- `scaling.method = "lundquist"` / internal alias `lundquist_bsa`

The decision rule is scientific defensibility, not lowest RMSE alone. Calibration should correct patient-specific mismatch, not compensate for a poor scaled baseline.

## Zhang Scaling

The Zhang option is weight-based. In the current implementation it scales:

| Quantity | Current scaling variable |
|---|---|
| heart rate | body weight ratio |
| chamber elastance | body weight ratio, different LV/LA and RV/RA exponents |
| unstressed volume | body weight ratio |
| vascular resistance | body weight ratio, separate systemic and pulmonary exponents |
| vascular compliance | body weight ratio |
| valve open resistance | body weight ratio |
| inertance | not scaled |
| pressure targets | not scaled |
| flow targets | not scaled |

This is motivated by Zhang et al.'s pediatric multiscale cardiovascular-model paper, which explicitly introduces allometric scaling laws for infants, children, and adolescents: [PubMed 31005012](https://pubmed.ncbi.nlm.nih.gov/31005012/).

## Lundquist Scaling

The Lundquist option is BSA-based. In the current implementation it scales:

| Quantity | Current scaling variable |
|---|---|
| heart rate | BSA ratio |
| chamber elastance | BSA ratio, different LV/LA and RV/RA exponents |
| unstressed volume | BSA ratio |
| vascular resistance | BSA ratio |
| vascular compliance | BSA ratio |
| vascular inertance | BSA ratio |
| valve open resistance | BSA ratio |
| pressure targets | not scaled |
| flow targets | not scaled |

This is motivated by Lundquist et al.'s 0D cardiovascular-model scaling work using patient anthropometrics, age, weight, height, and sex: [University of Twente record](https://research.utwente.nl/en/publications/patient-specific-size-and-age-scaling-in-a-zero-dimensional-cardi/).

## Same-Patient Scaling Factors

For Reyna (`weight = 13.4 kg`, `height = 95 cm`, `BSA = 0.588 m2`):

| Quantity | Zhang | Lundquist |
|---|---:|---:|
| HR baseline before clinical override | 123.16 bpm | 107.08 bpm |
| systemic resistance factor | 2.193 | 2.942 |
| pulmonary resistance factor | 3.181 | 2.942 |
| systemic compliance factor | 0.191 | 0.340 |
| LV elastance factor | 2.286 | 2.942 |
| RV elastance factor | 3.455 | 5.047 |
| LV/RV V0 factor | 0.266 | 0.340 |
| inertance factor | not scaled | 2.942 |

## Bounded Pipeline Comparison

Both methods were run through the same bounded pipeline:

```matlab
setenv('UNIFIED_VSD_DO_GSA','0')
setenv('UNIFIED_VSD_DO_PLOTS','0')
setenv('UNIFIED_VSD_DO_OVERLAY','0')
setenv('UNIFIED_VSD_FAST_CALIBRATION','1')
setenv('UNIFIED_VSD_MAX_FUN_EVALS','60')
setenv('UNIFIED_VSD_MAX_ITERATIONS','8')
setenv('UNIFIED_VSD_SCALING_MODE','lundquist') % or 'zhang'
clinical = patient_reyna();
main_run('pre_surgery', clinical)
```

Outputs:

- Lundquist run: `results/runs/20260514_155253_reyna_pre_surgery`
- Zhang run: `results/runs/20260514_160222_reyna_pre_surgery`
- Comparison bundle: `results/scaling_method_comparison/summary_20260514.csv`

## Results Summary

These values use the saved best/scientific candidate, not the rolled-back accepted artifact.

| Method | Baseline primary RMSE | Best primary RMSE | Best full RMSE | Hard RMSE | Soft RMSE | Param warnings | Primary gate failures |
|---|---:|---:|---:|---:|---:|---:|---|
| Lundquist | 0.2363 | 0.1065 | 0.1238 | 0.0932 | 0.0934 | 7 | QpQs, CO_Lmin |
| Zhang | 0.2317 | 0.1734 | 0.1891 | 0.1012 | 0.1146 | 5 | RAP_mean, PAP_mean, QpQs, CO_Lmin |

Selected calibrated residuals:

| Method | CO_Lmin | QpQs | PAP_mean | SAP_mean | RAP_mean | RVEDV check | RVESV |
|---|---:|---:|---:|---:|---:|---:|---:|
| Lundquist | -8.79% | -5.99% | -3.87% | -0.13% | +0.17% | +25.88% | -10.07% |
| Zhang | -7.24% | -11.53% | -5.02% | -0.44% | +6.56% | +33.14% | +17.32% |

## Interpretation

Lundquist is currently the better scientific default for Reyna because:

- It gives a much better calibrated primary RMSE under the same bounded optimization settings.
- It preserves SAP_mean, RAP_mean, and PAP_mean more tightly after calibration.
- It produces a lower full RMSE even while keeping RVEDV visible as a consistency-check-only target.
- Its main weakness is parameter distortion: 7 plausibility warnings versus 5 for Zhang.

Zhang is not rejected as invalid. It may be useful for broader pediatric sensitivity analysis because it produces fewer parameter-plausibility warnings in this bounded run. However, for this Reyna dataset it leaves larger residuals in QpQs/RAP/PAP and does not provide a better physiologic baseline after calibration.

## Default Selection

Keep `lundquist_bsa` as the default for now.

Reason: Zhang has slightly better parameter-warning count, but Lundquist gives the more competitive pressure-flow-volume fit under the new target governance while retaining explicit RVEDV transparency.

## Remaining Uncertainty

- The bounded comparison used reduced fmincon evaluations for reproducibility. A longer run may improve either method.
- The current Lundquist implementation uses a BSA law that is dimensionally simple and traceable, but the exact source exponent mapping should be reviewed against the Lundquist paper text before making a thesis-level claim.
- The Zhang implementation uses weight-based exponents from the project scaling logic. Its inertance is not scaled, so a direct method comparison is not perfectly symmetrical.
- Neither scaling method resolves the clinical RV stroke-volume inconsistency. That is a data-governance issue, not a scaling issue.
