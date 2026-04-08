# Baseline Accuracy and Parameter Ratio Report (2026-04-07)

## Summary
- Total baseline checks: 15
- Passed checks: 14
- Failed checks: 1
- Overall accuracy: 93.33%

## Accuracy Formula
\[
	ext{Accuracy} = \frac{\text{Passed}}{\text{Total}} \times 100\% = \frac{14}{15} \times 100\% = 93.33\%
\]

## Ratio Definition (Real vs Simulated)
For parameters with a reference range \([L, U]\):
\[
	ext{Reference Midpoint} = \frac{L+U}{2},\qquad
	ext{Ratio} = \frac{\text{Simulated}}{\text{Reference Midpoint}} \times 100\%
\]

Absolute percentage error to midpoint:
\[
	ext{Error \%} = \left|\frac{\text{Simulated} - \text{Reference Midpoint}}{\text{Reference Midpoint}}\right| \times 100\%
\]

For one-sided limits:
\[
	ext{Limit Usage Ratio} = \frac{\text{Simulated}}{\text{Allowed Limit}} \times 100\%
\]

## Per-Parameter Real vs Simulated Ratio

| Parameter | Simulated | Reference (Real/Physio) | Ratio (%) | Status |
|---|---:|---:|---:|---|
| SAP_max (mmHg) | 113.529 | [100, 140] (mid 120.0) | 94.61 | PASS |
| SAP_min (mmHg) | 55.561 | [60, 90] (mid 75.0) | 74.08 | FAIL |
| SAP_mean (mmHg) | 82.939 | [70, 100] (mid 85.0) | 97.58 | PASS |
| HR (bpm) | 75.000 | [65, 85] (mid 75.0) | 100.00 | PASS |
| QpQs (-) | 1.000 | [0.95, 1.05] (mid 1.0) | 100.00 | PASS |
| LVEF (fraction) | 0.666 | [0.55, 0.75] (mid 0.65) | 102.46 | PASS |
| PAP_mean (mmHg) | 13.492 | [10, 20] (mid 15.0) | 89.95 | PASS |
| RAP_mean (mmHg) | 4.811 | [0, 8] (mid 4.0) | 120.28 | PASS |
| Total blood-volume drift (mL) | 0.220 | < 1.0 (limit) | 22.00 of limit | PASS |
| Last-cycle drift (mL) | 0.0002 | < 0.10 (limit) | 0.20 of limit | PASS |
| Stroke work SW_LV (mmHg.mL) | 6797.815 | [3000, 12000] (mid 7500.0) | 90.64 | PASS |
| LVEDP (mmHg) | 4.666 | [0, 15] (mid 7.5) | 62.21 | PASS |
| LVESP (mmHg) | 109.358 | [80, 130] (mid 105.0) | 104.15 | PASS |
| Emax_derived (mmHg/mL) | 4.800 | [2.85, 6.65] (mid 4.75) | 101.05 | PASS |
| Steady-state flag | TRUE | Must be TRUE | 100.00 | PASS |

## Key Point
The only remaining mismatch is SAP_min, which is below the physiological lower bound. All other gates satisfy the baseline criteria.

## Strict <10% Error Check

Using midpoint-based error for range-defined metrics, the <10% target is **not yet fully met**.

| Parameter | Error (%) | <10% Target |
|---|---:|---|
| SAP_max | 5.39 | PASS |
| SAP_min | 25.92 | FAIL |
| SAP_mean | 2.42 | PASS |
| HR | 0.00 | PASS |
| QpQs | 0.00 | PASS |
| LVEF | 2.46 | PASS |
| PAP_mean | 10.05 | FAIL |
| RAP_mean | 20.28 | FAIL |
| Stroke work SW_LV | 9.36 | PASS |
| LVEDP | 37.79 | FAIL |
| LVESP | 4.15 | PASS |
| Emax_derived | 1.05 | PASS |

Summary for range-based metrics:
- Pass <10%: 8/12
- Fail <10%: 4/12 (SAP_min, PAP_mean, RAP_mean, LVEDP)

Note:
- This strict <10% criterion is different from gate-range pass/fail. A metric can be inside an acceptable physiological range but still be >10% from the midpoint of that range.
