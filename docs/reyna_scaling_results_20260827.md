# Reyna Scaling Experiment Results — 2026-08-27

> ## ⚠ Superseded — do not quote these numbers
>
> This document's results were calibrated against clinical inputs that the
> publication-readiness protocol-data reconciliation (2026-09-05) has since
> corrected: heart rate 136→119 bpm, VSD diameter 3.025→3.665 mm, and
> pre-surgery chamber volumes NaN→consistency-only 32/23.6/30.5/12 mL. See
> `docs/CHANGES_SINCE_PR22.md` §14 for the full before/after and
> `docs/publication_readiness_prd.md` for the re-run this triggered. This
> document is retained for its process findings (methodology, defects
> found, experiment-runner design), not for its numbers.

## Scope and controls

These results were generated on the isolated branch
`codex/luna-reyna-publishability-20260827` from the refreshed `origin/main`
base. Every reported run used:

- `patient_reyna()` and `pre_surgery`;
- frozen clinical input profile;
- governed primary metrics `RAP_mean`, `PAP_mean`, `SAP_mean`, `QpQs`, and
  `CO_Lmin`;
- full 14-parameter calibration set when GSA was disabled;
- serial MATLAB execution;
- the same higher-budget controls for the matched confirmations:
  `MaxFunctionEvaluations=60`, `MaxIterations=8`, and polish enabled.

## Matched higher-budget confirmations

| Arm | Seeds | Baseline primary RMSE | Calibrated primary RMSE | Full RMSE | Calibration status | Runner status |
|---|---:|---:|---:|---:|---|---|
| Fair Zhang | disabled | 0.267002 | **0.060445** | 0.063985 | **ACCEPT** | validated |
| Fair Lundquist BSA | disabled | 0.221945 | 0.221945 | 0.419975 | PHYSIOLOGICAL_BUT_POOR_FIT | validated |
| Operational Lundquist BSA | enabled | 0.515040 | 0.099028 | 0.120217 | PHYSIOLOGICAL_BUT_POOR_FIT | validated |

The fair Zhang candidate passed the 10% gate for all five primary metrics. Its
primary errors were approximately:

| Metric | Absolute error |
|---|---:|
| `RAP_mean` | 1.37% |
| `PAP_mean` | 3.45% |
| `SAP_mean` | 3.50% |
| `QpQs` | 1.52% |
| `CO_Lmin` | 6.18% |

The candidate still has five parameter-plausibility warnings, although all
parameters remain inside registry bounds. It should therefore be treated as a
strong publication candidate requiring scientific review, not as automatic
proof that Zhang is universally superior.

The fair Lundquist result did not improve with the higher budget. The
operational Lundquist result improved substantially from its seeded baseline,
but still failed two primary acceptance gates. This demonstrates that the
historical seed materially affects the observed Lundquist result.

## Screening repeatability

With polish disabled and the same screening budget, three fair Lundquist
repeats produced the identical primary RMSE `0.221945`. The fair Zhang
screening result was `0.267002`. These repeats are a deterministic
repeatability check, not an uncertainty interval.

## GSA status

The requested GSA execution was attempted with `GsaN=128`. MATLAB correctly
stopped at the GSA preflight because UQLab is not installed on the current
MATLAB path. The runner exported the failure with an actionable message:

> UQLab is required for DO_GSA=1 but was not found on the MATLAB path. Set
> UNIFIED_VSD_UQLAB_PATH to the external UQLab core folder.

Therefore the current numeric comparison is a calibration confirmation, not a
GSA-complete publication run. The code now supports external UQLab and SoBioS
paths through `UNIFIED_VSD_UQLAB_PATH` and `UNIFIED_VSD_SOBIOS_PATH`.

## Interpretation

The current evidence answers the practical question more precisely than a
single Zhang-vs-Lundquist number:

1. Zhang can improve Reyna strongly when given the same frozen clinical inputs
   and enough optimizer budget. Its best current governed primary RMSE is
   `0.060445`, with `ACCEPT` status.
2. Lundquist does not improve the fair-prior run under the tested controls.
3. The operational Lundquist result benefits from its historical seed, but the
   resulting candidate remains `PHYSIOLOGICAL_BUT_POOR_FIT` in this rerun.
4. The result does not prove Zhang is the generally correct pediatric scaling
   law. It supports reporting Zhang as the current Reyna primary candidate and
   Lundquist as a sensitivity/operational comparator, subject to UQLab-backed
   GSA and plausibility review.
5. The critical Reyna chamber-volume versus catheter-flow inconsistency is
   unchanged and remains a data-governance limitation rather than a scaling
   success.

## Evidence paths

- Fair Zhang experiment: `results/luna_experiments/reyna_p1_scaling_v1_20260827_223131/`
- Fair Zhang run: `results/runs/20260827_223135_reyna_pre_surgery/`
- Fair Lundquist experiment: `results/luna_experiments/reyna_p1_scaling_v1_20260827_222108/`
- Fair Lundquist run: `results/runs/20260827_222111_reyna_pre_surgery/`
- Operational Lundquist experiment: `results/luna_experiments/reyna_p1_scaling_v1_20260827_224007/`
- Operational Lundquist run: `results/runs/20260827_224011_reyna_pre_surgery/`
- GSA preflight experiment: `results/luna_experiments/reyna_p1_scaling_v1_20260827_222022/`
