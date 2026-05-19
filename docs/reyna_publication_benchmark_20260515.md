# Reyna Publication Benchmark - 2026-05-15

## Current Model Result

Accepted run:

`results/runs/20260515_173742_reyna_targeted_grid_polish`

| Metric group | Result |
|---|---:|
| Governed primary RMSE | 0.088354 |
| Governed primary MSE | 0.007806 |
| Full transparent RMSE | 0.109307 |
| Full transparent MSE | 0.011948 |
| Calibration-fit targets under 10 percent | 10 / 12 |
| Pressure-flow targets under 10 percent | 8 / 8 |
| All available clinical rows under 10 percent | 10 / 14 |

Remaining visible failures:

| Metric | Error | Governance |
|---|---:|---|
| `LAP_mean` | -12.274 percent | validation-only |
| `LVEDV` | +14.653 percent | calibration target |
| `LVEF` | +13.909 percent | calibration target |
| `RVEDV` | +25.649 percent | consistency-check only |

The pressure-flow block is strong. The residual failures are concentrated in
the echo-derived chamber volume/function block, consistent with the clinical
consistency audit showing disagreement between Fick/QpQs-derived stroke volume
and echo-derived stroke volume.

## Literature Benchmark

| Study | Disease / cohort | Calibration targets | Reported fit | Comparison with Reyna |
|---|---|---|---|---|
| Bozkurt et al. 2022 | Pediatric dilated cardiomyopathy, 3 patients | Pressures, ventricular volumes, cardiac output, ventricular dimensions | Hemodynamic differences after optimization were reported below 10 percent | Stronger strict-target closure than Reyna, but not VSD and used CMR-derived ventricular geometry |
| Valenti 2023 VSD | Adult/clinical VSD thesis case | Pressures and Qp/Qs-style flow outputs | Mean squared relative error 2.77e-2; pulmonary flow not reproduced well | Reyna governed primary MSE 7.81e-3 and pressure-flow block is stronger, but Reyna has volume-target conflict |
| Valenti 2023 VSD with vasodilator | Same VSD patient under vasodilators | Pressures and pulmonary/systemic flows | Mean squared relative error 1.10e-2; pulmonary pressure and flow still imperfect | Reyna full transparent MSE 1.19e-2 is similar; governed primary MSE is better |
| Valenti 2023 ASD/VSD | Combined septal defect | Pressures and flows | Mean squared relative error 4.81e-2; pressures/flows not all recovered | Reyna is stronger numerically |
| Zhang et al. 2019 | Generic infant/child/adolescent cardiovascular model | Age-related hemodynamics against literature/in vivo ranges | Broad age-scaling validation, not patient-specific VSD fitting | Supports pediatric allometric scaling, but is not a direct case-fit benchmark |
| Lundquist et al. 2026 | 0D model across age, sex, body size | Aggregate hemodynamic realism by Z-scores | Pediatric aggregate Z-score 0.69 | Supports our choice to prefer anthropometric/BSA scaling, but not a VSD calibration benchmark |
| Colebank et al. 2026 preprint | Pediatric isolated VSD simulations | VSD size, age-dependent scaling, vascular sensitivity | Foundation for digital shadows, not direct patient-specific validation | Supports novelty/need: relatively few VSD-specific pediatric models exist |

## Publishability Assessment

The result is credible enough for a thesis chapter, conference abstract, or
engineering-methods manuscript if framed carefully:

- A patient-specific pediatric VSD 0D model.
- A comparison of Zhang-style and Lundquist-style pediatric scaling priors.
- A target-governance method for inconsistent catheter and echo data.
- A pressure-flow calibration that places all direct hemodynamic targets within
  10 percent.

It is not yet strong enough to claim a clinically validated digital twin,
because:

- Only one patient scenario is currently convincing.
- Several visible volume/function targets remain above 10 percent.
- The RVEDV discrepancy is large unless explained as a clinical-data conflict.
- Uncertainty and identifiability should be reported formally.

## Required Before Manuscript Submission

1. Add a benchmark table in the paper that separates pressure-flow targets from
   echo-volume targets.
2. Report the clinical consistency audit before calibration.
3. Add a flow-implied ventricular volume calculation:
   `RVEDV_flow_implied = RVESV + Qp * 1000 / HR`.
4. Run an uncertainty ensemble using clinical measurement uncertainty and report
   whether clinical targets fall within model prediction intervals.
5. Add at least one second case, post-surgery case, or synthetic recovery test.
6. Include parameter identifiability/sensitivity rankings for fitted parameters.
7. Keep `RVEDV` visible, but explicitly label it consistency-check-only.

## Sources

- Bozkurt et al. 2022: https://link.springer.com/article/10.1007/s13239-022-00611-9
- Valenti 2023 thesis: https://www.politesi.polimi.it/handle/10589/215783
- Shimizu et al. 2017 review: https://doi.org/10.1007/s12576-017-0585-1
- Zhang et al. 2019: https://pubmed.ncbi.nlm.nih.gov/31005012/
- Lundquist et al. 2026: https://researchinformation.umcutrecht.nl/en/publications/patient-specific-size-and-age-scaling-in-a-zero-dimensional-cardi/
- Colebank et al. 2026 preprint: https://arxiv.org/abs/2602.04008
