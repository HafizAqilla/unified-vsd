# Response to Supervisor Comments — 2026-08-28

Source: `REV-SYLVI Hafiz Skripsi_Revisi22Juni_2 (2).docx`
Reviewer: Sylvi Febriana Rachmawati Irnadiastputri, 2026-06-24 (24 comments)
Branch: `codex/reyna-zhang-fullmetric-10pct`

This tracks only the comments that this branch's code work addresses. Comments
on citation placement, abbreviation use, and Indonesian punctuation are
manuscript edits and are not listed.

---

## Addressed by code in this branch

### [273] "Bagaimana estimasi nilai volume pra-tindakan, apakah akan lebih tinggi/rendah?"

*How are pre-procedure volume values estimated — will they be higher or lower?*

This was the sharpest comment in the set, and it exposed a real defect.

The chamber volumes in the manuscript are **H+1 post-operative echo**. They were
being used as pre-operative calibration targets, which is not an estimation
problem — it is a category error. Reyna's chambers after VSD closure are not
Reyna's chambers with an open VSD.

The block is now removed as a pre-surgery target entirely
(`config/calibration_recipes/reyna_pre_surgery.m`), and the pre-operative
chamber state is reported as what it actually is: a **model prediction with no
comparator**, via the new `predicted_chamber_state_report`.

That report answers the direction question explicitly. For a restrictive
left-to-right VSD:

| Quantity | Expected pre-operatively | Reason |
|---|---|---|
| `LVEDV` | **higher** than post-closure | shunt raises pulmonary blood flow, so pulmonary venous return volume-loads the LV; closure removes that load |
| `LVEF` | **higher or equal** | part of LV ejection goes into the lower-resistance RV/pulmonary bed, reducing effective afterload; closure restores full systemic afterload |
| `RVEDV`, `RVESV` | **not specified** | in a VSD the shunt is ejected in systole largely into the pulmonary artery, so the RV is typically pressure-loaded rather than volume-loaded; direction depends on defect size and pulmonary vascular state |

The report screens each prediction against the paediatric plausibility bands in
`clinical_reference_ranges.m` and flags any prediction that moves *against* the
expected direction. It states in its own output that the post-operative column
is a **direction check, not an accuracy target**.

### [295] "Ini membahas sistol, jadi perlu membahas diastol juga secara eksplisit"

*This discusses systole, so diastole needs explicit discussion too.*

Directly connected to a defect this branch fixes. `PAP_min` **is** the pulmonary
artery diastolic pressure, and `SAP_min` is the systemic diastolic pressure.
Both were declared calibration targets but assigned no tier, so they were
graded inside the reported RMSE while being invisible to the objective.

Measured consequence at baseline: `PAP_min` sat at **23.83%** error and
`SAP_min` at **13.98%**, while every systolic and mean pressure was within 10%.
The manuscript discussed systole better than diastole because the model was
*fitted* on systole and mean pressures and not on diastole.

Both diastolic pressures are now fitted (`soft` tier, weights 0.45 / 0.40), so
there is now a diastolic result to discuss rather than an artefact to explain.

### [221] "Apa itu primary governed? Full transparent? Bedanya apa, implikasinya apa?"

*What is primary governed? Full transparent? What is the difference, what are the implications?*

The distinction was real but under-explained, and one side of it was wrong.

- **primary governed** — RMSE over targets that are both measured for this
  patient and admitted as fitting targets. Excludes consistency-check-only,
  derived, and holdout rows.
- **full transparent** — RMSE over every target with a finite patient
  comparator, including rows deliberately excluded from fitting.

The implication is that `full_transparent` is normally the *worse* number, and
publishing both is what prevents an acceptance claim resting on a
favourably-chosen subset.

What was missing is that neither figure showed **how many individual metrics
passed**. A run could report a good RMSE while a single metric sat far outside
the acceptance band — and it did. `export_full_metric_gate` now writes a
per-metric table on every run and reports the count as *n of N* against a
disclosed denominator.

### [242] / [243] / [252] "Perbanyak grafik dibanding tabel... jadikan grafik"

*Use more graphs than tables; make the error comparison a line/bar graph; this table is hard to read, make it a graph.*

Added `plot_calibration_error_comparison`, which writes two figures per run:

1. **Grouped bar chart** — baseline versus calibrated absolute error per
   metric, sorted worst-first, with the 10% acceptance band drawn as a
   reference line so the threshold is visible rather than mental arithmetic.
2. **Paired slope chart** — the before/after direction for each metric, green
   where error fell and red where it rose. This is the "efek kalibrasi
   (before-after)" chart requested in [243], and is deliberately distinct from
   the sensitivity heatmap on the following manuscript page.

Both export at 200 dpi to the run's `figures/` folder, and are written even
when waveform plotting is disabled, because they are reporting evidence rather
than diagnostics.

### [217] "Bagaimana perbandingan zhang vs lundquist (dalam kata-kata)? Apakah sesuai rentang pediatrik yg umum?"

*How do Zhang and Lundquist compare, in words? Do they fall within common paediatric ranges?*

**Second part** — yes, screening ranges exist (`clinical_reference_ranges.m`)
and the predicted chamber report now screens every model prediction against
them, reporting "within paediatric screening range: n of N". These are broad
plausibility bands, not normative percentiles, and the report says so.

The screen is not decorative. Before it was enforced, the calibrated candidate
returned `LVEF = 0.891` against a paediatric band of 0.40–0.85, with a 4.0 mL
end-systolic volume — a near-empty ventricle. The current candidate is inside
the band on **6 of 6** chamber quantities. So the answer to "does it match
common paediatric ranges" is now yes, and it is yes because the ranges were
made a constraint rather than a report.

**First part — an important caveat before writing that comparison.** The
Zhang-versus-Lundquist numbers as currently published are not interpretable in
words, because the Lundquist arm's optimiser never moved: baseline and
calibrated primary RMSE were identical to six decimal places (`0.221945`). With
`MaxFunctionEvaluations = 60` against 14 free parameters — roughly four
gradient steps — neither arm converged. The honest statement is:

> Under a budget too small to converge either arm, the Zhang-scaled starting
> point was closer to the clinical targets than the Lundquist-scaled one.

That is a statement about **priors**, not about which scaling law is correct.
The branch adds an `OPTIMIZER_DID_NOT_MOVE` classification so this failure mode
can no longer be presented as a comparison result. A defensible written
comparison needs both arms re-run to convergence first.

### [229] "Belum ada pembahasan fisiologis untuk ini"

*No physiological discussion for this yet.*

The predicted chamber state report supplies the physiological reasoning for the
chamber results (shunt volume-loading of the LV, reduced effective afterload,
and why the RV direction is not specified), tied to explicit expectations
rather than post-hoc description.

### [291] "Bagaimana perbandingan dengan pediatrik sehat?"

*How does it compare with a healthy paediatric?*

Partially addressed: every chamber prediction is now screened against the
paediatric plausibility bands and the pass count is reported. A full
healthy-paediatric comparison would additionally need a matched healthy control
run (`config/patient_healthy_adult.m` is adult, not paediatric) — recorded as
outstanding below.

---

## Not addressed — recorded as outstanding

| Comment | Request | Status |
|---|---|---|
| [283] | "harusnya kita coba juga noise simulation" | Not implemented. Would need measurement noise injected on clinical targets and the calibration repeated to give error bars. Complements the multi-start machinery this branch adds. |
| [291] | Full healthy-paediatric comparison | Needs a matched paediatric control profile; only the plausibility screen is in place. |
| [249] | "Definisikan timing yang tepat. Di usia apa?" | Manuscript wording. Note the model does carry `age_years = 3.17` and an age-validity annotation. |
| [272] | "Pastikan ini sudah jelas di bagian metode pengambilan data" | Manuscript. The H+1 timing issue this branch surfaced should be stated explicitly in the data-collection method. |
| [138] | "Di alur ada 'simulink' tapi dalam desain tidak ada" | Manuscript. The repository is plain MATLAB — no Simulink model exists, so the flowchart should drop it. |
| [14], [15], [173], [225], [248], [286], [303], [306], [308] | Citation placement, abbreviations, punctuation, clinical relevance framing | Manuscript edits, outside this branch. |

---

## Recommended follow-up on the data itself

`clinical.post_surgery` in `patient_reyna()` currently has **every field set to
`NaN`** — the post-operative scenario has no clinical targets at all and cannot
be validated against anything.

The five H+1 chamber values removed from `pre_surgery` are exactly the
comparators `post_surgery` lacks. Relocating them is a clinical data-governance
decision and this branch does not make it, but it would convert the project's
most persistent liability into the post-surgery scenario's first validation
evidence. Two caveats: the 60% internal stroke-volume inconsistency travels
with the data, and whether H+1 represents a converged post-closure state is a
clinical judgement that should be stated either way.
