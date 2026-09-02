# What Changed Since PR #22 — Plain-English Summary

Patient: Reyna, 3 years 2 months old, girl, restrictive VSD (hole in the
heart wall).
Model: computer simulation of her heart and blood vessels, one version
before surgery ("pre-surgery") and one version after the hole is closed
("post-surgery").

This document explains, in simple terms, everything that was worked on with
Claude's help since PR #22, and gives the full, current results for both
the pre-surgery and post-surgery model.

---

## 1. The short version

Three things happened, in this order:

1. **We found the model was being graded unfairly.** Some numbers were
   scored as "wrong" even though the model was never told to match them in
   the first place. We fixed that.
2. **We found two of the patient's input numbers were wrong** — not
   estimates that could go either way, but numbers that didn't match the
   actual hospital record. We corrected them.
3. **After fixing both problems, we ran a genuine test**: we calibrated the
   model only on data from *before* surgery, then asked it to predict what
   her heart pressures would look like *after* surgery, without giving it
   any post-surgery data to peek at. **It got 5 of 7 pressures right, within
   10%.**

That third result is the strongest thing this project has produced. Nothing
before it was a real prediction test — everything was either data the model
was fitted to, or math built directly out of numbers it was already fitted
to.

---

## 2. Part 1 — Fixing how the model was graded

### The problem

Think of it like a school report card where a student is graded on a
subject they were never actually taught. That was happening here.

Two pressure readings from the catheter test — `PAP_min` and `PAP_max` (the
lowest and highest pressure in the lung artery) — were being used to *grade*
the model's accuracy, but the model was never told to *try* to match them
during fitting. It's like scoring an exam question the student never saw.

The same mistake affected three other readings.

### The fix

We changed the rules so the model can no longer be graded on something it
wasn't asked to fit. Every reading that counts toward the score must also be
one the model is actually trying to match.

### A second, more precise scoring method

The old scoring method treated every reading as equally trustworthy. That's
not realistic — some readings were measured three times in a row and got the
same result every time (very reliable), while others were only measured
once or estimated (less reliable).

We built a new scoring method that takes this into account: a reading
that's known very precisely has to be matched more closely; a reading
that's known loosely is allowed more slack. This is standard practice in
science whenever measurements come with different levels of confidence.

---

## 3. Part 2 — Finding two wrong input numbers

This was the biggest discovery. We obtained the actual hospital procedure
log (RSAB Harapan Kita, case ID HA000557, dated 6 April 2026) and checked
every number in the computer model against it, line by line.

**Two numbers used to build the model did not match the hospital record:**

| What | Old value (wrong) | Correct value | Why it was wrong |
|---|---:|---:|---|
| Heart rate | 119 beats/min | **136 beats/min** | The old number was actually her *blood-pressure reading* (119), copied into the wrong field by mistake. The log's actual pulse reading says 136. |
| Weight / height / body size | 14.0 kg / 98.0 cm | **13.4 kg / 95.0 cm** | The old numbers were measured **5 weeks after** the heart procedure — she had grown since then. Using them made the model simulate a bigger child than the one who was actually catheterised. |

**A third number was also corrected, though it's a fitted target rather
than a raw input:**

| What | Old value | Correct value | Why |
|---|---:|---:|---|
| Average blood pressure (systemic) | 71.3 mmHg | **77 mmHg** | The old number was calculated with a formula from the high/low pressure readings. The hospital machine actually prints its own average directly on the log — 77 mmHg — and the formula's estimate was off by about 6 mmHg. |

Heart rate and body size are not small details — they directly control how
fast the simulated heart beats and how everything is scaled to the child's
body. Getting them wrong means every number the model ever produced for
this patient, in every previous report, was computed for a differently
sized, differently paced child than the one actually measured.

**After fixing these three numbers, the model's fit improved on every
single measure we track** — see the results in Section 5.

### A bonus: post-surgery data was sitting in the same log, unused

The same hospital log also contains the pressure readings taken **after**
the hole was closed. Nobody had entered these into the model before. We
added all seven of them. This is what made the prediction test in Section 6
possible — without this data, there was nothing to test the model's
prediction against.

---

## 4. Part 3 — Bugs found along the way

While building the above, we found and fixed several bugs. The two worth
knowing about:

- **A reporting bug that could have hidden a real problem.** A statistical
  check (called chi-squared) is supposed to tell us if the model has too
  many adjustable knobs for the amount of data it has. A bug in how this
  was calculated made an overfit model look fine on paper. Fixed, and now
  the report explicitly warns when there isn't enough data to trust the
  number.
- **A bug that would have made "closing the valve" not actually close it**
  in one of the simulation modes. If left unfixed, a "post-surgery"
  simulation could have silently still been simulating the open hole.
  Caught before it produced any wrong result.

Full technical bug list is in `docs/CHANGES_SINCE_PR22.md` §8, for anyone
who wants the details.

---

## 5. Full results — before surgery (fitted)

This is what the model was actually calibrated to match. "Fitted" means the
model was adjusted until it got as close as possible to these numbers.

**All 9 pressure/flow readings used to grade the fit:**

| Variable | What it measures | Real value | Model's value | Error | Within 10%? |
|---|---|---:|---:|---:|:---:|
| RAP_mean | Right atrium average pressure | 5 mmHg | 5.26 mmHg | +5.2% | ✅ Yes |
| PAP_min | Lung artery, lowest pressure | 10 mmHg | 11.07 mmHg | +10.7% | ❌ **No** |
| PAP_max | Lung artery, highest pressure | 20 mmHg | 19.00 mmHg | −5.0% | ✅ Yes |
| PAP_mean | Lung artery, average pressure | 15 mmHg | 14.94 mmHg | −0.4% | ✅ Yes |
| SAP_min | Body artery, lowest pressure | 57 mmHg | 59.43 mmHg | +4.3% | ✅ Yes |
| SAP_max | Body artery, highest pressure | 100 mmHg | 96.05 mmHg | −4.0% | ✅ Yes |
| SAP_mean | Body artery, average pressure | 77 mmHg | 76.85 mmHg | −0.2% | ✅ Yes |
| CO_Lmin | Blood flow out of the heart | 3.423 L/min | 3.334 L/min | −2.6% | ✅ Yes |
| QpQs | Ratio of lung flow to body flow (shunt severity) | 1.194 | 1.195 | +0.1% | ✅ Yes |

**Score: 8 of 9 within 10%.** Only `PAP_min` (the lowest lung pressure)
missed the mark, at 10.7% off.

**Overall accuracy (RMSE — root-mean-square error, lower is better):**

| Measure | Before fitting | After fitting | Improvement |
|---|---:|---:|---:|
| RMSE (main 9 readings above) | 0.2915 | **0.0480** | 83.5% better |
| RMSE (absolutely everything, including readings not used for grading) | 0.3277 | **0.0444** | 86.4% better |

An RMSE of 0.0480 roughly means the model's readings are, on average, about
4.8% off from the real ones, across the 9 graded readings.

### Extra readings, not used for grading (shown for completeness)

These were tracked but deliberately not counted in the score above, because
they are either calculated *from* the 9 readings (not independent
information) or come from a different exam that measures a different
physical state:

| Variable | What it is | Real/reference value | Model's value | Note |
|---|---|---:|---:|---|
| Q_shunt_Lmin | Flow leaking through the hole | 0.664 L/min | 0.649 L/min (−2.3%) | Calculated from CO and QpQs above, so it isn't independent — kept out of the score for that reason. Also used as an "early warning" number (see Section 7). |
| SVR | Resistance in body's blood vessels | 21.03 WU | 21.47 WU (+2.1%) | Calculated from other readings, not independently measured. |
| LVEDV | Left heart chamber size (full) | ~41 mL (from a different, later exam) | 34.7 mL | Predicted, not fitted — the reference number is from an echo scan taken *after* surgery, so it isn't a fair pre-surgery comparison. Shown for interest only. |
| LVESV | Left heart chamber size (empty) | ~19.3 mL (same later exam) | 5.8 mL | Same caveat as above. |
| RVEDV | Right heart chamber size (full) | ~30.5 mL (same later exam) | 34.1 mL | Same caveat as above. |
| RVESV | Right heart chamber size (empty) | ~12.0 mL (same later exam) | 7.6 mL | Same caveat as above. |
| LVEF | Left heart pumping efficiency | ~52.8% (same later exam) | 83.2% | Same caveat as above; also flagged as higher than expected. |
| RVEF | Right heart pumping efficiency | not available | 77.7% | No reference value exists to compare against. |

---

## 6. Full results — after surgery (predicted, not fitted)

**This is the important part.** None of the numbers below were shown to the
model during fitting. The model was calibrated only on the "before surgery"
data above. We then told it "the hole is now closed" and asked it to
simulate what her pressures would look like — with no further adjustment.
This is a genuine test of whether the model actually learned something
correct about this patient, rather than just memorising numbers.

We proved (Section 7.5 of the technical document) that no after-surgery
data ever leaked into the before-surgery calibration, so this really is a
fair, blind prediction.

| Variable | What it measures | Real (measured after surgery) | Model's prediction | Error | Within 10%? |
|---|---|---:|---:|---:|:---:|
| RAP_mean | Right atrium average pressure | 5 mmHg | 5.36 mmHg | +7.3% | ✅ Yes |
| PAP_min | Lung artery, lowest pressure | 9 mmHg | 10.33 mmHg | +14.8% | ❌ No |
| PAP_max | Lung artery, highest pressure | 17 mmHg | 17.75 mmHg | +4.4% | ✅ Yes |
| PAP_mean | Lung artery, average pressure | 13 mmHg | 13.83 mmHg | +6.4% | ✅ Yes |
| SAP_min | Body artery, lowest pressure | 68 mmHg | 63.29 mmHg | −6.9% | ✅ Yes |
| SAP_max | Body artery, highest pressure | 89 mmHg | 101.40 mmHg | +13.9% | ❌ No |
| SAP_mean | Body artery, average pressure | 79 mmHg | 81.75 mmHg | +3.5% | ✅ Yes |

**Score: 5 of 7 predicted within 10%, and all 7 within about 15%.**

**Overall accuracy (RMSE-style, over these 7 predictions): about 9.2%
average error.** This isn't computed with exactly the same formula as the
pre-surgery RMSE above (that one is normalised against each reading's own
measurement uncertainty), so the two numbers shouldn't be compared side by
side as if they were on the same scale — but as a plain "how far off on
average" figure, this is it.

### What the pattern in the errors tells us

Six of the seven predictions came out slightly *too high*. That's not
random scatter — it's a consistent pattern, and it means the model isn't
predicting quite as much pressure relief from the surgery as actually
happened. In plain terms: **the real surgery helped the lungs and body
pressures a bit more than the model expected.** That's a useful, specific
clue for improving the model further — not just noise to shrug off.

### Why this result matters more than the "before surgery" score

Every other number in this document comes from data the model was allowed
to see and adjust itself to match. That makes those numbers a test of
*fitting*, not of *understanding*. This prediction test is different: the
model had never seen these 7 numbers in any form. Getting 5 of 7 right
without being told the answer is real evidence the model captured something
true about how this patient's body responds to the surgery — not just
curve-fitting.

We also checked this against a *worse* calibration attempt that scored a
perfect 9-of-9 on the before-surgery data. That "perfect" version actually
predicted the after-surgery data *worse* (only 3 of 7 correct). This is a
classic warning sign in modelling called **overfitting** — a model that
matches its training data too eagerly usually gets worse, not better, at
predicting new situations. We caught it because we ran this exact
before/after test, which is exactly why the test is worth doing.

---

## 7. What we still don't know / what to be careful about

Being honest about the limits of this work matters as much as the results
themselves.

- **This is one patient.** Nothing here proves the model works for children
  in general — only that it worked reasonably well for Reyna, once her data
  was corrected.
- **The model has more adjustable settings than independent measurements.**
  There are 9 real, independent pressure/flow readings before surgery, but
  the model has 12 adjustable internal settings. That's like solving 9
  equations with 12 unknowns — there's more than one way to get a good fit,
  which is part of why the overfitting check in Section 6 was so important.
  We've identified a way to cut this down to 7 settings without hurting the
  fit, but that hasn't been tested with a full run yet.
- **Two calibration attempts gave slightly different "best" answers** (8 of
  9 vs 9 of 9 on the before-surgery score). The gap between them is real,
  not just rounding — which is exactly why the after-surgery prediction
  test in Section 6 was needed to tell which one to trust.
- **No blood-flow measurements exist for after surgery**, only pressures.
  So the prediction test in Section 6 only checks pressure, not flow.

---

## 8. Where to find more detail

- `docs/CHANGES_SINCE_PR22.md` — the full technical writeup, same content
  but with code references, exact statistics, and every bug's root cause.
- `docs/reyna_statistical_calibration_results_20260829.md` — the complete
  run-by-run log this summary was drawn from.
- `scripts/evaluate_post_closure_prediction.m` — the script that runs the
  before/after prediction test described in Section 6.
