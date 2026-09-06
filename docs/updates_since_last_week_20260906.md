# What We Did This Past Week — Plain-English Update (as of 2026-09-06)

This is a simple summary of everything done on the model and the repo since
last week, written for anyone who doesn't want to read the technical
commit history. If you want the full technical detail, see
`docs/CHANGES_SINCE_PR22.md` (the complete change log) and
`docs/reyna_publication_readiness_results_20260906.md` (the actual
numbers).

## 1. The goal: make the repo ready for publication

The request was to clean up the repository, fix anything broken or
misleading, add proper citations, and then re-run the model on the most
complete and correct version of Reyna's data — so the results can actually
be trusted and quoted in the manuscript/thesis.

This was done in phases, all merged into one branch
(`codex/publication-readiness`), opened as
[PR #26](https://github.com/HafizAqilla/unified-vsd/pull/26).

## 2. Cleanup

- Removed real patient identifiers (hospital name, MRN, full names) from
  every tracked file in the repo. The real identifiers now live only in
  one file on this computer that is never uploaded to GitHub.
- Deleted 14 old/dead scripts and files that were no longer used, without
  changing any actual model behavior.
- Found and fixed 3 real bugs sitting quietly in the code:
  - One mis-routed how a certain patient case was categorized.
  - One meant a patient config file was missing a few fields it should
    have had.
  - One was a leftover placeholder value that didn't match how real
    patient files are written.

## 3. Fixing the clinical data itself

This is the part that actually changes numbers, so it's worth explaining
carefully. We got hold of the official IRB-approved protocol form for
Reyna's case and checked every single number in the code against it,
line by line. Two things turned out to be wrong, and one thing turned out
to be more complicated than we first thought:

1. **Heart rate.** The code had 136 bpm. The protocol form says 119 bpm is
   the correct *average* heart rate for that session. Fixed.
2. **VSD (hole) diameter.** The code had two different values in two
   different files that disagreed with each other (3.025 mm vs 3.665 mm) —
   a silent copy-paste-style bug. The protocol form says 3.665 mm. Fixed,
   and the duplicate is now deleted so this can't happen again.
3. **Chamber (heart chamber) volumes — this one took two tries.**
   - First pass: we thought these volumes were measured *before* surgery,
     at the same time as everything else, so we added them into the
     pre-surgery data.
   - Second pass (later the same day): it turned out the label on the
     protocol form ("pre-release occluder") actually means the closure
     device was already placed and blocking the hole — it just hadn't
     been mechanically detached yet. In other words, these volumes were
     actually measured **after** the hole was already closed, not before.
     We moved them to the post-surgery data instead.
   - Along the way we also found a genuinely separate, older bug: the
     code had no way to even record two of these volume numbers for
     after-surgery patients (a missing field in the underlying table). Fixed
     that too, since otherwise those two numbers would keep being
     silently invisible to every future calculation regardless of the
     timing question above.

## 4. Citations

- Added a proper reference list file (`docs/references.bib`).
- Double-checked the two main scientific papers this model's pediatric
  scaling is based on (Zhang 2019 and Lundquist 2025) against the actual
  publisher records — and caught that one paper's title had been
  mis-transcribed in an old doc for months. Fixed.
- Fixed an outdated/incorrect claim in the README about which scaling
  method the model actually prefers by default.

## 5. Running the model — what's done, what's next

We then ran the model itself against the corrected data.

**Done so far — one full pass:**
- Compared the two competing pediatric scaling methods (Zhang vs.
  Lundquist) head-to-head, on equal footing. **Zhang came out clearly
  ahead** for this patient's data.
- Took the winning method and ran it properly (a real 6-attempt search,
  not just one), landing on the best result the model has produced for
  Reyna so far — a good hemodynamic (pressure/flow) fit, though still
  short of the strictest "fully accepted" bar.
- Checked whether the model could predict what happens *after* surgery,
  using only the *before*-surgery fit, as an honest test. The pressure
  predictions were reasonable-to-fair. The heart-chamber-function
  predictions were poor — the model overestimated post-surgery heart
  function by a wide margin. This is reported honestly, not hidden,
  because a model that looks great on paper but fails this kind of check
  is exactly the kind of overclaiming this whole cleanup was meant to
  prevent.
- Also ran an analysis showing that far fewer parameters (7 instead of 12)
  could describe the data just as well, with much better statistical
  footing — a good lead for simplifying the model, but not yet confirmed
  by actually re-running the calibration with that smaller set.

**Still to do — the 3x repeat run:**
Everything above was run **one time**. To know whether these numbers are
stable or just a lucky (or unlucky) single attempt, the full pipeline
needs to be run **three independent times** and the spread of results
compared. This is planned for **this evening**. The manuscript draft
already has a placeholder table built for exactly this (see
`docs/manuscript/reyna_case_report_manuscript.tex`, Section
"Run-to-run uncertainty analysis") — once the three runs are done, that
table gets filled in for real, and that's what the manuscript should
ultimately quote, not tonight's single-run numbers.

## 6. Where to look

- Full technical change log: `docs/CHANGES_SINCE_PR22.md`
- Full results write-up (the real numbers): `docs/reyna_publication_readiness_results_20260906.md`
- Manuscript draft (outline stage, not complete): `docs/manuscript/reyna_case_report_manuscript.tex`
- The PR with every commit: https://github.com/HafizAqilla/unified-vsd/pull/26
