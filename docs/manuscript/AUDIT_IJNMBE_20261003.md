# Pre-submission audit — VSD P1 manuscript for IJNMBE

**Date:** 2026-10-03.
**What was audited:** `revision_20260927/main.tex` and `supplement.tex`, as compiled in `VSD_Paper_Keisya_Hafiz (7).pdf` (33 pp.). The three generated noise sections were also read, as were the repository and the run outputs where the paper makes claims about them.

**Benchmark papers:**

| Short name | Paper |
|---|---|
| **Owashi** | Owashi et al. 2020, *PLOS ONE*: model-based LV pressure in aortic stenosis, 12 patients |
| **Lechuga** | Lechuga et al. 2026, *Compr Physiol*: patient-specific LPM in pulmonary hypertension, 25 subjects |
| **Bozkurt** | Bozkurt et al. 2022, *Cardiovasc Eng Technol*: LPM optimisation in pediatric dilated cardiomyopathy, 3 patients |
| **Das** | Das et al. 2026, *Simul Model Pract Theory*: CFD plus machine learning, coronary arteries |

**Supervisor comments:** these are included verbatim in §0. They are folded into the plan (§36).

**Section numbers:** "§2.x" and "§3.x" refer to the compiled PDF. Line numbers (`L123`) refer to `main.tex`. Wherever information cannot be found in the manuscript, it is marked **NOT REPORTED**.

---

## 0. Supervisor comments and how this audit handles them

> 1. *"yang post closure prediction, bahas keterbatasannya cuman ngeliatin perubahan tepat sesaat setelah ditutup. pengembangan kedepannya itu, bisa dikembangin dengan ngeliatin perubahan terhadap waktunya gimana."*
> 2. *"yang post closure calibration blum dimasukkin ke papernya. post closure calibration dibuat untuk bikin simulasi anak yang ada VSD kalau udah ditutup, simulasinya gmn (bisa diubah ubah parameter kayak resistansi, kapasitansi, dll nya buat ngeliatin respon anaknya gimana)"*
> 3. *"post closure dan pre closure cek pakai box plot"*
> 4. *"3 data running -? buat bandingin post and pre closure"*

| # | Comment | Assessment | Action (details in §36) |
|---|---|---|---|
| S-1 | Discuss that the prediction covers only the acute state, and name time-resolved change as future work | Correct and important. The post-closure record mixes two times: catheter pressures at about 5 min and echo volumes at about 1 day. The model holds every chamber property fixed, so by construction it represents the instant after closure. This is currently only implied (Limitations bullets 4 and 6). It is also the most likely reason the predicted left-heart unloading is absent (Yuan 2025: LV dimensions normalise mostly within the first month). | **Text only (P1).** Add a limitation stating that the prediction is an acute, fixed-property prediction. Add a future-work paragraph: follow-up echo at about 1, 3 and 6 months, and time-varying chamber parameters (reverse remodelling of $E_{B,LV}$, $V_{0,LV}$) calibrated to each follow-up. |
| S-2 | Add a post-closure calibration so the closed-VSD child can be simulated, and run what-if changes in R, C, etc. | Scientifically valuable *if it is clearly separated from the prediction*. Your integrity rule 2 (post-closure data never used for calibration, priors or tuning of the *prediction*) still holds, because the prediction is already frozen. A post-closure calibration then answers a new question: **which parameters had to change at closure** for the model to match the post-closure record? That difference is a direct measure of what the fixed-property prediction misses (S-1), and it makes the chamber parameters identifiable, because post-closure volumes exist. The existing post-surgery run (`results/runs/20260906_212347_reyna_post_surgery`) **cannot be used**: it was made with the superseded `lundquist_bsa` exploratory scaling and an older target set (manifest: `ScalingRunMode: exploratory`, RMSE 0.146). | **New computation (P1).** Re-run the post-closure calibration with the current pipeline (weight scaling, same objective form), *after* freezing the prediction. Use the post-closure targets, and let chamber parameters ($E_{B,LV}$, $E_{A,LV}$, $V_{0,LV}$, $E_{B,RV}$) be free. Report it in its own subsection, "Post-closure calibration (fit, not prediction)", with a table of Δparameter (post-fit vs pre-fit). Then run one-at-a-time sweeps of R and C around the post-closure fit as the "what-if" demonstration (1 figure). Pre-register the parameter set before running it, and write it in the run protocol. |
| S-3 | Compare pre and post with box plots | Good idea, but a box plot needs a distribution. Box plots must come from the **50-replicate measurement-noise simulation**: one value per replicate for each pre-closure fitted output and each post-closure predicted output, with the measured values overlaid. This is the same design Lechuga uses for parameters across groups (their Fig. 7). | **Depends on the noise batch (P1).** Replace Fig. 8 with paired box plots, pre (calibrated) vs post (predicted), for the 9 pressures and flows and the 4 volumes, each with a measured-value marker. Add a second panel with box plots of the 7 calibrated parameters (ratio to reference). |
| S-4 | Use the "3 data running" to compare pre and post | **Do not use the three seed runs (101/202/303) for this.** They return *identical* parameters and outputs (Supplement S1: the selected start is the deterministic, clinically seeded one in all three runs). A box plot of three identical values is a line, and a reviewer would read it as a misunderstanding of variability. The 18 multi-start end-points (6 starts × 3 runs, RMSE 0.091–0.151) are a legitimate but secondary spread: they show optimiser variability, not uncertainty. | Explain this to the supervisor. Use the noise replicates (S-3). Optionally add the 18 multi-start end-points to the Supplement as an "optimiser spread" panel. |

---

## 1. Editorial triage

### A. Scope fit
The manuscript is **a patient-specific modelling study with a prospective-style prediction test**, built from established components:
- a closed-loop 0D model adapted from Ferrero 2026;
- Zhang 2019 allometric scaling;
- PCE-Sobol' screening (UQLab);
- multi-start `fmincon`;
- Monte Carlo measurement-noise propagation.

None of these components is new numerically. The only author-introduced mathematics is:
- the smoothed signed orifice law, Eq. (6);
- the composite objective, Eq. (10);
- the bookkeeping constraint, Eq. (8).

**Risk:** yes, it currently reads as *"a standard 0D model, standard GSA and standard optimiser applied to one child"*. Three reasons:
1. The Introduction frames the gap clinically (residual shunt, aortic regurgitation, operator dependence). The model cannot address residual shunt or aortic regurgitation.
2. The numerical-method content that IJNMBE values (verification, convergence, surrogate accuracy, identifiability, UQ) exists but is scattered or hidden:
   - identifiability is in Supplement S2;
   - the noise UQ is 1/50 complete;
   - no solver-tolerance or convergence evidence is given.
3. n = 1.

**What must be made clearer to avoid that reading:** the paper's real question is computational and methodological. *Given a sparse pre-closure catheter record, which post-closure quantities can a calibrated 0D model predict, how large is the prediction uncertainty induced by measurement error and non-identifiability, and which quantities are structurally unpredictable because no pre-closure measurement informs them?* The manuscript already contains a strong instrument for this: the GSA-informed "predictability map" (§3.1 → §4.2 chamber function → §3.4 volume errors). It is not stated as the contribution.

### B. Novelty
- **Apparent novelty, in one sentence:** a pediatric VSD 0D model calibrated to one child's pre-closure catheterization and tested, without refitting, against that child's withheld post-closure record.
- **Claimed novelty** (Intro, last sentence): the "first workflow that calibrates a pediatric VSD model to a child's pre-closure measurements and then tests its prediction against the same child's post-closure measurements, while reporting both the reproducibility and the identifiability of the calibrated parameters."
- **What the evidence demonstrates:**
  - (i) The pre→post forward test: genuinely demonstrated, and the strongest element.
  - (ii) "Reproducibility": **not demonstrated in a meaningful sense.** The three seeds give identical output *because the selected start is seed-independent* (Supp. S1). That is a property of the selection rule, not evidence of robustness.
  - (iii) Identifiability: demonstrated locally (cond(S) = 61.9), but it now lives in the Supplement. The single noise replicate (parameters shift by up to +141% under a −2.7σ $Q_s$ perturbation) contradicts any reading of the parameters as "identified".
  - (iv) The predictability split between pressures (predicted to within 1–14%, except PA systolic) and volumes (unconstrained): demonstrated, but buried in the Discussion.
- **Discrepancy:** the claim leans on reproducibility and identifiability. The evidence supports *prediction plus a quantified limit of predictability*.
- **Strong contributions that are currently buried:**
  - (a) shunt-phase decomposition (40% of shunt volume crosses after aortic-valve closure) and the loss of LV isovolumic phases. This is a mechanistic output no catheter gives.
  - (b) the predicted absence of left-heart unloading, explained by the simultaneous fall in heart rate (per-beat venous return −4% vs per-minute −13%). This is a model-derived, testable hypothesis.
  - (c) measurement-method disagreement (Fick → echo) as a sufficient explanation of the measured $Q_s$ fall and SVR rise.
- **Missing comparisons needed to establish novelty:**
  - against Colebank 2026 (pediatric VSD, uncalibrated) and Ferrero 2026 (0D shunt interpretation, no prediction). These are cited; state the difference in one table row each.
  - against a **naïve predictor**, e.g. "post = pre measured" or "post = pre model with Qp set to Qs". Without a baseline predictor a reviewer cannot judge whether 6.9 mmHg RMSE is informative. **This is cheap and high-leverage (P0).**

### C. Realistic desk-rejection risks
1. **Single patient**, with no statement of why the other enrolled patients (the config contains a cohort file) are not included. IJNMBE editors routinely desk-reject n = 1 application papers without a methodological contribution.
2. **No numerical verification** at all: no solver-tolerance study, no cycle-count convergence, no conservation check beyond the 0.03 mL closure statement, no PCE convergence. For IJNMBE this is a scope-level omission.
3. **The GSA training simulations are not at periodic steady state** (all 256 stopped at the 10-cycle cap, §3.1 last sentence). The parameter selection, which drives everything downstream, rests on transient outputs.
4. **Abstract is stale and contradicts the body:** 5.46%/6.65%, "3 of 13", R.SC/C.SVEN/vsd.Cd code names, "within 10%" acceptance language that §2.9 says is not used.
5. **Main-text UQ is placeholder** (n = 1 of 50; Tables 8–9 have "–" columns).
6. **Template is not the journal's:** `aaj2026.sty` ("AI Agent Journal").
7. **Red placeholders:** CRediT, funding, COI, acknowledgments.
8. **Code "available on request" only.** IJNMBE expects a data-availability statement and strongly favours archived code. *And see §26: the public GitHub repository already contains the patient's first name and full clinical record, which conflicts with the de-identification statement in §2.1. This must be resolved before submission.*

---

## 2. Benchmark comparison

### What the benchmarks do that we do not (practice-level, not topic-level)

| Practice | Owashi | Lechuga | Bozkurt | Das | Ours | Topic-driven or genuinely better reporting? |
|---|---|---|---|---|---|---|
| Problem → gap → aim in ≤ 2 paragraphs | yes (aims by end of p. 2) | yes (gap and aim in 2 paragraphs) | yes | no (long literature table) | 5 paragraphs, aim at end of P5 | Better reporting: tighter funnel |
| Explicit contribution list | Discussion opens with 4 numbered contributions | key findings 1–4 | – | numbered conclusions | "Five findings" in §4.1 (good) | We match |
| Every parameter with definition, formula/value, reference | Table, plus S1 File | **Table 2: definition / formula / reference for every parameter** | Tables 1, 2, 4, 5: bounds before and after, fixed values | Tables 2–4 | Seeding table (Table 3) is good; nominal table (App. B) has **no source column** ("locally harmonized … in the released code") | **Better reporting: we must add provenance** |
| Bounds of fitted parameters and their origin | – | bounds rule stated (1%–600% of baseline; timing fractions) | **explicit ub/lb table and update rule** | – | bounds given (Table 5), **rule NOT REPORTED** | Better reporting |
| Solver, tolerance, cycles | – | `ode15s`, 40 cycles, last cycle | `ode15s`, max step 1e-3 s, rtol 1e-3, Simulink R2017a | algorithm, Δt, residual 1e-6, hardware, runtime | `ode15s`, rtol 1e-5, atol 1e-6, MaxStep T/20, SS criterion (good) | We are *more* complete than three benchmarks on settings, but **no convergence evidence** |
| Convergence / independence test | – | – | – | **mesh-independence plot; cycle-independence < 0.01%** | none | Das shows the IJNMBE-expected practice. For 0D the analogue is tolerance and cycle independence |
| Hardware / runtime | – | – | – | yes (Ryzen 5600X, 32 GB, 4–5 h/case) | **NOT REPORTED** (84 min/noise replicate is known) | Better reporting, cheap |
| Waveform-level comparison | **12 measured vs simulated LV pressure traces, e% per patient** | **PAP and RVOT-flow traces, R² per subject** | waveforms shown, not compared | – | only horizontal lines for scalar readings | Topic/data-driven: we have no stored traces. If the cath lab stored pressure traces, using them would be a large gain |
| Test on data not used for fitting | **Monte-Carlo cross-validation, 6/6 split** | none (fit only) | none (fit only; "<10%" is fit error) | literature validation | **held-out post-closure state** | **Our strongest element; better than Lechuga and Bozkurt** |
| Uncertainty of fit | mean ± SD band of simulated curve (Fig. 7) | mean ± SD across subjects | – | 5-fold CV | 1/50 replicates | Must finish |
| Parameter distribution figure | – | **box plots per group (Fig. 7)** | – | – | none | Supervisor S-3 matches this |
| Nomenclature list | – | – | **yes** | **yes** | none (abbreviations in prose) | Better reporting (IJNMBE papers often include one) |
| Supporting data | **S1/S2 Data XLSX** | "on request" | – | Supp. Table S1 full dataset | "on request" | Better reporting |
| Limitations | short (n, filling pressure) | 3 numbered | none | 1 sentence | 10 bullets, specific (**better than all four**) | Keep |
| Discussion vs physiology | moderate | good, cites remodelling literature | parameter values vs literature values (useful) | thin | **very strong (Tables 10–11)** | Keep; trim repetition |

### Required comparison table (20 dimensions)

| # | Dimension | Our manuscript | Benchmarks | Gap | Required fix |
|---|---|---|---|---|---|
| 1 | Title | Copies the Lechuga pattern almost word for word ("Patient-Specific Lumped-Parameter Model for Quantifying … and Predicting …"); "Quantifying shunt flow" is not a result (the model under-predicts the measured shunt by 25%) | Owashi, Bozkurt: method plus condition, plain | Over-claim and near-duplication of the Lechuga title | §4 |
| 2 | Abstract | Stale numbers, code names, 10% acceptance language, no post-closure numbers that match the body | All give n, method, main metric | Contradicts body | §5 (user's text; flagged only) |
| 3 | Intro structure | 5 paragraphs; problem paragraph motivates residual shunt and AR, which the model cannot predict | Funnel problem → modelling → gap → aim | Motivation–capability mismatch | §6 |
| 4 | Gap definition | "few studies calibrated a VSD model to individual children"; "identifiability rarely reported" | Lechuga: "PV loop requires invasive equipment → model from routine data" | Gap is about *calibration existing*, not about *what it can predict* | §6 |
| 5 | Novelty statement | Last Intro sentence, includes reproducibility (trivial) and identifiability (Supplement) | Owashi: 4 contributions in the Discussion | Claim ≠ evidence | §1B |
| 6 | Mathematical formulation | Chambers, valves, shunt, appendix ODEs: good. Objective Eq. (10) has 3 **undefined** terms ($J_{sys}$, $J_{shape}$, $J_r$) | Lechuga Eq. 6 fully explicit | Objective not reproducible | §8 |
| 7 | Numerical method | Settings complete | Das adds hardware and runtime | Hardware and runtime missing; FD gradient settings missing | §10 |
| 8 | Parameter reporting | Scaled, seed, bounds, calibrated (Table 5): excellent. Nominal values: no source | Lechuga, Bozkurt: source per row | Provenance | §13 |
| 9 | Initial/boundary conditions | Table 3: good (closed loop, so no BC) | Lechuga Table 2 | Small gaps (LA defaults' source) | §8 |
| 10 | Numerical convergence | Steady-state criterion only | Das: mesh and cycle independence | **None demonstrated** | §11 |
| 11 | Validation | Held-out post-closure: strong | Owashi cross-validation | No baseline predictor; no uncertainty band on prediction yet | §12 |
| 12 | Sensitivity/UQ | GSA good; UQ 1/50 | Owashi Morris; Lechuga local; Das CV | GSA on truncated sims; PCE accuracy per output not shown; no CI on indices | §15 |
| 13 | Results structure | Ordered: GSA → pre → hemodynamics → post → PV → noise | Similar | Results contain some Methods (§3.1 last paragraph, §3.4 protocol-form) | §16 |
| 14 | Figure quality | Clean ggplot; legends waste space; Fig. 8 placeholder; parity panels with n = 1–2 | Lechuga: overlaid traces, boxes | §18 | §18 |
| 15 | Discussion depth | Excellent physiology mapping | – | Repeats Results numbers; one speculative closing sentence | §21 |
| 16 | Limitations | 10 specific bullets | short | Missing: acute-only (S-1), truncated GSA, FD-gradient noise, single operating point of noise, σ choices | §22 |
| 17 | Reproducibility | Scripts named; "on request" | Owashi data XLSX | Objective terms, bounds rule, stage composition, prior confidences NOT REPORTED | §25 |
| 18 | Data/code availability | On request | PLOS: in SI | Archive needed; **privacy conflict in public repo** | §26 |
| 19 | Conclusion | Number-dense, honest; ends with an uncomputed "next step" | Owashi: constrained | Contains a generated noise sentence based on n = 1 | §23 |
| 20 | Information density | High, sometimes too high (5 numbers per sentence in §3.4 paragraph 2) | Lechuga: 1–2 numbers per sentence, rest in tables | Move lists of numbers into tables | §28 |

---

## 3. Section-by-section audit

Format: **[severity] location — problem → why it matters → missing information → action → new computation? / text-only?**

### Title, Abstract, Keywords
- **[CRITICAL] Abstract (L80).**
  - **Problem:** 4 stale or contradicting statements (listed in the L66–77 comment). Also: "multi-start calibration … within 10% range in every run, with rejection and acceptance" contradicts §2.9 ("No acceptance threshold is applied").
  - **Why it matters:** desk rejection.
  - **Action:** the authors rewrite it (it is your text). The numbers to use are in §5. *Text-only.*
- **[MODERATE] Keywords.** "parameter identifiability" is now Supplement material. Replace it with "uncertainty quantification" or "model prediction". *Text-only.*

### 1 Introduction (user's text; recommendations only)
- **[MAJOR] P2 (L93).**
  - **Problem:** it motivates residual shunt (36.4%) and aortic regurgitation (3.1%). The model has no aortic-valve regurgitation and treats closure as complete ($A_{VSD}=0$).
  - **Why it matters:** a reviewer will ask why the model does not simulate residual shunt.
  - **Action:** either (a) re-aim P2 at "the hemodynamic state after closure is not predicted before the procedure" (pressures, flows, chamber loading), or (b) add one sentence acknowledging that residual shunt and AR are outside scope.
  - Cheap extension: a residual-shunt scenario is a one-parameter change ($A_{VSD}$ = fraction of pre value) and could go in the post-closure what-if section (S-2).
- **[MAJOR] P5 aims (ii)–(iii).**
  - **Problem:** "determine their reproducibility by different seeded runs" and "audit the identifiability" are now Supplement items, and the seed result is trivial (§1B).
  - **Action:** re-weight the aims toward (a) prediction and (b) the uncertainty of the prediction (noise simulation). Keep identifiability as supporting evidence.
- **[MODERATE] P4.** Literature is enumerated one study at a time (by your choice). Add *what each study could not do* relative to this paper. For example: "Colebank et al. … without calibration to an individual child"; "Ferrero et al. … interpreted catheterization data but did not predict a post-intervention state". The Discussion §4.6 already has this logic; one clause per study is enough.
- **[MINOR] Grammar.** Offered as a marked list on request (rule: not edited silently).
  - "The shunt itself depended on" → "depends on"
  - "which created patient-specific risk factors"
  - "These findings created a need"
  - "patient-specifics"

### 2.1–2.2 Study population and clinical data
- **[MAJOR] §2.1, L110.**
  - **Problem:** "P1 was the first enrolled patient with a complete pre- and post-closure record."
  - **Missing:** how many patients were enrolled by June 2026, and how many have pre-closure records only. **NOT REPORTED.**
  - **Why it matters:** single-patient selection is the first reviewer question. If other patients have pre-closure catheter data, a pre-closure calibration across them (no prediction) would multiply the evidence for the screen and the calibration.
  - **Action:** report the count and the reason each was not modelled. Text-only, or new computation if further patients are added.
- **[MAJOR] §2.2, L116.**
  - **Problem:** the Fick method details are **NOT REPORTED**: measured or assumed VO₂, saturation sites (SVC/mixed venous; PA; PV or assumed), hemoglobin.
  - **Why it matters:** the Discussion attributes the $Q_s$ fall to the assumed-VO₂ bias (Ralston 2023). That argument requires stating that VO₂ was assumed.
  - **Action:** state the VO₂ source and the saturation sites. *Text-only (needs the protocol form).*
- **[MODERATE] §2.2.** The SVR formula is not stated. The numbers match $(\mathrm{MAP}-\mathrm{RAP})/Q_s$: (71.33−5)/3.423 = 19.38 pre and (75−5)/2.54 = 27.56 post. State it once.
- **[MODERATE] §2.2.** The post-closure echo $Q_s$ method (LVOT VTI × area? Simpson SV × HR?) is **NOT REPORTED**. It matters because the Discussion's main measurement argument rests on it.
  - Internal check: the measured post LV SV of 22.33 mL × 108 bpm = 2.41 L/min, while the measured RV SV of 18.5 mL × 108 = 2.00 L/min. The two disagree with each other and with echo $Q_s$ = 2.54.
  - Report this. It **strengthens** the measurement-inconsistency argument in §4.4.
- **[MINOR] Table 1.** "Systolic/diastolic pressure, cuff — Use: reported". Use the defined code "C".

### 2.3 Mathematical model
- **[MODERATE] Fig. 2 vs text.**
  - The schematic shows one $C_{PC}$, but the equations (App. A) use $C_{PCOX}+C_{PCNO}$.
  - $R_{SC}/2$ is coloured as calibrated (correct, via the $R_{sys}$ scale), but the legend says "calibrated ($S_T \ge 0.10$)". The scale-group logic is not shown.
  - **Action:** add "(via resistance scale)" to the legend, and draw $C_{PCOX}$ and $C_{PCNO}$.
- **[MODERATE] Eq. (5).**
  - **Problem:** $R_{open}$ is given unscaled, but Table 2 scales it with α = −0.50.
  - **Action:** say "nominal" and give the scaled P1 value (4.0e-3 × (13.4/70)^−0.5 = 9.14e-3 mmHg·s/mL).
- **[MODERATE] App. A.** Pressures are clamped at −5 mmHg. The clamp is not smooth, which contradicts the stated reason for smoothing the valves. State the clamp form (hard max? smooth?) and whether it was active in any reported simulation.
- **[MINOR] Notation.** The appendix table uses $E_{A,LV}$, $T_{c,LV}$, $\epsilon_{valve}$; the text uses $E_{\mathrm{A},ch}$, $T_{vc}$, $\varepsilon_v$. Unify them (§27).

### 2.4 Scaling and initialization
- **[MAJOR] Table 2 / L271.**
  - **Problem:** the heart-rate exponent −0.30, the class assignments and the inertance exponent come from different sources. The validity regime is not stated, although your own run manifest says `AgeValidityRegime: preschool_extrapolated`: the 3.2-yr age lies between infant and school-age anchors.
  - **Why it matters:** a reviewer will ask whether Zhang's exponents were derived for this age range.
  - **Action:** one sentence stating the age/weight range of Zhang's derivation, and that P1 lies inside or outside it. *Text-only.*
- **[MODERATE] Eq. (8).** It is ambiguous whether $\sum_{i\neq SVEN} V_i(0)$ includes the four chambers. State "over all 4 chambers and 5 vascular compartments other than SVEN", and give the resulting $V_{0,SVEN}$ value for P1.
- **[MODERATE] Table 3.** The LA/LV filling defaults (6/6/8 mmHg) have no source. Cite one, or call them "assumed".

### 2.5 Numerical solution
- **[CRITICAL for IJNMBE] No verification evidence.** See §11. *New computation (cheap).*
- **[MODERATE] Missing items.**
  - Initial step size.
  - Jacobian (numerical? sparsity pattern supplied?).
  - Runtime per simulation and per calibration.
  - Hardware.
  - How "peak values of the states" are compared for steady state (max over the batch? per state?).

### 2.6 Global sensitivity analysis
- **[CRITICAL] Training simulations truncated at 10 cycles (L355).**
  - **Problem:** all 256 hit the cap (§3.1). The indices therefore describe a transient.
  - **Why it matters:** the retained-parameter set, and hence calibration and prediction, could differ if computed at steady state. The steady-state outputs and the cycle-10 outputs may differ most for slow compartments (venous volumes, RA pressure), which is exactly where $E_{B,RV}$ and $C_{SVEN}$ were found influential.
  - **Action (new computation, moderate):** for a random subset of ≥ 20 of the 128 design points, run to full steady state and report the relative difference between cycle-9–10 and steady-state outputs for the 9 targets. If it is small (< about 2%), state it; this closes the issue. If not, rerun the screen at steady state, which is expensive.
- **[MAJOR] PCE accuracy and index uncertainty.**
  - **Problem:** LOO errors 0.043–0.155 are reported only as a range. Which outputs have 0.155? **NOT REPORTED.**
  - **Action:** give a per-output LOO table (Supplement), and bootstrap confidence intervals on $S_T$ (UQLab supports bootstrap PCE). Retention decisions near the threshold ($C_{SAR}$ 0.12, $R_{SVEN}$ 0.10, $E_{A,RV}$ 0.080, $E_{B,LV}$ 0.073) are not robust without them.
  - A PCE convergence check (N = 64 vs 128) would also answer "is N = 128 enough". *New computation: cheap (re-uses existing training data, except N > 128).*
- **[MAJOR] Group index = max of member indices (L355).**
  - **Problem:** the members of a resistance group were sampled independently but are calibrated as one multiplier. The relevant quantity is the **group (closed) total index** of the group, which is ≥ the max of the members' individual indices. The max underestimates the group's influence.
  - **Action:** state this as a conservative choice, or compute the group index from the PCE (cheap; analytic from the coefficients).
- **[MODERATE] Input ranges are heterogeneous.**
  - **Problem:** resistances span [0.4, 2.5], $C_d$ [0.8, 1.2]. The ranking depends on these ranges, and their justification is **NOT REPORTED**.
  - **Action:** one sentence with the rationale (for example, $C_d$ is narrow because it is already seeded from the measured gradient and shunt flow).
- **[MODERATE] "Not in the calibration parameter set of the code" (L355, L682, L754).** This reads as a software artefact. The scientific statement is: "venous compliances and atrial unstressed volumes were excluded from the candidate set before the screen because …" (give the reason). If there was no scientific reason, say so plainly. That is also a limitation, because $C_{SVEN}$ exceeded the threshold.

### 2.7 Calibration
- **[CRITICAL for reproducibility] Objective function, Eq. (10).** The following are **NOT REPORTED**:
  - the functional forms of $J_{sys}$ (how are RA, MAP and $Q_s$ "coupled"?), $J_{shape}$ and $J_r$;
  - the prior confidences $c_j$;
  - the physiological ranges used by $J_r$;
  - the definition of the "hard band";
  - the bound-proximity penalty form $J_b$ (quadratic? log barrier?).

  A graduate student cannot rebuild the objective. **Action:** write all terms in Supplementary Section S3 with their numbers. *Text-only (read from code).*
- **[CRITICAL for reproducibility] Stage composition.** The following are **NOT REPORTED**:
  - which of the 7 parameters are free in stages A–F;
  - "up to five retained parameters jointly": which five;
  - "systemic-output refinement": which targets.

  **Action:** add a stage table (stage × free parameters × targets × algorithm × λ settings) to the Supplement.
- **[MAJOR] Bounds rule NOT REPORTED.** Two parameters end on their lower bounds (Table 5). A reviewer will ask (a) where the bounds came from and (b) whether widening them changes the fit and the prediction. **Action:** state the rule. Optionally run one calibration with the two lower bounds relaxed and report the effect on PA systolic pressure and the prediction. *New computation: 1 calibration, about 1.5 h.*
- **[MAJOR] Optimiser convergence.** `fmincon` exit flags, first-order optimality and the finite-difference step are **NOT REPORTED**. With rtol 1e-5 and a 0.1 mL steady-state tolerance, the objective is numerically noisy, and FD gradients can be dominated by that noise. The 300-evaluation / 40-iteration caps suggest stages end on the cap, not on convergence. **Action:** report the exit reason per stage for the selected start (from `calib_diagnostics`). Show that tightening ODE tolerances does not move the optimum (part of §11).
- **[MODERATE] Selection rule** uses the 10% setting. This is correctly declared as internal, but it means the *reported* optimum is not the lowest-RMSE optimum (Supp. S1: 0.0909 vs 0.1047). Keep the honest Supplement statement and add one sentence in the main text.

### 2.8–2.9 Post-closure prediction; comparison metrics
- **[MAJOR] σ choices (L387).**
  - **Problem:** $\sigma(Q_s)$ = 0.50 L/min (14.6%) pre (Fick), but **5%** post (echo). Echo cardiac output is not more precise than Fick, so the 5% inflates the post-closure $z$ (+4.97) and $\chi^2$.
  - **Action:** justify each σ by source (device or method literature, 2021+). Otherwise use one method-specific σ table and state it. If σ(post $Q_s$) changes, $\chi^2/N$ changes, so recompute it (no simulation needed).
- **[MAJOR] $\chi^2$ degrees of freedom.** $N-p = 9-7 = 2$ treats the 9 targets as independent. But $Q_s$ is derived from $Q_p$ and $Q_p/Q_s$, and the objective also fits the shunt and $J_{shape}$. Add one sentence acknowledging the dependence, or drop the reduced $\chi^2$ and keep $\chi^2/N$.
- **[MAJOR] L381.** "Differences … are therefore attributable to the model together with its pre-closure parameters" is wrong by omission: they are also attributable to measurement error and to the 5 min / 1 day timing (§4.4 says so). Fix the sentence.

### 2.10 Measurement-noise simulation
- **[MAJOR] Perturbation design.** Noise is applied to pre-closure data only, and the post-closure *reference* is unperturbed. The resulting spread is therefore "prediction uncertainty from pre-closure measurement error". Comparing it with measured post values without post-closure σ mixes two things. **Action:** in the post-closure comparison, show both the predictive spread (from replicates) and the post-measurement σ (error bars), as Owashi's Fig. 7 overlays mean ± SD on the measurement.
- **[MODERATE] "Could not start" replicates (about 3/50).** Pre-specified and honest, but the reporting rule (exclusion from statistics) must be stated before results exist. It is: keep it.

### 3 Results
- **[MAJOR] §3.4 paragraph 2.** 20+ numbers in one paragraph that duplicate Table 6. Keep 3–4 headline numbers (systolic RFA +1.2%, MAP −6.4%, PA systolic +40.7%, RVEDV +42.1%) and point to the table.
- **[MAJOR] §3.4 / §4.1 / §4.7.** "Every measured pressure change" is **wrong**:
  - systolic RFA was measured 100 → 89 (fell), model 88.1 → 90.1 (rose);
  - diastolic RFA was measured 57 → 68 (rose), model 52.8 → 52.0 (fell).

  **Action:** correct the three locations (L650, L742, Conclusion if present) and add both rows to Table 11 as "not consistent". *Text-only.*
- **[MODERATE] §3.2 last sentence (L436).** "26.0% and 36.4% above the model": ε is relative to the measurement, so the model is 26.0% and 36.4% **below** the cuff readings. The same direction error appears at L730 ("11.9% above the model").
- **[MODERATE] §3.2 vs Table 4.** Shunt cycle mean 0.495 (L508, Fig. 4E label) vs 0.496 (Table 4, L434). Use one value from the export.
- **[MODERATE] §3.3 vs §4.4.** The shunt carries 13.6% of pulmonary flow (model) vs "16% of pulmonary flow" (measured: 0.664/4.087). Label which is which.
- **[MODERATE] §3.2.** "Each fitted target was controlled by one or two parameters", yet $Q_s$ has three parameters ≥ 0.14 (Table 10). Soften to "one dominant parameter".
- **[MINOR] §3.2 Fig. 4 text.** The diastolic minimum of shunt flow is "above 0.10 L/min" here and "0.11 L/min" in §4.3. Use one number.

### 4 Discussion
- **[MAJOR] §4.1 first finding (L650).** "The dominant parameter of each target was the determinant that established physiology assigns" contradicts Table 10 (7 consistent, 2 partly). Fix it.
- **[MODERATE] Repetition.** The shunt phase split (6/37/40/17%) and per-beat venous return (30.6 → 29.4 mL) are each stated twice (Results and Discussion). Keep them in Results and refer to them in the Discussion.
- **[MODERATE] §4.6 last sentence.** "A dedicated post-closure calibration is the natural upper bound …" refers to an analysis that is not in the paper. Either delete it, or **make it true** by doing S-2. It then becomes a strong quantitative statement (prediction error vs post-fit error).
- **[MODERATE] §4.6.** "those pressures were never fitted" has an unclear referent. Say "the seven post-closure catheter pressures".
- **[MODERATE] Section title "Sources of the largest pre-closure errors".** It covers post errors too. Rename it "Sources of the largest errors".
- **Keep:** Tables 10–11 and the "where the measurement does not agree with physiology" paragraph. These are the most persuasive parts of the paper.

### Conclusion
- **[MAJOR] It includes `noise_conclusion.tex` based on n = 1.** Replace it once the batch is done. Until then, the conclusion should not carry a UQ claim.
- **[MODERATE] Last sentence** ("Adding a pre-closure echocardiographic volume measurement … is the concrete next step") is good. Add the time-course point (S-1).

---

## 4. Title audit

Current: *Patient-Specific Lumped-Parameter Model for Quantifying Shunt Flow and Predicting Post-Closure Hemodynamics in Pediatric Ventricular Septal Defect*

- **Near-identical to the Lechuga title** ("Patient-Specific Lumped-Parameter Model for Quantifying Vessel-Specific Remodeling and Predicting Right Ventricular Function in Pulmonary Hypertension"). An editor who knows that paper will notice it.
- **"Quantifying shunt flow"** over-claims: shunt flow is a measured input (Fick), and the model under-reproduces it by 25%.
- **"Pediatric ventricular septal defect"** implies a population; the paper has one child.
- Missing: the prediction-test character and the uncertainty element.

Proposals:
1. **Conservative/direct:** *Pre-closure calibration and post-closure prediction of hemodynamics with a lumped-parameter model in a child with ventricular septal defect*
2. **Method-focused:** *Sensitivity-guided calibration and measurement-uncertainty propagation in a lumped-parameter model of ventricular septal defect: a pre- to post-closure prediction test*
3. **Biomedical-impact-focused:** *How much of the circulation after ventricular septal defect closure can a pre-closure catheterization predict? A patient-specific lumped-parameter study*

---

## 5. Abstract audit (your text — not rewritten)

| Required element | Present? | Note |
|---|---|---|
| Problem/context | yes, too long (3 sentences) | Cut to 1 |
| Gap | partly ("few studies … patient-specific") | Make it "no study has tested a calibrated VSD model's prediction of the post-closure state" |
| Objective | implicit | State it |
| Method | yes, but with "zero dimensional … closure represented by eliminating shunt pathway" | Add: 7 of 14 parameters by PCE-Sobol', 6-start staged optimisation, 50-replicate noise MC |
| Setup | missing: n = 1, 3.2-yr, 13.4 kg, perimembranous restrictive, Qp/Qs 1.19 | Add |
| Validation | stale: "3 of 13 within 10%" | Replace |
| Main quantitative result | stale (5.46/6.65%) | Replace |
| Interpretation | "digital twins … pediatrics and adults" | Over-reach (adults not studied) |
| Undefined abbreviations / code names | R.SC, C.SVEN, vsd.Cd | Use symbols or words |

**Numbers that should appear:**
- RMSE of relative error, 9 targets: 0.274 → 0.105.
- 7 of 14 parameters retained; dominant ones: systemic resistance, $C_d$, $E_{B,RV}$.
- Post-closure, without refitting:
  - catheter-pressure RMSE 6.9 mmHg vs 5.2 mmHg pre;
  - systolic/mean arterial pressure +1.2% / −6.4%;
  - largest errors: PA systolic +40.7% (carried over) and RVEDV +42.1% (unconstrained).
- Noise result: CV or 95% interval of the predicted volumes (after the batch).
- One mechanistic finding: no predicted LV unloading, because HR fell at the same time.

**Suggested architecture (≤ 250 words):**
1. Background and gap, 2 sentences.
2. Aim, 1 sentence.
3. Patient and data, 1 sentence.
4. Methods, 2 sentences.
5. Pre-closure results, 1 sentence.
6. Post-closure prediction, 2 sentences.
7. Uncertainty, 1 sentence.
8. Conclusion, constrained to what is predictable and what is not, 1 sentence.

---

## 6. Introduction audit (functions per paragraph)

| Paragraph | Function | Verdict |
|---|---|---|
| P1 | Epidemiology + VSD pathophysiology + spontaneous closure | Useful; trim the regional statistics (Asia, Indonesia) to one sentence unless the journal audience needs them; the pathophysiology sentences are needed later (restrictive vs nonrestrictive) |
| P2 | Clinical outcomes of closure (residual shunt, AR, weight) | **Function mismatch:** motivates outcomes the model does not simulate. Re-aim to "the post-closure hemodynamic state is not known before closure" |
| P3 | Imaging limits → patient-specific models → digital twins → 0D | Good; the core motivation. Could absorb P2's last sentence |
| P4 | Literature list | Enumerative by design; add "what it could not do" clauses |
| P5 | Gap + aims + novelty | Gap is stated as "not combined into a single workflow". Weak, because combining established tools is not a gap in itself. Stronger gap: no published VSD model has been tested prospectively on a post-closure state of the same patient, and no study has quantified which post-closure quantities a pre-closure record can constrain |

By the end of the Introduction, a reader knows *what was done* but not *what previous studies could not do*. **Missing logic:**
1. Colebank 2026 and Ferrero 2026 show that 0D models reproduce VSD physiology.
2. Neither tested a prediction against post-intervention data in the same patient.
3. A prediction test is the only way to learn which post-closure quantities are predictable from pre-closure data.
4. That answer depends on which parameters the pre-closure data inform (GSA plus identifiability) and on the measurement error (UQ).
5. This study does all three in one child.

---

## 7. Literature and claim audit (selected; full classification in §31)

| Claim | Location | Category | Issue |
|---|---|---|---|
| "Fick output computed with an assumed oxygen consumption has a mean bias of +27.5%" | L725 | literature | Correct use only if P1's VO₂ was assumed; **NOT REPORTED** (see §2.2) |
| "LV dimensions in children normalize after transcatheter VSD closure, most of the change occurring in the first month" | L723 | literature (Yuan 2025) | Good; also supports S-1 |
| "LVEF falls early after VSD repair because of loss of low-pressure outlet" | L689 | literature (Safa 2023) | Safa is about *surgical repair*; transcatheter closure differs (no bypass). Say so |
| "RV systolic pressure 30% of LV is the pressure relation of a restrictive defect" | L678 | established fact | Needs a citation (restrictive definition: RVSP < ~ 50% LVSP, or gradient > ~ 64 mmHg); Rao 2024 likely covers it — verify |
| "Retained parameter set … should be re-screened for each patient" | L678 | **interpretation/speculation** | Presented as a recommendation; fine if labelled "we expect" |
| "C_d has a direct clinical counterpart" | L678 | interpretation | OK |
| "Closure does not reduce systemic output" | L725 | established physiology, no citation | Add a citation |
| "A proximal and distal PA compartment is the direct structural remedy" | L728 | hypothesis | Label it as such ("would be expected to") or test it |
| Colebank 2026 | refs | **arXiv preprint** | Check whether it has been published; IJNMBE reviewers dislike load-bearing preprints |
| Missing literature categories | – | – | (a) VVUQ framework for computational models (e.g., ASME V&V 40-2018; note it is outside your 2021–2026 rule — decide whether a standard is exempt); (b) profile-likelihood or Bayesian identifiability in cardiovascular 0D models after 2021; (c) 0D studies of post-intervention prediction in other CHD (Fontan, ASD closure) as the closest methodological precedent for a pre→post test; (d) Fick vs echo CO agreement in children (Ralston 2023 and Schramm 2023 are cited; good) |

---

## 8. Model definition and mathematical reproducibility

**Can a graduate student rebuild the model from the paper alone?** The ODE system, yes, with one gap. The calibration, no.

**Missing-information checklist**

| Item | Status |
|---|---|
| State vector, chamber ODEs, vascular ODEs | Reported (Eq. 1–2, App. A) |
| Elastance law, activation | Reported; atrial activation "same form" — parameters in App. B (OK) |
| Valve law + values | Reported; $R_{open}$ scaled value missing |
| Shunt law + units | Reported; dimensionally consistent (checked: $10^6 C_d A\sqrt{2\gamma\Delta P/\rho}$ in m²·m/s → mL/s ✓; the seed $C_d$ = 0.356 is reproduced from 0.664 L/min, $D$ = 3.665 mm, ΔP = 34.5 mmHg ✓) |
| Pressure clamp at −5 mmHg | Form **NOT REPORTED** |
| Nominal parameters | Values reported; **sources NOT REPORTED** |
| Scaling exponents | Reported with sources ✓ |
| Initial conditions | Reported (Table 3) ✓; $V_{0,SVEN}$ value for P1 **NOT REPORTED** |
| Blood volume | 70 mL/kg; source **NOT REPORTED** |
| Steady-state criterion | Reported ✓ |
| Output definitions (EDV, ESP at max P/V, SW) | Reported ✓ |
| Objective $J_{sys}$, $J_{shape}$, $J_r$, $J_b$ forms | **NOT REPORTED** |
| Prior confidences $c_j$ | **NOT REPORTED** |
| Bounds rule | **NOT REPORTED** |
| Stage A–F free-parameter sets | **NOT REPORTED** |
| FD gradient settings, exit criteria | **NOT REPORTED** |
| Sobol' start-point generation (scrambling type, skip) | Partly (seed 101); scramble method **NOT REPORTED** |
| GSA design (Halton; skip/leap?) | Partly |
| PCE truncation, LARS settings | Reported ✓ |

**Symbols:**
- $T_{vc}/T_{vr}$ (text) vs $T_{c,LV}/T_{r,LV}$ (App. B).
- $\varepsilon_v$ vs $\epsilon_{valve}$.
- $C_{PC}$ (Fig. 2) vs $C_{PCOX}+C_{PCNO}$.
- "$R_{\mathrm{sys}}$ scale" (Table 5) vs "systemic resistance scale" (text). Fine, but define the symbol once.
- "$\Delta P_{ref}$ 20 mmHg" (App. B) appears nowhere else except the footnote.
- **σ is overloaded:** "$\sigma$ declared uncertainty" vs "SD" in the noise tables. Acceptable.

---

## 9. Model assumptions

| Assumption | Justification given? | Probable effect | Breaks down when |
|---|---|---|---|
| Lumped 0D, single PA compartment | Partly | PA pulse pressure too wide (+5.1 / +6.9 mmHg) — **observed**, largest error | Narrow PA pulse with normal mean (this patient) |
| Time-varying elastance, linear ESPVR/EDPVR | Literature | Volumes depend on unscreened $V_0$; LV filling not well constrained | Diastolic dysfunction; large volume changes |
| Orifice law with constant $C_d$ through the cycle | Gorlin; clinical Bernoulli | Shunt timing correct in sign | Defect size changes during the cycle (perimembranous defects with aneurysmal tissue) |
| Chamber properties unchanged at closure | Stated (Limitations) | **No LV unloading predicted** | Any time after the acute state (S-1) |
| SVR unchanged at closure | Stated | Cannot reproduce the measured SVR rise | Sympathetic response, anaesthesia changes |
| Same anaesthesia/loading state pre and post | **Not stated** | Pressures pre and post were both under cath anaesthesia, volumes post were not (1 day later, likely awake) | Volume comparison |
| 5-min pressures and 1-day volumes = one state | Stated | Unknown | — |
| Blood volume 70 mL/kg, fixed | Literature | Stressed volume not screened → resistances absorb venous effects (stated) | Hypovolemia / fluids during cath |
| Atrial unstressed volumes and venous compliances not calibrated | **Code-based reason** | $C_{SVEN}$ influence carried by other params | — |
| Gaussian, independent measurement errors | Implicit | Fick-derived $Q_s$, $Q_p/Q_s$, shunt are correlated | Derived quantities |

No added complexity is demanded. The 0D choice is appropriate for the question. The two-compartment PA is the only structural change the results directly motivate.

---

## 10. Numerical implementation audit

| Numerical item | Reported? | Current information | Missing information | Required fix |
|---|---:|---|---|---|
| Spatial discretization | n/a | 0D | — | State "no spatial discretization; ODE system of 14 states" (one clause, pre-empts reviewer template questions) |
| Time integrator | ✓ | `ode15s` (variable-order NDF) | order limits, initial step | Add "MaxOrder default 5" or the value used |
| Tolerances | ✓ | rtol 1e-5, atol 1e-6 | atol per-state scaling? | State it |
| Max step | ✓ | $T_{HB}$/20 | — | — |
| Jacobian | ✗ | — | numerical or analytic, sparsity | Add |
| Steady state | ✓ | 5-cycle batches, 0.1 mL / 0.1 mmHg / 0.5% | Number of cycles actually needed for the 3 reported states | Add (e.g., "reached in N cycles") |
| Transient removal | ✓ | outputs from 2 further cycles | — | — |
| Event handling | ✗ | valves smooth; clamp at −5 mmHg | clamp form | Add |
| Software versions | ✓ | MATLAB R2025a, UQLab 2.2, R 4.5.2 | toolboxes (Optimization, Statistics), ggplot2 version | Add |
| Optimiser | partial | `fmincon` interior-point/SQP, 300 evals/40 it | FD type and step, TolFun/TolX/OptimalityTolerance, exit flags | Add (Supplement stage table) |
| Random seeds | ✓ | 101 (multi-start), 20260926 (noise) | Halton design seed/skip; Sobol' scramble type | Add |
| Hardware | ✗ | — | CPU, RAM, cores, parallel workers | Add |
| Runtime | ✗ | — (known: about 84 min per calibration) | single sim, GSA, calibration, noise batch | Add |
| Conservation | partial | post cycle closes < 0.03 mL | total blood-volume drift over the cycle for all 3 states | Add (cheap) |

---

## 11. Convergence (0D equivalents of mesh and time-step independence)

There is no mesh. The relevant tests are below, all **cheap** (minutes to a few hours). Present them as a 1-panel figure or a 1-table Supplement section.

1. **Solver-tolerance independence.** Re-simulate the calibrated pre-closure and the predicted post-closure states with rtol/atol = 1e-6/1e-7 and 1e-7/1e-8 and MaxStep $T_{HB}$/50 and /100. Report the max relative change over the 25 outputs. Target: < 0.1%, which is far below the measurement σ.
2. **Cycle independence / steady state.** Report the 25 outputs at cycles N, N+5, N+20 after the steady-state flag.
3. **Conservation.** Total volume $\sum V_i$ + $V_{PC}$-equivalent vs BV over the last cycle.
4. **GSA truncation check.** See §3 (20 design points to full steady state).
5. **PCE convergence.** LOO error and top-5 $S_T$ at N = 64, 96 and 128 (sub-designs of the existing Halton set; no new simulations).
6. **Objective smoothness.** Evaluate $J$ along one line through the optimum for two tolerance settings. This shows whether FD gradients are meaningful.

GCI or Richardson extrapolation **adds nothing here** (no grid). Do not include it.

---

## 12. Verification vs validation

- The manuscript **does not misuse "validation"** in the main text: a good point. It uses "comparison" and "prediction". The run manifests and old docs do use "validation" (`ValidationHoldoutMetrics`), so keep that word out of the paper.
- **Verification evidence present:** $Q_{VSD}$ = 0 exactly after closure (checked); post cycle closure < 0.03 mL.
- **Verification missing:** §11 items 1–3, 6.
- **Validation (prediction against withheld data):**

| Comparison | Reference data | n | Uncertainty | Metric | Threshold | Rationale |
|---|---|---|---|---|---|---|
| Pre fitted targets | cath, Fick | 9 | σ declared | ε, z, RMSE, $\chi^2$ | none (by design) | stated |
| Pre not-fitted | cath | 7 | no σ | ε | none | — |
| Post prediction | cath (5 min), echo (1 day) | 16 | σ declared (post $Q_s$ σ questionable) | ε, z, $\chi^2/N$, pressure RMSE | none | stated |
| Literature agreement | Tables 10–11 | qualitative | — | consistent/partly/not | — | Correctly labelled as consistency, not validation |

- **Strongest claim:** post-closure systemic and PA diastolic and mean pressures are predicted to within 1–14% without refitting.
- **Weakest claim:** "predicts the direction of every measured catheter pressure change" (false for 2/7).
- **Additional validation that would matter most:**
  - (a) a **naïve-baseline comparison**: prediction error vs "post = pre measured". If the model does not beat "nothing changes" for the pressures, the paper must say so. For a shunt with Qp/Qs 1.19 this is a real risk: the measured pressures barely change.
  - (b) the predictive interval from the noise MC, covering or not covering each post measurement.
  - (c) if any further patient has pre and post data, a second case.

---

## 13. Parameter provenance

| Group | Source given | Measured / fitted / assumed / literature | Issue |
|---|---|---|---|
| 7 calibrated params | Table 5 | fitted (seeded from clinical where available) | Bounds rule missing |
| $C_{SAR}$, $C_{PAR}$, $C_d$ seeds | Table 3 with formulas | clinically derived | ✓ (the "½ peak gradient" choice is declared as this work's) |
| Allometric exponents | Table 2 | literature (Zhang, Lundquist) | ✓ |
| Nominal adult values (App. B, about 70 numbers) | "adapted from Ferrero … locally harmonized … in the released code" | literature + assumed | **Not traceable.** Needs a source column per row (Ferrero, Colebank, Shi 2011, or "this work") |
| Timing fractions | App. B | literature? | Source missing (Colebank? Ferrero?) |
| BV 70/82 mL/kg | text | literature | Source missing |
| LA/LV filling defaults 6/6/8 mmHg | Table 3 | assumed | Say "assumed" + cite |
| Valve R_open/R_closed, ε | Eq. 5 | modelling choice | OK; say "chosen for numerical smoothness" |
| $\varepsilon_s$ = 0.1 mmHg | Eq. 6 | modelling choice | OK |
| σ values | §2.9 | assumed | Source per class missing |

Mixing populations: adult reference values (Ferrero: adult shunt patients) are scaled to a 13.4-kg child by exponents derived on a different reference (Zhang). State this once as a limitation of the scaling step.

---

## 14. Data provenance and preprocessing

| Dataset | Source | Acquisition | Preprocessing | Missing |
|---|---|---|---|---|
| Demographics | medical record | — | BSA formula **NOT REPORTED** (Mosteller? Haycock?) | formula |
| Cath pressures pre/post | cath protocol form | fluid-filled? catheter type, zero level, averaging over beats | MAP by protocol formula ✓ | **catheter/transducer, zeroing, number of beats averaged NOT REPORTED** |
| Fick flows | cath | — | — | **VO₂ (measured/assumed), saturation sites, Hb NOT REPORTED** |
| Post volumes | echo report | LV biplane, RV 3D | — | machine, operator blinded?, which software; **two conflicting LV records** (declared ✓) |
| Post $Q_s$ | echo | — | — | **method NOT REPORTED** |
| Defect diameter | echo | — | — | view, colour vs 2D, systolic/diastolic frame |

Ethics ✓ (IRB number, consent). Availability: on request. See §26 for the repository conflict.

---

## 15. Sensitivity, uncertainty and robustness

**What uncertainty could reverse or materially weaken the conclusions?**
1. **Systemic-flow measurement (Fick, assumed VO₂).** The single replicate shows −2.7σ on $Q_s$ moves $R_{sys}$ +92%, $E_{B,RV}$ +141%, and moves predicted RVEDV from 43.3 to 31.9 mL (meas. 30.5). *The headline volume error could be mostly measurement-driven.* The noise MC answers this; it is the most important pending computation.
2. **Truncated GSA.** Could change the retained set (§3).
3. **Lower-bound-pinned parameters.** If relaxing the bounds lowers $R_{pul}$ further, the PA pulse-pressure error may shrink. That would change the "structural remedy" argument.
4. **σ(post $Q_s$).** Changes the post $\chi^2$ and $z$ but not ε.

**Not needed:** more local sensitivity analysis, more seeds, a second GSA method.

**Cheap and useful:** Fisher-information approximate SDs from the existing scaled sensitivity matrix ($\mathrm{Cov} \approx (S^\top S)^{-1}$). This gives local parameter CIs to compare with the MC spread.

---

## 16. Results audit

| Subsection | Question | Evidence | Answers it? | Problems |
|---|---|---|---|---|
| 3.1 GSA | Which parameters does the pre-closure record inform? | Fig. 3, numbers | Yes | Last paragraph contains Methods-type statements (LOO, truncation) → move LOO table to Supplement; keep truncation as a stated caveat |
| 3.2 Pre calibration | How well does the model fit? | Table 4, Fig. 4, Table 5 | Yes | "above the model" direction; 0.495/0.496 |
| 3.3 Pre hemodynamics | What does the model add? | 1 paragraph | Partly | Too thin for the paper's strongest mechanistic output; merge the shunt-phase decomposition (now in 3.2 and 4.3) here |
| 3.4 Post prediction | How good is the forward prediction? | Table 6, Fig. 5, Fig. 6 | Yes | Overlong number list; no baseline predictor; "every pressure" error elsewhere |
| 3.5 PV loops | How do chamber loads change? | Fig. 7, Table 7 | Yes | — |
| 3.6 Noise | How uncertain is the prediction? | n = 1 | **No (pending)** | — |
| *(new)* 3.7 Post-closure calibration | Which parameters must change at closure? | S-2 | — | To add |

Chain check (evidence → number → interpretation): broken at "every measured pressure change" (§4.1, §4.7) and "dominant parameter … assigned by physiology" (§4.1).

---

## 17. Statistics

- **n = 1. No inferential statistics are used.** That is correct; keep it ("no between-patient statistics were computed" ✓).
- **$R^2$ on parity plots is pseudo-statistical here.** $R^2$ = 0.96–0.99 over pressures spanning 5–100 mmHg is driven by the range, not by agreement. With n = 7 points it carries little information. **Recommendation:** drop $R^2$ from Fig. 6 and the text, and keep RMSE in physical units plus signed ε per quantity. (Lechuga and Owashi report $R^2$ over waveform samples, which is a different, meaningful use.)
- Noise MC: mean, SD, CV and 2.5–97.5 percentiles are descriptive, which is appropriate. Rename "95% interval" to "2.5–97.5th percentile range" in the table header (it currently reads like a confidence interval).
- **Reporting precision is inconsistent:** pressures with 2 decimals in Table 4 (100.00) vs 1 decimal in Table 7. Cath pressures are integers, so use 1 decimal for model values and integers for measured ones.

---

## 18. Figure audit

| Figure | Decision | Main problem | Exact change |
|---|---|---|---|
| **Fig. 1** workflow | Modify | Good, but the dashed arrows to "Supplement" and "Noise" overlap and are ambiguous; inputs and outputs per box not shown | Add the data used at each step at the left margin (pre cath → fit; post record → compare only). Make the post-record box visibly "withheld" (lock icon or hatched). Add the new post-closure-calibration branch (S-2) |
| **Fig. 2** schematic | Modify | $C_{PC}$ vs two compliances; mixed italic/roman labels; legend does not explain scale groups; ground symbols inconsistent | Draw $C_{PCOX}$/$C_{PCNO}$; uniform label font; legend: "blue = calibrated directly or through a resistance scale"; add flow-direction arrows on both loops |
| **Fig. 3** GSA | Keep, small fixes | Panel A labels "0.08 0.07" collide with the dashed threshold; 45° labels are hard to read; Panel B colour scale is unlabeled in words | Horizontal bar chart (labels readable); move numbers outside the bars; B: column labels horizontal, colourbar title "$S_T$", outline cells ≥ 0.10; optionally a side panel with the calibrated-point index (now Fig. S2) for the 9 targets |
| **Fig. 4** pre waveforms | Modify | 3-row legends waste space; Panel E lacks the baseline; no phase markers; measured LV/RV lines are both dashed (colour only) | One shared legend; shade isovolumic contraction / ejection / relaxation / filling in all panels (this makes the 6/37/40/17% split visible); add the baseline to E |
| **Fig. 5** post waveforms | **Merge with Fig. 4** | Shown separately, the reader cannot see the change at closure, which is the paper's question | Pre (blue) and post (orange) overlaid in the same 4 panels, measured pre (black solid) and post (black dashed) lines; time axis normalised to cycle fraction (HR differs) |
| **Fig. 6** parity | Replace | Panels with n = 1–2 and RMSE values on 1 point; $R^2$ inflated by range; mixes $Q_p/Q_s$ and EF in one panel | **Signed relative-error dot plot**: one row per quantity (same order as Tables 4/6), x = ε (%), pre-fitted (filled), pre not-fitted (open), post-predicted (triangle), horizontal error bars ±σ/measured; optional shaded box for the noise-MC 2.5–97.5% range once available. This replaces both parity and part of Fig. 8B |
| **Fig. 7** PV loops | Keep | Legend for dashed/dotted lines sits on a separate row; RV y-axis cut at 27 mmHg | Single legend row; label EDV/ESV markers directly on the plot |
| **Fig. 8** noise | Replace (after batch) | n = 1, single points | Panel A: box plots of the 7 parameters (ratio to reference, log x-axis); Panel B: **paired pre vs post box plots** of outputs with measured markers (S-3) |
| Fig. S1/S2 tornado | Keep | — | Add bootstrap CIs |
| Fig. S3 before/after | Keep or delete | Duplicates Table 4 | Delete if space is needed |

**Missing figures that add information:**
1. **Numerical verification** (§11): tolerance and cycle convergence, 1 small panel or table (Supplement).
2. **Post-closure calibration** (S-2): Δparameter bar chart (post-fit / pre-fit) plus what-if sweeps.
3. **Pre/post box plots** (S-3), as part of new Fig. 8.

All figures: use one font family and size (~8 pt at print size); colour-blind-safe pre/post pair (blue/orange is fine). Check greyscale: blue vs grey baseline is too close in greyscale, so use a dashed baseline.

---

## 19. Caption audit

- **Fig. 4:** good. Add the phase shading definition and "HR 119 bpm".
- **Fig. 5:** add "HR 108 bpm; parameters unchanged from Fig. 4".
- **Fig. 6:** after replacement, define ε and σ in the caption.
- **Fig. 7:** add the HR of each state and "measured volumes from echo one day after closure; pressures not measured".
- **Fig. 8:** currently self-contained but reports n = 1. Regenerate it.
- **Fig. 3:** say "scaled baseline, before calibration, N = 128 PCE training simulations".

---

## 20. Table audit

| Table | Issue | Fix |
|---|---|---|
| 1 Clinical | "reported" use code; BSA formula | Use "C"; footnote the BSA formula |
| 2 Exponents | Footnote letters out of order (b before a) | Reorder |
| 3 Seeding | Good | Add $V_{0,SVEN}$ result |
| 4 Pre | 2-decimal measured pressures | Integers for cath values |
| 5 Calibrated params | No bounds source | Add a "bound rule" footnote |
| 6 Post | Good; identical rows to Table 4 ✓ | — |
| 7 PV | Good | — |
| 8/9 Noise | "95% interval" label | "2.5–97.5th percentile" |
| 10/11 Physiology | 11: missing the 2 arterial rows | Add systolic and diastolic RFA as "not consistent" |
| App. B | No source column | Add one |
| **New: solver/settings table** | Settings are spread over §2.5, §2.6, §2.7, App. B | One Supplementary table: integrator, tolerances, MaxStep, SS criteria, GSA design, PCE settings, optimiser per stage, seeds, hardware, runtime |
| **New: assumptions table** (optional) | Assumptions spread over text | §9 table above, in the Supplement |

---

## 21. Discussion audit

- **Result vs interpretation vs mechanism vs speculation** are well separated in §4.2–4.4 (strong).
- **Unsupported or over-stated statements:**
  - "exactly as clinical hemodynamics does" (L680): drop "exactly".
  - "guidelines … recommend exactly this restriction" (L684): drop "exactly".
  - "A proximal and distal pulmonary arterial compartment is the direct structural remedy": a hypothesis; say "is expected to reduce".
  - "the retained parameter set … should be re-screened for each patient": a recommendation from n = 1; say "we expect … to differ for nonrestrictive defects".
- **Repeats Results:** shunt phase split; per-beat venous return.
- **Missing numerical significance:**
  - what the conditioning (cond 61.9) and the noise result mean for predictive use;
  - whether the GSA truncation could change the parameter set.

---

## 22. Limitations audit

**Current list (10).** All are meaningful. **Missing:**

| Missing limitation | Cause | Affects | Direction / magnitude | Threatens main conclusion? | Fix |
|---|---|---|---|---|---|
| **Acute-only prediction (S-1)** | Fixed chamber properties; 5-min / 1-day data | All post volumes, LA pressure | Over-predicts LV size at later times (Yuan 2025) | No: it scopes it | Text + future work: follow-up echo, time-varying parameters |
| GSA on non-periodic training outputs | 10-cycle cap | Retained set | Unknown | Possibly | §3 check |
| FD gradients on a numerically noisy objective | ODE tolerances | Optimum location | Small if §11 passes | No if checked | §11 |
| Noise MC perturbs pre data only, at one operating point | Design | Prediction uncertainty | Under-states total uncertainty | No | State it |
| σ values are assumed, not measured | No device-specific data | z, $\chi^2$ | Directly scales $\chi^2$ | No | State it with sources |
| No naïve-baseline comparison | — | Value of the prediction | Unknown | **Yes, until done** | Compute |
| Anaesthesia/loading differs between cath and next-day echo | Clinical protocol | Volume comparison | Unknown | No | State it |
| Single-compartment PA (already implicit) | Structure | PA systolic | +5 to +7 mmHg | No | Already discussed; keep |

---

## 23. Conclusion audit

- **New information in the Conclusion:** none, except the generated noise sentence (n = 1). Remove it until the batch finishes.
- **Over-claiming:** none in the current Conclusion. It is constrained and quantitative (good).
- **Narrowest defensible conclusion now:** *In one child with a restrictive perimembranous VSD, a 0D model calibrated on 7 sensitivity-selected parameters predicted post-closure systemic and PA diastolic/mean pressures within 1–14% without refitting. It did not predict PA systolic pressure, chamber volumes or left-heart unloading, because no pre-closure measurement constrains the chamber parameters and chamber properties were held fixed at closure.*
- **Stronger conclusions possible after:**
  - (a) noise MC: "post-closure volume errors are within / outside the measurement-induced predictive range";
  - (b) post-closure calibration: "closure required a change of X% in LV elastance / unstressed volume";
  - (c) naïve baseline: "the model reduces pressure-prediction error by X% relative to assuming no change".

---

## 24. Ambiguity register

The manuscript is already unusually quantitative. Remaining items:

| Location | Statement | Why ambiguous | Required |
|---|---|---|---|
| L118 | "approximately five minutes after closure"; "approximately one day later" | Timing is a declared limitation | Exact times if in the record (e.g., "6 min", "22 h") |
| L110 | "measurements sufficient to initialize and evaluate" | Inclusion criterion not checkable | List the minimum required fields |
| L245 | "the continuous-switch approach used in closed-loop lumped-parameter models" | Which models | Cite or delete |
| L301 | "$\nu$ = 70 mL/kg above one year" | Source | Cite |
| L306 | "the code redistributes a measured SVR or PVR" | Procedure undefined; did not run | Delete (not used for P1) or define |
| L355 | "not in the calibration parameter set of the code, which was fixed before the runs" | Reason not given | Give a scientific reason or state "no reason; limitation" |
| L373 | "penalizes outputs outside pediatric physiological ranges" | Which ranges, which outputs | Table in the Supplement |
| L373 | "coupled term $J_{sys}$" | Formula missing | Equation |
| L376 | "a systemic-output refinement" | Targets, parameters | Stage table |
| L376 | "stronger prior and bound penalties" | Given numerically ✓ | — |
| L387 | "declared uncertainty" | Declared by whom, based on what | Source per class |
| L403 | "Every replicate is reported" | ✓ | — |
| L484 | "The objective therefore preferred values below the lower bound" | Not shown | "The optimum lies on the bound; the unconstrained optimum was not computed" |
| L508 | "matches the diastolic and mean readings" | How close | Give ε (−7.5%, −2.5%) or point to the table |
| L680 | "narrow sampling range of $C_{SAR}$ … explains part of this difference" | "Part": how much | Delete "part" or test (compare $S_T$ with a widened range) |
| L721 | "small hemodynamic changes", "slightly" | Relative to what | Give numbers (already present nearby; delete the adjectives) |
| L723 | "changes in loading conditions during and after the procedure" | Which | Name them (anaesthesia, fluids, HR) |
| L725 | "poor concordance" (Schramm) | Quantity | Give the reported metric |
| L742 | "should not be used for those purposes until …" | OK | — |
| Supp. S1 | "exactly reproducible under the seeds tested" | True but misleading | Add "because the selected start does not depend on the seed" (already explained in the next sentence; merge) |

---

## 25. Reproducibility challenge (paper and supplement only)

| Step | Label | Missing |
|---|---|---|
| 1. Acquire the data | Partially | Table 1 complete; acquisition details missing (catheter, Fick VO₂, echo Q_s method) |
| 2. Preprocess | Reproducible | MAP, Q_s, shunt and SVR formulas (SVR implicit) |
| 3. Build the model | Partially | Nominal values given; clamp form; $R_{open}$ scaled value |
| 4. Scale | Reproducible | — |
| 5. Initialise | Partially | Σ in Eq. 8 ambiguous; LA defaults' source |
| 6. Solve | Reproducible | Jacobian option; otherwise complete |
| 7. Convergence test | **Not reproducible** (not done) | §11 |
| 8. GSA | Partially | Halton skip/leap; group-index rule; per-output LOO |
| 9. Calibrate | **Not reproducible** | $J_{sys}$, $J_{shape}$, $J_r$, $J_b$, $c_j$, bounds rule, stage composition, FD settings, Sobol' scramble |
| 10. Select start | Reproducible | — |
| 11. Predict post | Reproducible | — |
| 12. Metrics | Reproducible | — |
| 13. Noise MC | Reproducible | — (draw rules given) |
| 14. Figures | Partially | Scripts named, not archived |
| 15. Identifiability (Supp.) | Reproducible | FD step for S |

---

## 26. Code and data reproducibility

### ⚠ Privacy conflict — resolve before submission (P0)
`github.com/HafizAqilla/unified-vsd` is **PUBLIC**. It tracks `config/patient_reyna.m`, `config/patient_profile_Razka.m`, `config/patient_cohort_cases.m`, run protocols and docs whose file names contain the patient's **first name**, together with her age, sex, weight, height and full catheterization values. The local (untracked) config also contains other children's first names.

This conflicts with:
- §2.1, "Data were de-identified and the patient is referred to as P1";
- the Data-availability statement ("on request, subject to IRB approval").

A first name plus age, sex, hospital and procedure window can identify a child.

**Recommended actions (for the authors and supervisor to decide; I have not changed anything):**
1. Make the repository private now.
2. Ask the IRB / data protection contact whether a disclosure report is needed.
3. Before any public release, rename every patient file to a code (`patient_P1.m`).
4. Rewrite history or create a fresh repository without the named files.
5. Remove free-text clinical notes.

### Ideal archived package (Zenodo DOI at submission)
```
/code         model, scaling, calibration, GSA, noise (MATLAB), figure scripts (R)
/config       patient_P1.m (de-identified, coded), run protocol, seeds
/data         P1 clinical record as CSV with a data dictionary (only if IRB permits; else a synthetic example with the same fields)
/results      exported tables (CSV) and figure_data used by the paper
/figures      the generated PDFs
/environment  MATLAB release + toolboxes, UQLab version, R sessionInfo(), renv.lock
README.md     one command per figure/table; expected runtimes; expected outputs with checksums
LICENSE       (the repo already has one)
CITATION.cff  (exists)
```
**Essential:** code, config, exported tables, figure scripts, environment, README.
**Optional:** raw MAT files, intermediate GSA designs.

**Figure checklist:** for each of Figs 1–8 and S1–S3, list the script, input file, command and runtime in the README.

---

## 27. Terminology and notation (canonical choices)

| Concept | Variants found | Canonical |
|---|---|---|
| systemic flow | "systemic flow", "cardiac output", "$Q_s$" | "systemic flow $Q_s$" (define "= cardiac output" once) |
| prediction state | "post-closure prediction", "predicted post-closure", "forward prediction" | "post-closure prediction" |
| baseline | "scaled baseline", "baseline" | "scaled baseline" |
| resistance groups | "systemic resistance scale", "$R_{\mathrm{sys}}$ scale", "R.SC" (abstract) | "$R_{\mathrm{sys}}$ (systemic resistance scale)" |
| timing | $T_{vc}$ / $T_{c,LV}$ | $T_{vc}$, $T_{vr}$ everywhere |
| smoothing width | $\varepsilon_v$ / $\epsilon_{valve}$ | $\varepsilon_v$, $\varepsilon_s$ |
| pulmonary capillary compliance | $C_{PC}$ / $C_{PCOX}+C_{PCNO}$ | the latter |
| "fitted targets" / "governed targets" (Supp.) | both | "fitted targets" |
| "the code" | L306, L355, L682, L754, Supp. | "the implementation", or state the scientific reason |
| patient | P1 ✓ | — |

---

## 28. Writing quality (targets, not a rewrite)

- **Split:** L527 (§3.4 paragraph 2: 9 claims, 20+ numbers) and L355 (the GSA paragraph covers design + grouping + exclusion + retention + truncation + re-screen + caveat). Split each into 2–3 paragraphs or move numbers to tables.
- **Shorten:**
  - L245 "the hyperbolic-tangent form below is the implementation adopted in this work" → "We use".
  - L253 "The square-root law is unmodified; the sign-preserving factor … is introduced here" → "We add a sign-preserving factor that makes the law differentiable at zero gradient".
- **Delete:** "exactly" (×2), "therefore" where no inference follows (L484), the last sentence of §4.6.
- **Agency:** "was checked numerically before every post-closure simulation" → "we checked …" (fine either way; be consistent).
- **Quantify:** see §24.

---

## 29. Logical continuity

Problem (post-closure state unknown) → gap (no tested prediction) → question (what is predictable) → method ✓ → **verification ✗** → validation (held-out) ✓ → results ✓ → interpretation ✓ → conclusion ✓.

Broken links:
1. Intro P2 problem (residual shunt/AR) ≠ question addressed.
2. Verification step absent.
3. Aims (ii)–(iii) promise reproducibility and identifiability; Results deliver UQ (main) and identifiability (Supp.).
4. §4.6 refers to a post-closure calibration that does not exist.

---

## 30. Internal consistency (global cross-check)

| # | Item | Location A | Location B | Status |
|---|---|---|---|---|
| 1 | Abstract numbers | Abstract | §3 | **Contradict** (stale) |
| 2 | Direction of pressure change | §4.1, §4.7 "every" | Table 6 vs Table 4 (RFA sys/dia) | **Contradict** |
| 3 | Shunt mean | 0.495 (§3.2, Fig. 4E) | 0.496 (Table 4) | Rounding conflict |
| 4 | Shunt fraction | 13.6% (§3.3, model) | 16% (§4.4, measured) | Unlabelled |
| 5 | Shunt diastolic minimum | ">0.10" (§3.2) | "0.11" (§4.3) | Minor |
| 6 | "above the model" | §3.2, §4.5 | ε definition | **Sign/direction wrong** |
| 7 | Physiology agreement | §4.1 "each target" | Table 10 (2 partly) | **Contradict** |
| 8 | Acceptance threshold | Abstract "acceptance" | §2.9 "none" | **Contradict** |
| 9 | Code names | Abstract (R.SC, C.SVEN, vsd.Cd) | body symbols; C_SVEN not retained | **Contradict** |
| 10 | Measured post SV | LV 22.3 mL, RV 18.5 mL | echo $Q_s$ 2.54 → 23.5 mL/beat | Clinical record internally inconsistent; report it |
| 11 | $R_{open}$ | Eq. 5 unscaled | Table 2 scaled | Clarify |
| 12 | $C_{PC}$ | Fig. 2 | App. A | Notation |
| 13 | Model LV SV pre | 30.0 mL (EDV−ESV) | Qp/HR = 30.6 mL | Consistent (diastolic shunt drains LV during filling); one sentence would pre-empt the question |
| 14 | Section reference | L376 "Supplementary Section S1" | supplement S1 ✓ | OK |
| 15 | Supplement refers to "main text Section 2.7", "Eq. 10" | Supp. L42, L101 | main §2.7, Eq. 10 ✓ | OK; re-check after edits |

Checked and consistent: $C_d$ seed (0.356), $C_{SAR}$ and $C_{PAR}$ seeds, $z$(PA sys) = 5.10, $\chi^2$ = 31.30, reduced $\chi^2$ = 15.65, SVR derivations, post LV/RV SV × HR vs model $Q_s$ (3.19/3.20 vs 3.17 L/min), noise replicate numbers (Table 8 vs Fig. 8A ratios).

---

## 31. Claim–evidence matrix

| Major claim | Evidence | Sufficient? | Missing evidence | Action |
|---|---|---:|---|---|
| 7 of 14 parameters govern the pre-closure targets | PCE-Sobol', N = 128 | Partly | Steady-state training outputs; bootstrap CI; per-output LOO | §3 checks |
| Calibration reduces RMSE 0.274 → 0.105 | Table 4 | Yes | — | — |
| Pipeline is reproducible across seeds | Supp. S1 | **No** (trivially identical) | — | Reframe as "deterministic selected start"; drop from novelty |
| Parameters are practically identifiable | cond(S) = 61.9 | **No** | Noise MC (1 replicate already shows large shifts) | Reframe as "locally well-conditioned, globally sensitive to measurement error" |
| The model predicts post-closure pressures | Table 6, RMSE 6.9 mmHg | Partly | Naïve baseline; predictive interval | Compute |
| The model predicts the direction of every measured pressure change | — | **False** | — | Correct |
| Shunt behaves as a restrictive defect | GSA hierarchy, gradient | Yes | — | — |
| The model does not predict LV unloading because HR fell | Per-beat return arithmetic | Yes (model-internal) | Sensitivity: simulate post with HR 119 | **Cheap check (1 sim)**: confirms the mechanism directly |
| The measured $Q_s$ fall is a method artefact | Ralston 2023, Schramm 2023 | Partly | VO₂ source; echo method | Report the methods |
| Post-closure volume errors are of the size measurement error produces | 1 replicate | **No** | 50 replicates | Finish the batch |
| PA systolic error is structural (single compartment) | Argument | Partly | Bounds relaxation test | Optional 1 calibration |

---

## 32. Main text vs Supplement vs repository

| Material | Location |
|---|---|
| Objective terms in full, stage table, bounds rule, prior confidences, physiological ranges | **Supplementary Methods** |
| Solver/settings table, hardware, runtime | Supplement (1 line of runtime in the main text) |
| Verification (§11) | Supplement, 1 sentence + numbers in the main text |
| Per-output PCE LOO, bootstrap CIs | Supplement |
| Seed repeats, identifiability | Supplement (already) |
| Post-closure calibration (S-2) | **Main text** (new short subsection + 1 figure) |
| Pre/post box plots (S-3) | **Main text** (new Fig. 8) |
| Multi-start 18-point spread (S-4) | Supplement |
| Run manifests, provenance CSVs | Repository |
| Fig. S3 (before/after bars) | Delete (duplicates Table 4) |

---

## 33. IJNMBE submission-readiness

Verify the items marked † against the current Wiley author guidelines; requirements change.

| Item | Status |
|---|---|
| Journal template (Wiley NJD LaTeX) | ✗ (`aaj2026.sty`) |
| Data availability statement | ✓ present; must be reconciled with the public repo |
| Code / supporting files | ✗ on request only → Zenodo archive |
| Numerical convergence/accuracy evidence | ✗ |
| Figures final | ✗ (Fig. 8 placeholder, Fig. 5/6 redesign) |
| Novelty/significance statement† | ✗ prepare (3–4 sentences) |
| Graphical abstract† | ✗ prepare (workflow + one pre/post panel) |
| Supporting information file | ✓ exists; extend |
| Ethics statement | ✓ IRB number present |
| CRediT, funding, COI, acknowledgments | ✗ red placeholders |
| Cover letter stating the n = 1 rationale and the methodological contribution | ✗ |
| Abstract matching results | ✗ |
| Preprint references (Colebank 2026 arXiv) | Check publication status |

---

## 34. What Reviewer 2 would attack (10)

| # | Probable comment | Valid? | Severity | Preventive revision |
|---|---|---|---|---|
| 1 | "A single patient. The workflow is assembled from established components. What is the methodological contribution?" | Yes | Critical | Reframe as a predictability/UQ study; add the naïve baseline; finish the noise MC; add any further patients with pre data |
| 2 | "No numerical verification: tolerance, cycle and conservation tests are absent." | Yes | Major | §11 (cheap) |
| 3 | "The GSA used non-periodic training simulations; the parameter selection is unreliable." | Yes | Major | 20-point steady-state check |
| 4 | "Two parameters sit on their lower bounds; the bounds are unexplained; the fit is constrained by arbitrary limits." | Yes | Major | State the bounds rule; relaxed-bound run |
| 5 | "The objective has undefined terms; it is not reproducible." | Yes | Major | Supplementary equations and table |
| 6 | "Parameters are not identifiable: one noise replicate moved them by up to 141%." | Yes | Major | Full MC; present the predictive distribution, not point estimates |
| 7 | "Post-closure data mix 5-min catheter and 1-day echo, different methods for $Q_s$, two conflicting LV records. The comparison is not clean." | Yes | Moderate | Already declared; add the σ justification and the measured-SV inconsistency |
| 8 | "Is the prediction better than assuming no change?" | Yes | Major | Naïve baseline |
| 9 | "Chamber properties fixed at closure guarantees no remodelling; the volume comparison is designed to fail." | Partly | Moderate | Post-closure calibration (S-2) quantifies the required change; acute-only limitation (S-1) |
| 10 | "$R^2$ of 0.96–0.99 on 7 points spanning 5–100 mmHg is not meaningful." | Yes | Minor | Drop $R^2$; use ε dot plot |

---

## 35. What should NOT be changed

- The **held-out post-closure design** and the rule that nothing post-closure enters calibration (keep it explicit in Table 1's "Use" column).
- **Tables 4 and 6 with identical rows**, ε and z: clear and auditable.
- **Table 3** (seeding with formulas): better than any benchmark.
- **Discussion Tables 10–11** and the consistent/partly/not-consistent wording.
- **Honest negative results:** PA systolic, RVEDV, no LV unloading, measurement disagreement.
- **The Limitations list** (add to it; do not shorten it).
- The **0D model structure**. Do not add complexity, except possibly a second PA compartment, and only if you decide to test the structural hypothesis.
- **No inferential statistics on n = 1.**
- The Appendix ODEs (complete and correct).

---

## 36. Prioritised revision plan

### P0 — must fix before submission
| # | Location | Problem | Action | New computation? | Benefit |
|---|---|---|---|---|---|
| P0-1 | Public GitHub repo | Patient name + clinical data public | Authors/supervisor: make private, consult IRB, rename to P1, fresh public archive later | No | Ethics; avoids rejection or retraction |
| P0-2 | Abstract | Stale, contradicts body | Authors rewrite (numbers in §5) | No | Desk survival |
| P0-3 | §4.1, §4.7, Table 11 | "Every measured pressure change" false | Correct; add RFA rows as "not consistent" | No | Integrity |
| P0-4 | §2.5 / Supp. | No verification | Tolerance, cycle and conservation tests (§11 items 1–3, 6) | **Yes, cheap (< 1 h)** | IJNMBE scope |
| P0-5 | §2.6 | GSA on truncated sims | 20-point steady-state check; per-output LOO; bootstrap CIs | **Yes, moderate (a few h)** | Defends the parameter selection |
| P0-6 | §2.7 / Supp. S3 | Objective and stages not reproducible | Write $J_{sys}$, $J_{shape}$, $J_r$, $J_b$, $c_j$, ranges, bounds rule, stage table | No (read code) | Reproducibility |
| P0-7 | §3.6 / Fig. 8 / Conclusion | UQ n = 1 | Finish the noise batch (about 22 h) and regenerate | **Yes (running)** | Core result |
| P0-8 | §3.4 | No reference for prediction quality | Naïve baseline (post = pre measured; post = pre model) error table | **Yes, trivial** | Answers Reviewer 2 #8 |
| P0-9 | Template, back matter | Wrong template; placeholders | Wiley NJD template; CRediT, funding, COI | No | Submission |
| P0-10 | Data/code statements | On request only | Zenodo archive with a coded patient file | No | IJNMBE expectation |

### P1 — strongly recommended
| # | Location | Problem | Action | New computation? | Benefit |
|---|---|---|---|---|---|
| P1-1 | New §2.11, §3.7 (S-2) | No post-closure calibration | Post-closure calibration with chamber params free, after freezing the prediction; Δparameter table; R/C what-if sweeps | **Yes, about 2–4 h** | Quantifies what the prediction misses; supervisor request |
| P1-2 | Fig. 8 (S-3, S-4) | No pre/post distribution view | Paired box plots from the noise replicates; do NOT use the 3 identical seed runs | After P0-7 | Supervisor request; visual UQ |
| P1-3 | Limitations + future work (S-1) | Acute-only not stated | Add the limitation and the time-course future work | No | Supervisor request; scopes the volume errors |
| P1-4 | Intro P2, P5; novelty | Motivation/claim mismatch | Authors revise (see §6) | No | Desk survival |
| P1-5 | §2.9 σ | Post $Q_s$ 5% unjustified | Justify or change; recompute $\chi^2$ | No (arithmetic) | Defensible $\chi^2$ |
| P1-6 | §2.2 | Fick VO₂, echo $Q_s$ method, catheter details missing | Add from the protocol form | No | Supports the measurement argument |
| P1-7 | Table 5 bounds | Two params on bounds | Bounds rule; optional relaxed-bound calibration | Optional (1.5 h) | Reviewer 2 #4 |
| P1-8 | §4.4 mechanism | HR-fall explanation untested | Simulate post with HR 119 | **Yes, 1 sim** | Converts the argument into evidence |
| P1-9 | App. B | Nominal values without source | Source column | No | Provenance |
| P1-10 | Fig. 5/6 | Pre/post not overlaid; parity panels weak | Merge Figs 4/5; replace parity with the ε dot plot | No (re-plot) | Clarity |
| P1-11 | §2.1 | Other enrolled patients unexplained | State the count and reasons | No | Reviewer 2 #1 |

### P2 — quality improvements
- Fix the direction wording ("above/below the model"), 0.495/0.496, 13.6% vs 16%, the "each target" claim, and the §4.6 last sentence.
- Remove the Results/Discussion duplicates (shunt phase split, per-beat return).
- Split the long paragraphs (L355, L527).
- Notation unification (§27); nomenclature table.
- Settings table, hardware, runtime.
- Drop $R^2$; rename "95% interval".
- Fig. 2 and Fig. 3 fixes; captions (§19).
- Title (§4).

### P3 — optional polish
- Table precision (integers for cath readings).
- Footnote order in Table 2.
- Delete Fig. S3.
- Grammar list for the Intro (on request).
- Greyscale check of all figures.

### Five highest-leverage changes
1. **Finish the noise MC and make it the centre of the Results** (with the supervisor's pre/post box plots). It turns "one fit" into "a prediction with quantified uncertainty".
2. **Numerical verification and the GSA steady-state check** (P0-4, P0-5). Cheap, and it removes the main IJNMBE-scope objection.
3. **Naïve-baseline comparison** (P0-8). It answers whether the prediction is informative at all, with trivial cost.
4. **Post-closure calibration as a separate, labelled analysis** (S-2). It quantifies what changes at closure, and it fixes the dangling §4.6 claim.
5. **Reproducible objective/stage specification plus a de-identified code archive** (P0-6, P0-10, P0-1).

---

## 37. Final verdict

**Current state:** *scientifically promising but major revision needed.*

**Scientific contribution:**
- The held-out pre→post prediction in a pediatric VSD is a real and uncommon contribution. Combined with the GSA-based explanation of what is and is not predictable, it is the paper's value.
- Novelty currently leans on reproducibility and identifiability claims the evidence does not carry.
- Reframed as a predictability and uncertainty study, the contribution is clear and defensible.

**Numerical rigor:**
- Solver settings and the steady-state criterion are reported better than in three of the four benchmarks.
- But there is no verification evidence (tolerance, cycle, conservation).
- The sensitivity screen rests on non-periodic training simulations.
- Optimiser convergence is undocumented.
- All of this is fixable with modest computation.

**Validation strength:**
- Strong in design: withheld post-closure data, no refitting, every measured quantity reported with ε and z.
- Weak in reference: no naïve baseline; post-closure σ partly unjustified; predictive uncertainty incomplete (n = 1 of 50).

**Reproducibility:**
- The forward model is reproducible from the paper.
- The calibration is not: four objective terms, prior confidences, bounds rule and stage composition are missing.
- Code is on request only.
- The public repository holds identifiable patient data, which must be fixed first.

**Writing/structure:**
- Direct, number-dense and honest.
- Discussion quality is high.
- Problems: a stale abstract, an Introduction whose clinical motivation does not match what the model can do, a few factual misstatements (direction of change, physiology agreement), and over-long number lists in §3.4.

**Figures:**
- Clean and consistent styling.
- Fig. 8 is a placeholder.
- Figs 4 and 5 should be merged to show the change at closure.
- The parity plot should be replaced by a signed-error plot.
- One verification figure or table and the supervisor's box plots are missing.

**IJNMBE scope fit:**
- Borderline today: an application of standard tools to n = 1.
- It fits if the paper foregrounds the computational questions (verification, sensitivity-informed predictability, UQ propagation from calibration to prediction) and documents them reproducibly.

**What would prevent me from recommending acceptance today:**
1. No numerical verification, and a GSA on non-converged simulations.
2. An incomplete uncertainty analysis (1/50) that already suggests the parameters, and hence the volume predictions, are not determined by the data.
3. A non-reproducible objective function and stage specification.
4. No comparison with a naïve predictor.
5. An abstract and two Discussion statements that contradict the results.
6. Identifiable patient data in a public repository, contrary to the de-identification statement.
