# Sense Foraging Study 2b: analysis pipeline (R)

The analysis the Study 2b preregistration specifies (`study2/preregistration_2b.md`, filed at https://osf.io/d4a5u from commit d3e60c9; OSF governs): the confirmatory factor structure of the 8 (H1), convergent and discriminant validity (H2–H3), criterion validity including the in-session pause (H4), incremental validity (H5), known groups (H6), and the exploratory analyses, H8 among them. The 10 is reported beside the 8 throughout. The registration says this code is the specification wherever its text is less precise; every place where the text allows more than one reading carries a comment starting `PREREG:`, listed at the end.

## Data (never in this repo)

- **Pull:** `extract.sql`, read-only, through the Supabase CLI from the radlab checkout:
  `npx supabase db query --linked -o json -f <path>/extract.sql > F:/radlab_data/sf_study2b/raw/extract_<date>.json`
- **What it pulls:** one row per fact, keyed by the random participant code. `source` says which: `enrollment` (consent, session status, step crashes), `questionnaire` (every submission with `step_index`, `received_at`, `resubmission_of`, `presentation_format` and the answers, which carry their own timing; the pause offer's answer is its whole record, `pause_1`), `vas` (the life-satisfaction face scale, an activity) and `step` (each step's entry and exit). Prolific IDs are never read. Staff test enrollments are pulled with `is_test` and excluded by `01`, so the flow counts them.
- **Where it lives:** `SF2B_DATA_DIR` (default `F:/radlab_data/sf_study2b`), local and unsynced. `raw/` holds extractions; the scripts write `derived/` (R objects, per mode) and `output/<mode>/`.
- **Definitions:** `study2b/materials/` holds the platform definitions as built: radlab `scripts/sense_foraging/study2b/` at c016b10d (the scale, the trial statements, the vignettes, the pause, Brief MAIA-2, BSSS-8, NCS-6, SPSQ-SF Sensory Discomfort, ASI-3 Physical Concerns, SAMD), and the measures shared with the old Study 2 (`study2/materials`, the same definitions as radlab `scripts/sense_foraging/study2/`). In this repository the wording of the third-party questionnaires is withheld (`study2b/materials/README.md`); their ids, keys, scoring and numeric response values are as built. `config.R` checks ids, keys and subscales against them.
- **Session steps** (`build_study2b.py`): 0 background, 1 vignettes, 2 the scale, 3 the trial statements, 4–12 the first half of the validity block (per-participant order), 13 the pause offer, 14–19 the second half, 20 debrief.

## Modes (`SF2B_MODE`)

| Mode | What the answers are | Use |
|---|---|---|
| `permuted` (default) | Real answers, each item shuffled across the participants who answered it: distributions kept, every correlation destroyed | Running the code on a real extraction without seeing structure |
| `simulated` | `raw/simulated_extract.json`, written by `00_simulate.R` from a known model in the real extraction's shape | Proving the code recovers structure; `99_check_recovery.R` reports it |
| `real` | The answers | **Only after the preregistration is filed** |

In `permuted` mode exclusions, attention checks, timing, background answers and the pause stay real; answer-derived values (longstring, scores) use the shuffled answers.

## The pilot run (before filing)

The first 60 places are run. Their extraction goes through the whole pipeline in `permuted` mode, and the preregistration is filed only after that run, at the commit it ends on. The remaining 640 places then open on the same posting.

**What the shuffle keeps, and what it destroys.** Each item's answers are shuffled across the participants who answered it, item by item and independently. What stays real:
- who enrolled, finished or stopped, and where;
- the attention checks and exclusions;
- all timing;
- the background answers, including practice;
- the pause: offer, choice, time and device;
- each item's own distribution, including its "Prefer not to answer" rate.

What goes is every link between two answers, within a person or across measures. Every correlation, loading and fit index in the run is therefore noise by construction. Nothing in it can show whether a hypothesis would hold, so nothing in it can tempt a change to one.

**What the run decides,** with the rule for each set before the pilot:

| Question | Read from | Rule |
|---|---|---|
| Does the code read the real export? | Pipeline counts against the platform's (admin completions) and Prolific's (approved submissions); every response key mapped to an item; "Prefer not to answer" read as missing; `01_pause_check.csv` | Any mismatch or unread key: the code is fixed and rerun. Fixes before filing are not deviations. |
| Does it run end to end on real records? | `run_times.csv`; the step log | A crash is fixed. A model that fails or is refused is expected: shuffled data at n = 60 have no structure to fit. Estimation was validated on simulated data. |
| Can the rest open unchanged? | Completion and breakoff by step; failure rate of each attention check; median time | One check failed by more than 10%, or more than 15% stopping at one screen: that screen is looked at before the rest open. If any material changes, the 60 are excluded and 60 places added (Section 3). |
| Pay | Median completion time without the pause | £12.71 an hour × that median (Section 3). |
| Power at the real rates | Pause take-up and engagement; practice-group sizes; projected complete-case n for H5 from each item's missing rate | Section 3 reports power at the observed rates, from `study2/power/`. No hypothesis, bound or measure changes. |

Before the 640 open, Norm sees only this run's report: counts, rates, times, run status and code changes. The platform's response views are not opened for the pilot.

**Filing, after the run.**
1. Fill the preregistration's remaining `[AT FILING]` entries from the run's report: the header (filing date), the pilot report (Section 3, Pilot) and power at the observed rates (Section 3, Sample size rationale). Commit.
2. Push that commit to `origin main`, the public repository (https://github.com/Normega/sense-foraging-scale), and check it is there: the registration cites the repository at that commit, and nothing is attached in its place.
3. `Rscript study2/filing/build_filing.R` from the repo root. It writes `F:/radlab_data/sf_study2b/filing/<commit>/`: the preregistration as PDF and HTML, `osf_fields.txt` (one block per OSF template field), a code archive kept as a private record, and a SHA-256 manifest. It fills the Code paragraph's commit itself. Any open marker or uncommitted change makes it a DRAFT build; file only a build without that mark.
4. On OSF: the Preregistration template, the fields pasted from `osf_fields.txt`, "Registration prior to analysis of the data". No archive is attached; the PDF may be.
5. Then open the 640 places.

## Steps (`run_all.R`, from the repo root)

0. **`00_simulate`** (simulated only). The scale from Study 2a's fitted one-factor model of the 13 statements with the agreement factor (`sim_population.json`, written once from 2a's real data by `make_sim_population.R`; aggregates only). Every validity score is a latent factor with a known correlation with the scale's factor and known within-battery correlations; the past-week items, the pause, practice and the vignettes are generated from latents with known values. One gender difference in a loading is planted (V10, women .20). Planted mess: staff tests, enrollments that never reach the scale, two and one failed checks, sessions stopping after the trial statements, before the pause and in the second half, speeders, longstrings, many "Prefer not to answer", flagged resubmissions and later re-entries.
1. **`01_prepare`:** one row per participant by the registered rules: first complete submission of each questionnaire (flagged resubmissions dropped first); exclusions (staff tests; no scale; two or more of the three checks failed, a check not reached is not a failure); flags (one check; < 2 s per position over the 14; longstring ≥ 10 of 14; "Prefer not to answer" on > 2 of 13); the 8 (7 of 8 answered) and the 10 (9 of 10); validity scores; the engaged pause; practice groups; separation; the trial's Access + Valuing score; positions within each half. The pause notes stay in `derived/`.
2. **`02_descriptives`:** item distributions, agreement rates, missing data per measure, alpha per score, the sample.
3. **`03_structure`:** H1 on the 8 and the 10, and the 8 without flagged participants; the models without the agreement factor; the 13 statements (the exploratory three's loadings); residual matrices; model-specific cutoffs.
4. **`04_validity`:** H2 (each of the four correlations reported on its own; H2 holds when at least three of the four lower bounds exceed .30, `h2_needed` in `config.R`), H3a–e; the latent matrix for the 8 and the 10; the agreement factor's correlation with every measure; the H3d bootstrap; the observed matrix; the confirmatory cells without flagged participants.
5. **`05_criterion`:** H4a and H4b (Spearman, Holm), the logistic regression with device class.
6. **`06_incremental`:** H5a (three hierarchical regressions, Holm) and H5b (proportional odds).
7. **`07_known_groups`:** H6 and its sensitivity analyses; the latent difference (MIMIC).
8. **`08_exploratory`:** H8 and the gender and country invariance; the trial's measure; the pause, further; order effects; the vignettes; BIDR dichotomous.
9. **`09_report`:** `09_decisions.csv`, one line per confirmatory test (H2: one per neighbour, then the at-least-three-of-four line with its count), and `09_report_study2b.json`, every number above (aggregates only).
99. **`99_check_recovery`** (simulated only): `99_recovery.csv`, pass or fail per check, and `99_validity_recovery.csv`.

Model settings in every fit (`sf2b_models.R`): ordinal, WLSMV, theta parameterization, `missing = "pairwise"`, `conditional.x = FALSE`; the agreement factor's loadings fixed at 1 on the unrecoded statements, variance free, uncorrelated with every other factor.

## Validated

On simulated data (seed 20261198, 9 November 2026, the registered 2,000 H3d resamples; 749 enrollments, 715 analysable): all 19 checks in `99_recovery.csv` pass.
- Every planted exclusion, flag and repeat submission was handled as registered.
- The 8's content loadings were within RMSE .015 of the population (over the seven statements without the planted difference), the agreement loadings within .020, and the H1 rule passed.
- 25 of 27 latent validity CIs covered the population value, and 27 of 27 covered the sample's own latent correlation; all 4 contrasts with an exact population value were covered.
- The planted gender loading difference (V10) was the one statement released at the loadings step.

The first run (seed 20261199, 200 resamples) also passed all 19. Rerun on 13 November with H2's at-least-three-of-four rule (seed 20261198, 2,000 resamples, 8 workers, 15 minutes): all 19 pass, every number as before; H2 holds 4 of 4.

H3d did not hold in that simulated sample (R .760, 95% CI .717–.812), and that is the simulated truth, not the code. In it, the four H2 neighbours correlate with each other only through the scale's factor, at .45–.50 each, which makes the population R .735. H3d's power has since been simulated (`study2/power/power_2b_h3d.R`). It is 42% at that truth, 83% at R .73 and 92–100% at R .71 or below.

One defect was found and fixed by the dry run. A validity factor with free loadings takes its sign from its first item, which for BFI-2 Aesthetic Sensitivity, SBI and several FFMQ and Brief MAIA-2 subscales is reverse-keyed, so their correlations came out negative. Every validity loading now starts at ±.6 by its keying, and `sf2b_orientation()` refuses any model whose factor is not oriented as keyed.

## Runtime

16 minutes end to end on 20 workers with the registered 2,000 H3d resamples (`run_times.csv`): 12 of them are `04_validity` (98 validity models, then the bootstrap, about 6 core-seconds per resample). With `SF2B_BOOT=200`, about 4 minutes.

## Readings of the preregistration (search the scripts for `PREREG:`)

Confirmed by Norm on 13 November (decision 41) and stated in the preregistration, Section 5, "Specifications fixed in the code". A change to any of them after filing is a deviation.

- **Validity scores (01).** "Scored per their manuals": mean or sum per manual, with Study 2's 80% rule for missing items. For a 3-item subscale that means all three answered. In the dry run, with 1% of items unanswered at random, the H5a complete-case sample fell from about 650 to 410–470, because H5 needs all 13 three-item FFMQ and Brief MAIA-2 subscales.
- **Practice groups (01).** Regular practice "yes" with frequency "none" or "less than weekly" is in neither H6 group. A missing regular-practice answer leaves every group missing.
- **Gender (01, 06).** H5: man, woman, other (non-binary or self-described), a factor; "prefer not to say" leaves the complete-case sample. Invariance: men and women.
- **Engaged pause (01).** A "pause" record without timing counts as not engaged; none should occur, and `01_pause_check.csv` counts them.
- **H1 (03).** "The largest" residual is the largest in absolute value; freeing it adds that residual covariance; "every loading" is compared item by item. Content ω is the power simulation's formula: agreement loadings and the residual (co)variances count as error.
- **Model-specific cutoffs (03, reported only).** The 5th percentile of CFI and the 95th of RMSEA and SRMR in 200 data sets generated from the fitted model, and from it with two residual correlations of .30.
- **BFI-2 Open-Mindedness total (04).** One factor of its 12 items, in its own model and in contrast 3's.
- **Single-item measures (04).** Life satisfaction, the past-week items and the engaged pause enter the latent matrix as single-indicator factors with no residual variance (the correlation with the item's latent response). H4's decisions are on the observed 8, as registered.
- **H3b (04).** Contrasts 1, 2 and 4 come from the measure's own model (BFI-2 facets, FFMQ facets, Brief MAIA-2 subscales); contrasts 3 and 5 from a model of just their two measures.
- **H3d (04).** The bootstrap CI is the percentile interval; a resample whose fit fails is left out and counted.
- **H4 (05).** The logistic regression's covariate is device class (desktop, phone, tablet) as recorded by the pause.
- **H5 (06).** Both steps use the cases complete on every variable of step 2; H5b's coefficient is per SD of the 8.
- **H6 (07).** Student's t; d's CI is the coefficient's CI over the pooled SD. The latent difference is a MIMIC coefficient with S's residual variance 1.
- **H8 (08).** Section 5 lists "configural, then loadings, then thresholds"; for ordinal items Wu & Estabrook (2016) require thresholds before loadings, which is the order used. The text should say so. Partial invariance releases, one at a time, the statement whose release gives the best CFI, at most 4 per step. Sparse categories (< 5 in a group) are merged toward the centre.
