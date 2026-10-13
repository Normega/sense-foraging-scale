# Sense Foraging Study 2: analysis pipeline (R)

The analysis the Study 2 preregistration specifies (`study2/preregistration_draft.md`, as of commit 51c41d5): confirmatory
structure of the 20-item short form, validity, incremental validity (H5a wellbeing, H5b practice), invariance, known
groups, two-week stability, and the candidate items (steps 7-10). The preregistration says this code is the specification wherever its text is less precise. Every
place where the text allowed more than one reading carries a comment beginning `PREREG:` in the script that implements
it; they are listed at the end of this file.

## Data (never in this repo)

- **Pull:** `extract.sql`, read-only, through the Supabase CLI from the main checkout:
  `npx supabase db query --linked -o json -f <path>/extract.sql > F:/radlab_data/sf_study2/raw/extract_<date>.json`
- **Run it after** `deidentify_external_enrollments` has run on Study 2 or its retest with `p_confirm = true`: that is what
  writes `external_identity_links`, the only link between a person's Study 2 code and their follow-up code. Before
  then every follow-up row is "unlinked" and H7 has no pairs.
- **What it pulls:** one row per fact, for both studies (Study 2 and its retest), keyed by random participant codes.
  `source` says which: `enrollment` (posting, consent, session status), `questionnaire` (every submission, with
  `step_index`, `received_at`, `resubmission_of`, `presentation_format` and the answers, which carry their own timing
  and, for the candidates, the presented order), `vas` (the life-satisfaction face scale, which is an activity and
  lands in `vas_responses`), and `step` (each step's entry and exit times). Prolific IDs are never read. Staff test
  enrollments are pulled with `is_test` and excluded by `01`, so the flow counts them.
- **Where it lives:** `SF2_DATA_DIR` (default `F:/radlab_data/sf_study2`), local and unsynced. `raw/` holds extractions
  and the follow-up invitation list; scripts write `derived/` (R objects, per mode) and `output/<mode>/`.

## Modes (`SF2_MODE`)

| Mode | What the answers are | Use |
|---|---|---|
| `permuted` (default) | Real answers, each item shuffled across the participants who answered it: distributions kept, every correlation destroyed | Running the code on the real extraction without seeing structure |
| `simulated` | `raw/simulated_extract.json`, written by `00_simulate.R` from a known model in the real extraction's shape | Proving the code recovers structure; `99_check_recovery.R` reports it |
| `real` | The answers | **Only after the preregistration is filed** |

In `permuted` mode, exclusions, attention checks, timing, posting, background answers and the "Prefer not to answer"
count stay real, as in Study 1; item-derived values (longstring, scores, agreement index) use the shuffled answers.

## Steps (`run_all.R`, from the repo root)

0. **`00_simulate`** (simulated mode only): a raw extraction from a known model, with planted mess (staff tests, partial
   sessions, attention-check failures including people who never reached a check, speeders, longstrings, many "Prefer
   not to answer", a flagged resubmission, a later re-entry, retests outside the window, an unlinked retest), and a
   follow-up of about 300 of 400 invited that repeats the short form and the past-week items. Writes the population
   truth to `derived/simulation_truth_simulated.rds`.
1. **`01_prepare`:** one row per participant by the registered rules: first complete submission of each questionnaire
   (flagged resubmissions dropped first); exclusions (staff tests; no short form; both checks answered and both failed);
   flags (one check failed or the only check reached; < 2 s per item; longstring >= 10; > 4 "Prefer not to answer");
   samples (structure = everyone analysable; validity = reached the end of the validity block; H6 = background and short
   form); practice groups; the country-group rule; observed scores (80% rule); the agreement index; validity-block
   positions; follow-up pairs (short form and past-week items) and the 14-21-day flag.
2. **`02_descriptives`:** item distributions, missing data per measure, internal consistency, the sample, seconds per
   questionnaire.
3. **`03_structure`:** H1. The primary ESEM with the random intercept and its decision; the per-factor rule and
   omega; Model A with its modification indices; the ESEM without the random intercept; a one-factor model; the
   negative-wording test for Letting be; the primary model without flagged participants.
4. **`04_validity`:** H2-H4. Latent-observed correlations, each validity score added in turn; the bootstrap tests
   H3b and H3d; the observed validity matrix; the random intercept's correlations with every measure; the confirmatory
   cells again without flagged participants.
5. **`05_incremental`:** H5a: hierarchical regressions (delta-R^2 >= .02 and Holm across the three outcomes), the
   version with spirituality and religiosity in Step 1, and the latent version. H5b: the proportional-odds regression of
   practice frequency (MASS::polr; 1-df likelihood-ratio test, odds ratio per SD, change in McFadden's R^2), Access and
   Valuing separately, and dedicated vs never by Firth logistic regression (logistf).
6. **`06_invariance`:** posting, gender and country group, for the four-factor and the Access + Valuing CFAs (ESEM-
   within-CFA when a CFA is not configurally acceptable); partial invariance; latent means where the intercepts step
   holds; longitudinal invariance of Access + Valuing.
7. **`07_known_groups`:** H6 and its three sensitivity analyses (and without flagged participants); the frequency
   trend; latent means across practice groups where scalar invariance holds.
8. **`08_stability`:** H7. ICC(2,1) with CIs, SEM, MDC95, with and without out-of-window pairs and flagged
   participants; responders vs non-responders; the agreement index's stability.
9. **`09_candidates`:** step 7. The split, the 39-item ESEM on the selection half, the four rules, at most two
   replacements per factor, bootstrap selection rates, and the wording test with the candidates.
10. **`10_revised_form`:** step 8. The revised form and the fixed form on the holdout half.
11. **`11_forage`:** step 9.
12. **`12_av14`:** step 10.
13. **`13_exploratory`:** order effects; the bifactor ESEM and its omega_h (from the model's loadings, not
    `psych::omega`); further outcomes; the following week (each follow-up past-week item on baseline Access + Valuing,
    controlling for the same item at baseline).
14. **`14_report_data`:** every number above in `14_report_study2.json` (aggregates only), and `14_decisions.csv`,
    one line per confirmatory decision.
99. **`99_check_recovery`** (simulated mode only): `99_recovery.csv`, pass or fail per check.

`98_seed_replication.R` (simulated mode, optional) reruns the simulation with other seeds and reports interval coverage
across them (`98_seed_replication*.csv`).

`draw_retest_sample.R` is not part of the run: it draws the follow-up invitations (a simple random sample of eligible
completers, recorded seed) from an extraction and writes `raw/retest_invited_<date>.csv` (Study 2 codes). Those codes
must then be turned into the Prolific allowlist by someone with access to the identity table.

## Output files (`output/<mode>/`)

| File | What it holds |
|---|---|
| `01_flow.csv`, `01_retest_flow.csv` | Participant flow, main session and follow-up |
| `01_flags.csv`, `01_submissions.csv` | Flag counts; submissions, flagged resubmissions and repeats per questionnaire |
| `02_item_descriptives.csv` | Short-form and candidate items: % per category, mean, SD, % missing |
| `02_missing_by_measure.csv`, `02_scores_available.csv` | Missing data per measure; n per score |
| `02_internal_consistency.csv`, `02_sample.json` | Alpha per score; the sample, flags, minutes per measure |
| `03_fit.csv` | Fit of the primary ESEM, the ESEM without RI, Model A, the one-factor model, the wording model |
| `03_h1_loadings.csv`, `03_h1_factor_correlations.csv`, `03_h1_factor_rule.csv` | Keyed, oriented loadings; factor correlations; the per-factor rule and omega |
| `03_h1_decision.json` | H1, the factor decisions, the wording test |
| `03_modelA_modification_indices.csv`, `03_wording_loadings.csv` | Model A's 40 largest MIs; the wording model's loadings |
| `04_latent_validity.csv` | Every factor x validity score: r, 95% CI, one-tailed p, the predicted cell and whether r falls in it |
| `04_observed_validity.csv`, `04_random_intercept_validity.csv` | The observed matrix; the RI's correlations |
| `04_validity_decisions.json`, `04_bootstrap_draws.csv` | H2, H3a-d, H4 (and without flagged participants); the bootstrap draws for H3b and H3d |
| `04_latent_validity_unflagged.csv` | The confirmatory cells without flagged participants |
| `05_incremental_observed.csv`, `05_incremental_latent.csv`, `05_h5_decision.json` | H5a (with its sensitivity version and the latent version) and H5b (with both sensitivity analyses) |
| `06_invariance.csv`, `06_invariance_summary.json` | Every invariance model, delta CFI/RMSEA, items released |
| `06_category_merges.csv`, `06_latent_means.csv` | Sparse-category merges; latent means where permitted |
| `07_known_groups.csv`, `07_h6_decision.json` | H6 and sensitivity (a)-(c) |
| `07_frequency_trend.csv`, `07_practice_invariance.csv`, `07_practice_latent_means.csv` | H6 secondary analyses |
| `08_stability.csv`, `08_responders.csv`, `08_h7_decision.json` | H7 |
| `09_selection_loadings.csv`, `09_candidates.csv`, `09_replacements.csv` | Step 7: the 39-item model, every candidate against the four rules, the replacements |
| `09_bootstrap_rates.csv`, `09_step7.json` | Selection rates; the revised form; the wording test with the candidates |
| `10_forms_holdout.csv`, `10_validity_by_form.csv`, `10_revised_form.json` | Step 8 |
| `11_forage.json`, `12_av14*.{json,csv}`, `13_*.{json,csv}` | Steps 9 and 10; exploratory analyses (`13_following_week.csv`: the following week) |
| `14_report_study2.json`, `14_decisions.csv` | Everything, aggregates only; the decisions |
| `99_recovery.csv` | Simulated mode: recovery and exclusion checks |
| `run_times.csv` | Minutes per step, resample counts and cores |

## Runtime

Measured on the simulated data (1,042 enrolled, 1,023 analysable, 299 follow-up pairs), 12 worker processes on this
machine: **44 minutes** end to end with the dry-run resample counts (`SF2_BOOT_H3=50`, `SF2_BOOT_STEP7=20`); see
`run_times.csv`. The parts that grow with the registered counts:

- **H3b and H3d** (`04`): about 42 core-seconds per resample (three ESEMs, each estimated from two starts). 2,000
  resamples: about 2 hours on 12 workers.
- **Step 7** (`09`): about 43 core-seconds per resample. 200 resamples: about 12 minutes on 12 workers.
- **Invariance** (`06`): about 25 minutes whatever the counts, most of it the ESEM-within-CFA models and their item-by-item
  partial-invariance search. In real data these run only when a CFA is not configurally acceptable (Study 1's data
  suggest they will).

A full run with the registered counts should take about 3 hours on 12 workers. Set `SF2_CORES` to use more or fewer;
resample index sets are drawn in the main process with recorded seeds, so results do not depend on the number of cores.

## Validated

On simulated data (2026-10-12, seed 20261099): `output/simulated/99_recovery.csv`. Every planted exclusion, flag and
repeat submission was handled as registered; loadings matched the estimand (Tucker congruence .99-1.00, RMSE .036),
factor correlations within .06, the random-intercept variance 0.204 (SE .016; true 0.22), 97% of the 128 latent and
96% of the 160 observed validity CIs covered the population values, and the planted non-invariant item (sf6_15,
posting) was the one item released. Of that sample's 17 single-sample intervals, 13 covered the population value: its
dedicated group happened to sit 1.5-2.5 SE low on the generating factors, so its three Access/Valuing d's missed.
`98_seed_replication.R` (optional, not in `run_all.R`) repeats the simulation with 8 more seeds: 108 of 112 intervals
(96%) covered, loadings RMSE .035-.040, random-intercept variance 0.194-0.242.

Recovery is judged against the estimand, not only the generating values: with cross-loadings of about .3, the
target-rotated ESEM converges to a rotation of the generating factors (`00_simulate.R`, "The ESTIMAND"), which raised
two factor correlations by about .2 (Access-Letting be .29 to .49; Understanding-Valuing .34 to .50) and shifts some
validity correlations by about .1. That is a property of ESEM, not of the code, and it applies to the real data too.

## Readings of the preregistration (search the scripts for `PREREG:`)

Where the text allowed more than one reading, the code takes the one below. Each should be either confirmed in the
preregistration's text or changed here before filing.

- **Parameterization.** Theta, as registered. With the random intercept's loadings fixed at 1 this is a different model
  from Study 1's ESEM (lavaan's default, delta): under delta every item's standardized agreement loading is equal; under
  theta the agreement factor's unstandardized loading equals each item's residual SD. Section 3's Study 1 fit values
  should say which they used.
- **H1.** "Partial support" is read as: not fully supported, scaled CFI >= .90 and RMSEA <= .08 (the text says CFI
  .90-.95). The random-intercept variance test is lavaan's two-sided Wald p.
- **Model-based omega** (H1 scoring rule, H7, step 7): (sum of keyed loadings)^2 / ((sum)^2 + sum(1 - loading^2)), each
  item's loading on its own factor only, so agreement and cross-loadings count as error; Study 1's formula with signed
  instead of absolute loadings.
- **Factor sign.** Each factor is oriented so the keyed loadings of its targeted items sum positive.
- **Observed scores.** The 80% rule is applied to every observed score, validity scales included (a three-item FFMQ-15
  facet therefore needs all three items); sum-scored scales are prorated (mean x k). BIDR-16 dichotomous scoring =
  number of 6-7 answers after reversal.
- **Validity sample.** "Reached the end of the validity block" = a submission at the 11th validity position (step 13),
  an exited step at or after it, or a completed session.
- **H3b.** SBI and BFI-2 Open-Mindedness are added to the same model (a difference between two correlations needs both).
  **H3d.** One model per factor with its four neighbours. **H4.** One-tailed tests on the latent correlations.
- **H5a.** Gender is a factor (man, woman, other = non-binary or self-described); "prefer not to say" leaves the
  complete-case sample; all steps of a regression use the cases complete on every variable of its last step. Holm is
  applied within each version (primary; + spirituality and religiosity). The latent version regresses the outcome on
  latent Access and Valuing and the covariates (single-indicator latents), with delta-R^2 against the paths fixed at 0.
- **H5b.** Gender as in H5a; posting is contemplative = 1. The odds ratio is per SD of the score in the analysed cases;
  McFadden's R^2 is against the thresholds-only model. "Access and Valuing entered separately" replaces the combined
  score in Step 2 (a 2-df test, each coefficient reported). The Firth model adds the Access + Valuing score to the same
  covariates; its p is logistf's penalized likelihood-ratio p and its CI the profile penalized-likelihood CI.
- **Invariance.** "Configurally acceptable" = the H1 fit criteria. ESEM invariance is fitted as ESEM-within-CFA. The
  Chen criterion is read as both conditions (delta CFI < -.010 AND delta RMSEA > .015); in the dry run the planted
  intercept difference first gave delta CFI -.011 with delta RMSEA +.008, invariant under this reading and not under
  "either". Partial invariance releases, one item at a time, the item whose release gives the best CFI, at most six per
  step. The random intercept's mean is freed with the content factor means at the intercepts step. The reference group is
  the alphabetically first. Sparse categories: the most extreme sparse category is merged first.
- **Country groups.** The six listed countries (bg_country us, uk, ca, au, ie, nz) vs "other"; a declined answer is
  missing; counted in the analysable sample.
- **H6.** Holm is two-tailed; d's CI is the coefficient's CI over the pooled SD. Someone who practises regularly but
  whose years, frequency or training is unknown is "other"; a declined regular-practice answer leaves the group missing.
  Sensitivity (c) is a MIMIC in the ESEM. The frequency trend codes the five categories 0-4.
- **H7.** The decision uses all pairs (out-of-window pairs are "kept and flagged"); 14-21 days = 14 <= days < 22 after the
  main short form. SD_baseline is the pairs' baseline SD (the full-sample SD version is reported too). The agreement
  index's stability uses the short-form pairs only.
- **Agreement index.** Only the designated pairs: Study 1's reversal partners with both members administered and the
  candidate block's opposite pairs (9 reversed, 9 forward). Two undesignated opposites (sf6_47/sfc_5, sf6_40/sfc_6) are
  not used.
- **Step 7.** When more candidates qualify than may be used, they are tried in order of keyed loading; a second
  replacement is judged against the already-revised set. The "Prefer not to answer" rate is the selection half's.
  Bootstrap resamples are stratified by posting; a replacement's selection rate counts it on its factor whatever item it
  replaced. The wording test "in step 7's model" is run on the selection half.
- **Step 8.** "Its validity correlations" = observed-score Pearson correlations in the holdout half.
- **Step 9.** The ESEM is the 20 short-form items + the 5 Forage items (the 44-item model beside it); "a factor of its
  own" = H1's per-factor rule. The Forage score's correlation with the agreement index is also given without the index's
  two Forage pairs. "The rest of BFI-2 Open-Mindedness" = the mean of its 8 Aesthetic Sensitivity and Creative
  Imagination items.
- **Flagged participants.** Section 5 flags cases "for sensitivity analyses" without naming them: H1, the H2-H4 cells,
  H6 and H7 are each repeated without flagged participants.
- **The following week.** Linear regression for all three items (time outdoors as numeric); every linked pair; a pair
  missing the item at either occasion drops out for that item.
- **Step 10.** Re-run items' Study 2 loadings come from the 14-item CFA and from step 7's ESEM.
- **Follow-up invitations.** "Eligible completers" = reached the end of the validity block, not a staff test, not excluded.
  `draw_retest_sample.R` draws 400 by default.
