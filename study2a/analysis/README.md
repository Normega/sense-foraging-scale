# Sense Foraging Study 2a: analysis pipeline (R)

The analysis registered at **osf.io/vt7gn** (filed 29 October 2026, before data collection; `study2a/preregistration_draft.md` is the text it was filed from, and OSF governs): the factor structure of the revised 58-statement set, a check that the two modes were understood, and the rules that choose the statements Study 2b will test. The registration says this code is written during the analysis, following its rules and Study 1's pipeline, and that its commit is cited with the results.

## Data (never in this repo)

- **Pull:** `extract.sql`, run read-only through the Supabase CLI from the radlab checkout:
  `npx supabase db query --linked -o json -f <path>/extract.sql > F:/radlab_data/sf_study2a/responses_<date>.json`
- **What it pulls:** one row per questionnaire submission by a real participant (plus one row for a session with no submissions), keyed by the random participant code. Prolific IDs and staff test runs are excluded.
- **Where it lives:** `SF2A_DATA_DIR` (default `F:/radlab_data/sf_study2a`), local and unsynced. The newest `responses_*.json` is read; scripts write `derived/` and `output/<mode>/` there.
- **Definitions:** `study2a/materials/` holds the live questionnaire definitions as built (radlab `scripts/sense_foraging/study2a/` at 1e23b565): statement wording, hypothesized factors, keying, look-alike pairs and ordering rules (`sf2a-scale.json`), the vignettes (`sf2a-vignettes.json`, `items.json`) and the Study 1 wording items (`sf2a-end.json`).

## Modes (`SF2A_MODE`)

| Mode | What the answers are | Use |
|---|---|---|
| `permuted` (default) | Real answers, each item shuffled across participants (statements, vignettes, IUS-12, end items): distributions kept, correlations destroyed | Running the code end to end without seeing structure |
| `simulated` | Statements generated from a known 5-factor model with agreement (variance .10), cut to the real statements' distributions; vignettes from a known per-person separation | Proving the code recovers structure; `99_check_recovery.R` reports congruence with the truth |
| `real` | The answers | The analysis |

Exclusions, attention checks, positions, timing and background answers are always the real values. Answer-derived indices (longstring, separation) are computed on whichever answers the mode uses.

## Steps (`run_all.R`, from the repo root)

`SF2A_K` sets the number of content factors in the primary model: default 4, decision K, made 30 October 2026 and recorded with its reasons in `config.R`. Simulation checks ran with 5 (`SF2A_K=5`). `SF2A_FROM=<step>` starts part-way.

1. **`01_prepare`:** one row per participant, by the registered rules. The first complete submission of each questionnaire is used; anyone failing two or more of the three checks is excluded. Flags: one check failed, < 2 s per statement (61 positions), longstring ≥ 15 (61 positions), > 11 of 58 declined; also the step-crash group (not registered). Separation and "tells the modes apart" are computed here.
2. **`02_items`:** distributions, floor and ceiling, missingness, agreement rates (4–6, unrecoded), and the 18 look-alike pairs' polychoric correlations.
3. **`03_vignettes`:** decision V. Sample mean sensing and doing ratings, the midpoint, each vignette against it, the share who tell the modes apart (Wilson CI) against 75%.
4. **`04_factor_number`:** decision K. Parallel analysis (polychoric, keyed, principal-axis, 500 datasets, 95th percentile); ordinal EFAs for k = 4, 5, 6 and the PA count; criterion 1, the hypothesized factor each factor matches, a method-factor check, fit.
5. **`05_primary_model`:** ESEM on the unrecoded statements with the random intercept (target: hypothesized for k = 5, else 04's geomin EFA), and without it; factor scores; decision F (Forage).
6. **`06_wording`:** decision L. The primary model plus a wording factor on the 19 reverse-keyed statements (or the wording factor in place of the random intercept if that fails or is improper).
7. **`07_ius`:** IUS-12 total with Letting be's forward and recoded reverse-keyed means; dependent-correlation test.
8. **`08_reliability`:** ordinal omega (model-based) and ordinal alpha per factor.
9. **`09_selection`:** decision S, a proposal for the PI: 5–6 statements per factor by the registered preferences.
10. **`10_sensitivity`:** congruence with the primary solution for participants who tell the modes apart, without flagged participants, without the step-crash group, without the random intercept, and a forward-only EFA; position effects on means and item-rest correlations.
11. **`11_exploratory`:** self-rating with factor scores; reworded trial statements with their Study 1 wording; the trial's Access + Valuing score with the Access and Valuing factors; current vs never practitioners (d) per statement and factor.
12. **`99_check_recovery`** (simulated, k = 5): congruence with the true loadings, the agreement variance, separation against the truth.

Model settings in every fit (`sf2a_models.R`, adapted from Study 2's): ordinal, WLSMV, theta parameterization, `missing = "pairwise"` (the registration says pairwise; Study 1 left lavaan's listwise default), and each ESEM fitted from lavaan's default starts and from its simple-structure CFA, keeping the lower minimum.

## Readings of the registration (`PREREG:` in the code)

Where the text does not settle a detail, the code takes one reading and marks it `PREREG:`:

- **Sample mean ratings (03):** the mean over participants of each person's mean sensing (doing) rating, among those with a separation score.
- **Primary loading (05, 06, 08, 09):** a statement's primary factor is the content factor with its largest absolute keyed loading; it counts toward "≥ .40" only in the keyed direction.
- **Factors when k ≠ 5 (05, 06, 11):** a factor is identified with the hypothesized factor most of its defining statements (keyed primary ≥ .40) were written for.
- **Factor scores (05):** lavaan's EBM scores; a case lavaan does not score is NA.
- **IUS-12 (07):** "more strongly" = larger in absolute value; Williams' T2 via `psych::r.test`, as Steiger (1980) recommends.
- **Reliability (08):** every content factor of the primary model, over the statements defining it; model-based omega with agreement and cross-loadings as error (Study 2's formula).
- **Selection (09):** .40 is the only hard minimum; "qualifies" = primary ≥ .40, rules 2–3 met. Candidates are ranked: qualifies, agreement unflagged, (if V failed) not naming a mode, ≥ .50, loading. Qualifying, unflagged short-form statements (source `sf6_*`) are kept first; the factor is filled to six with qualifying statements, then to five with any; a look-alike partner of a kept statement is skipped; reverse-keyed statements are swapped in for the lowest-ranked non-short-form forward ones until two are kept, where two qualify.
- **Trial score (11):** the mean of its 10 statements in Study 1 wording, of those answered.
