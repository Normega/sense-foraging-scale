# Sense Foraging Scale, Study 2b: confirmatory structure and construct validity of the 8-statement scale

**DRAFT, 9 November 2026; not filed.** Rewritten from `prereg_2b_plan.md` and decisions 1–36 in `decisions_after_2a.md`. The pre-2a draft (20 items, four factors, retest) is in the repository history.

- OSF Preregistration (standard template).
- University of Toronto REB Protocol #00051180; the amendment covering Study 2b is approved.
- **Before filing:** the 2b analysis code (`study2b/analysis/`) is written and has run end to end on data simulated from Study 2a's fitted model, with every recovery check passing (9 November). The code is the specification wherever this text is less precise; its readings of this text are listed in `study2b/analysis/README.md`. The run on real answers shuffled within each item needs real data. It is run on the pilot's extraction, and this registration is filed after it, before the remaining 640 places open (Section 3, Pilot).
- **Drafting choices for Norm to confirm** are listed in Section 6.

## 1. Study information

### Title

Development of the Sense Foraging Scale, Study 2b: confirmatory factor structure, convergent, discriminant and criterion validity of the 8-statement scale

### Authors

Norman Farb (principal investigator), University of Toronto. Co-investigators as listed on REB Protocol #00051180.

### Description

The Sense Foraging Scale measures how readily people move into, and value, the sensing mode: taking in what is around them and in their body with curiosity, without needing to do anything with it. Its counterpart is the doing mode: using what one knows to plan, solve, judge or get something done.

- **Study 1** (osf.io/48r2y; N = 500) found four correlated content factors and an agreement (acquiescence) factor in a 48-item pool.
- **Study 2a** (osf.io/vt7gn; N = 500) tested a revised set of 58 statements. It found one broad sensing factor, small reverse-keyed clusters and the agreement factor. The hypothesized parts (Access, Valuing, Forage, Understanding, Letting be) did not separate as subscales. From 2a we chose, in-sample:
  - **the 8**, the scale scored as one number: N4, sf6_36 (R), V10, sf6_27, sfc_36, sfc_x2 (R), sfc_35, U4 (R);
  - **the 10**, the 8 plus A9 and sfc_24, for studies of change.

Study 2b is the first test of the 8 in data it was not chosen on. It asks:
1. **Structure:** does one content factor plus the agreement factor fit the 8?
2. **Convergent validity:** does it relate as predicted to its closest established neighbours?
3. **Discriminant validity:** is it distinct from those neighbours, from neighbouring constructs that a reader might mistake it for (sensation seeking, need for cognition, sensory overwhelm, anxious monitoring of the body, finding meaning in events), and from socially desirable responding?
4. **Criterion validity:** does it predict a behaviour in the session, taking an optional pause to look and listen, and a reported behaviour, deliberately pausing in the past week?
5. **Incremental and known-groups validity:** does it add to mindfulness and interoceptive awareness in explaining wellbeing and practice, and do practitioners score higher?

### Hypotheses

Confirmatory hypotheses are H1–H6 on the 8. The 10 is reported beside it under the same rules, without its own hypotheses. H8 is exploratory (decision 30). Decision rules are in Section 5.

- **H1 (structure).** A one-factor model of the 8 with a random-intercept agreement factor fits.
- **H2 (convergent).** The 8 correlates positively and at least moderately with each of FFMQ-15 Observing, SBI Savoring the Moment, BFI-2 Aesthetic Sensitivity and Brief MAIA-2 Noticing.
- **H3 (discriminant).**
  - (a) Each H2 correlation is below .80: the scale is not any one of those measures.
  - (b) Five within-battery contrasts, each a larger correlation with the first measure than the second:
    1. BFI-2 Aesthetic Sensitivity > BFI-2 Intellectual Curiosity (not general curiosity about ideas);
    2. FFMQ Observing > FFMQ Describing (not mindfulness in general, or labelling);
    3. SBI Savoring the Moment > BFI-2 Open-Mindedness, total (not a general appetite for experience);
    4. Brief MAIA-2 Noticing > Brief MAIA-2 Not-Worrying (not low distress about the body);
    5. FFMQ Observing > SAMD Meaning-Detecting (not a spiritual or meaning-seeking openness).
  - (c) Its correlations with BIDR-16 Self-Deceptive Enhancement and Impression Management lie within ±.30.
  - (d) Its multiple correlation with the four H2 neighbours together is below .80: it is not those neighbours combined.
  - (e) Its correlations with the four new discriminant measures are bounded: BSSS-8 below .40; NCS-6, SPSQ-SF Sensory Discomfort and ASI-3 Physical Concerns within ±.30.
- **H4 (criterion).**
  - (a) Higher scores go with more past-week days of deliberately pausing to take in one's surroundings (`crit_pause_days`).
  - (b) Higher scores go with taking the optional pause in the session and staying with it (an engaged pause, Section 4).
- **H5 (incremental).**
  - (a) The 8 explains wellbeing beyond age, gender, the five FFMQ-15 facets and the eight Brief MAIA-2 subscales, for at least one of life satisfaction, PHQ-4 and PSS-4.
  - (b) The 8 predicts contemplative practice frequency beyond age, gender, spirituality, religiosity, the FFMQ-15 facets and the Brief MAIA-2 subscales.
- **H6 (known groups).** Current practitioners score higher than people who have never practised regularly, by at least d = 0.50. In 2a, in-sample, the 8's d was 0.67.

## 2. Design plan

### Study type

Observational study: a single-session online survey.

### Blinding

None. The optional pause is presented as a break; the debrief explains that whether it was taken, and for how long, was recorded, and why that was not said beforehand (approved by the REB).

### Study design

One session of about 22 minutes, plus any pause taken, in this order:

1. **Background questions:** the 13 questions from Study 1, unchanged.
2. **The two modes** (`sf2b-vignettes`): a three-paragraph introduction to the doing and sensing modes; nine brief vignettes, three a page in a per-participant order with each page holding at least one of each kind, rated 1 Doing mode to 6 Sensing mode, without feedback; then a self-rating on the same anchors. The introduction, V2 and V7 are reworded from 2a (decision 6), so separation scores are not directly comparable with 2a's.
3. **The scale** (`sf2b-scale`): the 10 (with the 8 inside it), three unscored reverse-keyed statements (N6, sfc_31, sf6_44) and one attention check, in one fixed order on two pages (Appendix A). Response scale 1 Strongly disagree to 6 Strongly agree, no midpoint, plus "Prefer not to answer".
4. **The PSY240 trial's other statements** (`sf2b-end`): eight statements in the trial's own wording and one attention check, fixed order. With sf6_36 and sf6_44 in step 3, they complete the trial's 10-statement Access + Valuing score (decision 13).
5. **Validity block, first half,** in a per-participant random order: FFMQ-15; BIDR-16; SBI Savoring the Moment; BSSS-8; NCS-6; PHQ-4; SPSQ-SF Sensory Discomfort; PSS-4; life satisfaction (one item, six-point face scale).
6. **The pause offer** (`sf2b-pause`): "If you'd like, you can take a short break here to look, listen and feel what's around you. There's nothing to do and nothing to report. It doesn't affect your payment." with two buttons, Take a pause and Continue the survey, their sides fixed per participant at random.
   - **Take a pause:** a quiet screen with "Take as long as you like. Tap 'I'm back' when you're ready." After 3 minutes on the page the line becomes "Whenever you're ready." After 10 minutes in all the pause ends itself. Then "Anything you noticed? (optional)".
   - **Continue the survey:** "Any reason? (optional)": No time right now / Not a good place / Not interested / Something else / Skip.
   - **Then, either way:** "Where are you right now?" (At home / At work or school / Outdoors / Travelling / Somewhere else) and "How tired are you right now?" (1 Not at all to 6 Very), each with "Prefer not to answer", asked after the decision so they cannot steer it.
7. **Validity block, second half,** in its own per-participant random order: Brief MAIA-2; BFI-2 Open-Mindedness (with its attention check); EQ Decentering; SAMD Meaning-Detecting; ASI-3 Physical Concerns; three past-week items.
8. **Debrief.**

Every question about the participant offers "Prefer not to answer"; the three attention checks do not. The halves of the validity block are fixed; past-week behaviour comes after the pause so its question about deliberately pausing cannot prompt the offer. The pause sits at least 5 minutes after the scale.

**Recruitment.** One Prolific posting with a neutral, general description (decision 24, amended by 34). No screening on practice: practice is measured in the background questions after consent.

### Randomization

- The order of each half of the validity block is a seeded random permutation per participant, stable on reload and recorded (`step_index` on every response).
- The vignettes' order is per participant within the page constraints, and recorded.
- The pause buttons' sides are set per participant from a hash of their link and recorded.
- The scale, the trial statements and everything else are in fixed order.

## 3. Sampling plan

### Existing data

**Registration prior to analysis of the data.** At filing, the only Study 2b data are the pilot's 60 sessions. Their answers have been seen only shuffled within each item (Pilot, below), which leaves no relation between any two answers. The 8 and the 10 were chosen on Study 2a's data, a separate sample. In 2a, in-sample (decision 14; WLSMV, theta, pairwise, the agreement factor as in Section 5):

| Form | CFI | RMSEA (90% upper) | SRMR | Largest residual | ω (content) | α | Practitioner d |
|---|---|---|---|---|---|---|---|
| The 8 | .985 | .061 (.081) | .031 | .079 | .838 | .80 | 0.67 |
| The 10 | .987 | .048 (.063) | .030 | .077 | .850 | .83 | 0.69 |

The same 8 in 2a: vignette self-rating r = .33; IUS-12 r = −.18. Study 1 and 2a had no other validity measures. The 8's statements are the PSY240 trial's measure only in part (sf6_36), and the trial's data are not used here.

### Data collection procedures

- **Participants:** Prolific participants aged 18+, fluent in English (Prolific's fluent-languages filter), in any country. Participants in Study 1 or Study 2a are excluded through Prolific's exclusion list.
- **Platform:** the RADlab research platform (radlab.zone). Prolific IDs are held in a separate table, replaced by a random code in all study data, and deleted after payment.
- **Compensation:** £12.71 an hour (the UK minimum wage), prorated, set from the pilot's median completion time (about £4.66 for 22 minutes), paid regardless of answers, of the attention checks and of the pause. Prolific's maximum time is set well beyond the median plus the 10-minute pause cap, so a long pause does not time out.
- **Pilot:** the first 60 places, identical materials, run before this registration is filed. Its extraction goes through the analysis code with every item's answers shuffled across participants, which keeps each item's distribution and destroys every relation between answers. Read from that run only:
  - completion and breakoff counts;
  - attention-check pass rates;
  - completion time;
  - each item's "Prefer not to answer" rate;
  - pause take-up;
  - practice-group sizes;
  - whether the code reads the export.

  No unshuffled answers are analysed. The rules for each reading are set in advance in `study2b/analysis/README.md` ("The pilot run"). If nothing changes after the pilot, pilot participants are included. If the materials change, they are excluded and 60 places are added.

### Sample size

Target: 700 completed sessions (60 pilot + 640), in one posting (decision 34).

### Sample size rationale

By simulation (`study2/power/`, 500 simulated studies per condition at N = 700, generated from 2a's fitted model of the 8 with its agreement factor and thresholds; validity measures as latent factors of 3–8 ordinal items):

| Test | Assumed truth | Power at 700 |
|---|---|---|
| H1, the decision rule | the model true / data like 2a's | 100% / 97% pass |
| | three statements sharing a second factor loading .60 (should fail) | 0.2% pass |
| H2 | true latent r .45 / .50 | 95% / 99% |
| H3b, a difference of .15 | 3-item measures correlating .30 / .50 | 76% / 86% |
| | 4- or 8-item measures correlating .30 | 84–86% |
| | 3-item measures barely correlated (contrast 4, Brief MAIA-2) | 64% |
| H3d, upper bound < .80 | true R ≤ .67 / .71 / .73 / .76 | 100% / 92% / 83% / 42% |
| H4b | latent r .20, take-up 15–50% | 88–97% |
| H5a | true ΔR² .03, at 700 / at the 410–470 complete cases of the 80% rule | 84–85% / 76–78% |
| | true ΔR² .04, at the 410–470 complete cases | 92–93% |
| H6 | true d .67 (382 vs 190) | 97% |

The smallest contrast of interest in H3b is a difference of .15 (decision 32). Contrast 4 compares two 3-item subscales that barely correlate, and has about 64% power; it is kept and its power reported (decision 36). H2 needs a true correlation of about .45 for 95% power (63% at .40).

H3d's power depends on the true multiple correlation R. That depends on two things no study has measured: the 8's correlation with each neighbour (.40–.50 simulated), and how much the neighbours share beyond the 8 (`study2/power/power_2b_h3d.R`; 300 simulated studies per condition).
- **Neighbours sharing at least .30 with each other:** R is at most .73, and power is 83–100%.
- **Neighbours related only through the 8, at .50 each:** R is .76, and power is 42%.

The simulation used the delta-method interval, because 2,000 bootstrap resamples for every simulated study would take days. In 16 studies that also had 300 bootstrap resamples, its upper bound was within .004 of the bootstrap's on average, and every decision was the same.

### Stopping rule

The number of places is fixed on Prolific in advance: 60, then 640. Data collection stops when these are filled. There is no interim analysis, and excluded participants are not replaced.

## 4. Variables

### Manipulated variables

None.

### Measured variables

- **The scale:** `sf2b_01` to `sf2b_13` (1–6 or "Prefer not to answer"), stored unrecoded. Appendix A gives each one's source id, role (8, 10 or exploratory) and keying.
- **Attention checks:** `sf2b_attn_agree` (correct = 5), `sf2b_attn_disagree` (correct = 2), `bfi2_attn_check` (BFI-2's).
- **The trial's statements:** `sf6_9`, `sf6_35`, `sf6_7`, `sf6_8`, `sf6_12`, `sf6_18`, `sf6_27`, `sf6_11` (1–6), with `sf2b_02` (sf6_36) and `sf2b_11` (sf6_44) in the scale.
- **Vignettes:** `vig_1` to `vig_9` and `vig_self` (1–6).
- **FFMQ-15:** Observing, Describing, Acting with awareness, Nonjudging, Nonreactivity (3 items each, 1–5).
- **Brief MAIA-2** (Rogowska, Tataruch & Klimowska, 2023): eight subscales of 3 items each (0–5), MAIA-2 items 2–4, 6, 9, 10, 11, 12, 15, 17, 19, 20, 25–27, 29–37; items 6, 9, 10, 11, 12, 15 reverse-scored.
- **SBI Savoring the Moment:** 8 items (1–7), items 2, 8, 14, 20 reverse-scored.
- **BFI-2 Open-Mindedness:** 12 items (1–5), facets Intellectual Curiosity, Aesthetic Sensitivity, Creative Imagination (4 each).
- **EQ Decentering:** 11 items (1–5).
- **BIDR-16:** Self-Deceptive Enhancement and Impression Management, 8 items each (1–7), scored continuously.
- **BSSS-8** (Hoyle et al., 2002): 8 items (1–5), mean.
- **NCS-6** (Lins de Holanda Coelho, Hanel & Wolf, 2020): 6 items (1–5), items 3 and 4 reverse-scored, mean.
- **SPSQ-SF Sensory Discomfort** (De Gucht & Woestenburg, 2025): 4 items (1–7), mean.
- **ASI-3 Physical Concerns** (Taylor et al., 2007): items 3, 4, 7, 8, 12, 15 (0–4), sum.
- **SAMD Meaning-Detecting** (Russo-Netzer & Icekson, 2023): 13 items (1–7), mean.
- **PHQ-4:** 0–12. **PSS-4:** 0–16, items 2 and 3 reverse-scored. **Life satisfaction:** one item, 1–6.
- **Past week:** `crit_pause_days` (0–7); `crit_notice_first` (1–5, "nothing unexpected happened" missing); `crit_outdoors` (0–4).
- **The pause** (`pause_1`): the choice; decision time; the decline reason; for a pause, seconds on the page (tab visible and focused), unfocused and hidden, how it ended (button or the 10-minute cap) and whether the 3-minute nudge appeared; the optional note; setting; tiredness; device class.
- **Background questions:** as in Study 1.
- **Paradata:** each response's step position; seconds per questionnaire.

### Indices

- **The 8 and the 10 (observed):** the mean of answered statements after reverse-scoring (7 − x), when at least 7 of the 8 (9 of the 10) are answered.
- **The 8 and the 10 (latent):** the content factor of the H1 model.
- **Engaged pause (H4b):** 1 when the participant chose the pause, stayed at least 20 seconds on the page, and was away (unfocused + hidden) for less than half the pause's total time; otherwise 0 (decision 22). A participant who did not reach the offer is missing.
- **Separation (vignettes):** mean rating of the sensing vignettes minus mean of the doing vignettes, for participants with at least three of each; **tells the modes apart** = separation of 1 or more (2a's rule).
- **Practice groups (H6):** current practitioner = regular practice "Yes" and current frequency 1–2 times a week or more; never = regular practice "No" (2a's definitions). Dedicated (3+ a week, a year or more, teacher-led training) is a sensitivity group.
- **Practice frequency (H5b):** five ordered levels: none in the past three months or no regular practice; less than weekly; 1–2 a week; 3–5 a week; daily or near-daily.
- **The trial's Access + Valuing score:** the mean of its 10 statements in its own wording (sf6_9, sf6_36 R, sf6_35 R, sf6_7, sf6_8, sf6_12, sf6_18, sf6_27, sf6_11, sf6_44 R).
- **Longstring:** the longest run of identical answers on the scale's 14 positions, in presented order.

## 5. Analysis plan

### Statistical models

R with lavaan, psych and MASS; versions reported. Scale and validity items are ordinal: WLSMV, theta parameterization, pairwise-available data (`missing = "pairwise"`, set explicitly). Fit indices are the scaled ones. The agreement factor is a random intercept with loadings fixed at 1 on every unrecoded statement of the scale in the model, a free variance, and no correlation with any content or validity factor (decision 5).

1. **H1, structure.** One content factor on the 8 plus the agreement factor. **H1 is supported** when all three hold (decisions 17, 33):
   - CFI ≥ .95;
   - local fit: no residual correlation (observed minus model-implied polychoric) exceeds .10 in absolute value; or, if one does, freeing the largest as a residual covariance changes the content ω by less than .03 and every standardized content loading by less than .05;
   - the unbiased SRMR close-fit test (Maydeu-Olivares, Shi & Rosseel, 2018; lavaan `lavResiduals`, H₀: SRMR ≤ .05) is not rejected at α = .05.

   Reported, not tested: RMSEA with its 90% CI, model-specific (dynamic) cutoffs, every loading, the agreement factor's variance and share of item variance, ordinal ω for the content factor, α, and the 10 under the same rule. The model without the agreement factor is reported for comparison.
2. **H2–H3, convergent and discriminant.** One model per validity measure: the H1 model plus that measure as a latent factor of its items (ordinal), correlated with the content factor and uncorrelated with the agreement factor. Subscales of one measure (FFMQ, Brief MAIA-2, BFI-2, BIDR) enter as correlated factors of that measure. Correlations have 95% Wald CIs.
   - **H2:** each of the four lower bounds exceeds .30.
   - **H3a:** each of the four upper bounds is below .80, the cutoff for "no problem" in Rönkkö and Cho's (2022) classification of the factor correlation's 95% upper limit.
   - **H3b:** in the model holding both measures of a contrast, the difference of the two latent correlations, defined as a parameter, has a 95% delta-method CI excluding 0 in the predicted direction. Each contrast is reported separately; H3b is supported when at least four of the five are. (Contrast 4 has about 64% power, Section 3.)
   - **H3c:** both BIDR CIs lie within ±.30.
   - **H3d:** in a model with the four H2 neighbours, the multiple correlation R = √(r′Σ⁻¹r) of the content factor with them, with a 2,000-resample bootstrap 95% CI; supported when its upper bound is below .80. Rönkkö and Cho's cutoff is for one correlation; applying it to R, which is at least as large as the largest of the four, is our extension and a stricter test.
   - **Reported, not tested:** each H3a upper bound and H3d's, classed by Rönkkö and Cho's (2022) Table 12: no problem (< .80), marginal (.80 to < .90), moderate (.90 to < 1), severe (≥ 1).
   - **H3e:** BSSS-8's upper bound below .40; NCS-6's, SPSQ-SF Sensory Discomfort's and ASI-3 Physical Concerns' CIs within ±.30.
   - **Secondary:** observed-score Pearson correlations with every measure, the full validity matrix below.
3. **H4, criterion.** On the observed 8:
   - (a) Spearman ρ with `crit_pause_days`;
   - (b) Spearman ρ with the engaged pause (0/1), and a logistic regression of the engaged pause on the 8 with device class as a covariate, reported beside it.
   - Each is supported when ρ is positive with p < .05 (two-sided), Holm-corrected across (a) and (b).
4. **H5, incremental.** Observed scores; complete cases.
   - **H5a:** for each of life satisfaction, PHQ-4 and PSS-4, a hierarchical regression: step 1 age, gender, the five FFMQ-15 facets and the eight Brief MAIA-2 subscales; step 2 the 8. Supported for an outcome when ΔR² ≥ .02 and the step-2 F-test is significant after Holm correction across the three; H5a holds when at least one outcome is supported.
   - **H5b:** a proportional-odds regression of practice frequency: step 1 age, gender, spirituality, religiosity, the FFMQ-15 facets and the Brief MAIA-2 subscales; step 2 the 8. Supported when the coefficient is positive with likelihood-ratio p < .05.
5. **H6, known groups.** Current practitioners against never-practitioners on the observed 8: Cohen's d (pooled SD) with 95% CI and a two-sided t-test. Supported when the difference favours practitioners, p < .05 and d ≥ 0.50. **Sensitivity:** adjusting for spirituality and religiosity; dedicated practitioners against never; the latent difference with agreement modelled.

### Predicted validity matrix

Latent correlations with the 8. "Near zero" is |r| < .20. **H** marks confirmatory cells.

| Measure | Prediction | H |
|---|---|---|
| FFMQ Observing | + moderate (.30–.60) | H2, H3a, H3b-2, H3b-5, H3d |
| FFMQ Describing | + small | H3b-2 |
| FFMQ Acting with awareness, Nonjudging | + small | |
| FFMQ Nonreactivity | + small to moderate | |
| SBI Savoring the Moment | + moderate | H2, H3a, H3b-3, H3d |
| BFI-2 Aesthetic Sensitivity | + moderate | H2, H3a, H3b-1, H3d |
| BFI-2 Intellectual Curiosity | + small | H3b-1 |
| BFI-2 Open-Mindedness, total | + small to moderate, below SBI | H3b-3 |
| BFI-2 Creative Imagination | + small | |
| Brief MAIA-2 Noticing | + moderate | H2, H3a, H3b-4, H3d |
| Brief MAIA-2 Not-Worrying | near zero | H3b-4 |
| Brief MAIA-2 Attention Regulation, Body Listening, Emotional Awareness, Self-Regulation | + small to moderate | |
| Brief MAIA-2 Not-Distracting, Trusting | + small | |
| EQ Decentering | + small to moderate | |
| SAMD Meaning-Detecting | + small to moderate, below FFMQ Observing | H3b-5 |
| BSSS-8 | + small | H3e (< .40) |
| NCS-6 | near zero | H3e (±.30) |
| SPSQ-SF Sensory Discomfort | near zero | H3e (±.30) |
| ASI-3 Physical Concerns | near zero | H3e (±.30) |
| BIDR SDE, IM | near zero to + small | H3c (±.30) |
| PHQ-4, PSS-4 | near zero to − small | |
| Life satisfaction | + small | |
| Past week: pause days | + small to moderate | H4a |
| Engaged pause | + small | H4b |
| Past week: notice first, time outdoors | + small | |

A large negative correlation with PHQ-4 or PSS-4 would mean the scale is partly a wellbeing measure; it is reported, not tested. The agreement factor's correlations with every measure are reported.

### Transformations

Reverse-keyed statements are recoded (7 − x) for observed scores and left unrecoded in models with the agreement factor. Validity scales are scored per their manuals, with a score computed when at least 80% of its items are answered (Study 2's rule), so a 3-item subscale needs all three. "Prefer not to answer" is missing.

### Inference criteria

Decision rules for H1–H6 are above. Two-tailed α = .05 unless stated, Holm-corrected within H4 and H5a. Everything else is exploratory and labelled so.

### Data exclusion

- **Excluded:** staff test enrollments; sessions without a submitted scale; participants who fail two or more of the three attention checks.
- **Partial sessions:** a check not reached is not a failure. A participant who submitted the scale is kept for H1; each later test uses those who reached its measures. The pause analyses use those who reached the offer.
- **Repeat submissions:** the first complete submission of each questionnaire; flagged re-submissions are dropped by the platform's resubmission rule.
- **Flagged for sensitivity analyses, not excluded:** failing exactly one check; fewer than 2 seconds per statement on the scale; a longstring of 10 or more of the scale's 14 positions; "Prefer not to answer" on more than 2 of the 13 statements.

### Missing data

"Prefer not to answer" and unanswered items are missing. Factor models use pairwise-available data. Observed scores need the item counts in Section 4. Regressions use complete cases. Missingness per measure is reported.

### Exploratory analysis

- **H8, understanding (decision 30).** Invariance of the H1 model across people who do and do not tell the modes apart: configural, then thresholds, then loadings (Wu & Estabrook, 2016), with the agreement factor's loadings fixed at 1 in both groups and its variance free. A step is non-invariant when ΔCFI < −.010 and ΔRMSEA > .015. In 2a, in-sample, the loadings differed. If they differ again, the score means something different for people who do not tell the modes apart.
- **The 10:** H1–H6's statistics for the 10, beside the 8's.
- **The three exploratory statements** (N6, sfc_31, sf6_44): their loadings when added to the 10, and their correlations with the agreement factor.
- **The trial's measure:** the correlation of the 8 with the trial's Access + Valuing score in its own wording, and of each with the validity measures, as a check that the 8 can stand in for the trial's measure (decision 13).
- **The pause, further:** seconds on the page, log(1 + x); acceptance alone; the 8 in decliners giving "not interested" against "no time" or "not a good place"; the pause against tiredness, setting and device; the 8's association with the pause compared with tiredness's; the notes, coded blind to scores.
- **Order effects:** a measure's position within its half against its mean and its correlation with the 8.
- **The vignettes:** separation, the self-rating's correlation with the 8, and each vignette's mean.
- **Posting-free invariance:** by gender (men and women) and by country group (the six countries in the REB application against everywhere else, if each has at least 200), as in H8.
- **BIDR dichotomous scoring** (Paulhus) for H3c.

## 6. Other

**Relation to the approved protocol.** The amendment covering Study 2b is approved (protocol #00051180). It covers the vignettes, the 8 and the 10, the trial statements, the battery with BSSS-8, NCS-6, SPSQ-SF Sensory Discomfort, ASI-3 Physical Concerns, SAMD Meaning-Detecting and the Brief MAIA-2, the pause and its debrief, one posting and the sample of 700.

**Measures used for non-commercial academic research.** SBI, PSS-4, BSSS-8, NCS-6, the SPSQ-SF and ASI-3 items, the Brief MAIA-2, BFI-2, EQ, BIDR-16, FFMQ-15 and PHQ-4 are used under their research-use terms; SAMD is CC BY, and its authors are told of the use.

**Drafting choices for Norm to confirm** (not in the decisions so far):
- H3b holds when at least four of the five contrasts hold; each is also reported on its own.
- H3e's bounds: BSSS-8 below .40, the other three within ±.30.
- Validity measures enter the latent models as latent factors of their items, as simulated, not as observed scores.
- H3b uses the delta-method CI, as simulated (the pre-2a draft used a 2,000-resample bootstrap); H3d keeps the bootstrap.
- H4 is two-sided, Holm across (a) and (b); H6 uses 2a's current-practitioner group, with dedicated as sensitivity.
- The pilot's 60 are included if nothing changes.
- Flags: a longstring of 10+ of 14; more than 2 of 13 statements unanswered.
- **Decided 10 November:** the pilot is run first. The shuffled-answers run of the code is on its extraction, and this registration is filed after it. No unshuffled pilot answers are analysed.
- **Decided 10 November:** validity scores keep Study 2's 80% rule (Section 4, Transformations), so a 3-item FFMQ-15 or Brief MAIA-2 subscale needs all three items. H5 needs all 13 such subscales. In the simulated run, with 1% of items unanswered at random, its complete-case sample fell from about 650 to 410–470; Section 3 gives H5's power at that n.
- **Decided 12 November:** H3a's and H3d's bound is .80, after Rönkkö and Cho (2022), whose upper limits are reported in their classes (Section 5). H3d's power is 92–100% when the true R is .71 or below, 83% at .73 and 42% at .76 (Section 3). The low-power case is the one in the dry run of the code: neighbours related to each other only through the 8.

**Data sharing.** A de-identified dataset will be shared as described in the consent form: no exact timestamps, age in bands, coarsened rare demographic combinations, free text (the pause notes) recoded. The analysis code will be shared alongside it.

### References

- Chen, F. F. (2007). Sensitivity of goodness of fit indexes to lack of measurement invariance. *Structural Equation Modeling, 14*(3), 464–504.
- De Gucht, V., & Woestenburg, D. H. A. (2025). Development and validation of the SPSQ-26: A short form of the Sensory Processing Sensitivity Questionnaire. *European Journal of Psychological Assessment, 41*(5), 348–356.
- Hoyle, R. H., Stephenson, M. T., Palmgreen, P., Lorch, E. P., & Donohew, R. L. (2002). Reliability and validity of a brief measure of sensation seeking. *Personality and Individual Differences, 32*(3), 401–414.
- Lins de Holanda Coelho, G., Hanel, P. H. P., & Wolf, L. J. (2020). The very efficient assessment of need for cognition: Developing a six-item version. *Assessment, 27*(8), 1870–1885.
- Maydeu-Olivares, A., & Coffman, D. L. (2006). Random intercept item factor analysis. *Psychological Methods, 11*(4), 344–362.
- Maydeu-Olivares, A., Shi, D., & Rosseel, Y. (2018). Assessing fit in structural equation models: A Monte-Carlo evaluation of RMSEA versus SRMR confidence intervals and tests of close fit. *Structural Equation Modeling, 25*(3), 389–402.
- Mehling, W. E., Acree, M., Stewart, A., Silas, J., & Jones, A. (2018). The Multidimensional Assessment of Interoceptive Awareness, Version 2 (MAIA-2). *PLoS ONE, 13*(12), e0208034.
- Rogowska, A. M., Tataruch, R., & Klimowska, K. (2023). Validation of the shortened 24-item multidimensional assessment of interoceptive awareness, version 2 (Brief MAIA-2). *Scientific Reports, 13*, 21270.
- Rönkkö, M., & Cho, E. (2022). An updated guideline for assessing discriminant validity. *Organizational Research Methods, 25*(1). https://doi.org/10.1177/1094428120968614
- Russo-Netzer, P., & Icekson, T. (2023). An underexplored pathway to life satisfaction: The development and validation of the synchronicity awareness and meaning-detecting scale. *Frontiers in Psychology, 13*, 1053296.
- Taylor, S., Zvolensky, M. J., Cox, B. J., Deacon, B., Heimberg, R. G., Ledley, D. R., et al. (2007). Robust dimensions of anxiety sensitivity: Development and initial validation of the Anxiety Sensitivity Index-3. *Psychological Assessment, 19*(2), 176–188.
- Wu, H., & Estabrook, R. (2016). Identification of confirmatory factor analysis models of different levels of invariance for ordered categorical outcomes. *Psychometrika, 81*(4), 1014–1045.
- Other measure citations: as in Study 2's *Study2_Measures.docx*.

## Appendix A. The scale, in presented order

Response scale: 1 Strongly disagree, 2 Disagree, 3 Mildly disagree, 4 Mildly agree, 5 Agree, 6 Strongly agree; "Prefer not to answer". Instruction: "In general, how much do you agree with each statement? There are no right or wrong answers."

| Page | Id | Source | Statement | Role | Reversed |
|---|---|---|---|---|---|
| 1 | sf2b_01 | N4 | As I go about my day, I notice the light, colours and sounds around me with interest | 8 | |
| 1 | sf2b_02 | sf6_36 | Moving into the sensing mode is a struggle for me | 8 | yes |
| 1 | sf2b_03 | V10 | I enjoy taking a moment just to notice how my body feels, like my breathing, the soles of my feet, or posture | 8 | |
| 1 | sf2b_04 | sf6_27 | I enjoy exploring the sights, sounds and textures of ordinary moments | 8 | |
| 1 | sf2b_05 | N6 | When I walk somewhere familiar, I notice very little of what I pass | exploratory | yes |
| 1 | sf2b_06 | sfc_36 | Even in chores I have done many times, I notice sounds, smells or textures I had missed | 8 | |
| 1 | sf2b_07 | U4 | Even when I pause to take in what is around me, I react to stress the same way | 8 | yes |
| 2 | sf2b_08 | A9 | I can pay attention to what I'm feeling in my body in new or unfamiliar places | 10 | |
| 2 | sf2b_09 | sfc_x2 | It never seems to be the right time to go into the sensing mode | 8 | yes |
| 2 | sf2b_attn_agree | | To show that you are reading each statement, please select "Agree" for this one. | check (5) | |
| 2 | sf2b_10 | sfc_35 | During my day, I take a moment to explore what I can see, hear and touch | 8 | |
| 2 | sf2b_11 | sf6_44 | When I am busy, returning to the sensing mode feels like a low priority | exploratory | yes |
| 2 | sf2b_12 | sfc_24 | Looking closely at ordinary things around me often gives me a sense of wonder | 10 | |
| 2 | sf2b_13 | sfc_31 | I react to things before I have really taken them in | exploratory | yes |
