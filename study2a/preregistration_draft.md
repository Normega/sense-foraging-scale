# Sense Foraging Scale, Study 2a: factor structure of a revised 58-statement set

OSF Preregistration (standard template), filed 29 October 2026 at [osf.io/vt7gn](https://osf.io/vt7gn), before data collection. University of Toronto REB Protocol #00051180, Study 2a amendment pending.

**This is the draft it was filed from; the OSF registration is the version of record.** Where the two differ, OSF governs.

**Settled after filing (29 October 2026, before any data were collected).** The filed version left two items marked "to settle". They are settled here as follows, and were filed as an update to the OSF registration the same day (submitted 20:54 UTC): (1) the vignette threshold is 75% (three in four), the value expectation 5 already states; (2) the analysis code will be written during the analysis, following these rules and Study 1's pipeline, and its commit will be cited when the results are reported. The materials it describes are built on the RADlab platform, inactive (study `33d012fb-e5c3-5e46-a147-e565c04666e6`; definitions in radlab `scripts/sense_foraging/study2a/`).

## 1. Study information

### Title

Development of the Sense Foraging Scale, Study 2a: factor structure of a revised 58-statement set, with a check that the two modes were understood

### Authors

Norman Farb (principal investigator), University of Toronto. Co-investigators as listed on REB Protocol #00051180.

### Description

The Sense Foraging Scale is a self-report measure of how people move between two ways of paying attention: a sensing mode (looking, listening, smelling, tasting or feeling out of curiosity, taking in whatever is there) and a doing mode (goal-directed attention that uses what one already knows, such as searching for something or checking that all is well). It is intended as an outcome measure for brief training in moving into the sensing mode, including a 4-week course evaluated in a randomized trial.

Study 1 (registered at osf.io/48r2y) found four content factors in a 48-statement pool, provisionally named Access, Understanding the modes, Valuing and savouring, and Letting be, plus a general tendency to agree. It also showed problems that Study 2a is designed to address:

- Understanding had low reliability (α .54), and four of its five short-form statements drew 87–94% agreement, so they barely separate people.
- Letting be (α .59) was defined mostly by reverse-keyed statements, so it may partly be a wording factor.
- A proposed fifth factor, Forage (habits of curious noticing in ordinary moments), had no statements of its own.

Study 2a gives a new sample a revised set of 58 statements: 7 in their Study 1 wording, 11 reworded short-form statements, 24 earlier candidates and 16 new statements, written for five hypothesized factors. Before the statements, participants read a short description of the two modes and rate nine brief vignettes on them, which checks that the distinction came across.

The aim is to choose, by the rules registered here, the statements for each factor and the structure that Study 2b will test. Study 2b is a separate validity study of 1,000 new participants; its confirmatory analysis will be registered after Study 2a and will cite these rules.

### Hypotheses

Study 2a is exploratory. We state expectations, not confirmatory tests. The decisions are made by the rules in Section 5.

1. **Structure.** The statements will form four or five correlated content factors resembling Access, Understanding, Valuing, Letting be and possibly Forage. The number is decided by the procedure in Section 5, not assumed.
2. **Agreement.** A random-intercept factor for the general tendency to agree will have non-zero variance.
3. **Forage.** No directional expectation. In Study 1, statements about valuing discovery loaded with Valuing; the Forage statements ask about habits instead, and may form a factor of their own or join Access or Valuing.
4. **Letting be.** No directional expectation. The set adds 7 forward statements to test whether the factor survives once reversed wording is modelled.
5. **The two modes.** Participants will rate the sensing vignettes higher than the doing vignettes, and at least three in four will separate the two by 1 point or more.
6. **Intolerance of uncertainty.** IUS-12 will correlate more strongly with Letting be's reverse-keyed statements than with its forward statements.

## 2. Design plan

### Study type

Observational study: a single-session online survey.

### Blinding

No blinding is involved. Participants are not told the hypotheses or which statements are reverse-keyed; the debrief explains both.

### Study design

Cross-sectional, one session of about 15 minutes, in this order:

1. **Consent.**
2. **Background**, as in Study 1 (13 questions): age, gender, race/ethnicity, country, education, contemplative-practice history (four follow-up questions shown only to those reporting regular practice), religiosity, spirituality and tradition.
3. **The two modes** (questionnaire `sf2a-vignettes`):
   - a two-paragraph introduction naming the doing mode and the sensing mode (Appendix A);
   - nine vignettes, each rated on "How much is this person in the doing or sensing mode?" from 1 Doing mode to 6 Sensing mode (2–5 unlabelled), with no feedback. Four are sensing vignettes and five are doing vignettes (Appendix A);
   - a self-rating on the same scale: "The last time you had a few minutes with nothing to do, how much were you in the doing or sensing mode?"
4. **The statements** (`sf2a-scale`): the 58 statements and three instructed-response attention checks (Appendix B). Instruction: "In general, how much do you agree with each statement? There are no right or wrong answers." Six points, every point labelled: 1 Strongly disagree, 2 Disagree, 3 Mildly disagree, 4 Mildly agree, 5 Agree, 6 Strongly agree, with no midpoint. Every statement except the checks also offers "Prefer not to answer". Pages of 8/8/8/8/8/7/7/7.
5. **IUS-12** (Carleton, Norton & Asmundson, 2007), 12 items on its 5-point scale. Besides its use in Section 5, it puts time between the reworded statements and step 6.
6. **The PSY240 trial's six other primary-score statements in their Study 1 wording** (`sf2a-end`; sf6_7, sf6_8, sf6_12, sf6_18, sf6_27, sf6_11), fixed order, same six points. With sf6_9, sf6_36, sf6_35 and sf6_44, which are in the scale in Study 1 wording, they let Study 2a score the trial's 10-statement Access + Valuing form exactly as the trial does.
7. **Debrief.**

Vignettes and statements also offer "Prefer not to answer". The platform shows several questions per page as labelled cards (layout v2), on phones and computers alike.

Recruitment is one Prolific posting with a general-interest description, as Study 1's general posting. Nobody is screened on practice.

### Randomization

Vignette and statement orders are drawn per participant: a seeded, constrained random order, seeded with the participant's session id, so it is stable on reload. The rules:

- **Vignettes:** three a page, each page holding at least one sensing and one doing vignette. The self-rating comes last, on its own page.
- **Statements:**
  - statements of the same hypothesized factor at least 2 positions apart and at most 3 a page;
  - reverse-keyed statements 1 to 4 a page, never first on a page;
  - the 9 statements that name a mode at least 3 positions apart and at most 2 a page;
  - 18 look-alike pairs on different pages (Appendix B);
  - one attention check in each third of the 61 positions (positions 3–19, 22–39 and 42–59), never first or last on a page.

Each question's presented position is recorded, so order effects can be examined. The rules found a valid order for 2,000 of 2,000 test seeds.

## 3. Sampling plan

### Existing data

**Registration prior to creation of data.** No Study 2a data exist at registration.

### Data collection procedures

- **Participants:** Prolific participants aged 18+, fluent in English (Prolific's fluent-languages filter), in any country. Everyone who took part in Study 1, on either posting, is excluded through Prolific's exclusion list. Study 2a's participants will in turn be excluded from Study 2b.
- **Platform:** the RADlab research platform (radlab.zone). Prolific IDs are held in a separate table, replaced by a random code in all study data, and deleted after payment.
- **Compensation:** USD 12 per hour, prorated (about USD 3.00 for 15 minutes), paid regardless of answers, including the attention checks. If the pilot's median time is longer, the posting's time and pay are raised before the main sample.
- **Pilot:** the first 20 places, with the identical materials. The pilot is inspected for quality only: completion counts, attention-check pass rates, completion time, "Prefer not to answer" counts, and per-person spread of answers. No correlations, factor solutions or vignette summaries are computed on it. If nothing is changed after the pilot, pilot participants are included. If the materials change, they are excluded and 20 places are added.

### Sample size

Target: 500 completed sessions (20 pilot + 480). Expected analysable N is about 470–480 after incomplete sessions and exclusions, as in Study 1.

### Sample size rationale

- **Participants per statement.** About 8.6 for 58 statements.
- **Overdetermination.** The five hypothesized factors have 10–14 statements each. With strongly overdetermined factors, several hundred participants recover the structure even at the moderate-to-low communalities typical of new self-report items (MacCallum, Widaman, Zhang & Hong, 1999; Mundfrom, Shaw & Ke, 2005).
- **The models.** Ordinal estimation with a random-intercept factor needs more data than a plain EFA. Simulations of the random-intercept approach covered N = 200–1,000 (de la Fuente & Abad, 2020). Study 1 recovered a four-factor solution at the same N with 48 statements.
- **Practical constraint.** About 15 minutes per participant at USD 12 per hour.

### Stopping rule

The number of places is fixed on Prolific in advance: 20 (pilot), then 480. Data collection stops when these places are filled. No interim analysis informs stopping, and excluded participants are not replaced.

## 4. Variables

### Manipulated variables

None.

### Measured variables

- **Statements:** `sf2a_01` to `sf2a_58` (1–6, or "Prefer not to answer"). The 19 reverse-keyed statements are stored unrecoded. Appendix B gives each one's review id, hypothesized factor and keying. Seven carry their Study 1 id as `study1_item` (sf6_9, sf6_36, sf6_35, sf6_31, sf6_44, sf6_43, sf6_21).
- **Attention checks:** `sf2a_attn_disagree` (correct = 2), `sf2a_attn_agree` (correct = 5), `sf2a_attn_strongly_disagree` (correct = 1).
- **Vignettes:** `vig_1` to `vig_9` (= V1 to V9) and `vig_self` (1–6, or "Prefer not to answer").
- **IUS-12:** `ius12_01` to `ius12_12` (1–5).
- **Study 1 wording:** `sf6_7`, `sf6_8`, `sf6_12`, `sf6_18`, `sf6_27`, `sf6_11` (1–6, or "Prefer not to answer").
- **Background questions:** as in Study 1.
- **Paradata:** each question's presented position; seconds per page and in total.

### Indices

- **Factor scores and reliability:** factor scores from the final model; ordinal omega per factor.
- **Separation (vignettes):** a participant's mean rating of the four sensing vignettes (V3, V4, V6, V8) minus their mean rating of the five doing vignettes (V1, V2, V5, V7, V9). It is computed for participants who rated at least three of each kind. A participant **tells the modes apart** when their separation is 1 point or more.
- **Sample midpoint (vignettes):** halfway between the sample's mean sensing rating and its mean doing rating.
- **Agreement rate:** for each statement, the share of answering participants choosing 4–6, on the unrecoded answer.
- **The trial's Access + Valuing score:** the mean of its 10 statements in Study 1 wording (sf6_9, sf6_36 and sf6_35 reversed, sf6_7, sf6_8, sf6_12, sf6_18, sf6_27, sf6_11, sf6_44 reversed), scored as the trial scores it.
- **IUS-12 total:** the sum of the 12 items (12–60).
- **Practice group:**
  - "current practitioner" = regular practice "Yes" and current frequency 1–2 times a week or more;
  - "never" = regular practice "No";
  - others are reported separately.
- **Longstring:** the longest run of identical answers in presented order, on the 61 scale positions.

## 5. Analysis plan

### Statistical models

R, using psych, lavaan, GPArotation and EFAtools; versions are reported. Statement answers are treated as ordinal (polychoric correlations, WLSMV estimation). The analysis code follows Study 1's pipeline. It will be written during the analysis, following the rules registered here, and its commit will be cited when the results are reported.

1. **Pre-checks.**
   - Distributions, floor and ceiling, "Prefer not to answer" rates, and agreement rates for every statement.
   - Polychoric correlations within the 18 look-alike pairs, screened before the agreement factor is fitted, since a random intercept can absorb doublets.
2. **The two modes** (decision V).
   - Mean sensing and doing ratings for the sample, and for each vignette.
   - The share of participants who tell the modes apart.
   - Each vignette's mean against the sample midpoint.
3. **Number of factors** (decision K). Parallel analysis on the polychoric correlations of the 58 statements, reverse-keyed statements recoded (principal-axis factors, 500 random datasets, 95th percentile). Ordinal EFAs with 4, 5 and 6 factors (WLSMV, geomin), and with the parallel-analysis count if it falls outside that range.
4. **Primary model.** An exploratory structural equation model (Asparouhov & Muthén, 2009) on the 58 unrecoded statements, with:
   - the chosen number of correlated content factors, partially specified target rotation: for five factors, toward the hypothesized membership (non-member loadings targeted to 0, members free); for any other number, toward the membership in the geomin EFA of step 3;
   - a random-intercept factor for agreement: loadings fixed to 1 on all 58 statements, variance estimated, uncorrelated with the content factors (Maydeu-Olivares & Coffman, 2006; Aichholzer, 2014).
5. **Letting be against wording** (decision L). The primary model plus a wording factor: loadings free on all 19 reverse-keyed statements, uncorrelated with the content and agreement factors. If that model does not converge or gives an improper solution, the wording factor replaces the random intercept and the reason is reported.
6. **Intolerance of uncertainty.** Correlations of IUS-12 with the mean of Letting be's forward statements and with the mean of its reverse-keyed statements (recoded), compared with a test for dependent correlations (Steiger, 1980).
7. **Reliability.** Ordinal omega for each retained factor.
8. **Statement selection** (decision S).

### Decision rules

- **V. Did the two modes come across?**
  - A participant tells the modes apart when their separation is 1 point or more. Because separation compares a participant's ratings with each other, it does not depend on where they centre the scale.
  - The distinction came across if at least 75% of participants with a separation score tell the modes apart (three in four).
  - A vignette whose mean falls on the wrong side of the sample midpoint is flagged and rewritten before any further use.
  - Nobody is excluded on their vignette ratings.
  - If the distinction did not come across, a statement that names a mode is retained in step 8 only where no statement without one qualifies for the same factor.
- **K. Number of factors.** By these criteria, in order:
  1. every factor has at least 3 statements with primary loadings ≥ .40;
  2. interpretability against the hypothesized factors;
  3. fit, reported descriptively (CFI, RMSEA, SRMR).

  The choice and its reasons are documented.
- **F. Forage.** Forage stands as a factor if the chosen solution has a factor on which at least three Forage statements have primary loadings ≥ .40 and more than half of the statements with primary loadings ≥ .40 are Forage statements. Otherwise the Forage statements are reported with the factor they join.
- **L. Letting be.** Letting be stands if, in the model of step 5, at least three of its 7 forward statements load ≥ .40 on the Letting be factor. Otherwise it is reported as a wording factor, and Study 2b keeps it only as an exploratory score.
- **S. Statement selection.** For each retained factor, five or six statements, at least two of them reverse-keyed where at least two qualify. Prefer statements with:
  1. primary loading ≥ .50 (minimum .40);
  2. no secondary loading ≥ .30, and primary minus secondary ≥ .20;
  3. "Prefer not to answer" ≤ 5%;
  4. an agreement rate between 15% and 85% (a statement outside it is flagged, and an unflagged statement is preferred to a flagged one);
  5. at most one statement of each look-alike pair.

  A short-form statement stays unless it fails these rules and another statement for the same factor does better. Decision V applies to statements that name a mode. The final selection is made by the PI applying these rules; every departure is reported with its reason.

### Transformations

Reverse-keyed statements are recoded (7 − x) for parallel analysis, reliability, selection summaries and scoring, and left unrecoded in models with the random-intercept factor. IUS-12 is scored as published. No other transformations.

### Inference criteria

The factor analyses are exploratory: decisions follow the rules above, and fit indices are descriptive. Correlations are reported with 95% confidence intervals. Practitioner comparisons are descriptive.

### Data exclusion

- **Excluded:**
  - staff test enrollments;
  - sessions without a submitted scale;
  - participants who fail two or more of the three attention checks.
- **Repeat submissions:** the participant's first complete submission of each questionnaire is used.
- **Flagged and examined in sensitivity analyses, not excluded:**
  - failing exactly one attention check;
  - fewer than 2 seconds per statement on the scale;
  - a longstring of 15 or more;
  - "Prefer not to answer" on more than 11 of the 58 statements.

### Missing data

"Prefer not to answer" and unanswered questions are treated as missing. Polychoric correlations and WLSMV use pairwise-available data (lavaan `missing = "pairwise"`), so a participant with a skipped statement is kept; the number of participants with any missing statement is reported. Statements with more than 5% missing are reported. Factor scores use available statements.

### Exploratory analysis

**Sensitivity analyses**, each compared with the primary solution:

1. the primary model in participants who tell the modes apart;
2. the primary model without flagged participants;
3. the primary model without the random-intercept factor;
4. an EFA of the forward statements only.

**Further exploratory analyses:**

- the self-rating's correlation with each factor score;
- each reworded trial statement's correlation with its Study 1 wording, and the trial's Access + Valuing score's correlations with the Access and Valuing factor scores;
- differences between current practitioners and never-practitioners for each statement and factor, standardized, as a guide to Study 2b, not as a selection rule;
- the effect of presented position on statement means and loadings.

## 6. Other

**Next steps.** The statements chosen here go to Study 2b (1,000 new participants with established measures of mindfulness, interoceptive awareness, curiosity and wellbeing, and a two-week retest), whose confirmatory factor model and validity hypotheses will be registered before its data are collected. In Study 2b, Forage is treated as a companion measure of practice, not as a predictor of how much people practise.

**Data sharing.** A de-identified dataset will be shared as described in the consent form: no exact timestamps or durations, age in bands, coarsened rare demographic combinations, free text recoded. The analysis code will be shared alongside it.

### References

- Aichholzer, J. (2014). Random intercept EFA of personality scales. *Journal of Research in Personality, 53*, 1–4.
- Asparouhov, T., & Muthén, B. (2009). Exploratory structural equation modeling. *Structural Equation Modeling, 16*(3), 397–438.
- Carleton, R. N., Norton, M. A. P. J., & Asmundson, G. J. G. (2007). Fearing the unknown: A short version of the Intolerance of Uncertainty Scale. *Journal of Anxiety Disorders, 21*(1), 105–117.
- de la Fuente, J., & Abad, F. J. (2020). Comparing methods for modeling acquiescence in multidimensional partially balanced scales. *Psicothema, 32*(4).
- MacCallum, R. C., Widaman, K. F., Zhang, S., & Hong, S. (1999). Sample size in factor analysis. *Psychological Methods, 4*(1), 84–99.
- Maydeu-Olivares, A., & Coffman, D. L. (2006). Random intercept item factor analysis. *Psychological Methods, 11*(4), 344–362.
- Mundfrom, D. J., Shaw, D. G., & Ke, T. L. (2005). Minimum sample size recommendations for conducting factor analyses. *International Journal of Testing, 5*(2), 159–168.
- Steiger, J. H. (1980). Tests for comparing elements of a correlation matrix. *Psychological Bulletin, 87*(2), 245–251.

**Appendices.** Appendices A and B follow. They are also in `study2a/appendices/` as separate files, with the statements as a CSV.

<!-- APPENDICES -->

## Appendix A. The two modes

**Introduction** (shown to everyone before the vignettes):

> People notice things in different ways. Our brains can toggle between two modes: a doing mode, which is goal-directed and uses knowledge we already have, and a sensing mode, which explores our sensory world to discover new things.
>
> In the doing mode, we look or listen to find something, like our keys, or to check that all is well, like a body ache or a strange sound outside. In the sensing mode, we look, listen, smell, taste or feel just out of curiosity, taking in whatever is there. Everyone does both at times, and neither is better. There are no right or wrong answers.

**Vignettes.** "How much is this person in the doing or sensing mode?" 1 Doing mode to 6 Sensing mode.

| Id | Kind | Vignette |
|---|---|---|
| V1 (`vig_1`) | Doing: finding something | Priya walks along the shelves in the drug store, reading labels until finding the brand on her list. |
| V2 (`vig_2`) | Doing: keeping watch | Walking home late, Riley keeps listening for footsteps and watching the shadows on the empty street. |
| V3 (`vig_3`) | Sensing: around me | Waiting for a bus, Sam watches the light move on the wet road and listens to the cars go by. |
| V4 (`vig_4`) | Sensing: letting go of a goal | Lee is photographing the light on a building, then stops trying for the perfect shot and just takes in the view. |
| V5 (`vig_5`) | Doing: planning | In line at the ATM, Alex barely sees the street, busy working out a budget for next summer’s holiday. |
| V6 (`vig_6`) | Sensing: in the body | Walking down the hallway, Maya notices how each breath feels, then the warmth spreading through both hands. |
| V7 (`vig_7`) | Doing: checking the body | Jamie feels a racing heart and keeps checking it, to make sure nothing is wrong. |
| V8 (`vig_8`) | Sensing: taste | Drinking a morning coffee, Kai notices its bitterness and warmth, and that it tastes different from usual. |
| V9 (`vig_9`) | Doing: finding something | In a noisy room, Noor strains to catch a phone number in a voicemail. |

**Self-rating** (`vig_self`): "The last time you had a few minutes with nothing to do, how much were you in the doing or sensing mode?"

## Appendix B. The statements

Ids in data order. R = reverse-keyed. "Study 1" = the Study 1 id of a statement still in its Study 1 wording. "Apart from" = its look-alike partners, never on the same page.

### Access

*Definition (not shown to participants):* Knowing that curious noticing is open to you, around you and in your body, and being able to turn to it.

| Id | Review id | Keying | Study 1 | Statement | Apart from |
|---|---|---|---|---|---|
| `sf2a_01` | sf6_9 |  | sf6_9 | Moving into the sensing mode comes easily to me | sf6_36, sf6_35 |
| `sf2a_02` | sf6_36 | R | sf6_36 | Moving into the sensing mode is a struggle for me | sf6_9, sf6_35 |
| `sf2a_03` | sf6_35 | R | sf6_35 | I question my ability to move into the sensing mode | sf6_36, sf6_9 |
| `sf2a_04` | sf6_7 |  |  | When I notice I am caught up in a task, I can switch to the sights and sounds around me | N2 |
| `sf2a_05` | sf6_8 |  |  | In a new place, I can take in what is there, not just look for what I need to do next | A9 |
| `sf2a_06` | sfc_21 |  |  | Even in a stressful moment, I can be curious about what I see and hear around me |  |
| `sf2a_07` | sfc_22 |  |  | Even in dull places, like a waiting room, I can find something interesting to look at or listen to |  |
| `sf2a_08` | sfc_23 |  |  | I can tell which moments in my day leave room for exploring what is around me |  |
| `sf2a_09` | N2 | R |  | I often slip into the doing mode without meaning to | sf6_7 |
| `sf2a_10` | A9 |  |  | I can pay attention to what I’m feeling in my body in new or unfamiliar places | sf6_8 |

### Understanding

*Definition (not shown to participants):* Learning by sensing: what noticing shows you about yourself and a situation, including how stress pulls you into rushing and fixing.

| Id | Review id | Keying | Study 1 | Statement | Apart from |
|---|---|---|---|---|---|
| `sf2a_11` | sf6_3 |  |  | For me, exploring what I see and hear is a different state of mind from planning and fixing |  |
| `sf2a_12` | sf6_15 |  |  | When I am stressed, I notice myself falling into my usual way of rushing to fix things | U1, U4 |
| `sf2a_13` | sf6_37 | R |  | The sensing and doing modes seem pretty similar to me |  |
| `sf2a_14` | sfc_8 |  |  | When I try to use the sensing mode, I notice how easily I shift back into thinking about what I have to do next | sfc_9 |
| `sf2a_15` | sfc_9 | R |  | I don’t notice that I have been rushing to get things done until much later | sfc_8, U3 |
| `sf2a_16` | sfc_2 | R |  | My attention works the same way whether I am stressed or relaxed |  |
| `sf2a_17` | U1 |  |  | When I am stressed, a moment of looking around out of curiosity makes me less likely to react as usual | sf6_15 |
| `sf2a_18` | U3 |  |  | When I pause to notice my body, I can feel how much I had been rushing | sfc_9 |
| `sf2a_19` | U4 | R |  | Even when I pause to take in what is around me, I react to stress the same way | sf6_15 |
| `sf2a_20` | U6 |  |  | When I feel unsafe, I lose my curiosity about what is around me |  |
| `sf2a_21` | U7 |  |  | Noticing what I can see, hear or feel often shows me something about a situation I had not thought of |  |

### Valuing

*Definition (not shown to participants):* Treating time spent noticing out of curiosity as worthwhile, and enjoying it.

| Id | Review id | Keying | Study 1 | Statement | Apart from |
|---|---|---|---|---|---|
| `sf2a_22` | sf6_31 |  | sf6_31 | Seeing something from a new or unexpected perspective feels valuable to me in its own right |  |
| `sf2a_23` | sf6_12 |  |  | Even when I am busy, it feels worth pausing to notice the sights and sounds around me | sf6_44 |
| `sf2a_24` | sf6_18 |  |  | When stress narrows my focus, looking around out of curiosity matters even more to me |  |
| `sf2a_25` | sf6_27 |  |  | I enjoy exploring the sights, sounds and textures of ordinary moments |  |
| `sf2a_26` | sf6_11 |  |  | I value time spent just looking and listening, even when there is work to do | sf6_43 |
| `sf2a_27` | sf6_44 | R | sf6_44 | When I am busy, returning to the sensing mode feels like a low priority | sf6_12 |
| `sf2a_28` | sfc_24 |  |  | Looking closely at ordinary things around me often gives me a sense of wonder |  |
| `sf2a_29` | sfc_25 |  |  | When I look around out of curiosity, I enjoy noticing things I did not expect |  |
| `sf2a_30` | sfc_34 | R |  | When I try to spend time in the sensing mode, I feel I should be doing something more useful |  |
| `sf2a_31` | V10 |  |  | I enjoy taking a moment just to notice how my body feels, like my breathing, the soles of my feet, or posture |  |
| `sf2a_32` | V11 |  |  | I enjoy noticing smells in ordinary places, like rain or coffee |  |
| `sf2a_33` | sf6_43 | R | sf6_43 | Only time spent getting things done feels worthwhile to me | sf6_11 |

### Forage

*Definition (not shown to participants):* Habits of curious noticing in ordinary moments, around you and in your body. No amounts or frequencies.

| Id | Review id | Keying | Study 1 | Statement | Apart from |
|---|---|---|---|---|---|
| `sf2a_34` | sfc_18 |  |  | Even on rushed days, I find moments to notice the sky or the sounds around me | N4 |
| `sf2a_35` | sfc_17 |  |  | In familiar places, I get curious about details I have not noticed before |  |
| `sf2a_36` | sfc_36 |  |  | Even in chores I have done many times, I notice sounds, smells or textures I had missed |  |
| `sf2a_37` | sfc_x2 | R |  | It never seems to be the right time to go into the sensing mode |  |
| `sf2a_38` | N7 |  |  | When I am waiting somewhere, I look around out of curiosity instead of reaching for my phone |  |
| `sf2a_39` | sfc_35 |  |  | During my day, I take a moment to explore what I can see, hear and touch | F10 |
| `sf2a_40` | N4 |  |  | As I go about my day, I notice the light, colours and sounds around me with interest | sfc_18 |
| `sf2a_41` | N5 |  |  | I have small habits of noticing what is around me, like watching the clouds on my way somewhere |  |
| `sf2a_42` | N6 | R |  | When I walk somewhere familiar, I notice very little of what I pass |  |
| `sf2a_43` | F10 |  |  | During my day, I notice how my breathing or body feels, just out of curiosity | sfc_35 |
| `sf2a_44` | F11 |  |  | When I eat something familiar, I notice whether I still like it as much as I expect |  |

### Letting be

*Definition (not shown to participants):* Staying with what is noticed before explaining, resolving or acting on it. The pause can be brief: letting something be for a moment, then acting when you choose, still counts. Its opposite is the pull to add meaning, predict what comes next, or respond straight away.

| Id | Review id | Keying | Study 1 | Statement | Apart from |
|---|---|---|---|---|---|
| `sf2a_45` | sfc_32 |  |  | I can pause when something starts to bother me, before I do anything about it |  |
| `sf2a_46` | sf6_40 | R |  | I need to understand why I feel fleeting sensations in my body |  |
| `sf2a_47` | sf6_47 | R |  | I don’t like going to events unless I know the full itinerary |  |
| `sf2a_48` | sf6_21 |  | sf6_21 | I can let a sensation simply be, without figuring out what it means or where it leads | sfc_14 |
| `sf2a_49` | sfc_10 |  |  | When something feels unpleasant, I first notice it without trying to change it |  |
| `sf2a_50` | sfc_12 | R |  | When I notice something, I start explaining it to myself straight away | sfc_14 |
| `sf2a_51` | sfc_14 | R |  | When I notice something, I have trouble just letting it be | sfc_12, sf6_21 |
| `sf2a_52` | sfc_26 |  |  | When I notice a sensation, I take it in for a moment before deciding what to do about it | sfc_13 |
| `sf2a_53` | sfc_28 |  |  | When something unexpected happens, I pause and take in what is happening before I react | sfc_31 |
| `sf2a_54` | sfc_29 |  |  | When I notice a problem, I let it be for a moment before trying to fix it | sfc_30 |
| `sf2a_55` | sfc_30 | R |  | As soon as I notice something, I feel I have to do something about it | sfc_29 |
| `sf2a_56` | sfc_31 | R |  | I react to things before I have really taken them in | sfc_28 |
| `sf2a_57` | L15 |  |  | When a sensation in my body is hard to name, I stay curious about it for a while |  |
| `sf2a_58` | sfc_13 | R |  | When I notice a sensation, my mind quickly moves on to thoughts about it | sfc_26 |

### Attention checks

| Id | Statement | Correct |
|---|---|---|
| `sf2a_attn_disagree` | To show that you are reading each statement, please select “Disagree” for this one. | 2 |
| `sf2a_attn_agree` | To show that you are reading each statement, please select “Agree” for this one. | 5 |
| `sf2a_attn_strongly_disagree` | To show that you are reading each statement, please select “Strongly disagree” for this one. | 1 |

### At the end, in Study 1 wording

| Id | Statement |
|---|---|
| `sf6_7` | Once I notice I have slipped into the doing mode without meaning to, I can shift back to sensing |
| `sf6_8` | I can move into the sensing mode in new or unfamiliar places |
| `sf6_12` | Returning to the sensing mode feels worth prioritizing, even when I am busy |
| `sf6_18` | When stress narrows my focus, returning to the sensing mode matters even more to me |
| `sf6_27` | I enjoy exploring the world through my senses |
| `sf6_11` | I value time spent in the sensing mode, rather than always focusing on getting things done |
