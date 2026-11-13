# Study 2b: decisions after Study 2a (for the registration)

Study 2a (osf.io/vt7gn) is analysed. The registered analysis is `study2a/analysis/01`–`11`; the exploratory follow-ups that led to these decisions are `12`–`21`. The results and this plan are on senseforaging.com/research/study2a/results (data from `22_report_data.R`).

`preregistration_draft.md` in this folder still describes the 20-item, four-factor short form. It is rewritten from this list before 2b is registered. (Since rewritten, and renamed `preregistration_2b.md`.)

## Decided (Norm, 31 October 2026)

1. **One score, not subscales.** 2a found one broad sensing factor plus small reverse-keyed clusters; Access, Valuing and Forage do not separate, and no hypothesized factor gives a 3-statement subscale beyond the general factor (`13`, `14`).
2. **The scale is 8 statements:** N4, sf6_36 (R), V10, sf6_27, sfc_36, sfc_x2 (R), sfc_35, U4 (R). *Revised 2 November (decision 14):* this was 12's best 8 with sf6_44 (R) in U4's place, all three reverse-keyed statements naming a mode. On 31 October Norm kept sfc_x2 rather than swapping in U4 for it: the mode-naming statements are the strongest reverse-keyed indicators (`17`, `19`, `20`).
3. **The 10** adds A9 and sfc_24 (`18`). The 8 sits inside it, so 2b scores both. The 10 is for studies of change, where its smaller detectable change helps.
4. **Exploratory in 2b, unscored:** N6 and sfc_31, reverse-keyed statements that do not name a mode, and sf6_44 (asked anyway as one of the trial's statements). *Revised 2 November:* U4 was the third exploratory statement until it moved into the 8 (decision 14).
5. **The confirmatory model is a one-factor CFA with an agreement factor:** a random intercept with loadings fixed at 1 on every unrecoded statement, a free variance, and no correlation with the content factor. WLSMV, theta parameterization, pairwise. Without it every version misfits (RMSEA .10–.13 in 2a), for a known reason.
6. **The vignettes stay in 2b**, before the scale as in 2a, with V2 ("Doing: keeping watch") rewritten: it fell on the sensing side of the midpoint. Norm's wording (2 November 2026), in `materials/sf2b-vignettes.json`:
   - **Introduction:** the doing mode is now "look or listen to accomplish something, like finding our keys, or checking out a body ache or a strange sound outside", and the sensing mode starts its own paragraph (three paragraphs, not two).
   - **V2:** "Walking home late, Riley listens for footsteps to make sure no one is following her." (2a: "…keeps listening for footsteps and watching the shadows on the empty street.") The goal is now stated.
   - **V7:** "Jamie feels a racing heart and keeps checking it to make sure his heart rate goes back to normal." (2a: "…keeps checking it, to make sure nothing is wrong.")
   - The other seven vignettes, the self-rating and the ordering rules are unchanged from 2a. Because V2, V7 and the introduction changed, separation scores are not directly comparable with 2a's.
7. **The understanding test:** measurement invariance of the CFA across people who do and do not tell the modes apart (separation ≥ 1, the 2a rule): configural, then loadings, then thresholds. In 2a (in-sample) loadings were not equal across the groups (`21`).
8. **IUS-12 is dropped from 2b.** It was there for Letting be, which is not in the scale; its correlation with the 8 in 2a (−.17) is reported.
9. **The validity battery is the focus of the freed time:** convergent and divergent measures for a single sensing score (`prereg_2b_plan.md`).

## Decided (Norm, 31 October 2026, second round)

10. **A behavioural criterion: a free-choice pause in the session.** Participants are offered an optional pause to look and listen around them, with Continue available at once; whether they take it and for how long is the criterion. It needs a platform step and goes in the ethics amendment.
11. **New discriminant measures: the Brief Sensation Seeking Scale (BSSS-8) and a short Highly Sensitive Person scale.** Licences to check.
12. **No longitudinal component in 2b.** The follow-up session goes, and with it H7 (test–retest stability), the over-time invariance model and the prospective past-week analysis.
13. **The trial's 10 statements in Study 1 wording are carried** (7 statements not otherwise asked), as a check that the 8 can stand in for the PSY240 trial's measure in independent data. With decision 14 the 8 shares two statements with the trial's 10 (sf6_36, sf6_27), not three.

## Decided (Norm, 2 November 2026)

14. **U4 replaces sf6_44 in the 8** (and so in the 10). sf6_44 loaded weakest (.37 in the one-factor CFA with the agreement factor), and its residual correlation with sfc_x2 (−.110) was the 8's largest misfit. In 2a, in-sample (`21`, rerun 2 November):

    | | CFI | RMSEA (upper CI) | SRMR | Largest residual | ω (content) | α | Practitioner d | r with trial's 10 |
    |---|---|---|---|---|---|---|---|---|
    | The 8 with sf6_44 | .972 | .084 (.103) | .038 | .110 | .845 | .80 | .70 | .84 |
    | The 8 with U4 | .985 | .061 (.081) | .031 | .079 | .838 | .80 | .67 | .80 |
    | The 10 with U4 | .987 | .048 (.063) | .030 | .077 | .850 | .83 | .69 | .81 |

    Also compared: U4 in place of sfc_x2 (CFI .988, RMSEA .053, largest residual .069, but two weak statements: U4 .32 and sf6_44 .39); and the 8 plus U4 as a ninth (no better than the 8: RMSEA .081, the sf6_44–sfc_x2 residual grows to −.129). The 8 still has three reverse-keyed statements, two naming a mode. U4 was written for Understanding ("Even when I pause to take in what is around me, I react to stress the same way") and loads .37. Chosen on 2a's data, so 2b is the test. The lower correlation with the trial's 10 is mostly the loss of a shared statement.
15. **The vignette introduction, V2 and V7 are reworded** (see decision 6, `materials/sf2b-vignettes.json`).
16. **The Meaning-Detecting subscale of the Synchronicity Awareness and Meaning-Detecting Scale (SAMD-MD; Russo-Netzer & Icekson, 2023) is added as a discriminant measure** (`materials/sf2-samd-md.json`, 13 items, 1–7, wording as published). It asks whether sense foraging is a spiritual or meaning-seeking openness: its core is reading meaning into events (for example, its item 5), where the scale is about taking in what is there. Predicted positive but small to moderate (a few statements share content: bodily sensations, awe and wonder), and tested as a contrast: the 8 correlates more with FFMQ Observing than with SAMD-MD. The Synchronicity Awareness subscale (frequency of coincidences) is not used. Caveats: developed in Hebrew with Israeli samples, and the English wording is a back-translation not validated in English speakers; two-factor CFA fit CFI .93, RMSEA .06, SRMR .05. The article is CC BY and prints the items; the authors are to be told.

## Decided (Norm, 9 November 2026, on the 2b decisions page)

17. **H1's fit rule.** The model passes when CFI ≥ .95; no residual correlation exceeds .10, or freeing the largest changes ω and the loadings by less than .02; and the unbiased SRMR close-fit test (SRMR ≤ .05) is not rejected. RMSEA and the model-specific cutoffs are reported, not tested. (The background and the simulation are under "Fit-rule background" below.)
18. **NCS-6 is added** as a discriminant measure (near zero). SBPS-8 and NR-6 are not.
19. **The Brief MAIA-2 (24 items; Mehling et al., 2018; Rogowska et al., 2023) replaces the full MAIA-2.** Noticing stays a confirmatory convergent (H2) and Noticing > Not-Worrying stays a contrast (H3b); H5a controls for the Brief MAIA-2.
20. **Past-week time outdoors stays descriptive.**
21. **One pause offer, not two.** Offer 1 is dropped: it came too early in the session to make sense. The single offer sits partway through, in the validity block at least 5 minutes after the scale, with the where-are-you and how-tired questions after it. The primary criterion becomes one engaged pause (0/1).
22. **A pause counts** when taken, at least 20 s on the page, with time away under half the total (`criterion_pause_design.md`).
23. **A pilot of about 60** runs first on the full session; it sets the time and pay and checks the export, and is not part of the confirmatory sample.
24. **N = 500 in one posting, as in 2a** (replacing 1,000 in two postings). 2a's single posting gave 273 current and 136 never-practised participants, and 154 who did not tell the modes apart.
25. **Pay is minimum wage for the pilot's median time,** with a Prolific maximum time long enough that a long pause doesn't time out.
26. **The HSP short form is dropped.** In its place: the 4 Sensory Discomfort items of the SPSQ-SF (De Gucht & Woestenburg, 2025, European Journal of Psychological Assessment, 41, 348–356; items in its supplement) for sensory overwhelm, and the 6 Physical Concerns items of the ASI-3 for anxious monitoring of the body. Both discriminant; bounds set in the preregistration. Copies in the Shared drive's `SenseForaging/Assessment/study2b`.
27. **Licences:** the measures are used for non-commercial academic research and are treated as permitted; PSS-4 and SBI stay.
28. **A new platform study, "Sense Foraging Study 2b".** The old Study 2 and its retest study are retired.
29. **The name stays Study 2b.**
30. **H8 is exploratory** (Norm, 9 November): reported, no pass/fail. At N = 500 about 150 fall in the group that doesn't tell the modes apart (154 in 2a), enough for loadings but thin for six-point thresholds.
31. **Pay: the UK minimum wage, £12.71 an hour** (Norm, 9 November), for the pilot's median time (decision 25).
32. **The smallest H3b contrast of interest is a difference of .15** (Norm, 9 November), to be backed by simulation before registering. The simulation (`power/README.md`) found .15 underpowered at N = 500 (54–76%); see "Still open".

33. **H1's rule is loosened after the simulation** (Norm, 9 November): the model passes when CFI ≥ .95; no residual correlation exceeds .10, or freeing the largest changes ω by less than .03 and every loading by less than .05; and the unbiased SRMR close-fit test (SRMR ≤ .05) is not rejected. RMSEA and the model-specific cutoffs are reported. Replaces decision 17's .02 / .02. Simulated at N = 500: passes 2a-like data 96%, fails a three-statement second factor loading .60 99% (`power/README.md`).
34. **N = 700, and the smallest H3b difference stays .15** (Norm, 9 November). Replaces decision 24's 500; one posting, as 2a. The simulation puts .15 at about 80% near N = 700.
35. **The debrief's account of the pause is covered by the approved amendment** (Norm, 9 November): it says the break was part of the study and was recorded, which was not said beforehand.
36. **H3b contrast 4 (Brief MAIA-2 Noticing > Not-Worrying) is kept, with its power reported** (Norm, 9 November): about 64% at N = 700, because it compares two 3-item subscales that barely correlate. The other four contrasts are 76–86%.
37. **The pilot of 60 runs before filing** (Norm, 10 November): its extraction goes through the code with answers shuffled within each item, and the registration is filed after that run, before the other 640 places open. No unshuffled pilot answers are analysed; the decision rules for the run are in `study2b/analysis/README.md`.
38. **Validity scores keep Study 2's 80% rule** (Norm, 10 November): a score needs at least 80% of its items, so a 3-item subscale needs all three. H5's power at the resulting n (410–470 in the simulated run) is in Section 3.
39. **H3a's and H3d's bound is .80, after Rönkkö & Cho (2022)** (Norm, 12 November): the upper limit of the 95% CI below .80 is their "no problem" class. Applying it to H3d's multiple R is our extension and a stricter test. Each upper limit is also reported in their four classes, not tested.
40. **The seven drafting choices are confirmed** (Norm, 12 November): H3b holds when at least four of five contrasts do; H3e's bounds (BSSS-8 below .40, the other three within ±.30); validity measures as latent factors of their items; delta-method CIs for H3b, the bootstrap for H3d; H4 two-sided with Holm, H6 on 2a's current-practitioner group with dedicated as sensitivity; the pilot's 60 included if nothing changes; the flags (longstring 10+ of 14, more than 2 of 13 unanswered).
41. **The code's 13 readings of the text are confirmed** (Norm, 13 November): each is stated in the preregistration's Section 5, "Specifications fixed in the code". The text is then complete except the entries marked [AT FILING], which the pilot run supplies.

Correction to decision 13 (9 November): the trial's 10 needs **eight** extra statements in 2b, not seven. 2b asks only sf6_36 and sf6_44 in the trial's wording; the 8's sf6_27 is the reworded version. And 2a's `sf2a-end` had two of its six statements (sf6_7, sf6_27) in wordings that are not the trial's (`psy240-sf-short-20`), so 2a's "trial's 10" is not exactly the trial's measure (r with the 8 = .80, decision 14). 2b's `sf2b-end` uses the trial's wording.

The REB amendment is approved under the same protocol (#00051180).

## Still open

- Nothing in the text. At filing: the [AT FILING] entries, from the pilot run.

## Fit-rule background

- **H1's fit rule.** The previous draft required CFI ≥ .95, RMSEA ≤ .06 (upper CI ≤ .08), SRMR ≤ .08. In 2a, in-sample, the 8 gives CFI .972, SRMR .038, RMSEA .084 (upper CI .103), so it would fail RMSEA; the 10 gives .977, .036, .066 (.080); the 8 with U4 in place of sfc_x2 gives .988, .030, .053.
  - **Simulation (31 October, from the fitted 2a models, N = 1,000, 200 runs):** the RMSEA of .084 is not a small-df artefact. A correctly specified 8 gives RMSEA ≈ .005 (95th percentile .026). The conventional cutoffs never detected a planted flaw (residual correlations of .30 on two item pairs; CFI ≈ .986, SRMR ≈ .030). Model-specific cutoffs that catch it 95% of the time are about CFI ≥ .992, RMSEA ≤ .042, SRMR ≤ .024. Both 8s miss those in 2a.
  - **The misfit sat on the sfc_x2–sf6_44 pair** (−.110), with sfc_x2's other residuals next (sf6_27 −.089, sfc_35 .075, sf6_36 .069). Removing either statement removes most of it, which is why decision 14 swaps out sf6_44. (Corrected 2 November: this said the misfit sits on sfc_x2.)
  - **Rule adopted (decision 17):** the local-fit rule and the SRMR close-fit test as proposed, plus CFI ≥ .95 (Norm, 9 November).
- Carried over: face-scale "Prefer not to answer" (now a study setting, `ratings_allow_pna`); pages marked `translate="no"` (live); 2a participants excluded on Prolific.
