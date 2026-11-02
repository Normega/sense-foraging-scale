# Study 2b: decisions after Study 2a (for the registration)

Study 2a (osf.io/vt7gn) is analysed. The registered analysis is `study2a/analysis/01`–`11`; the exploratory follow-ups that led to these decisions are `12`–`21`. The results and this plan are on senseforaging.com/research/study2a/results (data from `22_report_data.R`).

`preregistration_draft.md` in this folder still describes the 20-item, four-factor short form. It is rewritten from this list before 2b is registered.

## Decided (Norm, 31 October 2026)

1. **One score, not subscales.** 2a found one broad sensing factor plus small reverse-keyed clusters; Access, Valuing and Forage do not separate, and no hypothesized factor gives a 3-statement subscale beyond the general factor (`13`, `14`).
2. **The scale is 8 statements** (12's best 8): N4, sf6_36 (R), V10, sf6_27, sfc_36, sfc_x2 (R), sfc_35, sf6_44 (R). The three reverse-keyed statements all name a mode. Norm kept sfc_x2 rather than swapping in U4: the mode-naming statements are the strongest reverse-keyed indicators (`17`, `19`, `20`).
3. **The 10** adds A9 and sfc_24 (`18`). The 8 sits inside it, so 2b scores both. The 10 is for studies of change, where its smaller detectable change helps.
4. **Exploratory in 2b:** three plain-worded reverse-keyed statements, U4, N6 and sfc_31, unscored, to see whether a mode-free reverse-keyed statement can replace a mode-naming one.
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
13. **The trial's 10 statements in Study 1 wording are carried** (7 statements beyond the 8), as a check that the 8 can stand in for the PSY240 trial's measure in independent data.

## To decide before registering

- **H1's fit rule.** The previous draft required CFI ≥ .95, RMSEA ≤ .06 (upper CI ≤ .08), SRMR ≤ .08. In 2a, in-sample, the 8 gives CFI .972, SRMR .038, RMSEA .084 (upper CI .103), so it would fail RMSEA; the 10 gives .977, .036, .066 (.080); the 8 with U4 in place of sfc_x2 gives .988, .030, .053.
  - **Simulation (31 October, from the fitted 2a models, N = 1,000, 200 runs):** the RMSEA of .084 is not a small-df artefact. A correctly specified 8 gives RMSEA ≈ .005 (95th percentile .026). The conventional cutoffs never detected a planted flaw (residual correlations of .30 on two item pairs; CFI ≈ .986, SRMR ≈ .030). Model-specific cutoffs that catch it 95% of the time are about CFI ≥ .992, RMSEA ≤ .042, SRMR ≤ .024. Both 8s miss those in 2a.
  - **The misfit sits on sfc_x2:** the four largest residual correlations all involve it (with sf6_44 −.110, sf6_27 −.089, sfc_35 .075, sf6_36 .069). With U4 in its place the largest is .069.
- Whether NCS-6, SBPS-8 or NR-6 are added as well; whether MAIA-2 is given in full.
- Carried over: face-scale "Prefer not to answer"; SBI permission; PSS-4 licence; pages marked `translate="no"`; 2a participants excluded on Prolific.
