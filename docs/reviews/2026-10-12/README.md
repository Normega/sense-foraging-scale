# Independent review, 12 October 2026

`sf_review.md` is an outside-style psychometric review of Study 1 and the Study 2 plan, written by a
separate Claude agent briefed to check claims against the code and outputs. The `r1`–`r7` scripts are
its re-analyses: they read Study 1's derived data and write nothing.

Checked before acting:
- Model A on Study 1 gives CFI .850 and RMSEA .086.
- `psych::omega`'s default flipping inflated ω_h.
- The Study 2 draft's screening line contradicted the REB application.
- The 20-item ESEM fits in-sample (CFI .969).

What was done the same day (Norm's decisions, 12 October):
- **Scoring.** The PSY240 trial's primary score is Access + Valuing and savouring. Understanding the
  modes and Letting be are secondary, and there is no 20-item total.
- **Code.** `omega()` calls pass `flip = FALSE`, and 05, 06, 10 and 11 were re-run.
  `12_access_valuing.R` was added.
- **Study 2 preregistration.**
  - It now has an ESEM as the confirmatory model, N = 1,000, no screening, and per-subscale H7.
  - H3d (redundancy with neighbours combined) is new.
  - The IUS-12 convergent test was removed, and a wording-factor test for Letting be was added.
  - Missing data, attention checks, scores, fit indices and the life-satisfaction format are fixed.
- **REB amendment.** The documents now report the sample increase and the face-scale format.
- **Findings page.** It is corrected in place, with a correction note and a completed deviations list.

Still open, for Norm:
- the past-week items at the retest;
- non-jargon parallel items;
- a check item before the validity block;
- a curiosity measure for Forage;
- whether the approved protocol names the six countries;
- a decline option on the face-scale step;
- a minimum ΔR² for H5.
