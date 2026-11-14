# Sense Foraging Scale

Development and validation of a self-report measure of how people move between two ways of paying attention:
the **sensing mode** (taking in what is here through the senses, without needing to do anything with it) and the
**doing mode** (using what one knows to plan, solve, judge, fix or get something done). The scale is the outcome
measure for brief training in moving into the sensing mode, including the Sense Foraging Foundations course.

University of Toronto REB Protocol #00051180. PI: Norman Farb (RADlab, University of Toronto).

| Study | What | Status |
|---|---|---|
| [Study 1](study1/) | Exploratory factor analysis of a 48-item pool (item pool 6), Prolific, N = 500 | Preregistered at [osf.io/48r2y](https://osf.io/48r2y); data collected 2026-10-11; analysed |
| [Study 2a](study2a/) | Factor structure of a revised 58-statement set, with vignettes checking that the two modes are understood, Prolific, N = 500 | Preregistered at [osf.io/vt7gn](https://osf.io/vt7gn) (29 October 2026); data collected and analysed |
| [Study 2b](study2b/) | Confirmatory structure and construct validity of the 8-statement scale chosen in 2a: validity measures, an in-session pause as a behavioural criterion, known groups; Prolific, N = 700 planned | Preregistered at [osf.io/d4a5u](https://osf.io/d4a5u) (14 November 2026, from commit d3e60c9); 61-person pilot run; main sample collecting |

**Study 2a in brief.** The 58 statements showed one broad sensing factor, small reverse-keyed clusters and an
agreement (acquiescence) factor; the hypothesized parts did not separate as subscales. From 2a, in-sample, came the
scale Study 2b tests: **the 8**, scored as one number, and **the 10** (the 8 plus two statements) for studies of
change. The decisions and the 2a results behind them are in `study2/decisions_after_2a.md`.

**Study 2b.** The preregistration text is `study2/preregistration_2b.md` (filed on OSF at
[osf.io/d4a5u](https://osf.io/d4a5u) from commit d3e60c9, after the 61-person pilot; OSF governs). The analysis code is `study2b/analysis/`, validated
on simulated data; the platform definitions are `study2b/materials/`; the power simulations are `study2/power/`.

## Layout

- `study1/`
  - `preregistration.md`: the registration as filed (OSF template sections).
  - `materials/`: the questionnaires exactly as administered (platform definitions), the item key, and the consent and debrief forms.
  - `analysis/`: the R pipeline. Its README explains modes, steps and how to run it.
- `study2a/`
  - `preregistration_draft.md`: the text the registration was filed from (OSF governs).
  - `materials/`: the statements, vignettes and end statements as built (platform definitions and `items.json`).
  - `analysis/`: the R pipeline: the registered analysis (`01`–`11`), the exploratory follow-ups (`12`–`21`) and the report data (`22`).
  - `appendices/`: the registration's appendices; `make_appendices.py` regenerates them.
- `study2/`: the Study 2 design that Study 2b replaced, and the shared 2b records.
  - `preregistration_2b.md`: the Study 2b preregistration text.
  - `decisions_after_2a.md`, `prereg_2b_plan.md`, `criterion_pause_design.md`: how 2b was decided.
  - `power/`: the 2b power simulations.
  - `filing/`: builds the OSF filing text from one commit.
  - `materials/`: the Study 2 platform definitions, including the validity measures that Study 2b shares.
  - `analysis/`: the pre-2a Study 2 pipeline (retired 9 November 2026), kept for the record.
- `study2b/`
  - `materials/`: the platform definitions of every questionnaire in the session.
  - `analysis/`: the R pipeline. Its README explains modes, steps, the pilot run and how to run it.
- `docs/`
  - `literature/`: the literature reviews behind design decisions (6-point scale without midpoint; reverse-keyed items and acquiescence).
  - `reviews/`: re-analysis scripts from an independent review of Study 1 (12 October 2026).
- `item-development/`: an earlier item pool (pool 5).

## Third-party questionnaires

The validity measures (FFMQ-15, MAIA-2 and Brief MAIA-2, SBI, BFI-2, EQ, BIDR-16, IUS-12, PHQ-4, PSS-4, SAMD, ASI-3,
BSSS-8, NCS-6, SPSQ-SF) are other authors' work. Their wording is withheld from this repository for copyright: the
definitions keep item ids, keying, scoring and numeric response values, and cite their sources. See
`study2/materials/README.md` and `study2b/materials/README.md`. The Sense Foraging statements, vignettes, attention
checks, background and past-week questions, and the consent and debrief forms are the project's own and are given in full.

## Data

**No participant data are in this repository, and none should be committed.** Responses live on the RADlab platform
and are pulled for analysis with each study's `analysis/extract.sql` into a local folder outside the repository.
Prolific IDs are held only in a separate platform table and deleted after payment; the extractions use random
participant codes. A de-identified dataset will be shared as described in the consent forms. `.gitignore` excludes
data and model files.

## Running the analyses

From the repository root, with R ≥ 4.4 and the packages listed in each study's `analysis/config.R`. Each pipeline has
a `simulated` mode (data from a known model) and a `permuted` mode (real answers shuffled within each item) that
exercise the code without real structure; see each `analysis/README.md`. For example:

```
SF1_MODE=real SF1_K=4 Rscript study1/analysis/run_all.R
SF2B_MODE=simulated Rscript study2b/analysis/run_all.R
```

## Licence

Code is under the MIT licence (`LICENSE`). Text, preregistrations and the project's own materials are under CC BY 4.0
(`LICENSE-CC-BY-4.0.md`). The third-party questionnaires are covered by neither, and their wording is not in the
repository.
