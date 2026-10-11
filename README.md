# Sense Foraging Scale

Development and validation of a self-report measure of how people move between two ways of paying attention:
the **sensing mode** (taking in what is here through the senses, without needing to do anything with it) and the
**doing mode** (using what one knows to plan, solve, judge, fix or get something done). The scale is the outcome
measure for brief training in moving into the sensing mode, including the Sense Foraging Foundations course.

University of Toronto REB Protocol #00051180. PI: Norman Farb (RADlab, University of Toronto).

| Study | What | Status |
|---|---|---|
| [Study 1](study1/) | Exploratory factor analysis of a 48-item pool (item pool 6), Prolific, N = 500 | Data collected 2026-10-11; preregistered at [osf.io/48r2y](https://osf.io/48r2y); analysed |
| [Study 2](study2/) | Confirmatory factor analysis of the 20-item short form, validity measures, two-week retest | Built on the platform (inactive); REB amendment and preregistration drafted |

**Findings for collaborators:** https://senseforaging.com/research/study1 (unlisted).

## Layout

- `study1/`
  - `preregistration.md`: the registration as filed (OSF template sections).
  - `materials/`: the questionnaires exactly as administered (platform definitions), the item key, and the consent and debrief forms.
  - `analysis/`: the R pipeline. Its README explains modes, steps and how to run it.
- `study2/`: materials as built, and the preregistration draft.
- `docs/literature/`: the literature reviews behind design decisions (6-point scale without midpoint; reverse-keyed items and acquiescence).
- `item-development/`: earlier item pools.

## Data

**No participant data are stored here, and none should be committed.** Responses live on the RADlab platform and
are pulled for analysis with `study1/analysis/extract.sql` into a local folder (`SF1_DATA_DIR`). Prolific IDs are held only in
a separate platform table and deleted after payment; the extraction uses random participant codes. A de-identified
dataset will be shared as described in the consent form. `.gitignore` excludes data and model files.

## Running Study 1's analysis

From the repository root, with R ≥ 4.4 and the packages listed in `study1/analysis/config.R`:

```
SF1_MODE=real SF1_K=4 Rscript study1/analysis/run_all.R
```

`SF1_MODE=permuted` (item answers shuffled) and `SF1_MODE=simulated` (data from a known model) exercise the code without
real structure; see `study1/analysis/README.md`.
