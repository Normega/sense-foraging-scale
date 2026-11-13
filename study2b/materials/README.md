# Study 2b materials

Platform definitions (JSON) of every questionnaire in the Study 2b session, as built on the RADlab platform.
`study2b/analysis/config.R` checks the item ids, keys and subscales the analysis uses against these files.

## Third-party questionnaires: wording withheld

The wording of the third-party questionnaires is withheld from this repository for copyright. Administering a scale
is not a right to republish it. In each file below, the instructions, the item wording and the response-option labels
are replaced by `[wording withheld: see the cited source]`. Everything the analysis code reads is kept: the item ids,
the numeric response values, reverse keying, subscale membership and scoring, and the citation or source note. The
wording is in the sources cited in the preregistration (`study2/preregistration_2b.md`, Section 4 and the References)
and below.

| File | Measure | Source |
|---|---|---|
| `sf2-ffmq15.json` | FFMQ-15 | Gu et al. (2016) |
| `sf2b-maia2-brief.json` | Brief MAIA-2 | Rogowska, Tataruch & Klimowska (2023); MAIA-2 items (Mehling et al., 2018) |
| `sf2-sbi-moment.json` | SBI, Savoring the Moment subscale | Bryant (2003) |
| `sf2-bfi2-om.json` | BFI-2, Open-Mindedness domain | Soto & John (2017) |
| `sf2-eq-decentering.json` | EQ, Decentering | Fresco et al. (2007) |
| `sf2-bidr16.json` | BIDR-16 | Hart et al. (2015) |
| `sf2-phq4.json` | PHQ-4 | Kroenke et al. (2009) |
| `sf2-pss4.json` | PSS-4 | Cohen et al. (1983) |
| `sf2-samd-md.json` | SAMD, Meaning-Detecting subscale | Russo-Netzer & Icekson (2023) |
| `sf2b-asi3-pc.json` | ASI-3, Physical Concerns subscale | Taylor et al. (2007) |
| `sf2b-bsss8.json` | BSSS-8 | Hoyle et al. (2002) |
| `sf2b-ncs6.json` | NCS-6 | Lins de Holanda Coelho, Hanel & Wolf (2020) |
| `sf2b-spsq-sd.json` | SPSQ-SF, Sensory Discomfort subscale | De Gucht & Woestenburg (2025) |

The `sf2-*` files are the measures shared with the earlier Study 2 build (`study2/materials/`, identical files).
Purely numeric labels (such as "2" on a 1-7 scale) are kept. The attention check inside `sf2-bfi2-om.json`
(`bfi2_attn_check`) is the lab's own and keeps its wording.

The project's own materials keep their wording: the scale (`sf2b-scale.json`), the trial statements (`sf2b-end.json`),
the vignettes (`sf2b-vignettes.json`), the pause (`sf2b-pause.json`), the background questions (`sf-background.json`)
and the past-week items (`sf2-past-week.json`). The life-satisfaction face scale is a platform activity, not a file here.
