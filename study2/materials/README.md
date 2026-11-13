# Study 2 materials

Platform definitions (JSON) of the questionnaires in the Study 2 build, with the consent and debrief forms. The
definitions are the RADlab platform's, as built; `forms_changes.md` lists the forms' changes from the REB-approved text.

## Third-party questionnaires: wording withheld

The wording of the third-party questionnaires is withheld from this repository for copyright. Administering a scale
is not a right to republish it. In each file below, the instructions, the item wording and the response-option labels
are replaced by `[wording withheld: see the cited source]`. Everything the analysis code reads is kept: the item ids,
the numeric response values, reverse keying, subscale membership and scoring, and the citation or source note. The
wording is in the cited sources.

| File | Measure | Source |
|---|---|---|
| `sf2-ffmq15.json` | FFMQ-15 | Gu et al. (2016) |
| `sf2-maia2.json` | MAIA-2 | Mehling et al. (2018) |
| `sf2-sbi-moment.json` | SBI, Savoring the Moment subscale | Bryant (2003) |
| `sf2-bfi2-om.json` | BFI-2, Open-Mindedness domain | Soto & John (2017) |
| `sf2-eq-decentering.json` | EQ, Decentering | Fresco et al. (2007) |
| `sf2-bidr16.json` | BIDR-16 | Hart et al. (2015) |
| `sf2-ius12.json` | IUS-12 | Carleton, Norton & Asmundson (2007) |
| `sf2-phq4.json` | PHQ-4 | Kroenke et al. (2009) |
| `sf2-pss4.json` | PSS-4 | Cohen et al. (1983) |
| `sf2-samd-md.json` | SAMD, Meaning-Detecting subscale | Russo-Netzer & Icekson (2023) |

Purely numeric labels (such as "2" on a 1-7 scale) are kept. The attention checks inside `sf2-bfi2-om.json`
(`bfi2_attn_check`) and `sf2-maia2.json` (`maia_attn_check`) are the lab's own and keep their wording.

The project's own materials keep their wording: the candidate statements (`sf2-sfs-candidates.json`), the Study 2b
vignettes (`sf2b-vignettes.json`), the past-week items (`sf2-past-week.json`), and the consent and debrief forms.
