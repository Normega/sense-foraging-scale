# Study 2: confirmatory structure, validity and stability of the short form

**Status:** built on the RADlab platform, inactive, pending REB approval of the Study 2 amendment.
The preregistration is drafted, not filed.

- `preregistration_draft.md`: hypotheses H1–H7 and the predicted validity matrix. Points still to settle are marked [Norm].
- `materials/`: every questionnaire exactly as the platform administers it. The short form is
  `../study1/materials/sf-short-20.json`; background questions are `../study1/materials/sf-background.json`.
  The folder also holds the consent and debrief forms, for the main session and the follow-up, and
  `forms_changes.md`, which lists every change from the approved Rev2 forms.
- `analysis/`: the analysis pipeline, with simulated, permuted and real modes. Its README gives the
  run order, every output, and each reading of the preregistration the code had to choose.

## Design

One session of about 20 minutes, on Prolific, with a general posting (600) and a contemplative posting (400); 1,000 in all since
the 12 October revision. The postings differ in framing only: no one is screened on practice, as the approved protocol says.
There is no country filter (decided 12 October): adults fluent in English anywhere, as Study 1 recruited in practice.

1. Background questions (as in Study 1).
2. The 20-item short form, in one fixed order.
3. Twenty-four candidate items, `sf2-sfs-candidates`, unscored, in a per-participant order: 19 new (5 of them Forage, the seeking itself) and 5 Study 1 items with unchanged wording (`study1_item` maps them). Replacements are chosen on a random half of the sample; any revised form is checked on the other half.
4. The validity block, in a random order for each participant:

   | Measure | Purpose |
   |---|---|
   | FFMQ-15 | approved in the original protocol |
   | MAIA-2 | approved in the original protocol |
   | BIDR-16 | approved in the original protocol |
   | PHQ-4 | approved in the original protocol |
   | Life satisfaction | approved in the original protocol |
   | IUS-12 | new: discriminant for Letting be |
   | SBI Savoring the Moment | new: convergent for Valuing and savouring |
   | BFI-2 Open-Mindedness | new: discriminant, appetite for experience in general |
   | EQ Decentering | new |
   | PSS-4 | new: criterion |
   | Three past-week behaviour items | new: criterion |

**Follow-up.** About 14 days later, about 400 random completers are invited to repeat the short form and the three
past-week items, aiming for about 300 completions (raised from 200 and 150 on 12 October). This is a separate platform study, linked to Study 2 through `external_identity_links`.

## Licences

Sources and licence notes for each new measure are in the amendment's `Study2_Measures.docx`.

| Measure | Licence position |
|---|---|
| PSS-4 | Needs a Mapi ePROVIDE licence for electronic use |
| SBI | Permission requested from F. Bryant |
| BFI-2 | Free for research use |
| EQ | Free with citation |
| IUS-12 | No stated terms |
