# Study 1: exploratory factor structure of item pool 6

- **Design:** 48 agree–disagree statements (34 positively keyed, 14 written as logical opposites of named items), rated on a 6-point scale with no midpoint, every point labelled, plus "Prefer not to answer". The statements were shown in a constrained random order per participant, together with two instructed-response attention checks.
- **Sample:** Prolific, two postings with different framing (general and contemplative-leaning), 11 October 2026. Preregistered at [osf.io/48r2y](https://osf.io/48r2y), registered before any structural analysis.
- **Full write-up:** https://senseforaging.com/research/study1

## Results (summary)

- **Sample:** 516 consented, 502 submitted the scale, 2 failed both attention checks, **500 analysed**.
  - Median 6.4 minutes on the statements.
  - 53% practise at least weekly; 29% never have.
- **Factors:** four content factors, not the seven hypothesized; the seven-factor solution failed the registered criterion.
  1. **Access:** getting into the sensing mode; strongest factor. Reversed items 35 and 36 load highest.
  2. **Understanding the modes:** insight into the pull of doing, and valuing discovery.
  3. **Valuing and savouring:** valuing time in the sensing mode, and its felt reward.
  4. **Letting be:** versus needing to know or add meaning; carried largely by reversed items together with item 21.
- **General agreement (random-intercept factor):** about 10% of each item's variance. It explains the "reversed-items factor" that appears when agreement is not modelled.
- **General factor:** ω_h = .70, ω_total = .93.
- **Proposed short form:** 17 items. Ordinal α .85 (Access) and about .66 for the three 4-item subscales.
- **Practitioners vs never-practitioners:** no difference on any factor (|d| ≤ .16).
- **Robustness:** the structure holds without flagged respondents, without translating respondents, and on positive items alone (congruence .84–1.00).

## Short form (decided 2026-10-11)

**20 items:** four factors, five items each. Every factor mixes forward and reverse-keyed (R) items, and the selection was checked on agreement-free reliability as well as raw reliability:

| Factor | Items |
|---|---|
| Access | 35R, 36R, 9, 8, 7 |
| Understanding the modes | 22, 3, 15, 31, 37R |
| Valuing and savouring | 12, 18, 11, 27, 44R |
| Letting be | 40R, 46R, 43R, 47R, 21 |

- **How it departs from the registered rules:**
  - item 7 replaces 1, the near-duplicate of 9;
  - 37R and 21 are added for keying balance;
  - 27 is added for reliability.
- **How it is presented:** one fixed order (`materials/sf-short-20.json`), the same for everyone and at every timepoint. It is the PSY240 course trial's outcome from 17 October 2026.
- **Where new items would help:**
  - Understanding the modes, the weakest factor, which has no working reverse item;
  - a second Valuing reversal;
  - forward Letting-be items.

  Candidates are drafted in the findings page; if used, they would run as unscored extras in Study 2.

## Deviations from the registration

1. **Country:** the postings were not filtered by country, so 69% of the sample lived outside the six planned countries.
2. **Number of factors:** k = 4, not parallel analysis's 5. The fifth factor was agreement, and fails the registered criterion once agreement is modelled.
3. **Factor orientation:** each factor was oriented before the selection rules were applied.
4. **Short form:** 17 items with 2 qualifying reversed pairs. The registered targets (21–28 items, 4–6 pairs) assumed seven factors.

## Materials

`materials/` holds the platform definitions exactly as run:
- `sf-pool6.json`: the scale, including the interleave rules, factor key and reversal pairs;
- `sf-background.json`: the background questions;
- `item_key.csv`;
- the consent and debrief, with `forms_changes.md` listing their changes from the REB-approved Rev2 text.
