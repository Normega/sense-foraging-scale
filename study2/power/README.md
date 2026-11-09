# Study 2b power simulation (N = 500)

`power_2b.R`, run 9 November 2026 (seed 20261109; 500 simulated studies per SEM condition, 4,000 for the others). Populations come from Study 2a's fitted one-factor model of the 8 with the agreement factor, so loadings, agreement variance and thresholds are 2a's. Results: `power_2b_results.csv` (every test), `power_2b_h1_grid.csv` (fit-rule variants). `power_2b_h1_raw.rds` holds the H1 replications.

## H1: the fit rule (decision 17)

Pass rates, %. Populations: the model true (*fit*); 2b exactly like 2a (*observed*: 2a's polychoric correlations); two residual correlations of .30 (*small flaw*); three statements sharing a second factor loading .40 or .60.

| Rule (all with CFI ≥ .95 and the SRMR close-fit test) | fit | observed | small flaw | 2nd factor .40 | 2nd factor .60 |
|---|---|---|---|---|---|
| Registered: residual ≤ .10, or freeing it moves ω and loadings < .02 | 100 | **55** | 88 | 33 | 0 |
| Residual ≤ .10, or freeing it moves ω < .03 and loadings < .05 | 100 | **96** | 98 | 72 | 1 |
| Residual ≤ .12, or freeing it moves ω < .03 and loadings < .05 | 100 | 97 | 100 | 84 | 1 |
| Residual ≤ .15 | 100 | 99 | 100 | 97 | 1 |

- CFI ≥ .95 and the close-fit test pass in every population; the local rule decides.
- As registered, H1 is a coin flip if 2b looks like 2a: the largest sample residual of 2a-like data exceeds .10 in 45% of samples (median .098), and freeing it moves some loading by more than .02 (median .029).
- 2a's own misfit is about the size of a three-statement second factor loading .40; no rule separates them at N = 500. A second factor loading .60 fails every rule.
- Recommended: the second row (ω < .03, loadings < .05).

## H2–H6

| Test | Condition | Power % |
|---|---|---|
| H2, CI lower bound > .30 (3-item neighbour) | true r = .40 / .45 / .50 | 54 / 84 / 98 |
| H3b contrast, difference .15 | 3 items per measure, measures correlate 0 / .30 / .50 | 54 / 63 / 69 |
| | 4 items (loading .75) / 8 items, measures correlate .30 | 72 / 76 |
| | r .35 vs .20, 3 items | 60 |
| H3b contrast, difference .20 | 3 items, measures correlate .30 | 87 |
| H3b contrast, difference .10 | 3 items, measures correlate .50 | 34 |
| H4, the pause (0/1), Spearman | take-up 15 / 30 / 50%, latent r .20 | 75 / 87 / 89 |
| | latent r .25 | 91 / 97 / 98 |
| H5, ΔR² ≥ .02 and p < .05 | true ΔR² .02 / .03 / .04 | 50 / 80 / 95 |
| H6, d ≥ .50 and p < .05 (273 current vs 136 never) | true d .55 / .67 / .80 | 68 / 94 / 100 |

- A difference of .15 between two latent correlations is not well powered at N = 500 (54–76%, by the number of items in the measures compared). For 80% at .15, N ≈ 700; at N = 500, about 85–90% for .20.
- H2 needs a true correlation of about .45 or more for 80%.
- H6 uses 2a's groups (current practitioners vs never), not Study 1's "dedicated" definition (34 of 500 in Study 1).
- The contrasts used the delta-method CI from the latent model; the old draft registered a 2,000-resample bootstrap. Either can be registered; this is what was simulated.
