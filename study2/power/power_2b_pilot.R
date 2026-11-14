# Study 2b power for H4b, H5a and H6 at the pilot's observed rates, N = 700 (preregistration,
# Section 3, Sample size rationale: "Power for H4b, H5a and H6 at the pilot's observed pause
# take-up, missing-data rates and practice-group sizes").
#
# The pilot (extraction of 13 November 2026, 61 completed sessions, run through
# study2b/analysis/ in permuted mode) gives three kinds of rate. Each is an aggregate count from
# that run's outputs (01_pause_check.csv, 02_sample.json, 02_item_descriptives.csv); no answer
# enters here.
#   H4b  the engaged pause: 2 of the 61 who reached the offer (3 took the pause).
#   H6   practice groups among the 61 analysable: 26 current practitioners, 20 never, 15 neither.
#   H5a  missing data: 4 "Prefer not to answer" among the 61 x 131 validity items (one each in
#        Brief MAIA-2, BIDR-16, EQ Decentering and ASI-3; none in FFMQ-15, PHQ-4, PSS-4 or life
#        satisfaction); age missing for 1 of 61; gender "prefer not to say" for 0 of 61.
# Rates from about 60 people are uncertain, so every test is also run at the ends of the rate's
# exact (Clopper-Pearson) 95% interval.
#
# The models are the earlier scripts', unchanged; only the rates differ:
#   H4b  power_2b_review.R's H4 block: the observed 8 with reliability .84, the engaged pause a
#        latent normal correlating r with the 8's true score, cut at the take-up; Spearman,
#        two-sided, Holm across H4a and H4b, H4a (past-week pause days, latent r .20) generated
#        beside it. Supported when rho > 0 and the Holm-adjusted p < .05.
#   H5a  power_2b_h3d.R's H5a block: 16 step-1 columns explaining 30%, delta R^2 >= .02 and
#        p < .05 (and the Holm worst case, .05 / 3), at the projected complete-case n.
#   H6   power_2b.R's H6 block: Student's t, d >= .50 and p < .05, at the group sizes the
#        observed proportions give at 700.
#
# Projected complete-case n for H5a (each of its three regressions needs age, gender, the 15
# FFMQ-15 and 24 Brief MAIA-2 items, since every 3-item subscale needs all three under the 80%
# rule, the outcome's items (life satisfaction 1, PHQ-4 4, PSS-4 4, all needed) and 7 of the 8):
# missingness independent across items and people, so P(complete) is the product over required
# items of (1 - each item's missing rate), times the same for age and gender, and n = 700 x P.
#   per item     each item's own observed rate (0 for most items)
#   pooled       one rate for every validity item, 4 / (61 x 131), and its 95% upper end
# The 8 needs 7 of 8 statements; at the scale's observed rate (0 of 61 x 13) and its upper end
# that term is above .999 and is included.
#
# Run from the repo root:  Rscript study2/power/power_2b_pilot.R
# Writes study2/power/power_2b_pilot.csv.

# Set Up ---------
## Load libraries ---------
packages <- c("data.table")
new_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(new_packages)) install.packages(new_packages)
for (thispack in packages) {
  library(thispack, character.only = TRUE, quietly = TRUE, verbose = FALSE)
}

## Settings ---------
N <- as.integer(Sys.getenv("SF2B_N", "700"))
REPS <- as.integer(Sys.getenv("SF2B_REPS_FAST", "4000"))
SEED <- 20261114
OUT <- Sys.getenv("SF2B_POWER_OUT", "study2/power/power_2b_pilot.csv")
pct <- function(x) round(100 * mean(x, na.rm = TRUE), 1)
cp <- function(x, n) c(est = x / n, lo = if (x == 0) 0 else qbeta(.025, x, n - x + 1), hi = if (x == n) 1 else qbeta(.975, x + 1, n - x))

## The pilot's counts (aggregates from output/permuted/, extraction of 13 November 2026) ---------
PILOT <- list(offered = 61, engaged = 2, took = 3,
              analysable = 61, current = 26, never = 20,
              val_items = 131, val_pna = 4,            # PNA among 61 x 131 validity items
              scale_items = 13, scale_pna = 0,
              age_missing = 1, gender_missing = 0,
              item_missing = c(maia_09 = 1))           # the only H5a item with a missing answer
rates <- list(engaged = cp(PILOT$engaged, PILOT$offered), took = cp(PILOT$took, PILOT$offered),
              current = cp(PILOT$current, PILOT$analysable), never = cp(PILOT$never, PILOT$analysable),
              val_item = cp(PILOT$val_pna, PILOT$analysable * PILOT$val_items),
              scale_item = cp(PILOT$scale_pna, PILOT$analysable * PILOT$scale_items),
              age = cp(PILOT$age_missing, PILOT$analysable), gender = cp(PILOT$gender_missing, PILOT$analysable))
print(rbindlist(lapply(names(rates), function(k) data.table(rate = k, est = signif(rates[[k]][["est"]], 3),
                                                             lo = signif(rates[[k]][["lo"]], 3), hi = signif(rates[[k]][["hi"]], 3)))))
rows <- list()

# ---- H4b at the observed engaged-pause rate -------------------------------------------------------------
set.seed(SEED)
rel <- .84
P_DAYS <- c(.30, .15, .13, .11, .08, .07, .04, .12)
TH_DAYS <- qnorm(cumsum(P_DAYS)[-8])
h4_sim <- function(take, rho_b, rho_a = .20) {
  sim <- replicate(REPS, { t <- rnorm(N); x <- sqrt(rel) * t + sqrt(1 - rel) * rnorm(N)
    days <- findInterval(rho_a * t + sqrt(1 - rho_a^2) * rnorm(N), TH_DAYS)
    pause <- as.numeric(rho_b * t + sqrt(1 - rho_b^2) * rnorm(N) > qnorm(1 - take))
    a <- suppressWarnings(cor.test(x, days, method = "spearman", exact = FALSE))
    b <- suppressWarnings(cor.test(x, pause, method = "spearman", exact = FALSE))
    c(a$estimate, a$p.value, b$estimate, b$p.value, sum(pause)) })
  sim[is.na(sim)] <- 1   # no engaged pause in a sample: no test, not supported
  holm <- apply(sim[c(2, 4), ], 2, p.adjust, method = "holm")
  c(p05 = pct(sim[3, ] > 0 & sim[4, ] < .05), holm = pct(sim[3, ] > 0 & holm[2, ] < .05), engaged_median = median(sim[5, ]))
}
for (lab in c("est", "lo", "hi")) for (rho_b in c(.20, .30)) {
  take <- rates$engaged[[lab]]
  r <- h4_sim(take, rho_b)
  rows[[length(rows) + 1]] <- data.table(test = "H4b", rate = sprintf("engaged pause %.1f%% (%s)", 100 * take,
      c(est = "observed, 2 of 61", lo = "95% CI lower end", hi = "95% CI upper end")[[lab]]),
    condition = sprintf("latent r %.2f; Holm with H4a at latent r .20", rho_b), n = N,
    power = r[["holm"]], power_variant = r[["p05"]], detail = sprintf("power_variant = p < .05 without Holm; median %d engaged of %d", as.integer(r[["engaged_median"]]), N))
  print(rows[[length(rows)]])
}
# for reference: take-up 3 of 61 (any pause taken, the engaged rule not applied)
r <- h4_sim(rates$took[["est"]], .20)
rows[[length(rows) + 1]] <- data.table(test = "H4b", rate = sprintf("pause taken %.1f%% (3 of 61; reference only)", 100 * rates$took[["est"]]),
  condition = "latent r 0.20; Holm with H4a at latent r .20", n = N, power = r[["holm"]], power_variant = r[["p05"]],
  detail = sprintf("power_variant = p < .05 without Holm; median %d took of %d", as.integer(r[["engaged_median"]]), N))

# ---- H5a at the projected complete-case n ----------------------------------------------------------------
H5_ITEMS <- 15 + 24                       # FFMQ-15 and Brief MAIA-2 items, all needed
OUTCOME_ITEMS <- c(lifesat = 1, phq4_total = 4, pss4 = 4)
p8 <- function(q) 1 - pbinom(6, 8, 1 - q)  # P(at least 7 of the 8 answered)
proj <- function(q_item, q_scale, q_age, q_gender, outcome, per_item = NULL) {
  covs <- if (is.null(per_item)) (1 - q_item)^(H5_ITEMS + OUTCOME_ITEMS[[outcome]]) else per_item
  N * covs * p8(q_scale) * (1 - q_age) * (1 - q_gender)
}
per_item_term <- prod(1 - PILOT$item_missing / PILOT$analysable)   # every other H5a item: 0 missing in 61
H5N <- rbindlist(lapply(names(OUTCOME_ITEMS), function(o) data.table(outcome = o,
  per_item = proj(NA, rates$scale_item[["est"]], rates$age[["est"]], rates$gender[["est"]], o, per_item = per_item_term),
  pooled = proj(rates$val_item[["est"]], rates$scale_item[["est"]], rates$age[["est"]], rates$gender[["est"]], o),
  pooled_upper = proj(rates$val_item[["hi"]], rates$scale_item[["hi"]], rates$age[["hi"]], rates$gender[["hi"]], o))))
H5N[, (2:4) := lapply(.SD, round), .SDcols = 2:4]
cat("\n==== projected complete-case n for H5a at N =", N, "\n"); print(H5N)
set.seed(SEED + 1)
h5_sim <- function(n, dr2) {
  sim <- replicate(REPS, { x <- matrix(rnorm(n * 16), n); sc <- rnorm(n); bc <- sqrt(.30 / 16)
    y <- drop(x %*% rep(bc, 16)) + sqrt(dr2) * sc + sqrt(1 - .30 - dr2) * rnorm(n)
    m0 <- lm(y ~ x); m1 <- lm(y ~ x + sc)
    c(summary(m1)$r.squared - summary(m0)$r.squared, anova(m0, m1)$`Pr(>F)`[2]) })
  c(p05 = pct(sim[1, ] >= .02 & sim[2, ] < .05), holm_worst = pct(sim[1, ] >= .02 & sim[2, ] < .05 / 3))
}
ns <- c(per_item = min(H5N$per_item), pooled_upper = min(H5N$pooled_upper))
for (k in names(ns)) for (dr2 in c(.03, .04)) {
  r <- h5_sim(ns[[k]], dr2)
  rows[[length(rows) + 1]] <- data.table(test = "H5a", rate = c(per_item = "each item's observed missing rate",
      pooled_upper = "95% CI upper ends: validity items 0.13%, age 8.8%, gender 5.9%")[[k]],
    condition = sprintf("true delta R2 %.2f", dr2), n = ns[[k]], power = r[["p05"]], power_variant = r[["holm_worst"]],
    detail = "delta R2 >= .02 and p < .05; power_variant = Holm worst case (.05/3)")
  print(rows[[length(rows)]])
}

# ---- H6 at the observed group sizes -------------------------------------------------------------------------
set.seed(SEED + 2)
h6_sim <- function(n1, n0, d) pct(replicate(REPS, { a <- rnorm(n1, d); b <- rnorm(n0)
  dh <- (mean(a) - mean(b)) / sqrt((var(a) + var(b)) / 2); dh >= .50 && t.test(a, b, var.equal = TRUE)$p.value < .05 }))
# The observed proportions, and both at their 95% lower ends at once (the pessimistic corner;
# both upper ends together exceed 100% and are not a possible sample).
for (lab in c("est", "lo")) for (d in c(.55, .67)) {
  n1 <- round(N * rates$current[[lab]]); n0 <- round(N * rates$never[[lab]])
  rows[[length(rows) + 1]] <- data.table(test = "H6", rate = sprintf("current %.1f%%, never %.1f%% (%s)", 100 * rates$current[[lab]], 100 * rates$never[[lab]],
      c(est = "observed, 26 and 20 of 61", lo = "both at their 95% CI lower ends")[[lab]]),
    condition = sprintf("true d %.2f", d), n = n1 + n0, power = h6_sim(n1, n0, d), power_variant = NA_real_,
    detail = sprintf("%d current vs %d never; d >= .50 and p < .05", n1, n0))
  print(rows[[length(rows)]])
}

res <- rbindlist(rows)
fwrite(res, OUT)
fwrite(H5N, sub(".csv", "_h5_n.csv", OUT, fixed = TRUE))
print(res[, .(test, rate, condition, n, power)])
cat(sprintf("\nN = %d, %d replications per condition, seed %d\n", N, REPS, SEED))
