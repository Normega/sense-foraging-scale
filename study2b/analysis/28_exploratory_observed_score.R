# 28: the observed score users compute: norms, floor and ceiling, standard error of
# measurement, and observed beside latent correlations. EXPLORATORY, NOT REGISTERED.
#
# Not part of the filed registration (osf.io/d4a5u) and not in its list of exploratory
# analyses. Specified on 17 November 2026 (placement card S3-24, from the Fable review),
# after the registered analyses and 10-27 had been seen.
#
#   norms     the observed 8 (mean of the 8 after reverse-scoring, 7 of 8 answered, as
#             01 scored it): n, mean, SD, skew, excess kurtosis, percentiles 5-95,
#             overall, by gender (gender2) and by practice group (H6's groups and the rest)
#   ends      share at the floor (1.00), at the ceiling (6.00), and at or above 5.50
#   SEM       SD x sqrt(1 - reliability), with alpha (complete cases, as 03) and with the
#             registered model's content omega (.843); the 95% band is +/- 1.96 SEM
#   observed vs latent   the observed r of the 8 with every validity score beside the
#             latent r (04_latent_validity.csv, the 8, all participants)
# Writes aggregates only:
#   output/real/28_norms.csv, 28_ends.csv, 28_sem.csv, 28_observed_vs_latent.csv,
#   28_exploratory_observed_score.json
# Run from the repo root:  SF2B_MODE=real Rscript study2b/analysis/28_exploratory_observed_score.R

# Set Up ---------
## Load libraries ---------
packages <- c("jsonlite", "data.table", "psych")
new_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(new_packages)) install.packages(new_packages)
for (thispack in packages) {
  library(thispack, character.only = TRUE, quietly = TRUE, verbose = FALSE)
}

## Settings ---------
source("study2b/analysis/config.R")
source("study2b/analysis/sf2b_models.R")
stopifnot(SF2B$mode == "real")
T0 <- Sys.time()
A <- readRDS(sf2b_derived("sf2b")); D <- A$defs; M <- A$meta; SC <- A$scores; X <- A$scale_raw
stopifnot(identical(SC$participant, M$participant))
y <- SC$sf8
h1 <- fromJSON(sf2b_out("03_h1_decision.json"))
omega8 <- h1$h1$omega
stopifnot(is.numeric(omega8), length(omega8) == 1)

## Norms ---------
desc <- function(v, group) {
  v <- v[!is.na(v)]
  q <- quantile(v, c(.05, .10, .25, .50, .75, .90, .95), names = FALSE)
  data.table(group = group, n = length(v), mean = mean(v), sd = sd(v), skew = psych::skew(v), kurtosis = psych::kurtosi(v),
             p05 = q[1], p10 = q[2], p25 = q[3], p50 = q[4], p75 = q[5], p90 = q[6], p95 = q[7],
             floor = mean(v <= 1 + 1e-9), ceiling = mean(v >= 6 - 1e-9), at_or_above_5_5 = mean(v >= 5.5 - 1e-9))
}
pg <- fifelse(is.na(M$practice_group), "other", M$practice_group)
norms <- rbindlist(c(list(desc(y, "all")),
  lapply(sort(unique(na.omit(M$gender2))), function(g) desc(y[M$gender2 %in% g], paste0("gender: ", g))),
  lapply(sort(unique(pg)), function(g) desc(y[pg == g], paste0("practice: ", g)))))

## SEM ---------
items <- D$sf8; Xr <- X[, items]; storage.mode(Xr) <- "double"
rv <- intersect(items, D$reversed); Xr[, rv] <- 7 - Xr[, rv]
alpha8 <- psych::alpha(Xr[complete.cases(Xr), ], check.keys = FALSE, warnings = FALSE)$total$raw_alpha
sdy <- sd(y, na.rm = TRUE)
sem <- data.table(reliability = c("alpha", "omega (registered model)"), value = c(alpha8, omega8))
sem[, sem := sdy * sqrt(1 - value)][, band_95 := 1.96 * sem][, sd := sdy]

## Observed vs latent ---------
lat <- fread(sf2b_out("04_latent_validity.csv"))[form == "sf8" & sample == "all", .(score, latent_r = est, latent_lo = ci.lower, latent_hi = ci.upper)]
ovl <- rbindlist(lapply(lat$score, function(s) { r <- sf2b_r_ci(y, SC[[s]])
  data.table(score = s, n = r$n, observed_r = r$r, observed_lo = r$lower, observed_hi = r$upper) }))
ovl <- merge(lat, ovl, by = "score", sort = FALSE)
ovl[, ratio := observed_r / latent_r]

## Write ---------
num <- function(DT, k = 4) { cols <- names(DT)[vapply(DT, is.numeric, TRUE)]; DT[, (cols) := lapply(.SD, function(x) round(x, k)), .SDcols = cols]; DT }
fwrite(num(copy(norms)), sf2b_out("28_norms.csv")); fwrite(num(copy(sem)), sf2b_out("28_sem.csv"))
fwrite(num(copy(ovl)), sf2b_out("28_observed_vs_latent.csv"))
write_json(list(
  analysis = "28: the observed 8's norms, ends, SEM and observed vs latent correlations. EXPLORATORY, NOT REGISTERED (osf.io/d4a5u does not include it); specified 17 November 2026 after 10-27 were seen (card S3-24).",
  scoring = "Mean of the 8 after reverse-scoring (7 - x), 7 of 8 answered (01).", alpha_complete_cases = round(alpha8, 4), omega = omega8,
  norms = num(copy(norms)), sem = num(copy(sem)), observed_vs_latent = num(copy(ovl)),
  minutes = round(as.numeric(difftime(Sys.time(), T0, units = "mins")), 1)),
  sf2b_out("28_exploratory_observed_score.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
print(num(copy(norms), 3)); print(num(copy(sem), 3)); print(num(copy(ovl), 3)[, .(score, latent_r, observed_r, ratio)], nrows = 50)
sf2b_log("28 done")
