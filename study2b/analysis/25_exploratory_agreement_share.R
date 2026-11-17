# 25: how much of the observed score is agreement, for the 8 and for the 7 (the 8
# without U4). EXPLORATORY, NOT REGISTERED.
#
# Not part of the filed registration (osf.io/d4a5u) and not in its list of exploratory
# analyses. Specified on 17 November 2026, after the registered analyses and 10-24 had
# been seen, to inform whether to recommend dropping U4 (sf2b_07), which loads .18 on
# the content factor. Dropping it leaves 5 forward-keyed and 2 reverse-keyed statements,
# so agreement (the RI factor) no longer cancels as well in the keyed score, and 24
# found agreement differs between groups (race/ethnicity, practice).
#
# Three parts, each for the 8 and for the 7:
#   1. latent-response metric: the registered H1 model (sf2b_h1: S + RI, WLSMV, theta),
#      variance of the keyed sum of the latent responses split into content, agreement
#      and residual, from the standardized solution (content share = sf2b_omega).
#   2. observed metric: the same model fitted to the 1-6 answers as continuous (MLR,
#      FIML), and the model-implied variance of the keyed observed sum split the same
#      way. This is the metric of the score users compute.
#   3. group differences in the observed score (the mean score, 80% rule): d from 07's
#      sf2b_kg for 24's comparisons, and an observed-metric MIMIC (S ~ g, RI ~ g, MLR)
#      splitting the expected mean-score difference into its content part
#      (sum of keyed content loadings x beta_S) and its agreement part
#      (sum of keyed RI loadings x beta_RI; RI loadings are 1, so that sum is the number
#      of forward-keyed minus reverse-keyed statements), each over the pooled SD.
#      The MIMIC assumes equal intercepts; 24 found the intercepts step held with the RI
#      mean free in every comparison (ordinal metric).
#
# Race/ethnicity and practice groups are 24's, built by evaluating 24's top-level
# expressions up to its means models (option ids only; free text never read).
# Writes aggregates only:
#   output/real/25_variance_shares.csv
#   output/real/25_group_differences.csv
#   output/real/25_exploratory_agreement_share.json
# Run from the repo root:  SF2B_MODE=real Rscript study2b/analysis/25_exploratory_agreement_share.R

# Set Up ---------
## Load libraries ---------
packages <- c("jsonlite", "data.table", "lavaan")
new_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(new_packages)) install.packages(new_packages)
for (thispack in packages) {
  library(thispack, character.only = TRUE, quietly = TRUE, verbose = FALSE)
}

## Settings: 24's data and groups, nothing fitted ---------
ex24 <- parse("study2b/analysis/24_exploratory_invariance_race.R", keep.source = FALSE)
stop_at <- which(vapply(ex24, function(e) is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name("means_syntax")), TRUE))
stopifnot(length(stop_at) == 1)
for (i in seq_len(stop_at - 1)) eval(ex24[[i]], envir = globalenv())
stopifnot(exists("GROUPS"), exists("X"), exists("D"), !exists("INV"))
T0 <- Sys.time()

REV <- D$reversed
stopifnot(D$source[["sf2b_07"]] == "U4", "sf2b_07" %in% REV)
SETS <- list(eight = D$sf8, seven = setdiff(D$sf8, "sf2b_07"))
CHOICES <- list(
  seven = "The 7 is the registered 8 without U4 (sf2b_07), as in 15. No replacement statement.",
  keyed_sum = "Shares are of the keyed SUM (reverse-keyed statements recoded); the mean score has the same shares.",
  observed_model = "Observed metric: the H1 structure (S variance 1, loadings free; RI loadings 1 on every unrecoded statement, variance free, orthogonal to S) fitted to the answers as continuous, MLR, missing = ml, all rows.",
  scores = "Observed scores: mean of answered statements, 7 of 8 for the 8, 6 of 7 for the 7 (15's rule).",
  mimic = "Observed-metric MIMIC: the observed model plus S ~ g and RI ~ g (g = 0 reference, 1 comparison group), MLR, missing = ml, rows with a group. Parts of the expected mean-score difference = (keyed loading sum x beta) / k, over the pooled observed SD. The two parts need not add to the observed d exactly.")

keying <- function(items) ifelse(items %in% REV, -1, 1)
shares <- function(lam_s, lam_r, psi_r, theta, k) {
  cs <- sum(k * lam_s)^2; ag <- sum(k * lam_r)^2 * psi_r; er <- as.numeric(t(k) %*% theta %*% k)
  tot <- cs + ag + er
  list(content = cs / tot, agreement = ag / tot, residual = er / tot, agreement_to_content = ag / cs)
}

## 1 and 2: variance shares ---------
SH <- list(); FITS_OBS <- list()
for (s in names(SETS)) {
  it <- SETS[[s]]; k <- keying(it)
  h1 <- sf2b_h1(X, it, REV); stopifnot(h1$converged)
  st <- lavInspect(h1$fit, "std")
  a <- shares(st$lambda[it, "S"], st$lambda[it, "RI"], st$psi["RI", "RI"], st$theta[it, it], k)
  SH[[length(SH) + 1]] <- data.table(set = s, metric = "latent response (WLSMV)", n_items = length(it), n_reversed = sum(k < 0),
                                     content = a$content, agreement = a$agreement, residual = a$residual,
                                     agreement_to_content = a$agreement_to_content, omega_check = sf2b_omega(h1$fit, it, REV))
  fo <- suppressWarnings(lavaan::sem(paste(sf2b_scale_syntax(it), collapse = "\n"), data = as.data.frame(X[, it]),
                                     estimator = "MLR", missing = "ml"))
  stopifnot(lavInspect(fo, "converged"))
  e <- lavInspect(fo, "est")
  b <- shares(e$lambda[it, "S"], e$lambda[it, "RI"], e$psi["RI", "RI"], e$theta[it, it], k)
  SH[[length(SH) + 1]] <- data.table(set = s, metric = "observed (MLR)", n_items = length(it), n_reversed = sum(k < 0),
                                     content = b$content, agreement = b$agreement, residual = b$residual,
                                     agreement_to_content = b$agreement_to_content, omega_check = NA_real_)
  FITS_OBS[[s]] <- sf2b_fit_stats(fo)
  sf2b_log(s, ": agreement share latent ", round(a$agreement, 4), ", observed ", round(b$agreement, 4))
}
shares_dt <- rbindlist(SH)

## 3: group differences ---------
need <- c(eight = 7, seven = 6)
SC <- lapply(names(SETS), function(s) sf2b_score(X, SETS[[s]], REV, need = need[[s]])); names(SC) <- names(SETS)
PAIRS <- list(race2 = c("1_white", "2_not_white"), race_bw = c("1_white", "2_black"),
              race_other = c("1_white", "3_other_not_white"), practice = c("1_never", "2_current"))
PAIR_SRC <- c(race2 = "race2", race_bw = "race_bw", race_other = "race3", practice = "practice")
GD <- list()
for (p in names(PAIRS)) {
  gr <- GROUPS[[PAIR_SRC[[p]]]]$g
  g01 <- ifelse(gr == PAIRS[[p]][2], 1L, ifelse(gr == PAIRS[[p]][1], 0L, NA_integer_))
  for (s in names(SETS)) {
    it <- SETS[[s]]; k <- keying(it)
    kg <- sf2b_kg(SC[[s]], g01)
    dat <- as.data.frame(X[, it]); dat$g <- g01; dat <- dat[!is.na(dat$g), ]
    syn <- c(sf2b_scale_syntax(it), "S ~ bS*g", "RI ~ bR*g")
    fm <- suppressWarnings(lavaan::sem(paste(syn, collapse = "\n"), data = dat, estimator = "MLR", missing = "ml", fixed.x = TRUE))
    stopifnot(lavInspect(fm, "converged"))
    pe <- parameterEstimates(fm); lam <- lavInspect(fm, "est")$lambda
    bS <- pe$est[pe$label == "bS"]; bR <- pe$est[pe$label == "bR"]
    ks <- sum(k * lam[it, "S"]); kr <- sum(k * lam[it, "RI"])
    GD[[length(GD) + 1]] <- data.table(comparison = p, set = s, n_reference = kg$n_reference, n_group = kg$n_group,
      mean_reference = kg$mean_reference, mean_group = kg$mean_group, d = kg$d, d_lower = kg$d_lower, d_upper = kg$d_upper,
      keyed_ri_loading_sum = kr, beta_S = bS, beta_RI = bR,
      d_content_part = ks * bS / length(it) / kg$pooled_sd, d_agreement_part = kr * bR / length(it) / kg$pooled_sd)
  }
  sf2b_log(p, " done")
}
gd <- rbindlist(GD)
gd[, agreement_part_share := d_agreement_part / (d_content_part + d_agreement_part)]

## Write ---------
num <- function(DT, kk = 4) { cols <- names(DT)[vapply(DT, is.numeric, TRUE)]; DT[, (cols) := lapply(.SD, function(x) round(x, kk)), .SDcols = cols]; DT }
fwrite(num(copy(shares_dt)), sf2b_out("25_variance_shares.csv"))
fwrite(num(copy(gd)), sf2b_out("25_group_differences.csv"))
out <- list(
  analysis = "25: the agreement share of the observed score, the 8 vs the 7 (no U4). EXPLORATORY, NOT REGISTERED (osf.io/d4a5u does not include it); specified 17 November 2026 after 10-24 were seen.",
  choices = CHOICES, sets = SETS, observed_fit = FITS_OBS,
  shares = num(copy(shares_dt)), group_differences = num(copy(gd)),
  minutes = round(as.numeric(difftime(Sys.time(), T0, units = "mins")), 1))
write_json(out, sf2b_out("25_exploratory_agreement_share.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
print(num(copy(shares_dt)))
print(num(copy(gd))[, .(comparison, set, n_reference, n_group, d, d_lower, d_upper, keyed_ri_loading_sum, d_content_part, d_agreement_part, agreement_part_share)])
sf2b_log("25 done in ", round(as.numeric(difftime(Sys.time(), T0, units = "mins")), 1), " min")
