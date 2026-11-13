# Study 2b power for the tests that had none or only one row, N = 700 (Norm, 13 November 2026).
#
#   H2   joint: the four neighbours generated together, H2 fitted the registered way (one model
#        per measure: the 8 with the agreement factor plus that measure's latent factor; 95% Wald
#        CI lower bound > .30). Per-correlation power, P(all four), P(at least three of four).
#   H3c  BIDR-16 Self-Deceptive Enhancement and Impression Management (8 items each, 7 categories,
#        two correlated factors in one model): both 95% Wald CIs within +/-.30.
#   H3e  BSSS-8 (8 items, 5 categories) upper bound below .40; NCS-6 (6 items, 5 categories, items
#        3 and 4 reversed), SPSQ-SF Sensory Discomfort (4 items, 7 categories) and ASI-3 Physical
#        Concerns (6 items, 5 categories) each within +/-.30. The four are generated together and
#        each is fitted in its own model, as registered, so the joint power is simulated, not assumed.
#   H4a  Spearman rho of the observed 8 with past-week pause days (0-7, skewed), two-sided, Holm
#        across H4a and H4b (H4b generated beside it).
#   H5b  proportional-odds regression of practice frequency (five levels) on the 8 beyond 18
#        covariate columns; likelihood-ratio p < .05 and a positive coefficient.
#
# Populations, as in power_2b_h3d.R: the 8 from Study 2a's fitted 13-statement model
# (study2b/analysis/sim_population.json: content and agreement loadings, thresholds); each validity
# measure a latent factor of its real items (study2b/analysis/config.R), with their keying and
# category counts, items loading .70 (negative where reverse-keyed), mildly skewed thresholds; the
# agreement factor uncorrelated with every validity factor. Only the measures a test needs are
# generated: the other subscales of FFMQ-15, BFI-2 and Brief MAIA-2, which share those models in
# the pipeline, are left out.
#
# Every SEM condition is written to its CSV as soon as it finishes, and a rerun skips the
# conditions already written, so a stopped run resumes where it was.
#
# Run from the repo root:  SF2B_CORES=8 Rscript study2/power/power_2b_review.R
# Writes study2/power/power_2b_review_h2.csv, _h3c.csv, _h3e.csv, _h4a.csv, _h5b.csv.

# Set Up ---------
## Load libraries ---------
packages <- c("jsonlite", "data.table", "lavaan", "MASS")
new_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(new_packages)) install.packages(new_packages)
for (thispack in packages) {
  library(thispack, character.only = TRUE, quietly = TRUE, verbose = FALSE)
}
library(parallel)

## Settings ---------
source("study2b/analysis/config.R")
source("study2b/analysis/sf2b_models.R")

N <- as.integer(Sys.getenv("SF2B_N", "700"))
REPS <- as.integer(Sys.getenv("SF2B_REPS", "300"))          # per SEM condition
REPS_FAST <- as.integer(Sys.getenv("SF2B_REPS_FAST", "4000"))   # per H4a condition
REPS_H5 <- as.integer(Sys.getenv("SF2B_REPS_H5", "2000"))       # per H5b condition
SEED <- 20261113
OUT_DIR <- Sys.getenv("SF2B_POWER_OUT", "study2/power")   # a scratch folder for a smoke test
OUT <- function(x) file.path(OUT_DIR, sprintf("power_2b_review_%s.csv", x))

D <- sf2b_definitions(); f8 <- D$sf8
pop <- fromJSON("study2b/analysis/sim_population.json")
ix8 <- match(f8, pop$items)
lamS <- pop$loading_S[ix8]; lamRI <- pop$loading_RI[ix8]; TH8 <- pop$thresholds[f8]
its <- sf2b_score_items(D); rv <- unlist(lapply(D$validity, `[[`, "reverse"))
TH_K <- list(`5` = qnorm(c(.10, .30, .60, .85)), `6` = qnorm(c(.05, .15, .35, .65, .85)),
             `7` = qnorm(c(.04, .10, .22, .42, .68, .88)))
# Category counts, from each measure's range in config.R (checked against study2b/materials).
ncat <- vapply(names(its), function(s) {
  m <- Filter(function(v) s %in% names(v$scores), D$validity)[[1]]; as.integer(diff(m$range) + 1) }, 0L)
stopifnot(ncat[["bidr_sde"]] == 7, ncat[["bidr_im"]] == 7, ncat[["bsss8"]] == 5, ncat[["ncs6"]] == 5,
          ncat[["spsq_sd"]] == 7, ncat[["asi_pc"]] == 5, ncat[["ffmq_observing"]] == 5, ncat[["sbi_moment"]] == 7,
          ncat[["bfi_aesthetic"]] == 5, ncat[["maia_noticing"]] == 6,
          lengths(its[c("bidr_sde", "bidr_im", "bsss8", "ncs6", "spsq_sd", "asi_pc")]) == c(8, 8, 8, 6, 4, 6))

# The item correlation matrix implied by S, RI and the validity factors `sc`: r = each factor's
# correlation with S, rho = their correlations with each other, lam = item loading.
pop_gen <- function(sc, r, rho, lam = .70) {
  k <- length(sc); stopifnot(length(r) == k, all(dim(rho) == k))
  items <- c(f8, unlist(its[sc], use.names = FALSE))
  fac <- c("S", "RI", sc); Phi <- diag(k + 2); dimnames(Phi) <- list(fac, fac)
  Phi["S", sc] <- Phi[sc, "S"] <- r; Phi[sc, sc] <- rho; diag(Phi) <- 1
  L <- matrix(0, length(items), k + 2, dimnames = list(items, fac))
  L[f8, "S"] <- lamS; L[f8, "RI"] <- lamRI
  for (s in sc) L[its[[s]], s] <- lam * ifelse(its[[s]] %in% rv, -1, 1)
  P <- L %*% Phi %*% t(L); diag(P) <- 1
  stopifnot(min(eigen(Phi, only.values = TRUE)$values) > 0, min(eigen(P, only.values = TRUE)$values) > 0)
  TH <- c(TH8, unlist(lapply(sc, function(s) rep(list(TH_K[[as.character(ncat[[s]])]]), length(its[[s]]))), recursive = FALSE))
  names(TH) <- items
  list(P = P, TH = TH, items = items)
}
sym <- function(k, vals) { m <- diag(k); m[upper.tri(m)] <- vals; m[lower.tri(m)] <- t(m)[lower.tri(m)]; m }
equi <- function(k, x) { m <- matrix(x, k, k); diag(m) <- 1; m }

cut_ord <- function(z, t) findInterval(z, t) + 1L
gen <- function(pp) {
  Z <- MASS::mvrnorm(N, rep(0, ncol(pp$P)), pp$P)
  d <- as.data.frame(vapply(seq_len(ncol(pp$P)), function(j) cut_ord(Z[, j], pp$TH[[j]]), integer(N)))
  names(d) <- pp$items; d
}
# One registered validity model (the pipeline's 04_validity.R, ri_free = FALSE) for the
# factors `sc`: each correlation with S, its 95% Wald CI. A model that fails or whose factors
# are not oriented as keyed (the pipeline refuses those) gives NA, which counts as not supported.
fit_cor <- function(d, sc) {
  items <- c(f8, unlist(its[sc], use.names = FALSE))
  f <- sf2b_try(sf2b_sem(c(sf2b_scale_syntax(f8), sf2b_validity_syntax(sc, D, ri_free = FALSE)), d[, items], ordered = items))
  na <- setNames(rep(NA_real_, 3 * length(sc)), c(paste0("est_", sc), paste0("lo_", sc), paste0("up_", sc)))
  if (sf2b_failed(f) || length(sf2b_orientation(f, f8, sc, D))) return(na)
  pe <- parameterEstimates(f); x <- pe[match(paste0("r_", sc), pe$label), ]
  setNames(c(x$est, x$ci.lower, x$ci.upper), names(na))
}
# One simulated study: every measure in `models` (a list of factor sets, one model each).
one <- function(pp, models) { d <- gen(pp); unlist(lapply(models, function(sc) fit_cor(d, sc))) }

pct <- function(x) round(100 * mean(x %in% TRUE), 1)
cl <- makeCluster(SF2B$cores)
clusterSetRNGStream(cl, SEED)
clusterExport(cl, c("N", "f8", "its", "D", "cut_ord", "gen", "fit_cor", "one", "sf2b_sem", "sf2b_try", "sf2b_failed",
                    "sf2b_orientation", "sf2b_scale_syntax", "sf2b_validity_syntax", "sf2b_score_items"))
invisible(clusterEvalQ(cl, { suppressPackageStartupMessages({ library(lavaan); library(MASS) }); NULL }))

# Runs every condition not yet in the CSV; `summ` turns the replications into one row.
run_conds <- function(name, COND, models, summ) {
  out <- OUT(name)
  rows <- if (file.exists(out)) split(fread(out), by = "condition") else list()
  for (nm in setdiff(names(COND), names(rows))) {
    t0 <- Sys.time(); c0 <- COND[[nm]]
    pp <- pop_gen(c0$sc, c0$r, c0$rho, c0$lam %||% .70)
    r <- do.call(rbind, parLapply(cl, seq_len(REPS), function(i, pp, models) one(pp, models), pp = pp, models = models))
    rows[[nm]] <- cbind(data.table(condition = nm), summ(r, c0))
    fwrite(rbindlist(rows, fill = TRUE)[order(match(condition, names(COND)))], out)
    print(rows[[nm]]); cat(sprintf("  %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  rbindlist(rows, fill = TRUE)[order(match(condition, names(COND)))]
}

# ---- H2, the four neighbours together -----------------------------------------------------------
H2 <- SF2B_H2   # ffmq_observing, sbi_moment, bfi_aesthetic, maia_noticing
H2_COND <- list()
for (r in c(.40, .45, .50)) H2_COND[[sprintf("all four r %.2f, neighbours correlate .30", r)]] <- list(sc = H2, r = rep(r, 4), rho = equi(4, .30))
H2_COND[["mixed: r .45 / .45 / .45 / .38 (Obs / SBI / Aes / Not), neighbours correlate .30"]] <- list(sc = H2, r = c(.45, .45, .45, .38), rho = equi(4, .30))
h2_summ <- function(r, c0) {
  lo <- r[, paste0("lo_", H2), drop = FALSE]; up <- r[, paste0("up_", H2), drop = FALSE]
  pass <- lo > SF2B$h2_lower; pass[is.na(pass)] <- FALSE
  data.table(r_obs = c0$r[1], r_sbi = c0$r[2], r_aes = c0$r[3], r_not = c0$r[4], rho = .30,
             converged = pct(rowSums(is.na(lo)) == 0),
             power_obs = pct(pass[, 1]), power_sbi = pct(pass[, 2]), power_aes = pct(pass[, 3]), power_not = pct(pass[, 4]),
             power_all4 = pct(rowSums(pass) == 4), power_atleast3 = pct(rowSums(pass) >= SF2B$h2_needed),
             h3a_all4 = pct(rowSums(up < SF2B$h3a_upper, na.rm = TRUE) == 4))
}
cat("\n==== H2\n")
res_h2 <- run_conds("h2", H2_COND, as.list(H2), h2_summ)

# ---- H3c, BIDR-16 ----------------------------------------------------------------------------------
# SDE and IM correlate .40 with each other (an assumption; their scale scores typically correlate
# .3-.5). Both factors in one model, as registered.
BIDR <- SF2B_H3C
H3C_COND <- list()
for (r in c(0, .10, .15, .20, .25)) H3C_COND[[sprintf("both r %.2f", r)]] <- list(sc = BIDR, r = c(r, r), rho = equi(2, .40))
H3C_COND[["predicted (near zero to + small): SDE .20, IM .10"]] <- list(sc = BIDR, r = c(.20, .10), rho = equi(2, .40))
H3C_COND[["sensitivity: both r .10, item loadings .55"]] <- list(sc = BIDR, r = c(.10, .10), rho = equi(2, .40), lam = .55)
within <- function(lo, up, b) { x <- lo > -b & up < b; x[is.na(x)] <- FALSE; x }
h3c_summ <- function(r, c0) {
  a <- within(r[, "lo_bidr_sde"], r[, "up_bidr_sde"], SF2B$h3c_bound); b <- within(r[, "lo_bidr_im"], r[, "up_bidr_im"], SF2B$h3c_bound)
  data.table(r_sde = c0$r[1], r_im = c0$r[2], rho_sde_im = .40, loading = c0$lam %||% .70, converged = pct(!is.na(r[, "lo_bidr_sde"])),
             power_sde = pct(a), power_im = pct(b), power_both = pct(a & b),
             median_halfwidth = round(median(c(r[, "up_bidr_sde"] - r[, "lo_bidr_sde"], r[, "up_bidr_im"] - r[, "lo_bidr_im"]), na.rm = TRUE) / 2, 3))
}
cat("\n==== H3c\n")
res_h3c <- run_conds("h3c", H3C_COND, list(BIDR), h3c_summ)

# ---- H3e, the four new discriminant measures --------------------------------------------------------------
# Their correlations with each other (an assumption, none measured together): BSSS-NCS .20,
# BSSS-SPSQ -.10, BSSS-ASI 0, NCS-SPSQ -.10, NCS-ASI -.15, SPSQ-ASI .40.
E <- SF2B_H3E$score   # bsss8, ncs6, spsq_sd, asi_pc
rhoE <- sym(4, c(.20, -.10, -.10, 0, -.15, .40))   # upper triangle by column: 12, 13, 23, 14, 24, 34
H3E_COND <- list()
for (r in c(0, .10, .15, .20, .25)) H3E_COND[[sprintf("all four r %.2f", r)]] <- list(sc = E, r = rep(r, 4), rho = rhoE)
H3E_COND[["BSSS .30, others .25"]] <- list(sc = E, r = c(.30, .25, .25, .25), rho = rhoE)
H3E_COND[["predicted (BSSS + small, others near zero): BSSS .20, others .10"]] <- list(sc = E, r = c(.20, .10, .10, .10), rho = rhoE)
H3E_COND[["predicted, at their edge: BSSS .30, others .15"]] <- list(sc = E, r = c(.30, .15, .15, .15), rho = rhoE)
H3E_COND[["sensitivity: BSSS .20, others .10, item loadings .55"]] <- list(sc = E, r = c(.20, .10, .10, .10), rho = rhoE, lam = .55)
h3e_summ <- function(r, c0) {
  ok <- sapply(seq_along(E), function(k) { b <- SF2B_H3E[k]; s <- E[k]
    x <- r[, paste0("lo_", s)] > b$lower & r[, paste0("up_", s)] < b$upper; x[is.na(x)] <- FALSE; x })
  data.table(r_bsss = c0$r[1], r_ncs = c0$r[2], r_spsq = c0$r[3], r_asi = c0$r[4], loading = c0$lam %||% .70,
             converged = pct(rowSums(is.na(r[, paste0("lo_", E)])) == 0),
             power_bsss = pct(ok[, 1]), power_ncs = pct(ok[, 2]), power_spsq = pct(ok[, 3]), power_asi = pct(ok[, 4]),
             power_all4 = pct(rowSums(ok) == 4))
}
cat("\n==== H3e\n")
res_h3e <- run_conds("h3e", H3E_COND, as.list(E), h3e_summ)

# ---- H4a, past-week pause days, with H4b beside it for Holm -------------------------------------------------
# As power_2b.R: the observed 8 has reliability .84; each behaviour is a latent normal correlating
# r with the 8's true score (and with each other only through it). Pause days are cut from their
# latent so that 0-7 days have probabilities .30, .15, .13, .11, .08, .07, .04, .12 (a skewed count
# with a pile at 0 and a smaller one at "every day"). H4b: the engaged pause, take-up 30%, latent r .20.
set.seed(SEED)
rel <- .84
P_DAYS <- c(.30, .15, .13, .11, .08, .07, .04, .12)
TH_DAYS <- qnorm(cumsum(P_DAYS)[-8])
h4 <- list()
for (rho_a in c(.15, .20, .25)) {
  sim <- replicate(REPS_FAST, { t <- rnorm(N); x <- sqrt(rel) * t + sqrt(1 - rel) * rnorm(N)
    days <- findInterval(rho_a * t + sqrt(1 - rho_a^2) * rnorm(N), TH_DAYS)
    pause <- as.numeric(.20 * t + sqrt(1 - .20^2) * rnorm(N) > qnorm(.70))
    a <- suppressWarnings(cor.test(x, days, method = "spearman", exact = FALSE))
    b <- suppressWarnings(cor.test(x, pause, method = "spearman", exact = FALSE))
    c(a$estimate, a$p.value, b$estimate, b$p.value) })
  holm <- apply(sim[c(2, 4), ], 2, p.adjust, method = "holm")
  h4[[length(h4) + 1]] <- data.table(latent_r = rho_a, n = N, power_p05 = pct(sim[1, ] > 0 & sim[2, ] < .05),
    power_holm = pct(sim[1, ] > 0 & holm[1, ] < .05), power_holm_worst = pct(sim[1, ] > 0 & sim[2, ] < .025),
    h4b_power_holm = pct(sim[3, ] > 0 & holm[2, ] < .05), median_rho = round(median(sim[1, ]), 3),
    days_distribution = paste(P_DAYS, collapse = "/"), h4b = "engaged pause, take-up 30%, latent r .20")
}
h4 <- rbindlist(h4); fwrite(h4, OUT("h4a")); cat("\n==== H4a\n"); print(h4)

# ---- H5b, practice frequency ----------------------------------------------------------------------------
# Five levels with 2a's practice groups (273 current, 136 never and 91 others of 500; the current
# split .45 / .30 / .25 as in 00_simulate.R, the others split evenly between "none in the past
# three months" and "less than weekly"): none .36, less than weekly .09, 1-2 a week .25,
# 3-5 a week .16, daily .14.
# 18 covariate columns (age, gender as two dummies, spirituality, religiosity, five FFMQ-15 facets,
# eight Brief MAIA-2 subscales), simulated as independent standard normals. C is their equally
# weighted composite. The 8 = sqrt(R2x) C + sqrt(1 - R2x) e, so the covariates explain R2x of it
# (.40, and .25 as a check). The outcome's latent logit = 0.8 C + log(OR) x the 8 (per SD) +
# logistic error: the covariates' odds ratio is 2.2 per SD of their composite. Cut at the latent
# quantiles that give the five proportions.
P_FREQ <- c(.36, .09, .25, .16, .14)
h5_one <- function(n, b, r2x, th) {
  X <- matrix(rnorm(n * 18), n); C <- rowSums(X) / sqrt(18)
  x <- sqrt(r2x) * C + sqrt(1 - r2x) * rnorm(n)
  y <- factor(findInterval(.8 * C + b * x + rlogis(n), th), levels = 0:4, ordered = TRUE)
  x <- (x - mean(x)) / sd(x)
  m1 <- tryCatch(suppressWarnings(MASS::polr(y ~ X)), error = function(e) NULL)
  m2 <- tryCatch(suppressWarnings(MASS::polr(y ~ X + x)), error = function(e) NULL)
  if (is.null(m1) || is.null(m2)) return(c(ok = 0, coef = NA, p = NA))
  lr <- 2 * (as.numeric(logLik(m2)) - as.numeric(logLik(m1)))
  c(ok = 1, coef = unname(coef(m2)["x"]), p = pchisq(lr, 1, lower.tail = FALSE))
}
clusterExport(cl, "h5_one")
h5 <- list()
for (r2x in c(.40, .25)) for (or in c(1.2, 1.3, 1.4)) {
  b <- log(or)
  big <- 1e6; Cb <- rnorm(big); yb <- .8 * Cb + b * (sqrt(r2x) * Cb + sqrt(1 - r2x) * rnorm(big)) + rlogis(big)
  th <- unname(quantile(yb, cumsum(P_FREQ)[-5]))
  for (n in c(410, 470, 650)) {
    s <- do.call(rbind, parLapply(cl, seq_len(REPS_H5), function(i, n, b, r2x, th) h5_one(n, b, r2x, th), n = n, b = b, r2x = r2x, th = th))
    h5[[length(h5) + 1]] <- data.table(n = n, odds_ratio_per_sd = or, r2_with_covariates = r2x, converged = pct(s[, "ok"] == 1),
      power = pct(s[, "coef"] > 0 & s[, "p"] < .05), median_or = round(exp(median(s[, "coef"], na.rm = TRUE)), 3),
      levels = paste(P_FREQ, collapse = "/"))
    print(h5[[length(h5)]])
  }
}
h5 <- rbindlist(h5); fwrite(h5, OUT("h5b"))
stopCluster(cl)
cat(sprintf("\nN = %d, %d replications per SEM condition, %d per H4a condition, %d per H5b condition, seed %d\n",
            N, REPS, REPS_FAST, REPS_H5, SEED))
