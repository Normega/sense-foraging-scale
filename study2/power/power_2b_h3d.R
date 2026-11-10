# Study 2b power simulation for H3d, N = 700 (Norm, 10 November 2026: "Can you simulate power please").
#
# H3d: the multiple correlation R = sqrt(r' Sigma^-1 r) of the content factor with the four H2
# neighbours (FFMQ-15 Observing, SBI Savoring the Moment, BFI-2 Aesthetic Sensitivity, Brief
# MAIA-2 Noticing), from one model holding all four; supported when the 95% CI's upper bound is
# below .80. The fitted model is the pipeline's own (04_validity.R, task "h3d"): the 8 with the
# agreement factor, each neighbour a latent factor of its real items with their real keying and
# category counts (3, 8, 4 and 3 items).
#
# Population: the 8 from Study 2a's fitted 13-statement model (study2b/analysis/sim_population.json:
# content and agreement loadings, thresholds); neighbour items load .70 (negative where reverse-keyed),
# mildly skewed thresholds. What varies is what nobody has measured: the 8's latent correlation with
# each neighbour (r) and the neighbours' latent correlations with each other (rho). "Through the 8
# only" sets rho_jk = r_j r_k: the neighbours share nothing the 8 does not carry (the dry run's truth).
#
# The registered interval is a 2,000-resample percentile bootstrap, about 10 minutes per data set,
# too slow for 300 data sets per condition. Every replication therefore gets a delta-method interval
# (on R and on atanh R), and CALIB replications in each of two conditions also get the percentile
# bootstrap (BOOT resamples), so the approximation is checked against the registered interval.
#
# H5a at the complete-case n the 80% rule leaves (decided 10 November: the rule is kept) is
# simulated at the end; it needs no SEM.
#
# Run from the repo root:  Rscript study2/power/power_2b_h3d.R
# Writes study2/power/power_2b_h3d_results.csv, power_2b_h3d_calibration.csv, power_2b_h5_n.csv.

source("study2b/analysis/config.R")
source("study2b/analysis/sf2b_models.R")
library(parallel)

N <- as.integer(Sys.getenv("SF2B_N", "700"))
REPS <- as.integer(Sys.getenv("SF2B_REPS", "300"))
CALIB <- as.integer(Sys.getenv("SF2B_CALIB", "8"))
BOOT <- as.integer(Sys.getenv("SF2B_CALIB_BOOT", "300"))
SEED <- 20261110
BOUND <- SF2B$h3d_upper

D <- sf2b_definitions(); f8 <- D$sf8
pop <- fromJSON("study2b/analysis/sim_population.json")
ix8 <- match(f8, pop$items)
lamS <- pop$loading_S[ix8]; lamRI <- pop$loading_RI[ix8]; TH8 <- pop$thresholds[f8]
its <- sf2b_score_items(D); rv <- unlist(lapply(D$validity, `[[`, "reverse"))
H2 <- SF2B_H2
nb_items <- its[H2]
TH_K <- list(`5` = qnorm(c(.10, .30, .60, .85)), `6` = qnorm(c(.05, .15, .35, .65, .85)),
             `7` = qnorm(c(.04, .10, .22, .42, .68, .88)))
ncat <- c(ffmq_observing = 5, sbi_moment = 7, bfi_aesthetic = 5, maia_noticing = 6)
all_items <- c(f8, unlist(nb_items, use.names = FALSE))
TH <- c(TH8, unlist(lapply(H2, function(s) rep(list(TH_K[[as.character(ncat[[s]])]]), length(nb_items[[s]]))), recursive = FALSE))
names(TH) <- all_items

# The item correlation matrix implied by factors S, RI and the four neighbours.
pop_h3d <- function(r, rho, lam = .70) {
  stopifnot(length(r) == 4, all(dim(rho) == 4))
  fac <- c("S", "RI", H2); Phi <- diag(6); dimnames(Phi) <- list(fac, fac)
  Phi["S", H2] <- Phi[H2, "S"] <- r; Phi[H2, H2] <- rho; diag(Phi) <- 1
  L <- matrix(0, length(all_items), 6, dimnames = list(all_items, fac))
  L[f8, "S"] <- lamS; L[f8, "RI"] <- lamRI
  for (s in H2) L[nb_items[[s]], s] <- lam * ifelse(nb_items[[s]] %in% rv, -1, 1)
  P <- L %*% Phi %*% t(L); diag(P) <- 1
  stopifnot(min(eigen(Phi, only.values = TRUE)$values) > 0, min(eigen(P, only.values = TRUE)$values) > 0)
  Rtrue <- sqrt(drop(t(r) %*% solve(rho, r)))
  list(P = P, R = Rtrue)
}
equi <- function(x) { m <- matrix(x, 4, 4); diag(m) <- 1; m }

# The pipeline's H3d model, with the neighbours' (already free) covariances labelled so the
# delta method can find them.
pairs <- combn(4, 2)
cov_lab <- sprintf("s%d%d", pairs[1, ], pairs[2, ])
syn <- c(sf2b_scale_syntax(f8), sf2b_validity_syntax(H2, D, ri_free = FALSE),
         sprintf("%s ~~ %s*%s", H2[pairs[1, ]], cov_lab, H2[pairs[2, ]]))
r_lab <- paste0("r_", H2)
mult_R <- function(th) {
  r <- th[1:4]; S <- diag(4); S[t(pairs)] <- th[5:10]; S[t(pairs[2:1, ])] <- th[5:10]
  v <- tryCatch(drop(t(r) %*% solve(S, r)), error = function(e) NA_real_)
  if (is.na(v) || v < 0) NA_real_ else sqrt(v)
}
cut_ord <- function(z, t) findInterval(z, t) + 1L
gen <- function(P) {
  Z <- MASS::mvrnorm(N, rep(0, ncol(P)), P)
  as.data.frame(vapply(seq_len(ncol(P)), function(j) cut_ord(Z[, j], TH[[j]]), integer(N)))
}
fix_names <- function(d) { names(d) <- all_items; d }

one <- function(P, boot = 0L) {
  d <- fix_names(gen(P))
  f <- sf2b_try(sf2b_sem(syn, d, ordered = all_items))
  na <- c(ok = 0, R = NA, se = NA, up_delta = NA, up_atanh = NA, up_boot = NA, lo_boot = NA, boot_ok = NA)
  if (sf2b_failed(f)) return(na)
  pe <- parameterEstimates(f)
  th <- setNames(pe$est[match(c(r_lab, cov_lab), pe$label)], c(r_lab, cov_lab))
  V <- lavInspect(f, "vcov")[c(r_lab, cov_lab), c(r_lab, cov_lab)]
  R <- mult_R(th)
  g <- vapply(seq_along(th), function(k) { h <- 1e-5; e <- replace(numeric(length(th)), k, h)
    (mult_R(th + e) - mult_R(th - e)) / (2 * h) }, 0)
  se <- sqrt(drop(t(g) %*% V %*% g))
  out <- c(ok = 1, R = R, se = se, up_delta = R + qnorm(.975) * se,
           up_atanh = tanh(atanh(R) + qnorm(.975) * se / (1 - R^2)), up_boot = NA, lo_boot = NA, boot_ok = NA)
  if (boot > 0) {
    b <- vapply(seq_len(boot), function(i) {
      fb <- sf2b_try(sf2b_sem(syn, d[sample.int(N, N, TRUE), , drop = FALSE], ordered = all_items, fast = TRUE))
      if (sf2b_failed(fb)) return(NA_real_)
      C <- lavInspect(fb, "cor.lv"); r <- C["S", H2]; S <- C[H2, H2]
      v <- tryCatch(drop(t(r) %*% solve(S, r)), error = function(e) NA_real_)
      if (is.na(v) || v < 0) NA_real_ else sqrt(v)
    }, 0)
    q <- quantile(b, c(.025, .975), na.rm = TRUE)
    out[c("lo_boot", "up_boot", "boot_ok")] <- c(q[1], q[2], sum(!is.na(b)))
  }
  out
}

# ---- conditions ---------------------------------------------------------------------------------
COND <- list()
for (r in c(.40, .45, .50)) {
  COND[[sprintf("r %.2f, rho through the 8 only (%.2f)", r, r^2)]] <- list(r = rep(r, 4), rho = equi(r^2))
  for (rho in c(.30, .40, .50)) COND[[sprintf("r %.2f, rho %.2f", r, rho)]] <- list(r = rep(r, 4), rho = equi(rho))
}
# One unequal pattern, an assumption (no study has these four together): Observing and Noticing
# closest to the 8 and to each other.
mix_rho <- equi(0); mix_rho[t(pairs)] <- c(.35, .40, .55, .35, .30, .30)   # Obs-SBI, Obs-Aes, Obs-Not, SBI-Aes, SBI-Not, Aes-Not
mix_rho[lower.tri(mix_rho)] <- t(mix_rho)[lower.tri(mix_rho)]
COND[["unequal: r .50 / .40 / .45 / .50 (Obs / SBI / Aes / Not), rho .30-.55"]] <- list(r = c(.50, .40, .45, .50), rho = mix_rho)
POPS <- lapply(COND, function(c) pop_h3d(c$r, c$rho))
cat(sprintf("%-70s R = %.3f\n", names(POPS), vapply(POPS, `[[`, 0, "R")), sep = "")
CALIB_CONDS <- c("r 0.45, rho through the 8 only (0.20)", "r 0.45, rho 0.30")   # R .74 and .67: near the bound

cl <- makeCluster(SF2B$cores)

clusterSetRNGStream(cl, SEED)
clusterExport(cl, c("N", "syn", "all_items", "TH", "H2", "r_lab", "cov_lab", "pairs", "mult_R", "cut_ord", "gen",
                    "fix_names", "one", "sf2b_sem", "sf2b_try", "sf2b_failed"))
invisible(clusterEvalQ(cl, { suppressPackageStartupMessages({ library(lavaan); library(MASS) }); NULL }))
pct <- function(x) round(100 * mean(x, na.rm = TRUE), 1)

# Results are written after every condition, and a rerun skips the conditions already written,
# so a run that is stopped part-way resumes where it was.
OUT <- "study2/power/power_2b_h3d_results.csv"
rows <- if (file.exists(OUT)) split(fread(OUT), by = "condition") else list(); RAW <- list()
for (nm in setdiff(names(POPS), names(rows))) {
  t0 <- Sys.time()
  r <- do.call(rbind, parLapply(cl, seq_len(REPS), function(i, P) one(P), P = POPS[[nm]]$P))
  RAW[[nm]] <- r
  rows[[nm]] <- data.table(condition = nm, true_R = round(POPS[[nm]]$R, 3), converged = pct(r[, "ok"]),
                           power_delta = pct(r[, "up_delta"] < BOUND), power_atanh = pct(r[, "up_atanh"] < BOUND),
                           median_R = round(median(r[, "R"], na.rm = TRUE), 3), median_upper = round(median(r[, "up_delta"], na.rm = TRUE), 3))
  fwrite(rbindlist(rows)[order(match(condition, names(POPS)))], OUT)
  print(rows[[nm]]); cat(sprintf("  %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
res <- rbindlist(rows)[order(match(condition, names(POPS)))]

# ---- calibration: the delta-method upper bound against the registered percentile bootstrap ------
cal <- list()
for (nm in CALIB_CONDS) {
  t0 <- Sys.time()
  r <- do.call(rbind, parLapply(cl, seq_len(CALIB), function(i, P, B) one(P, B), P = POPS[[nm]]$P, B = BOOT))
  cal[[nm]] <- cbind(data.table(condition = nm, true_R = round(POPS[[nm]]$R, 3)), as.data.table(r))
  cat(sprintf("calibration %s: %.1f min\n", nm, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
cal <- rbindlist(cal)
fwrite(cal[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 4) else x)], "study2/power/power_2b_h3d_calibration.csv")
print(cal[, .(n = .N, mean_up_boot = mean(up_boot), mean_up_delta = mean(up_delta), mean_up_atanh = mean(up_atanh),
              mean_abs_diff_delta = mean(abs(up_delta - up_boot)), mean_abs_diff_atanh = mean(abs(up_atanh - up_boot)),
              same_decision_delta = pct((up_delta < BOUND) == (up_boot < BOUND)),
              same_decision_atanh = pct((up_atanh < BOUND) == (up_boot < BOUND))), by = condition])

# ---- H5a at the complete-case n the 80% rule leaves -------------------------------------------
# As in power_2b.R (controls explain 30%, true delta R^2 .03), but with the registered step 1: age,
# gender (two dummies), the five FFMQ-15 facets and the eight Brief MAIA-2 subscales, 16 columns.
# "Holm, worst case": only one of the three outcomes has an effect, so it must reach .05 / 3.
set.seed(SEED)
h5 <- list()
for (n in c(410, 470, 550, 700)) for (dr2 in c(.03, .04)) {
  sim <- replicate(4000, { x <- matrix(rnorm(n * 16), n); sc <- rnorm(n); bc <- sqrt(.30 / 16)
    y <- drop(x %*% rep(bc, 16)) + sqrt(dr2) * sc + sqrt(1 - .30 - dr2) * rnorm(n)
    m0 <- lm(y ~ x); m1 <- lm(y ~ x + sc)
    c(summary(m1)$r.squared - summary(m0)$r.squared, anova(m0, m1)$`Pr(>F)`[2]) })
  h5[[length(h5) + 1]] <- data.table(n = n, true_delta_R2 = dr2, power_p05 = pct(sim[1, ] >= .02 & sim[2, ] < .05),
                                     power_holm_worst = pct(sim[1, ] >= .02 & sim[2, ] < .05 / 3))
}
h5 <- rbindlist(h5); fwrite(h5, "study2/power/power_2b_h5_n.csv"); print(h5)
stopCluster(cl)
cat(sprintf("\nN = %d, %d replications per condition, calibration %d x %d resamples, seed %d\n", N, REPS, CALIB, BOOT, SEED))
