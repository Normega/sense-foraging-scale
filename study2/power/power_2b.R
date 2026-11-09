# Study 2b power simulation, N = 500 (Norm's decisions of 9 November 2026: one posting of 500;
# the smallest H3b contrast is a difference of .15 between two latent correlations; H8 exploratory).
#
# Populations are built from Study 2a's fitted one-factor model of the 8 (with the random-intercept
# agreement factor), so the scale's loadings, agreement variance and response thresholds are 2a's.
# Only the population is read from 2a's derived data; every output here is a simulation summary.
#
#   H1  the registered fit rule (decision 17) on three populations:
#         fit       2a's model-implied correlations (the model is true)
#         observed  2a's observed polychoric correlations (2b looks exactly like 2a)
#         flawed    two planted residual correlations of .30 (the model is wrong)
#   H2  latent r of the 8 with a 3-indicator neighbour, 95% CI lower bound > .30
#   H3b difference of two dependent latent correlations, 95% CI excluding 0 (delta method)
#   H4  the single pause (0/1) against the 8's score, Spearman, two-sided .05
#   H6  current practitioners vs never-practised, d >= .50 and p < .05 (2a's group sizes)
#   H5  incremental delta R^2 >= .02 and p < .05
#
# Run from the repo root:  SF2A_MODE=real Rscript study2/power/power_2b.R
# Writes study2/power/power_2b_results.csv (one row per condition) and prints the same.

source("study2a/analysis/config.R")
source("study2a/analysis/sf2a_models.R")
stopifnot(SF2A$mode == "real")
library(parallel)

N <- 500
REPS <- as.integer(Sys.getenv("SF2B_REPS", "500"))
SEED <- 20261109
OUT <- "study2/power/power_2b_results.csv"

A <- readRDS(sf2a_derived("sf2a")); def <- A$definition; X <- A$items_raw; items <- def$items; rv <- def$reversed
id <- function(src) items[match(src, def$source[items])]
f8 <- id(c("N4", "sf6_36", "V10", "sf6_27", "sfc_36", "sfc_x2", "sfc_35", "U4"))
keyed <- ifelse(f8 %in% rv, -1, 1)

syn8 <- c(paste0("S =~ NA*", paste(f8, collapse = " + ")), "S ~~ 1*S",
          paste0("RI =~ ", paste0("1*", f8, collapse = " + ")), "RI ~~ RI", "S ~~ 0*RI")
fit2a <- sf2a_sem(syn8, X[, f8], ordered = f8)
std <- lavInspect(fit2a, "std")
lamS <- std$lambda[f8, "S"]; lamRI <- std$lambda[f8, "RI"]
P_fit <- cov2cor(lavInspect(fit2a, "implied")$cov[f8, f8])
P_obs <- lavCor(as.data.frame(X[, f8]), ordered = f8, missing = "pairwise")[f8, f8]
D <- sqrt(diag(lavInspect(fit2a, "implied")$cov[f8, f8]))
th <- lavInspect(fit2a, "th")
TH <- lapply(f8, function(v) sort(th[startsWith(names(th), paste0(v, "|"))] / D[v]))
res_var <- 1 - lamS^2 - lamRI^2
P_flaw <- P_fit
for (p in list(c(1, 2), c(3, 4))) {   # two non-overlapping pairs, residual correlation .30
  i <- p[1]; j <- p[2]; P_flaw[i, j] <- P_flaw[j, i] <- P_flaw[i, j] + .30 * sqrt(res_var[i] * res_var[j])
}
# a consequential flaw: three statements (N4, sf6_36, V10) share a second factor loading .40
m <- c(rep(.40, 3), rep(0, 5)); P_dim <- cov2cor(P_fit + outer(m, m))
m2 <- c(rep(.60, 3), rep(0, 5)); P_dim2 <- cov2cor(P_fit + outer(m2, m2))   # a larger one, loading .60
stopifnot(min(eigen(P_flaw)$values) > 0, min(eigen(P_obs)$values) > 0, min(eigen(P_dim)$values) > 0, min(eigen(P_dim2)$values) > 0)
PART <- Sys.getenv("SF2B_PART", "all")
cat(sprintf("2a population: loadings %.2f-%.2f, agreement loading %.2f\n", min(abs(lamS)), max(abs(lamS)), mean(lamRI)))

# ---- shared helpers (exported to workers) ----------------------------------------------------
cut_ord <- function(z, t) findInterval(z, t) + 1L
TH5 <- qnorm(c(.10, .30, .60, .85))   # neighbour items: 5 categories, mildly skewed
gen_items <- function(n, P, TH, names) {
  Z <- MASS::mvrnorm(n, rep(0, ncol(P)), P)
  Y <- vapply(seq_len(ncol(P)), function(j) cut_ord(Z[, j], TH[[j]]), integer(n))
  colnames(Y) <- names; as.data.frame(Y)
}
omega_keyed <- function(fit) {
  s <- lavInspect(fit, "std"); l <- s$lambda[f8, "S"] * keyed; r <- s$lambda[f8, "RI"] * keyed
  Th <- s$theta[f8, f8] * outer(keyed, keyed)
  num <- sum(l)^2; num / (num + sum(r)^2 + sum(Th))
}
h1_one <- function(P) {
  d <- gen_items(N, P, TH, f8)
  f <- tryCatch(sf2a_sem(syn8, d, ordered = f8), error = function(e) NULL)
  na <- c(ok = 0, cfi = NA, rmsea = NA, maxres = NA, usrmr = NA, closefit_p = NA, d_omega = NA, d_lambda = NA, omega = NA)
  if (is.null(f) || !lavInspect(f, "converged")) return(na)
  fm <- fitMeasures(f, c("cfi.scaled", "rmsea.scaled"))
  rr <- lavResiduals(f, type = "cor.bentler", summary = TRUE)
  R <- rr$cov; R[upper.tri(R, TRUE)] <- NA
  sm <- rr$summary
  k <- which.max(abs(R)); ij <- arrayInd(k, dim(R))
  f2 <- tryCatch(sf2a_sem(c(syn8, paste(f8[ij[1]], "~~", f8[ij[2]])), d, ordered = f8), error = function(e) NULL)
  om <- omega_keyed(f); dom <- NA; dlam <- NA
  if (!is.null(f2) && lavInspect(f2, "converged")) {
    dom <- abs(omega_keyed(f2) - om)
    dlam <- max(abs(lavInspect(f2, "std")$lambda[f8, "S"] - lavInspect(f, "std")$lambda[f8, "S"]))
  }
  c(ok = 1, cfi = unname(fm[1]), rmsea = unname(fm[2]), maxres = max(abs(R), na.rm = TRUE), usrmr = sm["usrmr", "cor"],
    closefit_p = sm["usrmr.closefit.pvalue", "cor"], d_omega = dom, d_lambda = dlam, omega = om)
}
# rule variants, scored from the raw quantities
rule_pass <- function(r, res_max = .10, d_om = .02, d_lam = .02, cfi_min = .95) {
  local <- r[, "maxres"] <= res_max | (!is.na(r[, "d_omega"]) & r[, "d_omega"] < d_om & r[, "d_lambda"] < d_lam)
  r[, "cfi"] >= cfi_min & local & r[, "closefit_p"] > .05
}

# The 8 plus two neighbours A and B, each with k indicators loading lam: latent r_SA, r_SB, r_AB.
# RI is uncorrelated with both.
pop_three <- function(rSA, rSB, rAB, k = 3, lam = .70) {
  p8 <- length(f8); nn <- 2 * k; P <- diag(p8 + nn)
  P[1:p8, 1:p8] <- P_fit
  la <- c(rep(lam, k), rep(0, k)); lb <- c(rep(0, k), rep(lam, k))
  Pn <- outer(la, la) + outer(lb, lb) + rAB * (outer(la, lb) + outer(lb, la)); diag(Pn) <- 1
  P[p8 + 1:nn, p8 + 1:nn] <- Pn
  C <- outer(lamS, rSA * la + rSB * lb)
  P[1:p8, p8 + 1:nn] <- C; P[p8 + 1:nn, 1:p8] <- t(C)
  stopifnot(min(eigen(P)$values) > 0)
  list(P = P, k = k)
}
cor_one <- function(Pk) {
  k <- Pk$k; nb <- c(paste0("a", 1:k), paste0("b", 1:k))
  syn <- c(paste0("S =~ NA*", paste(f8, collapse = " + ")), "S ~~ 1*S",
           paste0("RI =~ ", paste0("1*", f8, collapse = " + ")), "RI ~~ RI", "S ~~ 0*RI",
           paste0("A =~ NA*", paste(nb[1:k], collapse = " + ")), "A ~~ 1*A",
           paste0("B =~ NA*", paste(nb[k + 1:k], collapse = " + ")), "B ~~ 1*B",
           "RI ~~ 0*A", "RI ~~ 0*B", "S ~~ rsa*A", "S ~~ rsb*B", "A ~~ B", "dif := rsa - rsb")
  d <- gen_items(N, Pk$P, c(TH, rep(list(TH5), 2 * k)), c(f8, nb))
  f <- tryCatch(sf2a_sem(syn, d, ordered = c(f8, nb)), error = function(e) NULL)
  if (is.null(f) || !lavInspect(f, "converged")) return(c(ok = 0, rsa = NA, rsa_lo = NA, rsa_hi = NA, dif = NA, dif_lo = NA))
  pe <- parameterEstimates(f)
  a <- pe[pe$label == "rsa", ]; df <- pe[pe$label == "dif" & pe$op == ":=", ]
  c(ok = 1, rsa = a$est, rsa_lo = a$ci.lower, rsa_hi = a$ci.upper, dif = df$est, dif_lo = df$ci.lower)
}

# ---- run ------------------------------------------------------------------------------------
cl <- makeCluster(min(20, detectCores() - 2))
on.exit(stopCluster(cl), add = TRUE)
clusterSetRNGStream(cl, SEED)
clusterExport(cl, c("N", "f8", "keyed", "syn8", "TH", "TH5", "cut_ord", "gen_items", "omega_keyed",
                    "h1_one", "cor_one", "sf2a_sem", "SF2A"))
invisible(clusterEvalQ(cl, { suppressPackageStartupMessages({ library(lavaan); library(MASS) }); NULL }))
runs <- function(fn, P) do.call(rbind, parLapply(cl, seq_len(REPS), function(i, g, Pm) g(Pm), g = fn, Pm = P))
pct <- function(x) round(100 * mean(x, na.rm = TRUE), 1)
rows <- list()

VARIANTS <- list(`registered (res .10, dOmega/dLambda .02)` = list(),
                 `dLambda .05` = list(d_lam = .05),
                 `res .15` = list(res_max = .15),
                 `dOmega .02 only` = list(d_lam = Inf))
H1RAW <- list()
for (pop in c("fit", "observed", "flawed", "second factor", "second factor .60")) {
  P <- switch(pop, fit = P_fit, observed = P_obs, flawed = P_flaw, `second factor` = P_dim, `second factor .60` = P_dim2)
  r <- runs(h1_one, P); H1RAW[[pop]] <- r
  for (v in names(VARIANTS)) {
    ps <- do.call(rule_pass, c(list(r = r), VARIANTS[[v]]))
    rows[[length(rows) + 1]] <- data.table(test = paste("H1 fit rule:", v), condition = pop, converged = pct(r[, "ok"]),
      power_or_pass = pct(ps), detail = sprintf("CFI>=.95 %.1f%%; max residual <=.10 %.1f%%; close fit not rejected %.1f%%; medians: CFI %.3f, RMSEA %.3f, max residual %.3f, USRMR %.3f, dOmega %.3f, dLambda %.3f",
        pct(r[, "cfi"] >= .95), pct(r[, "maxres"] <= .10), pct(r[, "closefit_p"] > .05), median(r[, "cfi"], na.rm = TRUE), median(r[, "rmsea"], na.rm = TRUE),
        median(r[, "maxres"], na.rm = TRUE), median(r[, "usrmr"], na.rm = TRUE), median(r[, "d_omega"], na.rm = TRUE), median(r[, "d_lambda"], na.rm = TRUE)))
  }
  print(rbindlist(tail(rows, length(VARIANTS)))[, .(test, condition, power_or_pass)])
}
saveRDS(H1RAW, "study2/power/power_2b_h1_raw.rds")
GRID <- CJ(res_max = c(.10, .12, .15), d_om = c(.02, .03), d_lam = c(.02, .05, .10, Inf))
grid <- GRID[, as.list(sapply(H1RAW, function(r) pct(rule_pass(r, res_max, d_om, d_lam)))), by = .(res_max, d_om, d_lam)]
fwrite(grid, "study2/power/power_2b_h1_grid.csv"); print(grid)
if (PART == "h1") { fwrite(rbindlist(rows), sub(".csv", "_h1.csv", OUT, fixed = TRUE)); quit(save = "no") }

for (rsa in c(.40, .45, .50)) {
  r <- runs(cor_one, pop_three(rsa, .20, .30, k = 3))
  rows[[length(rows) + 1]] <- data.table(test = "H2 lower bound > .30", condition = sprintf("true r = %.2f, 3-item neighbour", rsa), converged = pct(r[, "ok"]),
    power_or_pass = pct(r[, "rsa_lo"] > .30), detail = sprintf("median CI width %.3f", median(r[, "rsa_hi"] - r[, "rsa_lo"], na.rm = TRUE)))
  print(rows[[length(rows)]])
}

CONTRASTS <- list(c(.45, .30, 0, 3, .70), c(.45, .30, .30, 3, .70), c(.45, .30, .50, 3, .70), c(.35, .20, .30, 3, .70),
                  c(.45, .30, .30, 4, .75), c(.45, .30, .30, 8, .70), c(.40, .30, .50, 3, .70), c(.50, .30, .30, 3, .70))
for (cc in CONTRASTS) {
  r <- runs(cor_one, pop_three(cc[1], cc[2], cc[3], k = cc[4], lam = cc[5]))
  rows[[length(rows) + 1]] <- data.table(test = "H3b contrast", condition = sprintf("r %.2f vs %.2f (difference %.2f), measures correlate %.2f, %d items per measure loading %.2f",
      cc[1], cc[2], cc[1] - cc[2], cc[3], cc[4], cc[5]),
    converged = pct(r[, "ok"]), power_or_pass = pct(r[, "dif_lo"] > 0), detail = sprintf("median estimate %.3f", median(r[, "dif"], na.rm = TRUE)))
  print(rows[[length(rows)]][, .(condition, power_or_pass)])
}

# H4, H5, H6 on observed scores: fast, no SEM. The 8's score has 2a's reliability.
set.seed(SEED)
rel <- .84
score_from <- function(t) sqrt(rel) * t + sqrt(1 - rel) * rnorm(length(t))
for (take in c(.15, .30, .50)) for (rho in c(.15, .20, .25)) {
  p <- replicate(4000, { t <- rnorm(N); z <- rho * t + sqrt(1 - rho^2) * rnorm(N)
    y <- as.numeric(z > qnorm(1 - take)); suppressWarnings(cor.test(score_from(t), y, method = "spearman", exact = FALSE)$p.value) })
  rows[[length(rows) + 1]] <- data.table(test = "H4 pause (0/1)", condition = sprintf("take-up %.0f%%, latent r %.2f", 100 * take, rho),
    converged = 100, power_or_pass = pct(p < .05), detail = "Spearman, two-sided .05, score reliability .84")
}
n1 <- round(N * 273 / 500); n0 <- round(N * 136 / 500)
for (d in c(.55, .67, .80)) {
  ok <- replicate(4000, { a <- rnorm(n1, d); b <- rnorm(n0)
    dh <- (mean(a) - mean(b)) / sqrt((var(a) + var(b)) / 2); dh >= .50 && t.test(a, b)$p.value < .05 })
  rows[[length(rows) + 1]] <- data.table(test = "H6 known groups", condition = sprintf("true d = %.2f (%d current vs %d never)", d, n1, n0),
    converged = 100, power_or_pass = pct(ok), detail = "d >= .50 and p < .05")
}
for (dr2 in c(.02, .03, .04)) {
  ok <- replicate(4000, { R2b <- .30; x <- matrix(rnorm(N * 2), N)   # two controls explaining 30%
    sc <- rnorm(N); b <- sqrt(dr2); bc <- sqrt(R2b / 2)
    y <- bc * x[, 1] + bc * x[, 2] + b * sc + sqrt(1 - R2b - dr2) * rnorm(N)
    m0 <- lm(y ~ x); m1 <- lm(y ~ x + sc); a <- anova(m0, m1)
    (summary(m1)$r.squared - summary(m0)$r.squared) >= .02 && a$`Pr(>F)`[2] < .05 })
  rows[[length(rows) + 1]] <- data.table(test = "H5 incremental", condition = sprintf("true delta R2 = %.2f", dr2),
    converged = 100, power_or_pass = pct(ok), detail = "delta R2 >= .02 and p < .05; controls explain 30%")
}

res <- rbindlist(rows)
fwrite(res, OUT)
print(res[, .(test, condition, power_or_pass)])
cat(sprintf("\nN = %d, %d replications per SEM condition, seed %d\n", N, REPS, SEED))
