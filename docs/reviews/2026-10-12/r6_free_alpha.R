# Review check 6: can the page's hard-coded "agreement-free" alphas (.73 .41 .44 .57) be reproduced
# from siren's residual matrix? And omega_h without auto-flipping for the 15- and 10-item totals.
suppressPackageStartupMessages({library(siren); library(psych)})
A <- readRDS("F:/radlab_data/sf_study1/derived/sf1_real.rds")
s3 <- acquiextraction(as.data.frame(A$items_raw), n_factors = 4, corr = "Polychoric", raw_data = TRUE)
Res <- s3$Res; rev <- A$definition$reversed
S <- diag(ifelse(colnames(Res) %in% rev, -1, 1)); Rk <- S %*% Res %*% S; dimnames(Rk) <- dimnames(Res)
SF <- list(Access = c(35, 36, 9, 8, 7), Und = c(22, 3, 15, 31, 37), Val = c(12, 18, 11, 27, 44), LB = c(40, 46, 43, 47, 21))
ids <- function(v) paste0("sf6_", v)
a_cov <- function(C) { p <- ncol(C); (p / (p - 1)) * (1 - sum(diag(C)) / sum(C)) }
a_cor <- function(C) { C <- cov2cor(C); p <- ncol(C); (p / (p - 1)) * (1 - p / sum(C)) }
cat("alpha on keyed Res as covariance:", round(sapply(SF, function(v) a_cov(Rk[ids(v), ids(v)])), 2), "\n")
cat("alpha on keyed Res rescaled to correlations:", round(sapply(SF, function(v) a_cor(Rk[ids(v), ids(v)])), 2), "\n")
K <- A$items_keyed; P <- polychoric(K, correct = 0)$rho
for (set in list(T15 = unlist(SF[1:3]), AV10 = unlist(SF[c(1, 3)]), T20 = unlist(SF))) {
  o <- suppressWarnings(omega(P[ids(set), ids(set)], nfactors = if (length(set) > 10) 3 else 2, n.obs = 500, plot = FALSE, flip = FALSE))
  cat(length(set), "items: omega_h (no flip)", round(o$omega_h, 3), " omega_t", round(o$omega.tot, 3), "\n")
}
