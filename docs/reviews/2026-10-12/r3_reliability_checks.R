# Review check 3: (a) reproduce the hard-coded "decided" ordinal alphas on the page; (b) does psych::omega
# auto-flip items in 05 and 11 (flip = TRUE default)?; (c) what siren::acquiextraction returns.
suppressPackageStartupMessages({library(psych); library(siren); library(data.table)})
A <- readRDS("F:/radlab_data/sf_study1/derived/sf1_real.rds")
K <- A$items_keyed; R <- A$items_raw
SF <- list(Access = c(35, 36, 9, 8, 7), Und = c(22, 3, 15, 31, 37), Val = c(12, 18, 11, 27, 44), LB = c(40, 46, 43, 47, 21))
ids <- function(v) paste0("sf6_", v)
Rk <- polychoric(K, correct = 0)$rho
oa <- function(C) { p <- ncol(C); (p / (p - 1)) * (1 - p / sum(C)) }
cat("(a) ordinal alpha from the full 48-item polychoric (as 06 does):\n")
print(sapply(SF, function(v) round(oa(Rk[ids(v), ids(v)]), 3)))
cat("    from a polychoric of just the 5 items (as 11 does):\n")
print(sapply(SF, function(v) round(oa(polychoric(K[, ids(v)], correct = 0)$rho), 3)))
cat("\n(b) psych::omega on the 48 keyed items (05): items flipped by default\n")
om <- omega(readRDS("F:/radlab_data/sf_study1/derived/polychoric_keyed_real.rds"), nfactors = 4, n.obs = 500, fm = "minres", plot = FALSE)
fl <- om$key; cat("  flipped:", if (is.null(fl)) "none (key NULL)" else paste(names(fl)[fl < 0], collapse = " "), "\n")
cat("  omega_h", round(om$omega_h, 3), " omega_t", round(om$omega.tot, 3), "\n")
om2 <- omega(readRDS("F:/radlab_data/sf_study1/derived/polychoric_keyed_real.rds"), nfactors = 4, n.obs = 500, fm = "minres", plot = FALSE, flip = FALSE)
cat("  flip = FALSE: omega_h", round(om2$omega_h, 3), " omega_t", round(om2$omega.tot, 3), "\n")
print(round(om$omega.group, 3))
cat("  SL group factor sizes (items with |g-free loading| >= .2 on each group):\n")
sl <- om$schmid$sl; gcols <- grep("^F", colnames(sl), value = TRUE); print(colSums(abs(sl[, gcols]) >= .2))
cat("\n  omega on the 20 short-form items (11's call):\n")
P20 <- polychoric(K[, ids(unlist(SF))])$rho
o20 <- omega(P20, nfactors = 4, n.obs = 500, plot = FALSE)
cat("  flipped:", paste(names(o20$key)[o20$key < 0], collapse = " "), "| omega_h", round(o20$omega_h, 3), "omega_t", round(o20$omega.tot, 3), "\n")
o20b <- omega(P20, nfactors = 4, n.obs = 500, plot = FALSE, flip = FALSE)
cat("  flip = FALSE: omega_h", round(o20b$omega_h, 3), "omega_t", round(o20b$omega.tot, 3), "\n")
print(round(o20$schmid$sl[, 1:5], 2))
cat("\n(c) siren::acquiextraction on raw items: names of the returned object\n")
s3 <- acquiextraction(as.data.frame(R), n_factors = 4, corr = "Polychoric", raw_data = TRUE, display = FALSE)
print(names(s3)); print(str(s3, max.level = 1))
