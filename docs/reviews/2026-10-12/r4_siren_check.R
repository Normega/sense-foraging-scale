# Review check 4: what does siren::acquiextraction's Res contain, and is it R minus the
# acquiescence outer product? Also: is the balanced-core assumption met?
suppressPackageStartupMessages({library(siren); library(psych)})
A <- readRDS("F:/radlab_data/sf_study1/derived/sf1_real.rds")
R <- A$items_raw
s3 <- acquiextraction(as.data.frame(R), n_factors = 4, corr = "Polychoric", raw_data = TRUE)
P <- psych::polychoric(as.data.frame(R), smooth = FALSE, correct = FALSE)$rho
acq <- as.numeric(s3$acq)
cat("acq loadings: length", length(acq), " range", round(range(acq), 3), " mean", round(mean(acq), 3), "\n")
cat("acq loadings, positive items mean", round(mean(acq[1:34]), 3), "; reversed items mean", round(mean(acq[35:48]), 3), "\n")
D <- P - s3$Res
cat("P - Res: unique values (first 5):", head(unique(round(as.vector(D), 6)), 5), "\n")
cat("max |(P - Res) - outer(acq,acq)| =", signif(max(abs(D - outer(acq, acq))), 3), "\n")
cat("diag of Res (first 5):", round(diag(s3$Res)[1:5], 3), "\n")
cat("AQ_var reported:", round(s3$AQ_var, 3), " sum(acq)/... vs sum(acq^2):", round(sum(acq^2), 3), "\n")
