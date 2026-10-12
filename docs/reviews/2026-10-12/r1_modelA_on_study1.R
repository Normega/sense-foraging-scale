# Review check 1: fit the Study 2 preregistered Model A (4-factor simple-structure CFA + random
# intercept, WLSMV) to Study 1's data on the 20 short-form items. In-sample, so optimistic.
suppressPackageStartupMessages({library(lavaan); library(data.table)})
A <- readRDS("F:/radlab_data/sf_study1/derived/sf1_real.rds")
X <- as.data.frame(A$items_raw)
SF <- list(Access = c(35, 36, 9, 8, 7), Und = c(22, 3, 15, 31, 37), Val = c(12, 18, 11, 27, 44), LB = c(40, 46, 43, 47, 21))
rev <- paste0("sf6_", c(35, 36, 37, 44, 40, 46, 43, 47))
items <- paste0("sf6_", unlist(SF))
X <- X[, items]
cat("N rows:", nrow(X), " complete cases on 20 items:", sum(complete.cases(X)), "\n")
# lavaan does sign via free loadings; reversed items just get negative loadings (unrecoded data)
cf <- paste(sapply(names(SF), function(f) paste0(f, " =~ ", paste0("sf6_", SF[[f]], collapse = " + "))), collapse = "\n")
ri <- paste0("RI =~ ", paste0("1*", items, collapse = " + "), "\nRI ~~ NA*RI\n", paste0("RI ~~ 0*", names(SF), collapse = "\n"))
fitm <- function(model, miss) {
  m <- tryCatch(cfa(model, data = X, ordered = items, estimator = "WLSMV", parameterization = "theta", missing = miss, std.lv = TRUE),
                error = function(e) e)
  if (inherits(m, "error")) return(data.table(err = conditionMessage(m)))
  fm <- fitMeasures(m, c("ntotal", "chisq.scaled", "df", "cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.upper.scaled", "srmr",
                         "cfi.robust", "rmsea.robust", "rmsea.ci.upper.robust"))
  list(fit = round(fm, 3), m = m)
}
cat("\n== Model A (4 CFA + RI), listwise (lavaan default) ==\n"); a1 <- fitm(paste(cf, ri, sep = "\n"), "listwise"); print(a1$fit)
cat("\n== Model A, pairwise ==\n"); a2 <- fitm(paste(cf, ri, sep = "\n"), "pairwise"); print(a2$fit)
cat("\n== 4 CFA, no RI ==\n"); a3 <- fitm(cf, "listwise"); print(a3$fit)
s <- standardizedSolution(a1$m)
cat("\nStandardized loadings (Model A, listwise):\n"); print(s[s$op == "=~" & s$lhs != "RI", c("lhs", "rhs", "est.std", "se")], row.names = FALSE)
cat("\nFactor correlations:\n"); print(round(lavInspect(a1$m, "std")$psi, 2))
cat("\nRI variance:\n"); pe <- parameterEstimates(a1$m); print(pe[pe$lhs == "RI" & pe$op == "~~" & pe$rhs == "RI", c("est", "se", "pvalue")])
mi <- modindices(a1$m, sort. = TRUE); cat("\nTop 12 modification indices:\n"); print(head(mi[, c("lhs", "op", "rhs", "mi", "sepc.all")], 12), row.names = FALSE)
# Ordinal omega per factor from Model A (theta param): omega = (sum lambda)^2 / ((sum lambda)^2 + sum theta) on latent responses
cat("\nOmega per factor (latent-response, content only):\n")
for (f in names(SF)) {
  l <- s[s$op == "=~" & s$lhs == f, "est.std"]
  cat(f, round(sum(abs(l))^2 / (sum(abs(l))^2 + sum(1 - l^2)), 3), "\n")
}
