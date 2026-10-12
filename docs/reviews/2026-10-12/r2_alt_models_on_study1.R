# Review check 2: in-sample fit on Study 1 of (a) the 14-item Access+Valuing form (2 CFA + RI),
# (b) 15 items without Letting be (3 CFA + RI), (c) the 20 items as a 4-factor ESEM + RI (target
# rotation), to locate the misfit of Model A. All in-sample (selection sample): optimistic.
suppressPackageStartupMessages({library(lavaan); library(data.table)})
A <- readRDS("F:/radlab_data/sf_study1/derived/sf1_real.rds")
X0 <- as.data.frame(A$items_raw)
fit_cfa <- function(SF, label) {
  items <- paste0("sf6_", unlist(SF))
  cf <- paste(sapply(names(SF), function(f) paste0(f, " =~ ", paste0("sf6_", SF[[f]], collapse = " + "))), collapse = "\n")
  ri <- paste0("RI =~ ", paste0("1*", items, collapse = " + "), "\nRI ~~ NA*RI\n", paste0("RI ~~ 0*", names(SF), collapse = "\n"))
  m <- cfa(paste(cf, ri, sep = "\n"), data = X0[, items], ordered = items, estimator = "WLSMV", parameterization = "theta", std.lv = TRUE)
  fm <- fitMeasures(m, c("ntotal", "chisq.scaled", "df", "cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.upper.scaled", "srmr"))
  cat("\n==", label, "==\n"); print(round(fm, 3))
  s <- standardizedSolution(m); print(s[s$op == "=~" & s$lhs != "RI", c("lhs", "rhs", "est.std")], row.names = FALSE)
  print(round(lavInspect(m, "std")$psi, 2))
  invisible(m)
}
fit_cfa(list(Access = c(9, 7, 6, 2, 8, 14, 36, 35), Val = c(18, 12, 25, 27, 28, 44)), "14-item Access 8 + Valuing 6, CFA + RI")
fit_cfa(list(Access = c(35, 36, 9, 8, 7), Und = c(22, 3, 15, 31, 37), Val = c(12, 18, 11, 27, 44)), "15 items, 3 CFA + RI")
fit_cfa(list(Access = c(35, 36, 9, 8, 7), Val = c(12, 18, 11, 27, 44)), "10 items Access + Valuing (from the 20), CFA + RI")
# (c) ESEM on the 20 with target rotation
SF <- list(Access = c(35, 36, 9, 8, 7), Und = c(22, 3, 15, 31, 37), Val = c(12, 18, 11, 27, 44), LB = c(40, 46, 43, 47, 21))
items <- paste0("sf6_", unlist(SF)); fn <- names(SF)
member <- sapply(fn, function(f) items %in% paste0("sf6_", SF[[f]])); rownames(member) <- items
target <- matrix(0, 20, 4, dimnames = dimnames(member)); mask <- ifelse(member, 0, 1)
model <- paste0(paste0('efa("c")*', fn, collapse = " + "), " =~ ", paste(items, collapse = " + "), "\n",
                "RI =~ ", paste0("1*", items, collapse = " + "), "\nRI ~~ NA*RI\n", paste0("RI ~~ 0*", fn, collapse = "\n"))
m <- sem(model, data = X0[, items], ordered = items, estimator = "WLSMV", rotation = "pst",
         rotation.args = list(target = target, target_mask = mask, orthogonal = FALSE, rstarts = 30))
cat("\n== 20 items, 4-factor ESEM (target) + RI ==\n")
print(round(fitMeasures(m, c("ntotal", "chisq.scaled", "df", "cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.upper.scaled", "srmr")), 3))
s <- standardizedSolution(m); L <- s[s$op == "=~" & s$lhs %in% fn, ]
W <- dcast(as.data.table(L)[, .(lhs, rhs, est.std = round(est.std, 2))], rhs ~ lhs, value.var = "est.std"); print(W[match(items, rhs)])
