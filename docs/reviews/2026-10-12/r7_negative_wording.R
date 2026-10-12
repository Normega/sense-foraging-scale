# Review check 7: the primary k=3 and k=4 ESEM + RI, plus an orthogonal negative-wording method
# factor on the 14 reversed items (free loadings). Does "Letting be" survive once wording is modelled?
suppressPackageStartupMessages({library(lavaan); library(data.table)})
A <- readRDS("F:/radlab_data/sf_study1/derived/sf1_real.rds")
X <- as.data.frame(A$items_raw); def <- A$definition; items <- def$items; rev <- def$reversed
P <- readRDS("F:/radlab_data/sf_study1/derived/primary_real_k04.rds"); member <- P$member
run <- function(fn, member) {
  target <- matrix(0, 48, length(fn), dimnames = list(items, fn)); mask <- ifelse(member, 0, 1)
  model <- paste0(paste0('efa("c")*', fn, collapse = " + "), " =~ ", paste(items, collapse = " + "), "\n",
    "RI =~ ", paste0("1*", items, collapse = " + "), "\nRI ~~ NA*RI\n", paste0("RI ~~ 0*", fn, collapse = "\n"), "\n",
    "NEG =~ ", paste(rev, collapse = " + "), "\nNEG ~~ 0*RI\n", paste0("NEG ~~ 0*", fn, collapse = "\n"))
  m <- sem(model, data = X, ordered = items, estimator = "WLSMV", rotation = "pst",
           rotation.args = list(target = target, target_mask = mask, orthogonal = FALSE, rstarts = 20))
  print(round(fitMeasures(m, c("chisq.scaled", "df", "cfi.scaled", "tli.scaled", "rmsea.scaled", "srmr")), 3))
  s <- standardizedSolution(m)
  L <- s[s$op == "=~" & s$lhs %in% fn, ]; W <- dcast(as.data.table(L)[, .(lhs, rhs, v = round(est.std, 2))], rhs ~ lhs, value.var = "v")
  W <- W[match(items, rhs)]
  neg <- s[s$op == "=~" & s$lhs == "NEG", c("rhs", "est.std")]; W[, NEG := round(neg$est.std[match(rhs, neg$rhs)], 2)]
  print(W[rhs %in% paste0("sf6_", c(20, 21, 28, 29, 38:48))])
  cat("factor correlations:\n"); print(round(lavInspect(m, "std")$psi[fn, fn], 2))
}
cat("== k = 4 ESEM + RI + NEG (targets as primary model) ==\n"); run(colnames(member), member)
