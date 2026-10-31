# 04: number of factors (Section 5, step 3; decision K).
#
# Parallel analysis on the polychoric correlations of the 58 statements,
# reverse-keyed statements recoded: principal-axis factors, 500 random datasets,
# 95th percentile. Then ordinal EFAs with 4, 5 and 6 factors (WLSMV, geomin,
# pairwise), and with the parallel-analysis count if it falls outside that range.
# For each: criterion 1 (every factor has >= 3 statements with primary loadings
# >= .40), the hypothesized factor each factor matches (criterion 2,
# interpretability, is the PI's), and fit, descriptively (criterion 3).
# The choice is recorded as SF2A_K (config.R) with its reasons.
# Also, as Study 1: a method-factor check (a factor made mostly of reverse-keyed
# statements, or on which forward and recoded reversed statements load with
# opposite signs). In simulated data with agreement, parallel analysis can count
# agreement as one extra factor.

source("study2a/analysis/config.R")
source("study2a/analysis/sf2a_models.R")
A <- readRDS(sf2a_derived("sf2a"))
def <- A$definition; K <- A$items_keyed
set.seed(SF2A$seed)

R <- polychoric(K, correct = 0)$rho
saveRDS(R, sf2a_derived("polychoric_keyed"))
pa <- PARALLEL(R, N = nrow(K), n_datasets = SF2A$pa_datasets, percent = 95, eigen_type = "EFA",
               decision_rule = "percentile")
k_pa <- max(1L, as.integer(pa$n_factors[["EFA"]]))
sf2a_log("parallel analysis suggests ", k_pa, " factor(s)")
res <- pa$results[[1]]
fwrite(data.table(factor = res$x, real = round(res$y, 3), random_95 = round(res$references[["95 Percentile"]], 3),
                  random_mean = round(res$references[["Means"]], 3)), sf2a_out("04_parallel_analysis.csv"))

ks <- sort(unique(c(4:6, k_pa)))
Kdf <- as.data.frame(K)
fits <- rbindlist(lapply(ks, function(k) {
  sf2a_log("ordinal EFA, k = ", k)
  m <- tryCatch(efa(data = Kdf, nfactors = k, ordered = TRUE, estimator = "WLSMV", rotation = "geomin",
                    missing = "pairwise", output = "lavaan"), error = function(e) e)
  if (inherits(m, "error")) return(data.table(k = k, error = conditionMessage(m)))
  L <- inspect(m, "std")$lambda
  colnames(L) <- paste0("F", seq_len(k))
  prim <- apply(abs(L), 1, which.max); load <- apply(abs(L), 1, max)
  per_factor <- tabulate(prim[load >= SF2A$load_min], nbins = k)
  fm <- fitMeasures(m, c("cfi.scaled", "rmsea.scaled", "srmr"))
  fwrite(data.table(item = rownames(L), round(L, 3)), sf2a_out(sprintf("04_efa_k%02d_loadings.csv", k)))
  rev <- rownames(L) %in% def$reversed
  lab <- sf2a_label(L, def)
  facs <- rbindlist(lapply(seq_len(k), function(j) {
    on <- abs(L[, j]) >= .30; on40 <- prim == j & load >= SF2A$load_min
    comp <- table(factor(def$factor[rownames(L)[on40]], levels = def$factors))
    data.table(k = k, factor = colnames(L)[j], n_primary_ge_40 = sum(on40), matches = lab[j],
               composition = paste(sprintf("%s %d", names(comp), comp)[comp > 0], collapse = ", "),
               n_ge_30 = sum(on), pct_reversed = round(100 * mean(rev[on]), 0),
               method_suspect = sum(on) > 0 && (mean(rev[on]) >= 2/3 ||
                 (sum(on & rev) >= 3 && sum(on & !rev) >= 3 && all(sign(L[on & rev, j]) != sign(median(L[on & !rev, j]))))))
  }))
  fwrite(facs, sf2a_out(sprintf("04_efa_k%02d_factors.csv", k)))
  data.table(k = k, from_pa = k == k_pa, criterion1_met = all(per_factor >= 3),
             min_statements_ge_40 = min(per_factor), factors_with_3plus = sum(per_factor >= 3),
             distinct_matches = length(unique(na.omit(lab))), method_suspect_factors = sum(facs$method_suspect),
             cfi = round(fm[[1]], 3), rmsea = round(fm[[2]], 3), srmr = round(fm[[3]], 3))
}), fill = TRUE)
fwrite(fits, sf2a_out("04_factor_number_models.csv"))
print(fits)
