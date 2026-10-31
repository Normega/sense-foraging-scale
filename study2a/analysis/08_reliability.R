# 08: reliability (Section 5, step 7): ordinal omega for each retained factor.
#
# Model-based, from the primary model's keyed standardized loadings (latent-
# response scale): omega = (sum l)^2 / ((sum l)^2 + sum(1 - l^2)) over the
# statements whose primary loading (keyed, >= .40) is on that factor; agreement
# and cross-loadings count as error (Study 2's formula, sf2a_omega). Also ordinal
# alpha of the same statements from the polychoric matrix of keyed statements.
# PREREG: "retained factor" is read as every content factor of the primary model;
# the omega of the selected 5-6 statements is in 09.

source("study2a/analysis/config.R")
source("study2a/analysis/sf2a_models.R")
A <- readRDS(sf2a_derived("sf2a"))
def <- A$definition; k <- SF2A$k_primary
P <- readRDS(sf2a_derived(sprintf("primary_k%02d", k)))
R <- readRDS(sf2a_derived("polychoric_keyed"))
Kor <- P$keyed; fn <- P$fnames
prim <- colnames(Kor)[apply(abs(Kor), 1, which.max)]; pload <- Kor[cbind(rownames(Kor), prim)]

res <- rbindlist(lapply(fn, function(f) {
  ids <- rownames(Kor)[prim == f & pload >= SF2A$load_min]
  data.table(factor = f, matches = P$labels[[f]], n_statements = length(ids),
             ordinal_omega = if (length(ids)) round(sf2a_omega(Kor[ids, f]), 3) else NA_real_,
             ordinal_alpha = if (length(ids) >= 2) round(psych::alpha(R[ids, ids], warnings = FALSE)$total$raw_alpha, 3) else NA_real_,
             statements = paste(ids, collapse = " "))
}))
fwrite(res, sf2a_out(sprintf("08_reliability_k%02d.csv", k)))
print(res[, !"statements"])
