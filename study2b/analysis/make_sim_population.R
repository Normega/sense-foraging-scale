# The population the 2b simulation generates the scale from: Study 2a's fitted
# one-factor model of the 13 statements Study 2b administers (the 10 and the three
# exploratory ones), with the random-intercept agreement factor. Writes aggregates only
# (standardized loadings and thresholds), to study2b/analysis/sim_population.json, which
# is committed; 00_simulate.R reads it. Not part of run_all.R: run once, from the repo
# root, with Study 2a's derived data present:
#   SF2A_MODE=real Rscript study2b/analysis/make_sim_population.R

source("study2a/analysis/config.R")
source("study2a/analysis/sf2a_models.R")
stopifnot(SF2A$mode == "real")

A <- readRDS(sf2a_derived("sf2a")); X <- A$items_raw
s <- fromJSON("study2b/materials/sf2b-scale.json", simplifyVector = FALSE)
map <- setNames(vapply(s$items, `[[`, "", "sf2a_id"), vapply(s$items, `[[`, "", "id"))
X <- X[, unname(map)]; colnames(X) <- names(map)
it <- names(map)

syn <- c(paste0("S =~ NA*", paste(it, collapse = " + ")), "S ~~ 1*S",
         paste0("RI =~ ", paste0("1*", it, collapse = " + ")), "RI ~~ RI", "S ~~ 0*RI")
fit <- sf2a_sem(syn, X, ordered = it)
stopifnot(lavInspect(fit, "converged"))
std <- lavInspect(fit, "std")
D <- sqrt(diag(lavInspect(fit, "implied")$cov[it, it]))
th <- lavInspect(fit, "th")
pop <- list(
  source = "Study 2a (osf.io/vt7gn), real data, one content factor + random intercept on the 13 statements Study 2b administers; WLSMV, theta, pairwise",
  n = nrow(X),
  items = it,
  loading_S = round(unname(std$lambda[it, "S"]), 4),
  loading_RI = round(unname(std$lambda[it, "RI"]), 4),
  thresholds = lapply(setNames(it, it), function(v) round(unname(sort(th[startsWith(names(th), paste0(v, "|"))] / D[v])), 4)),
  fit = as.list(round(fitMeasures(fit, c("cfi.scaled", "rmsea.scaled", "srmr")), 3)))
write_json(pop, "study2b/analysis/sim_population.json", auto_unbox = TRUE, pretty = TRUE, digits = NA)
print(data.frame(item = it, S = pop$loading_S, RI = pop$loading_RI)); print(unlist(pop$fit))
