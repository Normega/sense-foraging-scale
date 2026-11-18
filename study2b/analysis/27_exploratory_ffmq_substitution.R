# 27: the 8 in place of Observing in the five-facet FFMQ structure. EXPLORATORY, NOT
# REGISTERED.
#
# Not part of the filed registration (osf.io/d4a5u) and not in its list of exploratory
# analyses. Specified on 17 November 2026 (placement card S3-23, from the Fable review),
# after the registered analyses and 10-26 had been seen.
#
# 26 already compares the 8's and Observing's correlations with the other four facets
# (from 04's FFMQ model). What it cannot show is how each sits in the facet structure as
# a whole. Here, second-order models of the FFMQ-15:
#   original     M (overall mindfulness) =~ Observing, Describing, Acting, Nonjudging,
#                Nonreactivity; each facet on its three items
#   substituted  the same with the 8's content factor S (with its agreement factor, as
#                registered) in place of Observing
#   both         the original plus S as a sixth first-order factor under M
# Reported: standardized second-order loadings, fit, and the variance of each first-order
# factor that M explains. A first-order factor's residual variance is fixed at 1 and the
# second-order loadings are free (M's variance 1), so the model is identified without a
# marker. Starting values as 04 (+/-.6 by keying), orientation checked as 04.
#
# Settings as registered: ordinal items, WLSMV, theta, missing = "pairwise",
# conditional.x = FALSE. Writes aggregates only:
#   output/real/27_second_order.csv, 27_fit.csv, 27_exploratory_ffmq_substitution.json
# Run from the repo root:  SF2B_MODE=real Rscript study2b/analysis/27_exploratory_ffmq_substitution.R

# Set Up ---------
## Load libraries ---------
packages <- c("jsonlite", "data.table", "lavaan")
new_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(new_packages)) install.packages(new_packages)
for (thispack in packages) {
  library(thispack, character.only = TRUE, quietly = TRUE, verbose = FALSE)
}

## Settings: 04's data, nothing fitted ---------
source("study2b/analysis/config.R")
source("study2b/analysis/sf2b_models.R")
stopifnot(SF2B$mode == "real")
T0 <- Sys.time()
ex04 <- parse("study2b/analysis/04_validity.R", keep.source = FALSE)
first_fit <- which(vapply(ex04, function(e) grepl("^sf2b_log\\(\"fitting", deparse(e)[1]), TRUE))
stopifnot(length(first_fit) == 1)
for (i in seq_len(first_fit - 1)) eval(ex04[[i]], envir = globalenv())
stopifnot(exists("DAT"), exists("score_items"), !exists("RES"))

I8 <- D$sf8; FACETS <- SF2B_MEASURES$ffmq; OBS <- "ffmq_observing"
RV <- c(D$reversed, unlist(lapply(D$validity, `[[`, "reverse")))
facet_lines <- function(f) {
  it <- score_items[[f]]
  c(paste0(f, " =~ NA*", it[1], " + ", paste0("start(", ifelse(it %in% RV, -.6, .6), ")*", it, collapse = " + ")),
    paste0(f, " ~~ 1*", f), paste0("RI ~~ 0*", f))
}
s_lines <- c(paste0("S =~ NA*", paste(I8, collapse = " + ")), "S ~~ 1*S",
             paste0("RI =~ ", paste0("1*", I8, collapse = " + ")), "RI ~~ RI", "S ~~ 0*RI")
second <- function(firsts) c(paste0("M =~ NA*", firsts[1], " + ", paste(paste0("start(.5)*", firsts), collapse = " + ")),
                             "M ~~ 1*M", "RI ~~ 0*M")
MODELS <- list(
  original    = list(firsts = FACETS),
  substituted = list(firsts = c("S", setdiff(FACETS, OBS))),
  both        = list(firsts = c("S", FACETS)))
build <- function(firsts) {
  syn <- c(s_lines, unlist(lapply(setdiff(firsts, "S"), facet_lines)), second(firsts))
  ord <- c(I8, unlist(score_items[setdiff(firsts, "S")]))
  if (!"S" %in% firsts) { syn <- syn[!grepl("^(S|RI) ", syn) & !grepl("^RI ~~ 0", syn)]; ord <- setdiff(ord, I8) }
  list(syn = syn, ord = unique(ord))
}

## Fit ---------
OUT <- list(); FIT <- list()
for (nm in names(MODELS)) {
  fs <- MODELS[[nm]]$firsts; b <- build(fs)
  f <- sf2b_try(sf2b_sem(b$syn, DAT, ordered = b$ord))
  if (sf2b_failed(f)) { sf2b_log(nm, " FAILED: ", sf2b_msg(f)); next }
  ss <- as.data.table(standardizedSolution(f))
  so <- ss[lhs == "M" & op == "=~", .(model = nm, factor = rhs, loading = est.std, lower = ci.lower, upper = ci.upper)]
  # orientation: M should load positively on Describing (keyed positive); flip all if not
  if (so[factor == "ffmq_describing", loading] < 0) so[, c("loading", "lower", "upper") := .(-loading, -upper, -lower)]
  so[, r2 := loading^2]
  OUT[[nm]] <- so; FIT[[nm]] <- c(list(model = nm), sf2b_fit_stats(f))
  sf2b_log(nm, ": ", paste(so$factor, round(so$loading, 3), collapse = ", "))
}
so <- rbindlist(OUT); fit <- rbindlist(FIT, fill = TRUE)

## Write ---------
num <- function(DT, k = 4) { cols <- names(DT)[vapply(DT, is.numeric, TRUE)]; DT[, (cols) := lapply(.SD, function(x) round(x, k)), .SDcols = cols]; DT }
fwrite(num(copy(so)), sf2b_out("27_second_order.csv")); fwrite(fit, sf2b_out("27_fit.csv"))
write_json(list(
  analysis = "27: the 8 in place of Observing in a second-order FFMQ-15 model. EXPLORATORY, NOT REGISTERED (osf.io/d4a5u does not include it); specified 17 November 2026 after 10-26 were seen (card S3-23).",
  models = lapply(MODELS, `[[`, "firsts"), second_order = num(copy(so)), fit = fit,
  minutes = round(as.numeric(difftime(Sys.time(), T0, units = "mins")), 1)),
  sf2b_out("27_exploratory_ffmq_substitution.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
print(num(copy(so))); print(fit[, .(model, cfi.scaled, rmsea.scaled, srmr)])
sf2b_log("27 done in ", round(as.numeric(difftime(Sys.time(), T0, units = "mins")), 1), " min")
