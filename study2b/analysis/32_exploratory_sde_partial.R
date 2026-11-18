# 32: the 8's validity correlations with BIDR self-deceptive enhancement partialled out.
# EXPLORATORY, NOT REGISTERED.
#
# Not part of the filed registration (osf.io/d4a5u) and not in its list of exploratory
# analyses. Specified on 17 November 2026 (placement card S3-28, from the Fable review),
# after the registered analyses and 10-31 had been seen. H3c failed: the 8's latent r
# with SDE is .48, so this asks how much of the validity pattern survives it.
#
# One model per measure: 04's model for that measure (Brief MAIA-2 with decision 46's
# bound) plus the BIDR SDE factor on its items (04's validity syntax; correlated with S
# and with the measure's factors, uncorrelated with the agreement factor). Labels:
# r_<x> = cor(S, x), r_bidr_sde = cor(S, SDE), q_<x> = cor(SDE, x); the partial
# correlation is a defined parameter
#   p_<x> := (r_<x> - r_bidr_sde * q_<x>) / sqrt((1 - r_bidr_sde^2) * (1 - q_<x>^2))
# with its delta-method 95% CI. Single items as single-indicator factors (04). The BIDR
# block itself is left out (SDE would be partialled from itself); IM is reported from a
# model of the 8, SDE and IM.
# Settings as registered: ordinal items, WLSMV, theta, missing = "pairwise",
# conditional.x = FALSE. Writes aggregates only:
#   output/real/32_partial.csv, 32_fit.csv, 32_exploratory_sde_partial.json
# Run from the repo root:  SF2B_MODE=real Rscript study2b/analysis/32_exploratory_sde_partial.R

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

I8 <- D$sf8; SDE <- "bidr_sde"; MAIA_BOUND <- 9
SINGLE_CRIT <- c("lifesat", "crit_pause_days", "crit_notice_first", "crit_outdoors", "engaged_pause")
BLOCKS <- SF2B_MEASURES[setdiff(names(SF2B_MEASURES), "bidr")]
BLOCKS$bidr_im <- "bidr_im"
partial_line <- function(x) sprintf("p_%s := (r_%s - r_%s * q_%s) / sqrt((1 - r_%s^2) * (1 - q_%s^2))", x, x, SDE, x, SDE, x)

model_of <- function(t) {
  syn <- sf2b_scale_syntax(I8); ord <- I8
  if (!is.null(t$single)) {
    v <- t$single
    syn <- c(syn, sf2b_validity_syntax(SDE, D),
             paste0("L_", v, " =~ 1*", v), paste0(v, " ~~ 0*", v), paste0("L_", v, " ~~ 1*L_", v),
             paste0("S ~~ r_", v, "*L_", v), paste0("RI ~~ 0*L_", v), paste0(SDE, " ~~ q_", v, "*L_", v), partial_line(v))
    return(list(syn = syn, ord = unique(c(ord, score_items[[SDE]], v)), crit = v))
  }
  scores <- unique(c(SDE, t$scores)); crit <- setdiff(scores, SDE)
  syn <- c(syn, sf2b_validity_syntax(scores, D))
  if ("maia_notworrying" %in% scores) syn <- c(syn, "maia_notworrying =~ l_maia12*maia_12", sprintf("l_maia12^2 < %s", MAIA_BOUND))
  syn <- c(syn, paste0(SDE, " ~~ q_", crit, "*", crit), vapply(crit, partial_line, ""))
  list(syn = syn, ord = unique(c(ord, unlist(score_items[scores]))), crit = crit)
}
TASKS <- c(lapply(names(BLOCKS), function(b) list(block = b, scores = BLOCKS[[b]], single = NULL)),
           lapply(SINGLE_CRIT, function(v) list(block = v, scores = character(), single = v)))
run_task <- function(t) {
  m <- model_of(t)
  f <- sf2b_try(sf2b_sem(m$syn, DAT, ordered = m$ord))
  if (sf2b_failed(f)) return(list(task = t, ok = FALSE, error = sf2b_msg(f)))
  bad <- sf2b_orientation(f, I8, if (is.null(t$single)) unique(c(SDE, t$scores)) else SDE, D)
  if (length(bad)) return(list(task = t, ok = FALSE, error = paste("not oriented as keyed:", paste(bad, collapse = ", "))))
  pe <- as.data.table(parameterEstimates(f))
  g <- function(l) pe[label == l]
  out <- rbindlist(lapply(m$crit, function(x) data.table(block = t$block, criterion = x,
    r = g(paste0("r_", x))$est, r_lo = g(paste0("r_", x))$ci.lower, r_hi = g(paste0("r_", x))$ci.upper,
    r_sde_x = g(paste0("q_", x))$est, r_8_sde = g(paste0("r_", SDE))$est,
    partial = g(paste0("p_", x))$est, partial_lo = g(paste0("p_", x))$ci.lower, partial_hi = g(paste0("p_", x))$ci.upper)))
  list(task = t, ok = TRUE, r = out, fit = c(list(block = t$block), sf2b_fit_stats(f)))
}
sf2b_log("fitting ", length(TASKS), " models")
RES <- sf2b_par_lapply(TASKS, run_task, data = list(DAT = DAT, D = D, I8 = I8, SDE = SDE, MAIA_BOUND = MAIA_BOUND,
                                                    score_items = score_items, model_of = model_of, partial_line = partial_line))
for (r in Filter(function(r) !r$ok, RES)) sf2b_log("FAILED: ", r$task$block, ": ", r$error)
pr <- rbindlist(lapply(Filter(function(r) r$ok, RES), `[[`, "r"))
pr[, change := partial - r]
fits <- rbindlist(lapply(Filter(function(r) r$ok, RES), `[[`, "fit"), fill = TRUE)

## Write ---------
num <- function(DT, k = 4) { cols <- names(DT)[vapply(DT, is.numeric, TRUE)]; DT[, (cols) := lapply(.SD, function(x) round(x, k)), .SDcols = cols]; DT }
fwrite(num(copy(pr)), sf2b_out("32_partial.csv")); fwrite(fits, sf2b_out("32_fit.csv"))
write_json(list(
  analysis = "32: the 8's validity correlations with BIDR SDE partialled out. EXPLORATORY, NOT REGISTERED (osf.io/d4a5u does not include it); specified 17 November 2026 after 10-31 were seen (card S3-28).",
  partial = num(copy(pr)), failed = lapply(Filter(function(r) !r$ok, RES), function(r) list(block = r$task$block, error = r$error)),
  minutes = round(as.numeric(difftime(Sys.time(), T0, units = "mins")), 1)),
  sf2b_out("32_exploratory_sde_partial.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
print(num(copy(pr), 3)[, .(criterion, r, r_sde_x, r_8_sde, partial, partial_lo, partial_hi, change)], nrows = 60)
sf2b_log("32 done in ", round(as.numeric(difftime(Sys.time(), T0, units = "mins")), 1), " min")
