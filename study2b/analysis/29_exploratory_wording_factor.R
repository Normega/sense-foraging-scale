# 29: a wording (reverse-keying) factor in place of, and beside, the agreement factor.
# EXPLORATORY, NOT REGISTERED.
#
# Not part of the filed registration (osf.io/d4a5u) and not in its list of exploratory
# analyses. Specified on 17 November 2026 (placement card S3-25, from the Fable review),
# after the registered analyses and 10-28 had been seen.
#
# Models of the 8 (ordinal, WLSMV, theta, missing = "pairwise", conditional.x = FALSE):
#   content_only   S alone
#   registered     S + the agreement factor RI (loadings 1 on every unrecoded statement,
#                  variance free, orthogonal to S): the H1 model
#   wording        S + W, a method factor on the reverse-keyed statements only (loadings
#                  free, variance 1, orthogonal to S)
#   both           S + RI + W (W orthogonal to S and RI)
# Reported: fit, S's keyed standardized loadings, content omega (method factors counted
# as error: (sum kl)^2 / ((sum kl)^2 + sum over method factors (sum k m)^2 psi + k'Theta k)),
# and S's latent r with FFMQ Observing, SBI Savoring and PSS-4 under each model (04's
# validity syntax; validity factors orthogonal to RI and W).
# Writes aggregates only:
#   output/real/29_fit.csv, 29_loadings.csv, 29_validity.csv, 29_exploratory_wording_factor.json
# Run from the repo root:  SF2B_MODE=real Rscript study2b/analysis/29_exploratory_wording_factor.R

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

I8 <- D$sf8; REV8 <- intersect(I8, D$reversed); k <- ifelse(I8 %in% REV8, -1, 1)
stopifnot(length(REV8) == 3)
S_LINES <- c(paste0("S =~ NA*", paste(I8, collapse = " + ")), "S ~~ 1*S")
RI_LINES <- c(paste0("RI =~ ", paste0("1*", I8, collapse = " + ")), "RI ~~ RI", "S ~~ 0*RI")
W_LINES <- c(paste0("W =~ NA*", REV8[1], " + ", paste0("start(.4)*", REV8, collapse = " + ")), "W ~~ 1*W", "S ~~ 0*W")
MODELS <- list(content_only = S_LINES, registered = c(S_LINES, RI_LINES), wording = c(S_LINES, W_LINES),
               both = c(S_LINES, RI_LINES, W_LINES, "RI ~~ 0*W"))
CRIT <- c("ffmq_observing", "sbi_moment", "pss4")

omega_of <- function(f) {
  st <- lavInspect(f, "std"); L <- st$lambda[I8, , drop = FALSE]; P <- st$psi
  cs <- sum(k * L[, "S"])^2
  meth <- sum(vapply(setdiff(colnames(L), "S"), function(m) sum(k * L[, m])^2 * P[m, m], 0))
  er <- as.numeric(t(k) %*% st$theta[I8, I8] %*% k)
  cs / (cs + meth + er)
}

## Fit ---------
FIT <- list(); LOAD <- list(); VAL <- list()
for (nm in names(MODELS)) {
  f <- sf2b_try(sf2b_sem(MODELS[[nm]], DAT, ordered = I8))
  if (sf2b_failed(f)) { sf2b_log(nm, " FAILED: ", sf2b_msg(f)); FIT[[nm]] <- list(model = nm, error = sf2b_msg(f)); next }
  st <- lavInspect(f, "std")
  FIT[[nm]] <- c(list(model = nm, omega = omega_of(f)), sf2b_fit_stats(f))
  LOAD[[nm]] <- data.table(model = nm, item = I8, source = D$source[I8], reversed = I8 %in% REV8,
                           keyed_content = st$lambda[I8, "S"] * k,
                           agreement = if ("RI" %in% colnames(st$lambda)) st$lambda[I8, "RI"] else NA_real_,
                           wording = if ("W" %in% colnames(st$lambda)) st$lambda[I8, "W"] else NA_real_,
                           residual = diag(st$theta)[I8])
  for (cr in CRIT) {
    syn <- c(MODELS[[nm]], sf2b_validity_syntax(cr, D, ri = "RI" %in% colnames(st$lambda)),
             if ("W" %in% colnames(st$lambda)) paste0("W ~~ 0*", cr))
    fv <- sf2b_try(sf2b_sem(syn, DAT, ordered = unique(c(I8, score_items[[cr]]))))
    if (sf2b_failed(fv)) { sf2b_log(nm, " ", cr, " FAILED: ", sf2b_msg(fv)); next }
    pe <- as.data.table(parameterEstimates(fv))[label == paste0("r_", cr)]
    VAL[[length(VAL) + 1]] <- data.table(model = nm, criterion = cr, r = pe$est, lower = pe$ci.lower, upper = pe$ci.upper)
  }
  sf2b_log(nm, " done")
}
fit <- rbindlist(FIT, fill = TRUE); load <- rbindlist(LOAD); val <- rbindlist(VAL)

## Write ---------
num <- function(DT, kk = 4) { cols <- names(DT)[vapply(DT, is.numeric, TRUE)]; DT[, (cols) := lapply(.SD, function(x) round(x, kk)), .SDcols = cols]; DT }
fwrite(num(copy(fit)), sf2b_out("29_fit.csv")); fwrite(num(copy(load)), sf2b_out("29_loadings.csv")); fwrite(num(copy(val)), sf2b_out("29_validity.csv"))
write_json(list(
  analysis = "29: a wording factor in place of, and beside, the agreement factor. EXPLORATORY, NOT REGISTERED (osf.io/d4a5u does not include it); specified 17 November 2026 after 10-28 were seen (card S3-25).",
  reversed = REV8, fit = num(copy(fit)), loadings = num(copy(load)), validity = num(copy(val)),
  minutes = round(as.numeric(difftime(Sys.time(), T0, units = "mins")), 1)),
  sf2b_out("29_exploratory_wording_factor.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
print(num(copy(fit), 3)[, intersect(c("model", "omega", "cfi.scaled", "rmsea.scaled", "srmr", "chisq.scaled", "df", "error"), names(fit)), with = FALSE])
print(dcast(num(copy(load), 3), item + reversed ~ model, value.var = "keyed_content"))
print(load[model %in% c("wording", "both") & reversed, .(model, item, agreement = round(agreement, 3), wording = round(wording, 3))])
print(dcast(num(copy(val), 3), criterion ~ model, value.var = "r"))
sf2b_log("29 done in ", round(as.numeric(difftime(Sys.time(), T0, units = "mins")), 1), " min")
