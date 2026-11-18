# 31: the 8 on phones vs computers. EXPLORATORY, NOT REGISTERED.
#
# Not part of the filed registration (osf.io/d4a5u) and not in its list of exploratory
# analyses. Specified on 17 November 2026 (placement card S3-27), after the registered
# analyses and 10-30 had been seen.
#
# Device is the one recorded at the pause task (M$pause_device), in the same session as
# the scale: desktop 552, phone 136, tablet 22, missing 10. Tablets and missing are left
# out. Phone (n = 136) is below SF2B$country_min (200), the registered minimum for a
# group comparison, so the invariance test is under-powered and flagged as such.
#   invariance   24's procedure unchanged (H8's steps via sf2b_run_mg, then 24's
#                intercepts step and latent means with the agreement mean free and fixed),
#                by evaluating 24's function definitions; its race groups are replaced
#                by the device groups. Reference: computer.
#   observed     07's sf2b_kg on the observed 8 (phone vs computer)
#   phones       the registered H1 model and rule (sf2b_h1) on phones alone, and on
#                computers alone
# Writes aggregates only:
#   output/real/31_invariance.csv, 31_latent_means.csv, 31_observed.csv, 31_h1_by_device.csv,
#   31_exploratory_device.json
# Run from the repo root:  SF2B_MODE=real Rscript study2b/analysis/31_exploratory_device.R

# Set Up ---------
## Load libraries ---------
packages <- c("jsonlite", "data.table", "lavaan")
new_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(new_packages)) install.packages(new_packages)
for (thispack in packages) {
  library(thispack, character.only = TRUE, quietly = TRUE, verbose = FALSE)
}

## Settings: 24's data and functions, nothing fitted ---------
ex24 <- parse("study2b/analysis/24_exploratory_invariance_race.R", keep.source = FALSE)
run_at <- which(vapply(ex24, function(e) is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name("INV")), TRUE))
stopifnot(length(run_at) == 1)
for (i in seq_len(run_at - 1)) eval(ex24[[i]], envir = globalenv())
stopifnot(exists("run_means"), exists("means_syntax"), !exists("INV"))
T0 <- Sys.time()

dev <- M$pause_device
GROUPS <- list(device = list(label = "device: computer vs phone",
                             g = fifelse(dev == "desktop", "1_computer", fifelse(dev == "phone", "2_phone", NA_character_))))
n_dev <- table(GROUPS$device$g); GROUPS$device$n <- as.list(n_dev)
GROUPS$device$below_min <- names(n_dev)[n_dev < SF2B$country_min]
sf2b_log("device: ", paste(names(n_dev), n_dev, sep = " = ", collapse = ", "))

## Invariance and latent means (24's procedure) ---------
inv <- sf2b_run_mg(X, I8, GROUPS$device$g, GROUPS$device$label)
if (!"noninvariant" %in% names(inv$table)) inv$table[, noninvariant := NA]
steps <- cbind(comparison = "device", inv$table); means <- data.table()
last <- tail(inv$table, 1)
if (grepl("^loadings", last$step) && !isTRUE(last$noninvariant)) {
  mm <- run_means("device", inv); steps <- rbind(steps, cbind(comparison = "device", mm$table), fill = TRUE); means <- mm$means
} else sf2b_log("device: loadings step not held; means not compared")

## Observed and H1 by device ---------
g01 <- ifelse(GROUPS$device$g == "2_phone", 1L, ifelse(GROUPS$device$g == "1_computer", 0L, NA_integer_))
obs <- sf2b_kg(A$scores$sf8, g01)
H1 <- rbindlist(lapply(c("1_computer", "2_phone"), function(gname) {
  rows <- which(GROUPS$device$g == gname)
  h <- sf2b_h1(X[rows, , drop = FALSE], I8, D$reversed)
  data.table(device = gname, n = length(rows), converged = h$converged, supported = h$supported,
             cfi = h$stats$cfi.scaled %||% NA_real_, rmsea = h$stats$rmsea.scaled %||% NA_real_, srmr = h$stats$srmr %||% NA_real_,
             max_residual = h$max_residual %||% NA_real_, omega = h$omega %||% NA_real_,
             min_keyed_loading_excl_U4 = if (h$converged) min(h$loadings[item != "sf2b_07", keyed_content]) else NA_real_,
             U4_keyed_loading = if (h$converged) h$loadings[item == "sf2b_07", keyed_content] else NA_real_)
}))

## Write ---------
num <- function(DT, k = 4) { cols <- names(DT)[vapply(DT, is.numeric, TRUE)]; DT[, (cols) := lapply(.SD, function(x) round(x, k)), .SDcols = cols]; DT }
fwrite(num(copy(steps)), sf2b_out("31_invariance.csv")); fwrite(num(copy(means)), sf2b_out("31_latent_means.csv"))
fwrite(num(copy(obs)), sf2b_out("31_observed.csv")); fwrite(num(copy(H1)), sf2b_out("31_h1_by_device.csv"))
write_json(list(
  analysis = "31: the 8 on phones vs computers. EXPLORATORY, NOT REGISTERED (osf.io/d4a5u does not include it); specified 17 November 2026 after 10-30 were seen (card S3-27). Phone n is below the registered 200 minimum: under-powered.",
  device = "M$pause_device, recorded at the pause task in the scale's session; tablets and missing left out.",
  n = as.list(n_dev), released = inv$free, steps = num(copy(steps)), latent_means = num(copy(means)), observed = num(copy(obs)), h1_by_device = num(copy(H1)),
  minutes = round(as.numeric(difftime(Sys.time(), T0, units = "mins")), 1)),
  sf2b_out("31_exploratory_device.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
print(num(copy(steps))[, .(step, released, cfi, rmsea, d_cfi, d_rmsea, noninvariant)])
if (nrow(means)) print(num(copy(means))[, .(model, group, n, S_d, S_d_lo, S_d_hi, RI_d, RI_d_lo, RI_d_hi)])
print(num(copy(obs))); print(num(copy(H1)))
sf2b_log("31 done in ", round(as.numeric(difftime(Sys.time(), T0, units = "mins")), 1), " min")
