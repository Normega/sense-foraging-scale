# 34: is the 8 the same measure in Study 2a and Study 2b? EXPLORATORY, NOT REGISTERED.
#
# Not part of either filed registration (osf.io/vt7gn, osf.io/d4a5u). Specified on 17
# November 2026 (placement card X-12, from the Fable review), after both studies'
# registered analyses had been seen.
#
# The 8 statements were answered in 2a (n = 500, among 58 statements, in the order 2a
# randomised) and in 2b (n = 720, on their own after the trial statements and the
# vignettes). Both use the same 1-6 agreement scale. Statements are matched by their
# source id (2a's definition$source, 2b's defs$source); their wording is compared and the
# result recorded (equal or not), never printed.
#   invariance   24's procedure unchanged (H8's steps via sf2b_run_mg, then 24's
#                intercepts step and latent means with the agreement mean free and fixed),
#                groups 2a (reference) vs 2b. Sample, context and position all differ, so a
#                difference cannot be assigned to any one of them.
#   H1 in 2a     the registered 2b H1 model and rule (sf2b_h1) on 2a's answers alone
# Reads 2a's derived object (F:/radlab_data/sf_study2a/derived/sf2a_real.rds), 2a's own
# pipeline output; writes aggregates only:
#   output/real/34_invariance.csv, 34_latent_means.csv, 34_h1_2a.csv, 34_exploratory_cross_study.json
# Run from the repo root:  SF2B_MODE=real Rscript study2b/analysis/34_exploratory_cross_study.R

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
stopifnot(exists("run_means"), !exists("INV"))
T0 <- Sys.time()

## Study 2a's answers to the 8, matched by source id ---------
A2 <- readRDS("F:/radlab_data/sf_study2a/derived/sf2a_real.rds")
stopifnot(identical(A2$mode, "real"))
d2 <- A2$definition
src8 <- D$source[I8]
id2a <- setNames(d2$items[match(src8, d2$source)], I8)
stopifnot(!anyNA(id2a))
rev_match <- identical(sort(names(src8)[src8 %in% d2$source[d2$items %in% d2$reversed]]), sort(intersect(I8, D$reversed)))
# Wording and response scale from the two studies' questionnaire definitions in the repo.
qmap <- function(f) {
  d <- fromJSON(f, simplifyVector = FALSE)
  comps <- Filter(function(c) !is.null(c$question), unlist(lapply(d$pages, `[[`, "components"), recursive = FALSE))
  ids <- vapply(comps, `[[`, "", "id")
  list(q = setNames(vapply(comps, function(c) trimws(c$question), ""), ids),
       s = setNames(vapply(comps, function(c) as.character(toJSON(c$scale, auto_unbox = TRUE)), ""), ids))
}
qb <- qmap("study2b/materials/sf2b-scale.json"); qa <- qmap("study2a/materials/sf2a-scale.json")
stopifnot(all(I8 %in% names(qb$q)), all(id2a %in% names(qa$q)))
text_2b <- qb$q[I8]
wording_equal <- setNames(qb$q[I8] == qa$q[unname(id2a)], I8)
scale_equal <- setNames(qb$s[I8] == qa$s[unname(id2a)], I8)
stopifnot(all(wording_equal), all(scale_equal))
X2 <- as.matrix(A2$items_raw[, unname(id2a)]); colnames(X2) <- I8; storage.mode(X2) <- "double"
stopifnot(all(X2 %in% c(1:6, NA)), all(X[, I8] %in% c(1:6, NA)))

Xall <- rbind(X2, as.matrix(X[, I8]))
gall <- c(rep("1_study2a", nrow(X2)), rep("2_study2b", nrow(X)))
X <- Xall   # run_means() reads the global X
GROUPS <- list(study = list(label = "study: 2a vs 2b", g = gall))
sf2b_log("2a n = ", nrow(X2), ", 2b n = ", sum(gall == "2_study2b"), "; keying matches: ", rev_match,
         "; wording equal: ", if (is.null(text_2b)) "not checked (2b defs carry no text)" else sum(wording_equal), " of 8")

## Invariance and latent means ---------
inv <- sf2b_run_mg(Xall, I8, gall, GROUPS$study$label)
if (!"noninvariant" %in% names(inv$table)) inv$table[, noninvariant := NA]
steps <- cbind(comparison = "study", inv$table); means <- data.table()
last <- tail(inv$table, 1)
if (grepl("^loadings", last$step) && !isTRUE(last$noninvariant)) {
  mm <- run_means("study", inv); steps <- rbind(steps, cbind(comparison = "study", mm$table), fill = TRUE); means <- mm$means
} else sf2b_log("study: loadings step not held; means not compared")

## H1 in 2a alone ---------
h <- sf2b_h1(X2, I8, D$reversed)
h1_2a <- data.table(study = "2a", n = nrow(X2), converged = h$converged, supported = h$supported,
                    cfi = h$stats$cfi.scaled, rmsea = h$stats$rmsea.scaled, srmr = h$stats$srmr, max_residual = h$max_residual,
                    max_residual_pair = h$max_residual_pair, omega = h$omega)
load_2a <- h$loadings[, .(item, source = src8[item], keyed_content, agreement, residual_variance)]

## Write ---------
num <- function(DT, k = 4) { cols <- names(DT)[vapply(DT, is.numeric, TRUE)]; DT[, (cols) := lapply(.SD, function(x) round(x, k)), .SDcols = cols]; DT }
fwrite(num(copy(steps)), sf2b_out("34_invariance.csv")); fwrite(num(copy(means)), sf2b_out("34_latent_means.csv"))
fwrite(num(copy(h1_2a)), sf2b_out("34_h1_2a.csv")); fwrite(num(copy(load_2a)), sf2b_out("34_loadings_2a.csv"))
write_json(list(
  analysis = "34: the 8 in Study 2a vs Study 2b. EXPLORATORY, NOT REGISTERED (neither osf.io/vt7gn nor osf.io/d4a5u includes it); specified 17 November 2026 after both studies' analyses were seen (card X-12).",
  matched_by = "source id", id_2a = as.list(id2a), keying_matches = rev_match,
  wording_equal = if (is.null(text_2b)) "not checked: 2b's definitions carry no statement text" else as.list(wording_equal),
  n = list(study2a = nrow(X2), study2b = sum(gall == "2_study2b")), released = inv$free,
  steps = num(copy(steps)), latent_means = num(copy(means)), h1_2a = num(copy(h1_2a)), loadings_2a = num(copy(load_2a)),
  minutes = round(as.numeric(difftime(Sys.time(), T0, units = "mins")), 1)),
  sf2b_out("34_exploratory_cross_study.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
print(num(copy(steps))[, .(step, released, cfi, rmsea, d_cfi, d_rmsea, noninvariant)])
if (nrow(means)) print(num(copy(means))[, .(model, group, n, S_d, S_d_lo, S_d_hi, RI_d, RI_d_lo, RI_d_hi, S_var)])
print(num(copy(h1_2a), 3)); print(num(copy(load_2a), 3))
sf2b_log("34 done in ", round(as.numeric(difftime(Sys.time(), T0, units = "mins")), 1), " min")
