# 02: descriptives. Item distributions (the 13 statements, the trial statements, every
# validity item), missing data per measure, internal consistency per observed score
# (alpha, items keyed), the sample, and minutes per questionnaire.

source("study2b/analysis/config.R")
source("study2b/analysis/sf2b_models.R")
A <- readRDS(sf2b_derived("sf2b")); D <- A$defs; M <- A$meta

item_desc <- function(X, measure, range) {
  cats <- seq(range[1], range[2])
  rbindlist(lapply(colnames(X), function(j) { x <- X[, j]; ok <- !is.na(x)
    pct <- if (any(ok)) round(100 * tabulate(match(x[ok], cats), length(cats)) / sum(ok), 1) else rep(NA_real_, length(cats))
    c(list(measure = measure, item = j, n = sum(ok), pct_missing = round(100 * mean(!ok), 1),
           mean = round(mean(x, na.rm = TRUE), 2), sd = round(sd(x, na.rm = TRUE), 2)),
      setNames(as.list(pct), paste0("pct_", cats))) }), fill = TRUE)
}
desc <- rbindlist(c(list(item_desc(A$scale_raw, "scale", c(1, 6)), item_desc(A$end_raw, "trial statements", c(1, 6))),
                    lapply(names(D$validity), function(k) item_desc(A$val_items[[k]], k, D$validity[[k]]$range))), fill = TRUE)
desc[, role := ifelse(measure == "scale", D$role[item], NA_character_)]
desc[, reversed := ifelse(measure == "scale", item %in% D$reversed, NA)]
fwrite(desc, sf2b_out("02_item_descriptives.csv"))

# Agreement rates on the scale (unrecoded 4-6), as reported for 2a.
agree <- data.table(item = D$items, source = D$source[D$items], role = D$role[D$items],
                    pct_agree = round(100 * colMeans(A$scale_raw[, D$items] >= 4, na.rm = TRUE), 1))
fwrite(agree, sf2b_out("02_agreement_rates.csv"))

# Missing data per measure: participants who reached it, and items unanswered among them.
miss <- rbindlist(c(list(
  data.table(measure = "scale", reached = nrow(A$scale_raw), pct_items_missing = round(100 * mean(is.na(A$scale_raw)), 2)),
  data.table(measure = "trial statements", reached = sum(rowSums(!is.na(A$end_raw)) > 0),
             pct_items_missing = round(100 * mean(is.na(A$end_raw[rowSums(!is.na(A$end_raw)) > 0, ])), 2))),
  lapply(names(A$val_items), function(k) { X <- A$val_items[[k]]; r <- rowSums(!is.na(X)) > 0
    data.table(measure = k, reached = sum(r), pct_items_missing = round(100 * mean(is.na(X[r, , drop = FALSE])), 2)) })))
fwrite(miss, sf2b_out("02_missing_by_measure.csv"))

# Internal consistency of every observed score (alpha on keyed items, complete cases).
alpha_of <- function(X, items, rev, range) {
  Z <- X[, items, drop = FALSE]; rv <- intersect(items, rev); if (length(rv)) Z[, rv] <- sum(range) - Z[, rv]
  Z <- Z[complete.cases(Z), , drop = FALSE]
  if (nrow(Z) < 10 || ncol(Z) < 2) return(NA_real_)
  suppressWarnings(psych::alpha(Z, check.keys = FALSE, warnings = FALSE)$total$raw_alpha)
}
ic <- rbindlist(c(
  list(data.table(score = "sf8", items = 8, alpha = alpha_of(A$scale_raw, D$sf8, D$reversed, c(1, 6))),
       data.table(score = "sf10", items = 10, alpha = alpha_of(A$scale_raw, D$sf10, D$reversed, c(1, 6)))),
  unlist(lapply(names(D$validity), function(k) { v <- D$validity[[k]]
    lapply(names(v$scores), function(s) data.table(score = s, items = length(v$scores[[s]]),
                                                   alpha = alpha_of(A$val_items[[k]], v$scores[[s]], v$reverse, v$range))) }), recursive = FALSE)))
ic[, alpha := r3(alpha)]
S <- A$scores
ic[, n_scored := vapply(score, function(s) sum(!is.na(S[[s]])), 0L)]
fwrite(ic, sf2b_out("02_internal_consistency.csv"))

tab <- function(x) as.list(table(x, useNA = "ifany"))
mins <- A$validity_secs[, lapply(.SD, function(x) round(median(x, na.rm = TRUE) / 60, 2)), .SDcols = -"participant"]
sample_info <- list(
  n_analysable = nrow(M), age_mean = round(mean(M$age, na.rm = TRUE), 1), age_sd = round(sd(M$age, na.rm = TRUE), 1),
  gender = tab(M$gender), country = tab(M$country), practice_group = tab(M$practice_group),
  dedicated = sum(M$dedicated), freq_code = tab(M$freq_code), tells_apart = tab(M$tells_apart),
  presentation_format = tab(M$presentation_format), pause_device = tab(M$pause_device),
  flags = list(one_check = sum(M$flag_one_check), speeder = sum(M$flag_speeder), longstring = sum(M$flag_longstring),
               pna = sum(M$flag_pna), any = sum(M$flagged)),
  minutes_median = list(scale = round(median(M$scale_secs, na.rm = TRUE) / 60, 2), validity = as.list(mins)))
write_json(sample_info, sf2b_out("02_sample.json"), auto_unbox = TRUE, pretty = TRUE)
print(ic); print(miss)
