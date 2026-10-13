# 02: descriptives. Item distributions for the short form and the candidates,
# "Prefer not to answer" and missing data per measure (Section 5, "Missing data":
# "The amount of missing data per measure is reported"), the sample, internal
# consistency of every scored measure, and seconds per questionnaire (paradata).

source("study2/analysis/config.R")
source("study2/analysis/sf2_models.R")
A <- readRDS(sf2_derived("sf2"))
D <- A$defs; M <- A$meta; S <- A$scores

item_desc <- function(X, set, factor_of, reversed) {
  d <- data.table(item = colnames(X), set = set, factor = unname(factor_of[colnames(X)]),
                  keying = ifelse(colnames(X) %in% reversed, "reverse", "forward"), n_answered = colSums(!is.na(X)))
  for (k in 1:6) set(d, j = paste0("pct_", k), value = round(100 * colMeans(X == k, na.rm = TRUE), 1))
  set(d, j = "mean", value = round(colMeans(X, na.rm = TRUE), 2)); set(d, j = "sd", value = round(apply(X, 2, sd, na.rm = TRUE), 2))
  set(d, j = "pct_missing", value = round(100 * colMeans(is.na(X)), 1))
  d
}
has_cand <- rowSums(!is.na(A$cand_raw)) > 0
items <- rbind(item_desc(A$sf_raw, "short form", D$sf_factor, D$sf_rev),
               item_desc(A$cand_raw[has_cand, , drop = FALSE], "candidate", setNames(D$cand$factor, D$cand$id), D$cand[reversed == TRUE, id]))
items[, text := c(D$sf_text, setNames(D$cand$text, D$cand$id))[item]]
fwrite(items, sf2_out("02_item_descriptives.csv"))

# Missing data per measure, among those who submitted it (analysable sample).
miss <- rbindlist(c(
  list(data.table(measure = "short form", n_submitted = nrow(A$sf_raw), pct_items_missing = round(100 * mean(is.na(A$sf_raw)), 2),
                  pct_people_any_missing = round(100 * mean(rowSums(is.na(A$sf_raw)) > 0), 1))),
  list(data.table(measure = "candidates", n_submitted = sum(has_cand), pct_items_missing = round(100 * mean(is.na(A$cand_raw[has_cand, ])), 2),
                  pct_people_any_missing = round(100 * mean(rowSums(is.na(A$cand_raw[has_cand, ])) > 0), 1))),
  lapply(names(A$val_items), function(k) { X <- A$val_items[[k]]; h <- rowSums(!is.na(X)) > 0
    data.table(measure = k, n_submitted = sum(h), pct_items_missing = round(100 * mean(is.na(X[h, ])), 2),
               pct_people_any_missing = round(100 * mean(rowSums(is.na(X[h, , drop = FALSE])) > 0), 1)) }),
  list(data.table(measure = "life satisfaction", n_submitted = sum(!is.na(S$lifesat)), pct_items_missing = 0, pct_people_any_missing = 0))))
fwrite(miss, sf2_out("02_missing_by_measure.csv"))
score_cols <- setdiff(names(S), "participant")
miss_scores <- data.table(score = score_cols, n_scored = sapply(score_cols, function(s) sum(!is.na(S[[s]]))),
                          n_validity_sample = sapply(score_cols, function(s) sum(!is.na(S[[s]][M$reached_end_validity]))))
fwrite(miss_scores, sf2_out("02_scores_available.csv"))

# Internal consistency (descriptive): alpha of every multi-item observed score, keyed.
alpha_of <- function(X) { C <- cor(X, use = "pair"); p <- ncol(C); (p / (p - 1)) * (1 - p / sum(C)) }
keyed <- function(X, rev, rg) { X[, intersect(colnames(X), rev)] <- (rg[1] + rg[2]) - X[, intersect(colnames(X), rev)]; X }
rel <- rbindlist(c(
  lapply(c("Access", "Understanding", "Valuing", "LettingBe"), function(f) {
    it <- D$sf_items[D$sf_factor == f]; data.table(score = f, n_items = length(it), alpha = r3(alpha_of(keyed(A$sf_raw[, it], D$sf_rev, c(1, 6))))) }),
  list(data.table(score = "AccessValuing", n_items = 10, alpha = r3(alpha_of(keyed(A$sf_raw[, D$sf_items[D$sf_factor %in% c("Access", "Valuing")]], D$sf_rev, c(1, 6)))))),
  unlist(lapply(names(D$validity), function(k) { v <- D$validity[[k]]
    lapply(names(v$scores), function(s) data.table(score = s, n_items = length(v$scores[[s]]),
      alpha = if (length(v$scores[[s]]) > 1) r3(alpha_of(keyed(A$val_items[[k]][, v$scores[[s]], drop = FALSE], v$reverse, v$range))) else NA_real_)) }), recursive = FALSE)))
fwrite(rel, sf2_out("02_internal_consistency.csv"))

pct <- function(x) { t <- table(x, useNA = "ifany"); as.list(setNames(as.integer(t), ifelse(is.na(names(t)), "missing", names(t)))) }
sample_desc <- list(
  n = nrow(M), n_validity = sum(M$reached_end_validity), posting = pct(M$posting), practice_group = pct(M$practice_group),
  practice_group_by_posting = M[, .N, by = .(posting, practice_group)][order(posting, practice_group)],
  gender = pct(M$gender), country = pct(M$country), country_rule = A$country_rule, education = pct(M$education),
  age = list(median = median(M$age, na.rm = TRUE), mean = r3(mean(M$age, na.rm = TRUE)), sd = r3(sd(M$age, na.rm = TRUE))),
  flags = M[, .(one_check = sum(flag_one_check), speeder = sum(flag_speeder), longstring = sum(flag_longstring),
                pna_over_4 = sum(flag_pna), any = sum(flagged))],
  minutes_short_form = r3(median(M$sf_secs, na.rm = TRUE) / 60), minutes_candidates = r3(median(M$cand_secs, na.rm = TRUE) / 60),
  seconds_per_measure = A$validity_secs[, lapply(.SD, function(v) round(median(v, na.rm = TRUE), 1)), .SDcols = -"participant"])
write_json(sample_desc, sf2_out("02_sample.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
sf2_log("items: ", sum(items$pct_missing > 5), " with > 5% missing; measures: ", nrow(miss))
print(miss)
