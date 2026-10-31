# 09: statement selection (Section 5, step 8; decision S). A PROPOSAL for the PI,
# who makes the final selection by these rules and reports every departure.
#
# For each retained factor, five or six statements, at least two of them
# reverse-keyed where at least two qualify. Preferred, in this order:
#   1. primary loading >= .50 (minimum .40)
#   2. no secondary loading >= .30, and primary minus secondary >= .20
#   3. "Prefer not to answer" <= 5%
#   4. agreement rate 15-85% (outside it: flagged; an unflagged statement is
#      preferred to a flagged one)
#   5. at most one statement of each look-alike pair
# A short-form statement stays unless it fails these rules and another statement
# for the same factor does better. Decision V: if the distinction did not come
# across, a statement that names a mode is retained only where no statement
# without one qualifies for the same factor.
#
# PREREG readings, used to make the proposal mechanical:
#   - .40 is the only hard minimum; "qualifies" = primary >= .40 and rules 2-3 met
#   - candidates are ranked: qualifies, then agreement unflagged, then (if V
#     failed) not naming a mode, then primary >= .50, then primary loading
#   - a short-form statement (one carried from Study 1's short form, in its Study 1
#     wording or reworded) that qualifies with an unflagged agreement rate is kept
#     first; up to six are kept by rank
#   - the factor is filled to six with qualifying statements, then to five with
#     any candidate; a look-alike partner of a kept statement is skipped
#   - if fewer than two are reverse-keyed and two or more reverse-keyed statements
#     qualify, the lowest-ranked kept forward statements that are not short-form
#     are swapped for the best qualifying reverse-keyed ones
#   - a retained factor is every content factor of the primary model

source("study2a/analysis/config.R")
source("study2a/analysis/sf2a_models.R")
A <- readRDS(sf2a_derived("sf2a"))
def <- A$definition; k <- SF2A$k_primary
P <- readRDS(sf2a_derived(sprintf("primary_k%02d", k)))
desc <- fread(sf2a_out("02_item_descriptives.csv"))
V <- fread(sf2a_out("03_vignettes_decision.csv"))
came_across <- isTRUE(V$came_across)
Kor <- P$keyed; fn <- P$fnames

tab <- data.table(item = def$items, source = def$source[def$items], hypothesized = def$factor[def$items],
                  keying = ifelse(def$items %in% def$reversed, "reverse", "forward"),
                  short_form = grepl("^sf6_", def$source[def$items]), names_mode = def$items %in% def$names_mode,
                  text = def$text)
tab[, primary_factor := colnames(Kor)[apply(abs(Kor), 1, which.max)]]
tab[, matches := unname(P$labels[primary_factor])]
tab[, primary := round(Kor[cbind(item, primary_factor)], 3)]
tab[, secondary := round(apply(abs(Kor), 1, function(r) sort(r, decreasing = TRUE)[2]), 3)]
tab[, gap := round(primary - secondary, 3)]
tab[, `:=`(missing_pct = desc$missing_pct[match(item, desc$item)], agree_pct = desc$agree_pct[match(item, desc$item)])]
tab[, `:=`(min_ok = primary >= SF2A$load_min, pref_50 = primary >= SF2A$load_pref,
           cross_ok = secondary < SF2A$cross_max & gap >= SF2A$gap_min, pna_ok = missing_pct <= SF2A$pna_max_pct,
           agree_flag = agree_pct < SF2A$agree_lo | agree_pct > SF2A$agree_hi)]
tab[, qualifies := min_ok & cross_ok & pna_ok]
partner <- function(i) { p <- def$pairs[def$pairs[, "a"] == i | def$pairs[, "b"] == i, , drop = FALSE]; setdiff(c(p), i) }
tab[, lookalike := vapply(item, function(i) paste(partner(i), collapse = " "), "")]

pick_factor <- function(f) {
  c0 <- tab[primary_factor == f & min_ok == TRUE]
  if (!nrow(c0)) return(character())
  c0[, mode_ok := came_across | !names_mode]
  setorder(c0, -qualifies, agree_flag, -mode_ok, -pref_50, -primary)
  keep <- character()
  add <- function(i) if (!(i %in% keep) && !any(partner(i) %in% keep)) keep <<- c(keep, i)
  for (i in c0[short_form == TRUE & qualifies == TRUE & agree_flag == FALSE & mode_ok == TRUE, item]) if (length(keep) < 6) add(i)
  for (i in c0[qualifies == TRUE, item]) if (length(keep) < 6) add(i)
  for (i in c0$item) if (length(keep) < 5) add(i)
  rq <- c0[qualifies == TRUE & keying == "reverse", item]
  if (length(rq) >= SF2A$min_reversed) {
    for (r in rq) {
      if (sum(keep %in% def$reversed) >= SF2A$min_reversed) break
      if (r %in% keep) next
      out <- rev(c0[item %in% keep & keying == "forward" & short_form == FALSE, item])
      out <- setdiff(out, partner(r))
      if (any(partner(r) %in% keep) || !length(out)) next
      keep[keep == out[1]] <- r
    }
  }
  c0[item %in% keep, item]
}
sel <- unlist(lapply(fn, pick_factor))
tab[, proposed := item %in% sel]
tab[, note := fifelse(short_form & !proposed & min_ok, "short-form statement not kept", "")]
setorder(tab, primary_factor, -proposed, -primary)
fwrite(tab, sf2a_out(sprintf("09_selection_k%02d.csv", k)))

R <- readRDS(sf2a_derived("polychoric_keyed"))
summ <- rbindlist(lapply(fn, function(f) {
  ids <- tab[proposed == TRUE & primary_factor == f, item]
  n_rq <- tab[primary_factor == f & qualifies == TRUE & keying == "reverse", .N]
  data.table(factor = f, matches = P$labels[[f]], n_selected = length(ids), n_reversed = sum(ids %in% def$reversed),
             reversed_qualifying = n_rq, n_short_form = sum(ids %in% tab[short_form == TRUE, item]),
             n_agree_flagged = tab[item %in% ids, sum(agree_flag)], n_names_mode = sum(ids %in% def$names_mode),
             ordinal_omega = if (length(ids)) round(sf2a_omega(Kor[ids, f]), 3) else NA_real_,
             ordinal_alpha = if (length(ids) >= 2) round(psych::alpha(R[ids, ids], warnings = FALSE)$total$raw_alpha, 3) else NA_real_,
             within_5_6 = length(ids) %in% SF2A$per_factor,
             reversed_rule_met = sum(ids %in% def$reversed) >= min(SF2A$min_reversed, n_rq),
             statements = paste(ids, collapse = " "))
}))
fwrite(summ, sf2a_out(sprintf("09_selection_summary_k%02d.csv", k)))
sf2a_log("proposed: ", length(sel), " statements; decision V came across = ", came_across)
print(summ[, !"statements"])
print(tab[short_form == TRUE, .(item, source, hypothesized, primary_factor, primary, qualifies, agree_flag, proposed)])
