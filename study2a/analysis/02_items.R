# 02: pre-checks (Section 5, step 1): distributions, floor and ceiling, "Prefer
# not to answer" rates and agreement rates for every statement, and the
# polychoric correlation within each of the 18 look-alike pairs, screened before
# the agreement factor is fitted (a random intercept can absorb doublets).
# Agreement rate = the share of answering participants choosing 4-6, on the
# unrecoded answer (Section 4, "Indices").

source("study2a/analysis/config.R")
A <- readRDS(sf2a_derived("sf2a"))
def <- A$definition; X <- A$items_raw

desc <- data.table(item = def$items, source = def$source[def$items], factor = def$factor[def$items],
                   keying = ifelse(def$items %in% def$reversed, "reverse", "forward"),
                   names_mode = def$items %in% def$names_mode, n = colSums(!is.na(X)), text = def$text)
for (k in 1:6) set(desc, j = paste0("pct_", k), value = round(100 * colMeans(X == k, na.rm = TRUE), 1))
desc[, `:=`(mean = round(colMeans(X, na.rm = TRUE), 2), sd = round(apply(X, 2, sd, na.rm = TRUE), 2),
            missing_pct = round(100 * colMeans(is.na(X)), 1),
            agree_pct = round(100 * colMeans(X >= 4, na.rm = TRUE), 1))]
desc[, `:=`(floor_pct = pct_1, ceiling_pct = pct_6, flag_missing_over_5pct = missing_pct > SF2A$pna_max_pct,
            flag_agreement = agree_pct < SF2A$agree_lo | agree_pct > SF2A$agree_hi)]
fwrite(desc, sf2a_out("02_item_descriptives.csv"))

pairs <- rbindlist(lapply(seq_len(nrow(def$pairs)), function(i) {
  a <- def$pairs[i, "a"]; b <- def$pairs[i, "b"]
  r <- polychoric(X[, c(a, b)], correct = 0)$rho[1, 2]
  data.table(a = a, b = b, source_a = def$source[a], source_b = def$source[b],
             same_keying = (a %in% def$reversed) == (b %in% def$reversed), polychoric_r = round(r, 3))
}))
fwrite(pairs, sf2a_out("02_lookalike_pairs.csv"))

sf2a_log("statements: missing > 5% on ", sum(desc$flag_missing_over_5pct), "; agreement outside 15-85% on ",
         sum(desc$flag_agreement), "; participants with any statement missing: ", sum(rowSums(is.na(X)) > 0))
print(pairs)
