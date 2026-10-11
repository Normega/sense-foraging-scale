# 10: every number the Study 1 findings page reports, in one JSON file
# (aggregates only: no participant-level rows, no participant codes). The page on
# senseforaging.com (/research/study1) is built from a copy of this file.

source("study1/analysis/config.R")
source("study1/analysis/sf1_models.R")
A <- readRDS(sf1_path("derived", sprintf("sf1_%s.rds", SF1$mode)))
def <- A$definition; k <- SF1$k_primary
P <- readRDS(sf1_path("derived", sprintf("primary_%s_k%02d.rds", SF1$mode, k)))
rd <- function(f) fread(sf1_out(f))
r2 <- function(x) round(x, 2)

# Loadings, keyed and oriented as in 06.
Lk <- sf1_loadings(P$with_ri, P$fnames)[def$items, ]; Lk[def$reversed, ] <- -Lk[def$reversed, ]
orient <- sf1_orientation(Lk); Lk <- sweep(Lk, 2, orient, `*`)
Ln <- sf1_loadings(P$without_ri, P$fnames)[def$items, ]; Ln[def$reversed, ] <- -Ln[def$reversed, ]
Ln <- sweep(Ln, 2, sf1_orientation(Ln), `*`)
s <- standardizedSolution(P$with_ri); ri_l <- s[s$op == "=~" & s$lhs == "RI", "est.std"][1]
pe <- parameterEstimates(P$with_ri); ri_v <- pe[pe$lhs == "RI" & pe$op == "~~" & pe$rhs == "RI", ]
phi <- function(m) { p <- lavInspect(m, "std")$psi[P$fnames, P$fnames]; p * outer(orient, orient) }
sel <- rd(sprintf("06_item_selection_k%02d.csv", k))
M <- A$meta

pct <- function(x) { t <- sort(table(x, useNA = "no"), decreasing = TRUE); as.list(round(100 * t / sum(t), 1)) }
age_band <- cut(M$age, c(17, 24, 34, 44, 54, 64, Inf), c("18-24", "25-34", "35-44", "45-54", "55-64", "65+"))

out <- list(
  generated = format(Sys.time(), "%Y-%m-%d %H:%M %Z"), mode = SF1$mode, k = k,
  registration = "https://osf.io/48r2y",
  flow = rd("01_flow.csv"),
  flags = rd("01_flags.csv")[excluded_attention == FALSE],
  sample = list(n = nrow(M), posting = as.list(table(M$posting)), practice = as.list(table(M$practice_group)),
                age_median = median(M$age, na.rm = TRUE), age_bands = pct(age_band), gender = pct(M$gender),
                country = pct(M$country), education = pct(M$education),
                median_minutes_scale = r2(median(M$scale_secs, na.rm = TRUE) / 60)),
  parallel = rd("03_parallel_analysis.csv")[factor <= 10],
  models = rd("03_factor_number_models.csv"),
  primary_fit = rd(sprintf("04_primary_fit_k%02d.csv", k)),
  acquiescence = list(variance = r2(ri_v$est), se = round(ri_v$se, 3), std_loading = r2(ri_l),
                      pct_item_variance = round(100 * ri_l^2, 1)),
  items = lapply(seq_along(def$items), function(i) {
    id <- def$items[i]; row <- sel[item == id]
    list(id = id, n = i, text = unname(def$text[i]), keying = if (id %in% def$reversed) "reverse" else "positive",
         hypothesized = unname(def$factor[i]), loadings = unname(r2(Lk[id, ])), loadings_no_ri = unname(r2(Ln[id, ])),
         primary = row$primary_factor, eligible = isTRUE(row$eligible), proposed = isTRUE(row$proposed),
         note = if (is.null(row$note) || is.na(row$note)) "" else row$note)
  }),
  factors = P$fnames,
  factor_correlations = list(with_ri = r2(phi(P$with_ri)), without_ri = r2(phi(P$without_ri))),
  omega_pool = rd(sprintf("05_omega_pool_k%02d.csv", k)),
  short_form = rd(sprintf("06_short_form_reliability_k%02d.csv", k)),
  reversal_pairs = rd(sprintf("06_reversal_pairs_k%02d.csv", k)),
  known_groups = rd(sprintf("07_known_groups_k%02d.csv", k)),
  posting_means = rd(sprintf("07_posting_means_k%02d.csv", k)),
  sensitivity = rd(sprintf("08_sensitivity_k%02d.csv", k)),
  position_median_abs_r = round(median(abs(rd("08_position_effects.csv")$r_position), na.rm = TRUE), 3),
  near_duplicates = rd("02_near_duplicates.csv")
)
# Proposed items in the short form also include partners added by the pair rule.
pairs_prop <- out$reversal_pairs[proposed == TRUE]
short <- unique(c(sel[proposed == TRUE, item], pairs_prop$reversed, pairs_prop$partner))
out$short_form_items <- short
write_json(out, sf1_out("10_report_study1.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
sf1_log("report data: ", sf1_out("10_report_study1.json"), " (", length(short), " short-form items)")
