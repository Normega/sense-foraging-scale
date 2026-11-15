# 11: what predicts taking the in-session pause? EXPLORATORY. NOT REGISTERED.
#
# Not part of the filed registration (osf.io/d4a5u, filed 14 November 2026) and not one
# of its listed exploratory analyses. Specified on 15 November 2026, after H4b's null
# result was known (05_h4_decision.json: rho of the 8 with the engaged pause -.022), so
# every number here is hypothesis-generating only. No registered cell, test or output
# is changed, and nothing here may be reported as a test of H4.
#
# Sample: the 720 analysable participants; the pause analyses use the 710 who reached the
# offer (step 13). Outcomes: took the pause (0/1, 73 events) and the engaged pause
# (0/1, the registered definition, 44 events). With so few events every estimate is
# imprecise and the multivariable model is kept small.
#
#   1. Who takes the pause, and who has an engaged one?
#      (a) Bivariate screen of every available predictor against each outcome: Spearman
#          rho for continuous and ordinal predictors (with an odds ratio per SD from a
#          logistic regression, log-transformed first where the predictor is a duration),
#          Fisher's exact test and its conditional odds ratio for binary ones, Fisher's
#          exact test (Monte Carlo p, 20,000 tables) for categorical ones, whose
#          per-level rates go to 11_pause_levels.csv. The 8's latent correlation with
#          each outcome comes from 04's single-indicator model (the outcome's latent
#          response). Holm and Benjamini-Hochberg corrections within each outcome,
#          across the whole screen.
#      (b) One Firth-penalized logistic regression per outcome (logistf). Its predictors
#          were fixed before the screen was run: the 8 (z), tiredness (z), mobile device
#          (phone or tablet vs desktop), log minutes from session start to the offer (z),
#          current practice (yes vs no) and age (z): 6 df, about 12 events per df for
#          took and 7 for engaged. Complete cases. Discrimination: the apparent AUC and
#          the AUC under stratified 10-fold cross-validation repeated 50 times (seeded),
#          reported as the mean and the 2.5th-97.5th percentile over repeats (the spread
#          over fold splits, not a confidence interval).
#   2. Among those who paused: on-page seconds, how the pause ended, and whether a note
#      was written. Descriptives; Spearman rho of on-page seconds and of the note (0/1)
#      with the 8 and the main scores (Holm within each family).
#   3. Decline reasons among decliners: the distribution, and the 8 by reason
#      (Kruskal-Wallis). The registered exploratory contrast of "not interested" with
#      "no time / not a good place" is 08's and is not redone here.
#   4. Is H4b's null a measurement problem? The 8 across declined / took but not engaged /
#      engaged; dose (the 8 against "paused at least X seconds on the page" for several
#      X, and against on-page seconds with decliners at 0); and four interaction checks,
#      the 8 x mobile device and the 8 x tiredness, on each outcome (Firth).
#
# Predictors, and when each was measured relative to the offer (the "timing" column):
#   session context  device class (recorded at the offer); local time of day and day of
#                    week of session start; minutes from session start to the offer;
#                    seconds per item on the questionnaires before it (the scale and the
#                    first half of the validity block); the side the pause button was on
#                    (randomized); the four flags. Of the three attention checks, two
#                    come before the offer and the BFI-2 check after it, so "failed one
#                    check" can be measured after.
#   self-reports     setting and tiredness, asked after the decision (by design), so a
#                    decision could colour them.
#   background       age, gender, country group (the six countries vs other), current
#                    practice, practice frequency (0-4).
#   scores           the 8 (observed and latent); past-week pause days (after the offer);
#                    every validity score. The first half of the validity block (FFMQ,
#                    BIDR, SBI, BSSS-8, NCS-6, PHQ-4, SPSQ-SF, PSS-4, life satisfaction)
#                    precedes the offer; the second half (Brief MAIA-2, BFI-2, EQ, SAMD,
#                    ASI-3, the past-week items) follows it, so those could be affected
#                    by having paused.
#   order            the halves are fixed (Section 4); the order within each is random
#                    per participant. The order predictor is the measure shown last
#                    before the offer (position 9 of the first half).
#
# Data. The pipeline's saved real-mode objects (derived/sf2b_real.rds). Three facts the
# derived object does not hold are read from the same extraction file that 01_prepare
# read (named in the object; not a new extraction, no database query), and only these:
# each step's entry time (session start = the earliest entry, as in decision 45; the
# offer = step 13's entry) and the pause record's local_hour (the browser's local hour).
# The local UTC offset is inferred to the whole hour from local_hour against the
# offer's server time, and session start's local time from it (a half-hour time zone is
# off by 30 minutes). Start day and hour mostly mark the recruitment batch, not a
# participant's choice: the pilot ran on a weekday and the main wave opened at the
# weekend, each filling within hours. Nothing from the extraction is written out. Prolific IDs are not in
# the extraction (extract.sql) and are never read.
#
# Output (aggregates only; no participant rows or codes): output/real/11_pause_*.csv and
# 11_pause.json.
#
# Run from the repo root:  SF2B_MODE=real Rscript study2b/analysis/11_exploratory_pause.R

# Set Up ---------
## Load libraries ---------
packages <- c("jsonlite", "data.table", "lavaan", "logistf")
new_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(new_packages)) install.packages(new_packages)
for (thispack in packages) {
  library(thispack, character.only = TRUE, quietly = TRUE, verbose = FALSE)
}

## Settings ---------
source("study2b/analysis/config.R")
source("study2b/analysis/sf2b_models.R")
stopifnot(SF2B$mode == "real")
SEED <- 20261115
CV_FOLDS <- 10; CV_REPEATS <- 50; MC_TABLES <- 20000
DOSE_SECS <- c(10, 20, 60, 120)
A <- readRDS(sf2b_derived("sf2b")); D <- A$defs; M <- A$meta; S <- A$scores
stopifnot(identical(M$participant, S$participant), identical(M$participant, A$validity_positions$participant),
          identical(M$participant, A$validity_secs$participant))
out <- list(status = "EXPLORATORY, NOT REGISTERED (specified 15 November 2026, after H4b's result was known)",
            mode = SF2B$tag, extraction = A$extraction)

# Session timing from the extraction (three fields only) ---------
ts_parse <- function(x) {   # as 01_prepare.R
  x <- unlist(lapply(x, function(v) if (is.null(v)) NA_character_ else v))
  off <- ifelse(grepl("[+-][0-9]{2}(:?[0-9]{2})?$", x), sub(".*([+-][0-9]{2})(:?([0-9]{2}))?$", "\\1\\3", x), "+00")
  base <- sub("[+-][0-9]{2}(:?[0-9]{2})?$", "", sub("Z$", "", sub("T", " ", x)))
  t <- as.POSIXct(base, tz = "UTC", format = "%Y-%m-%d %H:%M:%OS")
  h <- as.numeric(substr(off, 1, 3)); m <- suppressWarnings(as.numeric(substr(off, 4, 5))); m[is.na(m)] <- 0
  t - sign(h) * (abs(h) * 3600 + m * 60)
}
raw_file <- sf2b_path("raw", A$extraction)
stopifnot(file.exists(raw_file))
TM <- local({
  rows <- fromJSON(raw_file, simplifyVector = FALSE)$rows
  rows <- Filter(function(r) isTRUE(r$participant %in% M$participant) &&
                   (identical(r$source, "step") || (identical(r$source, "questionnaire") && identical(r$slug, "sf2b-pause"))), rows)
  st <- rbindlist(lapply(Filter(function(r) identical(r$source, "step"), rows), function(r)
    data.table(participant = r$participant, step = as.integer(r$step_index %||% NA), entered = r$entered_at %||% NA_character_)))
  st[, entered := ts_parse(as.list(entered))]
  pz <- Filter(function(r) is.null(r$resubmission_of) && !is.null(r$responses$pause_1), rows)
  pz <- rbindlist(lapply(pz, function(r) data.table(participant = r$participant, received = r$received_at %||% NA_character_,
                                                     local_hour = as.numeric(r$responses$pause_1$local_hour %||% NA))))
  pz[, received := ts_parse(as.list(received))]
  pz <- pz[order(participant, received)][!duplicated(participant)]
  tm <- merge(st[, .(start = min(entered, na.rm = TRUE)), by = participant],
              st[step == SF2B_STEPS$pause, .(offer = min(entered, na.rm = TRUE)), by = participant], by = "participant", all = TRUE)
  merge(tm, pz[, .(participant, local_hour)], by = "participant", all = TRUE)
})
TM <- TM[match(M$participant, TM$participant)]
utc_frac <- as.numeric(format(TM$offer, "%H", tz = "UTC")) + as.numeric(format(TM$offer, "%M", tz = "UTC")) / 60
utc_off <- round(TM$local_hour + .5 - utc_frac)
utc_off <- ((utc_off + 10) %% 24) - 10                      # whole-hour offsets, UTC-10 to UTC+13
start_local <- as.POSIXlt(TM$start + utc_off * 3600, tz = "UTC")
start_hour <- start_local$hour
out$timing_source <- list(with_start_and_offer = sum(!is.na(TM$start) & !is.na(TM$offer) & M$has_pause),
                          with_local_hour = sum(!is.na(TM$local_hour) & M$has_pause),
                          utc_offset_regions = as.list(table(cut(utc_off[M$has_pause], c(-Inf, -1.5, 3.5, Inf),
                                labels = c("UTC-10 to -2", "UTC-1 to +3", "UTC+4 to +13")), useNA = "ifany")))

# The analysis frame: the 710 offered ---------
half_a_slugs <- c(ffmq15 = "sf2-ffmq15", bidr16 = "sf2-bidr16", sbi = "sf2-sbi-moment", bsss = "sf2b-bsss8",
                  ncs = "sf2b-ncs6", phq4 = "sf2-phq4", spsq = "sf2b-spsq-sd", pss4 = "sf2-pss4")
half_a_items <- 14 + sum(vapply(names(half_a_slugs), function(k) length(D$validity[[k]]$items), 0))
pre_secs <- M$scale_secs + rowSums(as.matrix(A$validity_secs[, ..half_a_slugs]))
VP <- A$validity_positions
half_a_all <- c(half_a_slugs, lifesat = "lifesat")
last_before <- apply(as.matrix(VP[, ..half_a_all]), 1, function(p) { k <- which(p == 9); if (length(k) == 1) names(half_a_all)[k] else NA_character_ })
last_before <- c(ffmq15 = "FFMQ-15", bidr16 = "BIDR-16", sbi = "SBI", bsss = "BSSS-8", ncs = "NCS-6", phq4 = "PHQ-4",
                 spsq = "SPSQ-SF", pss4 = "PSS-4", lifesat = "life satisfaction")[last_before]
period <- cut(start_hour, c(-1, 5, 11, 17, 23), labels = c("night 0-5", "morning 6-11", "afternoon 12-17", "evening 18-23"))

FR <- data.table(
  took = M$pause_taken, engaged = M$engaged_pause, offered = M$has_pause,
  device = M$pause_device, mobile = as.integer(M$pause_device %in% c("phone", "tablet")),
  start_period = period, weekend = as.integer(start_local$wday %in% c(0, 6)),
  weekday = factor(c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")[start_local$wday + 1],
                   levels = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")),
  minutes_before = as.numeric(difftime(TM$offer, TM$start, units = "mins")),
  secs_per_item_before = pre_secs / half_a_items,
  button_left = as.integer(M$pause_button_side == "left"),
  flag_one_check = as.integer(M$flag_one_check), flag_speeder = as.integer(M$flag_speeder),
  flag_longstring = as.integer(M$flag_longstring), flag_pna = as.integer(M$flag_pna),
  tired = M$pause_tired, setting = M$pause_setting,
  setting_home = fifelse(M$pause_setting %in% c("home", "work", "outdoors", "transit", "other"), as.integer(M$pause_setting == "home"), NA_integer_),
  age = M$age, gender3 = factor(M$gender3, levels = c("man", "woman", "other")),
  six_countries = fifelse(is.na(M$country2), NA_integer_, as.integer(M$country2 == "six_countries")),
  practice_current = fifelse(is.na(M$practice_group), NA_integer_, as.integer(M$practice_group == "current")),
  freq_code = M$freq_code, practice_group3 = M$practice_group, last_before = factor(last_before),
  reason = M$pause_decline_reason, on_page_s = M$pause_on_page_s, total_s = M$pause_total_s,
  away_s = M$pause_unfocused_s + M$pause_hidden_s, hidden_s = M$pause_hidden_s, ended_by = M$pause_ended_by,
  nudged = M$pause_nudged, note = as.integer(M$pause_has_note), decision_s = M$pause_decision_s)
SC <- c("sf8", "crit_pause_days", "ffmq_observing", "ffmq_describing", "ffmq_acting", "ffmq_nonjudging", "ffmq_nonreactivity",
        "sbi_moment", "bsss8", "ncs6", "spsq_sd", "phq4_total", "pss4", "bidr_sde", "bidr_im", "lifesat",
        "maia_noticing", "maia_notdistracting", "maia_notworrying", "maia_attnregulation", "maia_emoawareness",
        "maia_selfregulation", "maia_bodylistening", "maia_trusting", "bfi_aesthetic", "bfi_intellectual", "bfi_creative",
        "eq_decentering", "samd_md", "asi_pc")
for (s in SC) FR[[s]] <- S[[s]]
O <- FR[offered == TRUE]
stopifnot(nrow(O) == 710, sum(O$took) == 73, sum(O$engaged) == 44)
sf2b_log("offered ", nrow(O), ", took ", sum(O$took), ", engaged ", sum(O$engaged))

# Helpers ---------
z <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
auc <- function(p, y) { r <- rank(p); n1 <- sum(y == 1); n0 <- sum(y == 0); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
rnd <- function(x) if (is.data.frame(x)) x else if (is.list(x)) lapply(x, rnd) else if (is.numeric(x)) signif(x, 4) else x
firth <- function(form, d, pl = TRUE) suppressWarnings(logistf(form, data = d, pl = pl, plconf = NULL))
firth_tab <- function(f) {
  data.table(term = names(coef(f)), odds_ratio = exp(coef(f)), lower = exp(f$ci.lower), upper = exp(f$ci.upper), p = f$prob)
}

# 1a. Base rates ---------
base <- rbind(
  data.table(group = "all offered", n = nrow(O), took = sum(O$took), engaged = sum(O$engaged)),
  O[, .(n = .N, took = sum(took), engaged = sum(engaged)), by = .(group = paste("device:", device))][order(group)])
base[, `:=`(took_pct = round(100 * took / n, 1), engaged_pct = round(100 * engaged / n, 1),
            engaged_pct_of_takers = round(100 * engaged / pmax(took, 1), 1))]
base <- rbind(data.table(group = "analysable", n = nrow(FR), took = NA, engaged = NA), base, fill = TRUE)
fwrite(base, sf2b_out("11_pause_base_rates.csv"))
wil <- function(k, n) prop.test(k, n, correct = FALSE)$conf.int
out$base_rates <- list(analysable = nrow(FR), offered = nrow(O), took = sum(O$took), engaged = sum(O$engaged),
                       took_rate = mean(O$took), took_rate_ci_wilson = wil(sum(O$took), nrow(O)),
                       engaged_rate = mean(O$engaged), engaged_rate_ci_wilson = wil(sum(O$engaged), nrow(O)),
                       engaged_share_of_takers = sum(O$engaged) / sum(O$took))

# 1a. Bivariate screen ---------
PRED <- rbindlist(list(
  list("device", "session", "at offer", "cat", FALSE),
  list("start_period", "session", "before", "cat", FALSE),
  list("weekend", "session (in practice: main wave vs pilot)", "before", "bin", FALSE),
  list("minutes_before", "session", "before", "cont", TRUE),
  list("secs_per_item_before", "session", "before", "cont", TRUE),
  list("button_left", "session", "at offer (randomized)", "bin", FALSE),
  list("flag_one_check", "session", "before or after (BFI-2 check after)", "bin", FALSE),
  list("flag_speeder", "session", "before", "bin", FALSE),
  list("flag_longstring", "session", "before", "bin", FALSE),
  list("flag_pna", "session", "before", "bin", FALSE),
  list("tired", "self-report at pause", "after decision", "cont", FALSE),
  list("setting_home", "self-report at pause", "after decision", "bin", FALSE),
  list("age", "background", "before", "cont", FALSE),
  list("gender3", "background", "before", "cat", FALSE),
  list("six_countries", "background", "before", "bin", FALSE),
  list("practice_current", "background", "before", "bin", FALSE),
  list("freq_code", "background", "before", "cont", FALSE),
  list("sf8", "scores", "before", "cont", FALSE),
  list("sf8_latent", "scores", "before", "latent", FALSE),
  list("crit_pause_days", "scores", "after", "cont", FALSE),
  list("ffmq_observing", "scores", "before", "cont", FALSE), list("ffmq_describing", "scores", "before", "cont", FALSE),
  list("ffmq_acting", "scores", "before", "cont", FALSE), list("ffmq_nonjudging", "scores", "before", "cont", FALSE),
  list("ffmq_nonreactivity", "scores", "before", "cont", FALSE), list("sbi_moment", "scores", "before", "cont", FALSE),
  list("bsss8", "scores", "before", "cont", FALSE), list("ncs6", "scores", "before", "cont", FALSE),
  list("spsq_sd", "scores", "before", "cont", FALSE), list("phq4_total", "scores", "before", "cont", FALSE),
  list("pss4", "scores", "before", "cont", FALSE), list("bidr_sde", "scores", "before", "cont", FALSE),
  list("bidr_im", "scores", "before", "cont", FALSE), list("lifesat", "scores", "before", "cont", FALSE),
  list("maia_noticing", "scores", "after", "cont", FALSE), list("maia_notdistracting", "scores", "after", "cont", FALSE),
  list("maia_notworrying", "scores", "after", "cont", FALSE), list("maia_attnregulation", "scores", "after", "cont", FALSE),
  list("maia_emoawareness", "scores", "after", "cont", FALSE), list("maia_selfregulation", "scores", "after", "cont", FALSE),
  list("maia_bodylistening", "scores", "after", "cont", FALSE), list("maia_trusting", "scores", "after", "cont", FALSE),
  list("bfi_aesthetic", "scores", "after", "cont", FALSE), list("bfi_intellectual", "scores", "after", "cont", FALSE),
  list("bfi_creative", "scores", "after", "cont", FALSE), list("eq_decentering", "scores", "after", "cont", FALSE),
  list("samd_md", "scores", "after", "cont", FALSE), list("asi_pc", "scores", "after", "cont", FALSE),
  list("last_before", "order", "before", "cat", FALSE)))
setnames(PRED, c("predictor", "domain", "timing", "type", "log_for_or"))

# The 8's latent correlation with each outcome: 04's single-indicator model, refitted.
latent_r <- function(v) {
  DAT <- as.data.frame(A$scale_raw); DAT[[v]] <- FR[[v]]
  syn <- c(sf2b_scale_syntax(D$sf8), paste0("L_", v, " =~ 1*", v), paste0(v, " ~~ 0*", v), paste0("L_", v, " ~~ 1*L_", v),
           paste0("S ~~ r_", v, "*L_", v), paste0("RI ~~ 0*L_", v))
  f <- sf2b_try(sf2b_sem(syn, DAT, ordered = c(D$sf8, v)))
  if (sf2b_failed(f)) return(list(n = sum(!is.na(FR[[v]])), r = NA, lower = NA, upper = NA, p = NA, note = sf2b_msg(f)))
  pe <- parameterEstimates(f); e <- pe[pe$label == paste0("r_", v), ]
  list(n = sum(!is.na(FR[[v]])), r = e$est, lower = e$ci.lower, upper = e$ci.upper, p = e$pvalue, note = "")
}
LAT <- list(took = latent_r("took"), engaged = latent_r("engaged"))
out$sf8_latent <- LAT

LEVELS <- list()
screen_one <- function(pr, outcome) {
  y <- O[[outcome]]; pp <- PRED[PRED$predictor == pr]
  row <- data.table(outcome = outcome, predictor = pr, domain = pp$domain, timing = pp$timing, type = pp$type,
                    n = NA_integer_, rho = NA_real_, rho_lower = NA_real_, rho_upper = NA_real_,
                    odds_ratio = NA_real_, or_lower = NA_real_, or_upper = NA_real_, or_unit = NA_character_,
                    rate_if_1 = NA_real_, rate_if_0 = NA_real_, n_x1 = NA_integer_, p = NA_real_, test = NA_character_)
  if (pp$type == "latent") {
    l <- LAT[[outcome]]
    row[, `:=`(n = l$n, rho = l$r, rho_lower = l$lower, rho_upper = l$upper, p = l$p,
               test = "latent r (WLSMV, single-indicator outcome), 04's model")]
    return(row)
  }
  x <- O[[pr]]; ok <- !is.na(x) & !is.na(y); row$n <- sum(ok)
  if (pp$type == "cont") {
    r <- sf2b_r_ci(x, y, "spearman")
    xx <- if (pp$log_for_or) log(x) else x
    g <- glm(y ~ zx, family = binomial, data = data.frame(y = y[ok], zx = z(xx[ok])))
    co <- summary(g)$coefficients["zx", ]
    row[, `:=`(rho = r$r, rho_lower = r$lower, rho_upper = r$upper, p = r$p, test = "Spearman",
               odds_ratio = exp(co[["Estimate"]]), or_lower = exp(co[["Estimate"]] - 1.96 * co[["Std. Error"]]),
               or_upper = exp(co[["Estimate"]] + 1.96 * co[["Std. Error"]]),
               or_unit = if (pp$log_for_or) "per SD of log" else "per SD")]
  } else if (pp$type == "bin") {
    ft <- fisher.test(table(factor(x[ok], 0:1), factor(y[ok], 0:1)))
    row[, `:=`(odds_ratio = unname(ft$estimate), or_lower = ft$conf.int[1], or_upper = ft$conf.int[2], or_unit = "1 vs 0 (conditional MLE)",
               rate_if_1 = mean(y[ok][x[ok] == 1]), rate_if_0 = mean(y[ok][x[ok] == 0]), n_x1 = sum(x[ok] == 1), p = ft$p.value,
               test = "Fisher exact")]
  } else {
    xf <- droplevels(factor(x[ok])); set.seed(SEED)
    tb <- table(xf, factor(y[ok], 0:1))
    ft <- fisher.test(tb, simulate.p.value = TRUE, B = MC_TABLES)
    row[, `:=`(p = ft$p.value, test = sprintf("Fisher exact, Monte Carlo (%d tables), %d levels", MC_TABLES, nlevels(xf)))]
    ref <- levels(xf)[which.max(table(xf))]
    LEVELS[[length(LEVELS) + 1L]] <<- rbindlist(lapply(levels(xf), function(l) {
      k <- xf == l; e <- sum(y[ok][k])
      o <- if (l == ref) list(estimate = NA, conf.int = c(NA, NA)) else
        fisher.test(table(factor(as.integer(xf[xf %in% c(l, ref)] == l), 0:1), factor(y[ok][xf %in% c(l, ref)], 0:1)))
      data.table(outcome = outcome, predictor = pr, level = l, n = sum(k), events = e, rate = e / sum(k),
                 reference = ref, odds_ratio_vs_reference = unname(o$estimate), lower = o$conf.int[1], upper = o$conf.int[2])
    }))
  }
  row
}
BIV <- rbindlist(lapply(c("took", "engaged"), function(o) rbindlist(lapply(PRED$predictor, screen_one, outcome = o))))
BIV[, `:=`(p_holm = p.adjust(p, "holm"), p_fdr = p.adjust(p, "BH")), by = outcome]
BIV[, `:=`(survives_holm = p_holm < .05, survives_fdr = p_fdr < .05)]
num_cols <- c("rho", "rho_lower", "rho_upper", "odds_ratio", "or_lower", "or_upper", "rate_if_1", "rate_if_0")
BIV[, (num_cols) := lapply(.SD, r3), .SDcols = num_cols]
BIV[, (c("p", "p_holm", "p_fdr")) := lapply(.SD, signif, 3), .SDcols = c("p", "p_holm", "p_fdr")]
fwrite(BIV, sf2b_out("11_pause_bivariate.csv"))
LEVELS <- rbindlist(LEVELS)
# Descriptive only: day of week (the screen carries weekend vs weekday), and practice
# as three groups and by frequency (the screen carries current vs not, and frequency 0-4).
desc_rates <- function(v, label) rbindlist(lapply(c("took", "engaged"), function(o)
  O[!is.na(get(v)), .(outcome = o, predictor = label, n = .N, events = sum(get(o)), rate = mean(get(o))), by = .(level = as.character(get(v)))][order(level)]))
LEVELS <- rbind(LEVELS, desc_rates("weekday", "weekday (descriptive)"), desc_rates("practice_group3", "practice group (descriptive)"),
                desc_rates("freq_code", "practice frequency 0-4 (descriptive)"), fill = TRUE)
LEVELS[, (c("rate", "odds_ratio_vs_reference", "lower", "upper")) := lapply(.SD, r3), .SDcols = c("rate", "odds_ratio_vs_reference", "lower", "upper")]
fwrite(LEVELS, sf2b_out("11_pause_levels.csv"))
out$bivariate <- list(tests_per_outcome = nrow(PRED),
                      nominal_p05 = BIV[p < .05, .(outcome, predictor, rho, odds_ratio, p, p_holm, p_fdr)],
                      survive_holm = BIV[survives_holm == TRUE, .(outcome, predictor)],
                      survive_fdr = BIV[survives_fdr == TRUE, .(outcome, predictor)])
# The one predictor that survives (current practice): within the main wave alone, so
# the recruitment batch cannot carry it.
out$practice_within_main_wave <- rbindlist(lapply(c("took", "engaged"), function(o) {
  d <- O[weekend == 1 & !is.na(practice_current)]; ft <- fisher.test(table(factor(d$practice_current, 0:1), factor(d[[o]], 0:1)))
  data.table(outcome = o, n = nrow(d), rate_current = mean(d[[o]][d$practice_current == 1]), rate_not = mean(d[[o]][d$practice_current == 0]),
             odds_ratio = unname(ft$estimate), lower = ft$conf.int[1], upper = ft$conf.int[2], p = ft$p.value) }))
sf2b_log("screen: ", BIV[, sum(p < .05)], " nominal p < .05 of ", nrow(BIV), "; Holm survivors ", BIV[, sum(survives_holm)],
         "; FDR survivors ", BIV[, sum(survives_fdr)])

# 1b. The multivariable model (predictors fixed in advance) ---------
MV <- c("z_sf8", "z_tired", "mobile", "z_log_minutes", "practice_current", "z_age")
MD <- O[, .(took, engaged, z_sf8 = z(sf8), z_tired = z(tired), mobile, z_log_minutes = z(log(minutes_before)),
            practice_current, z_age = z(age))]
MD <- MD[complete.cases(MD)]
cv_auc <- function(outcome) {
  d <- as.data.frame(MD); y <- d[[outcome]]; form <- reformulate(MV, outcome)
  set.seed(SEED)
  reps <- vapply(seq_len(CV_REPEATS), function(r) {
    fold <- integer(length(y))
    for (cl in 0:1) { i <- which(y == cl); fold[i] <- sample(rep_len(seq_len(CV_FOLDS), length(i))) }
    pred <- numeric(length(y))
    for (k in seq_len(CV_FOLDS)) {
      f <- firth(form, d[fold != k, ], pl = FALSE)
      pred[fold == k] <- as.numeric(as.matrix(cbind(1, d[fold == k, MV])) %*% coef(f))
    }
    auc(pred, y)
  }, 0)
  f <- firth(form, d)
  list(fit = f, apparent = auc(f$linear.predictors, y), cv = reps)
}
MOD <- list(took = cv_auc("took"), engaged = cv_auc("engaged"))
MODT <- rbindlist(lapply(names(MOD), function(o) cbind(outcome = o, firth_tab(MOD[[o]]$fit))))
MODT[, (c("odds_ratio", "lower", "upper")) := lapply(.SD, r3), .SDcols = c("odds_ratio", "lower", "upper")][, p := signif(p, 3)]
fwrite(MODT, sf2b_out("11_pause_model.csv"))
AUCT <- rbindlist(lapply(names(MOD), function(o) data.table(outcome = o, n = nrow(MD), events = sum(MD[[o]]), predictors_df = length(MV),
  events_per_df = round(sum(MD[[o]]) / length(MV), 1), auc_apparent = r3(MOD[[o]]$apparent), auc_cv_mean = r3(mean(MOD[[o]]$cv)),
  auc_cv_p025 = r3(quantile(MOD[[o]]$cv, .025)), auc_cv_p975 = r3(quantile(MOD[[o]]$cv, .975)),
  cv = sprintf("stratified %d-fold, %d repeats, Firth refitted in every fold", CV_FOLDS, CV_REPEATS))))
fwrite(AUCT, sf2b_out("11_pause_model_auc.csv"))
out$model <- list(predictors = MV, estimator = "Firth-penalized logistic regression (logistf), profile-likelihood CIs",
                  complete_cases = nrow(MD), coefficients = MODT, auc = AUCT)
sf2b_log("model AUC (CV mean): took ", AUCT$auc_cv_mean[1], ", engaged ", AUCT$auc_cv_mean[2])

# 2. Among those who paused ---------
TK <- O[took == 1]
q <- function(x) as.list(setNames(round(quantile(x, c(0, .25, .5, .75, 1), na.rm = TRUE), 1), c("min", "q25", "median", "q75", "max")))
DUR <- rbindlist(lapply(c("total_s", "on_page_s", "away_s", "hidden_s", "decision_s"), function(v)
  c(list(variable = v, n = sum(!is.na(TK[[v]])), mean = round(mean(TK[[v]], na.rm = TRUE), 1)), q(TK[[v]]))))
fwrite(DUR, sf2b_out("11_pause_duration.csv"))
short <- TK$on_page_s < SF2B$pause_min_secs; away <- TK$away_s >= SF2B$pause_max_away * TK$total_s
out$paused <- list(
  n = nrow(TK), duration = DUR,
  engaged = sum(TK$engaged), not_engaged = sum(TK$engaged == 0),
  not_engaged_why = list(under_20s_on_page_only = sum(short & !away), away_half_or_more_only = sum(!short & away), both = sum(short & away)),
  any_time_away = sum(TK$away_s > 0), on_page_ge_60s = sum(TK$on_page_s >= 60), on_page_ge_180s = sum(TK$on_page_s >= 180),
  ended_by = as.list(table(TK$ended_by, useNA = "ifany")), nudged_after_3_min = sum(TK$nudged %in% TRUE),
  wrote_note = sum(TK$note), no_note = sum(TK$note == 0),
  on_page_by_device = TK[, .(n = .N, median_on_page_s = round(median(on_page_s), 1), engaged = sum(engaged)), by = device][order(device)],
  on_page_device_kruskal_p = kruskal.test(on_page_s ~ droplevels(device), data = TK)$p.value,
  note_by_device = TK[, .(n = .N, notes = sum(note)), by = device][order(device)],
  note_by_engaged = TK[, .(n = .N, notes = sum(note), median_on_page_s = round(median(on_page_s), 1)), by = engaged][order(engaged)],
  ended_by_note = "72 of 73 ended with 'I'm back' and one at the 10-minute cap: too few to relate to anything; descriptive only")
Q2V <- c("sf8", "ffmq_observing", "sbi_moment", "bfi_aesthetic", "maia_noticing", "crit_pause_days", "tired", "age",
         "minutes_before", "decision_s")
q2 <- function(yv, label) {
  t <- rbindlist(lapply(Q2V, function(v) { r <- sf2b_r_ci(TK[[v]], TK[[yv]], "spearman")
    data.table(family = label, with = v, n = r$n, rho = r$r, lower = r$lower, upper = r$upper, p = r$p) }))
  t[, p_holm := p.adjust(p, "holm")]
}
Q2 <- rbind(q2("on_page_s", "on-page seconds (paused)"), q2("note", "wrote a note (0/1, paused)"))
Q2 <- rbind(Q2, data.table(family = "wrote a note (0/1, paused)", with = "on_page_s",
                           n = nrow(TK), rho = sf2b_r_ci(TK$on_page_s, TK$note, "spearman")$r,
                           lower = sf2b_r_ci(TK$on_page_s, TK$note, "spearman")$lower, upper = sf2b_r_ci(TK$on_page_s, TK$note, "spearman")$upper,
                           p = sf2b_r_ci(TK$on_page_s, TK$note, "spearman")$p, p_holm = NA_real_))
Q2[, (c("rho", "lower", "upper")) := lapply(.SD, r3), .SDcols = c("rho", "lower", "upper")][, (c("p", "p_holm")) := lapply(.SD, signif, 3), .SDcols = c("p", "p_holm")]
fwrite(Q2, sf2b_out("11_pause_duration_correlations.csv"))
NOTE <- rbindlist(lapply(c("sf8", "on_page_s", "tired", "age"), function(v)
  TK[, .(variable = v, n = .N, mean = round(mean(get(v), na.rm = TRUE), 2), median = round(median(get(v), na.rm = TRUE), 2)), by = .(wrote_note = note)][order(wrote_note)]))
fwrite(NOTE, sf2b_out("11_pause_note.csv"))
out$paused$correlations <- Q2; out$paused$note_groups <- NOTE

# 3. Decline reasons ---------
DC <- O[took == 0]
REA <- DC[, .(n = .N, pct = round(100 * .N / nrow(DC), 1), sf8_mean = round(mean(sf8, na.rm = TRUE), 2), sf8_sd = round(sd(sf8, na.rm = TRUE), 2),
              sf8_median = round(median(sf8, na.rm = TRUE), 2), tired_mean = round(mean(tired, na.rm = TRUE), 2),
              pct_mobile = round(100 * mean(mobile), 1), minutes_before_median = round(median(minutes_before, na.rm = TRUE), 1)), by = reason][order(-n)]
REA <- rbind(REA, O[took == 1, .(reason = "(took the pause)", n = .N, pct = NA, sf8_mean = round(mean(sf8, na.rm = TRUE), 2),
                                 sf8_sd = round(sd(sf8, na.rm = TRUE), 2), sf8_median = round(median(sf8, na.rm = TRUE), 2),
                                 tired_mean = round(mean(tired, na.rm = TRUE), 2), pct_mobile = round(100 * mean(mobile), 1),
                                 minutes_before_median = round(median(minutes_before, na.rm = TRUE), 1))])
fwrite(REA, sf2b_out("11_pause_decline_reasons.csv"))
given <- DC[reason != "skipped"]
out$decline_reasons <- list(table = REA, decliners = nrow(DC),
  sf8_by_reason_kruskal_all5 = kruskal.test(sf8 ~ factor(reason), data = DC)$p.value,
  sf8_by_reason_kruskal_given4 = kruskal.test(sf8 ~ factor(reason), data = given)$p.value,
  tired_by_reason_kruskal_all5 = kruskal.test(tired ~ factor(reason), data = DC)$p.value,
  not_interested_vs_cant = "08_exploratory.json (registered exploratory): d = -.21, 95% CI -.43 to .02, p = .072")

# 4. Measurement checks ---------
O[, pause_group := factor(fifelse(took == 0, "declined", fifelse(engaged == 1, "engaged", "took, not engaged")),
                          levels = c("declined", "took, not engaged", "engaged"))]
GRP <- O[, .(n = .N, sf8_mean = round(mean(sf8, na.rm = TRUE), 3), sf8_sd = round(sd(sf8, na.rm = TRUE), 3),
             pause_days_mean = round(mean(crit_pause_days, na.rm = TRUE), 2)), by = pause_group][order(pause_group)]
dose <- rbindlist(c(
  list(data.table(check = "took (any pause)", events = sum(O$took), as.data.table(sf2b_r_ci(O$sf8, O$took, "spearman")))),
  lapply(DOSE_SECS, function(s) { y <- as.integer(O$took == 1 & O$on_page_s >= s)
    data.table(check = sprintf("paused >= %d s on page", s), events = sum(y), as.data.table(sf2b_r_ci(O$sf8, y, "spearman"))) }),
  list(data.table(check = "engaged (registered)", events = sum(O$engaged), as.data.table(sf2b_r_ci(O$sf8, O$engaged, "spearman"))),
       data.table(check = "on-page seconds, decliners at 0", events = sum(O$took), as.data.table(sf2b_r_ci(O$sf8, fifelse(O$took == 1, O$on_page_s, 0), "spearman"))),
       data.table(check = "on-page seconds, paused only", events = sum(O$took), as.data.table(sf2b_r_ci(TK$sf8, TK$on_page_s, "spearman"))),
       data.table(check = "pause group 0/1/2 (declined / took not engaged / engaged)", events = sum(O$took),
                  as.data.table(sf2b_r_ci(O$sf8, as.integer(O$pause_group) - 1, "spearman"))))))
setnames(dose, "r", "rho")
dose[, (c("rho", "lower", "upper")) := lapply(.SD, r3), .SDcols = c("rho", "lower", "upper")][, p := signif(p, 3)]
INT <- rbindlist(lapply(c("took", "engaged"), function(o) {
  d1 <- O[, .(y = get(o), zs = z(sf8), mobile, zt = z(tired))]
  fm <- firth(y ~ zs * mobile, as.data.frame(d1[complete.cases(d1[, .(y, zs, mobile)])]))
  ft <- firth(y ~ zs * zt, as.data.frame(d1[complete.cases(d1[, .(y, zs, zt)])]))
  # simple slopes of the 8: desktop / mobile; tiredness at -1 SD / +1 SD
  sl <- function(f, w) { b <- coef(f); V <- vcov(f); g <- setNames(rep(0, length(b)), names(b)); g["zs"] <- 1; g[grep(":", names(b))] <- w
    e <- sum(g * b); se <- sqrt(drop(t(g) %*% V %*% g)); c(or = exp(e), lower = exp(e - 1.96 * se), upper = exp(e + 1.96 * se)) }
  mk <- function(f, lab, tm, w) { tb <- firth_tab(f); i <- tb[tb$term == tm]; stopifnot(nrow(i) == 1)
    s0 <- sl(f, w[1]); s1 <- sl(f, w[2])
    data.table(outcome = o, interaction = lab, n = f$n, interaction_or = i$odds_ratio, lower = i$lower, upper = i$upper, p = i$p,
               sf8_or_at_a = s0[["or"]], sf8_or_at_a_lower = s0[["lower"]], sf8_or_at_a_upper = s0[["upper"]],
               sf8_or_at_b = s1[["or"]], sf8_or_at_b_lower = s1[["lower"]], sf8_or_at_b_upper = s1[["upper"]]) }
  rbind(cbind(mk(fm, "sf8 x mobile device", "zs:mobile", c(0, 1)), a = "desktop", b = "phone or tablet"),
        cbind(mk(ft, "sf8 x tiredness (z)", "zs:zt", c(-1, 1)), a = "tiredness -1 SD", b = "tiredness +1 SD"))
}))
nc <- c("interaction_or", "lower", "upper", grep("^sf8_or", names(INT), value = TRUE))
INT[, (nc) := lapply(.SD, r3), .SDcols = nc][, p := signif(p, 3)]
INT[, p_holm := signif(p.adjust(p, "holm"), 3)]
MEAS <- rbind(cbind(part = "pause groups", GRP), cbind(part = "dose", dose), fill = TRUE)
fwrite(MEAS, sf2b_out("11_pause_measurement.csv"))
fwrite(INT, sf2b_out("11_pause_interactions.csv"))
out$measurement <- list(groups = GRP, groups_kruskal_p = kruskal.test(sf8 ~ pause_group, data = O)$p.value,
                        dose = dose, interactions = INT)

write_json(rnd(out), sf2b_out("11_pause.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
sf2b_log("done: ", sf2b_out("11_pause.json"))
