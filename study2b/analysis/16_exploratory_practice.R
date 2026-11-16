# 16: practice as a dose, with spirituality held constant. EXPLORATORY, NOT REGISTERED.
#
# Not part of the filed registration (osf.io/d4a5u) and not in its list of exploratory
# analyses. Specified on 15 November 2026, after H6 (d .70 current vs never; .45 adjusted
# for spirituality and religiosity; dedicated vs never .88), 11 (practitioners paused
# less) and 13 (a general self-description factor) had been seen. Every number here is
# hypothesis-generating only. No registered cell, test or output is changed.
#
#   1. Dose. The 8 by practice frequency: the registered 0-4 code (H5b: 0 = no regular
#      practice or none in the past three months, 1 = less than weekly, 2 = 1-2 a week,
#      3 = 3-5 a week, 4 = daily or near-daily) and the background item's own six
#      categories (no regular practice kept apart from regular practice with none in the
#      past three months). Means and SDs by level; the linear trend (OLS slope, Spearman
#      rho), departure from linearity (the level factor against the linear term, F), and
#      monotonicity (are the level means non-decreasing; each adjacent step with its
#      Welch CI). Observed 8 and latent (07's MIMIC with frequency as the covariate, and
#      with the levels as dummies), all participants and among regular practitioners.
#   2. Spirituality and religiosity. The same trend with both as covariates (07 adjusts
#      H6 for them); the 8 by each alone, Spearman rho (and the latent r with each as a
#      single ordinal indicator).
#   3. Practice types (bg_practice_types, a multi-select, read from the extraction): the
#      8's mean for every type reported by 30 or more regular practitioners, and its
#      difference from regular practitioners not reporting that type, adjusted for
#      frequency (and, second, also for spirituality and religiosity). Descriptive.
#   4. Dedicated vs current non-dedicated vs never.
#   5. Net of the general factor and of Observing: the dose trend with 13's observed
#      general factor G_obs (first principal component of the 28 other observed scores,
#      13's construction) and with FFMQ Observing (observed; and latent in the MIMIC) as
#      covariates.
#   6. The pause (11: practitioners paused less): take-up and engaged-pause rates by
#      frequency, by group and by type. Descriptive (73 took, 44 engaged).
#
# Data. The pipeline's saved real-mode objects (04's top-level code, as 10-15 do). One
# fact the derived object does not hold, the practice types, is read from the same
# extraction 01 read (named in the object; no database query), from the same submission
# 01 uses (each participant's first complete sf-background, flagged resubmissions dropped),
# option ids only: the free text of "other" is never read into the analysis or written.
# Nothing from the extraction is written out except counts. Prolific IDs are not in it.
#
# Settings as registered (Section 5) for every latent model: ordinal items, WLSMV, theta,
# missing = "pairwise", conditional.x = FALSE, the 8's agreement factor as registered
# (uncorrelated with every covariate: RI ~ 0*x, as 07). CHOICES below are in the JSON.
#
# Writes output/<tag>/16_*.csv and 16_exploratory_practice.json (aggregates only).
# Run from the repo root:  SF2B_MODE=real Rscript study2b/analysis/16_exploratory_practice.R

# Set Up ---------
## Load libraries ---------
packages <- c("jsonlite", "data.table", "lavaan", "psych")
new_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(new_packages)) install.packages(new_packages)
for (thispack in packages) {
  library(thispack, character.only = TRUE, quietly = TRUE, verbose = FALSE)
}

## Settings ---------
source("study2b/analysis/config.R")
source("study2b/analysis/sf2b_models.R")
stopifnot(SF2B$mode == "real")
T0 <- Sys.time()
ex04 <- parse("study2b/analysis/04_validity.R", keep.source = FALSE)
first_fit <- which(vapply(ex04, function(e) grepl("^sf2b_log\\(\"fitting", deparse(e)[1]), TRUE))
stopifnot(length(first_fit) == 1)
for (i in seq_len(first_fit - 1)) eval(ex04[[i]], envir = globalenv())
stopifnot(exists("DAT"), exists("score_items"), !exists("RES"))
ex10 <- parse("study2b/analysis/10_deviation_maia.R", keep.source = FALSE)
is_def <- function(e, nm) is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name(nm))
eval(Filter(function(e) is_def(e, "improper_of"), as.list(ex10))[[1]], envir = globalenv())
# lavaan's phantom factor for an observed covariate (x =~ x) also has a standardized
# loading of 1 by construction: not improper.
improper_list <- function(f) { x <- improper_of(f)$list
  x[!grepl("^std loading (L_[[:alnum:]_]+|grp) =~ ", x) & !grepl("^std loading ([[:alnum:]_]+) =~ \\1 = ", x)] }
stopifnot(identical(M$participant, S$participant), nrow(DAT) == nrow(M))

I8 <- D$sf8; OBS <- "ffmq_observing"; IO <- score_items[[OBS]]
TYPE_MIN_N <- 30
FREQ6 <- c("no regular practice", "regular, none in past 3 months", "less than weekly", "1-2 a week", "3-5 a week", "daily or near-daily")
TYPE_LABEL <- c(meditation = "meditation", yoga = "yoga", taichi_qigong = "tai chi / qigong", contemplative_prayer = "contemplative prayer",
                martial_arts = "martial arts", other = "other")

CHOICES <- list(
  frequency_codes = "freq5 = the registered H5b code 0-4 (01's freq_code; 0 joins no regular practice and regular practice with none in the past three months). freq6 = the background item's categories with those two kept apart (0-5). One participant without a regular-practice answer is missing on both.",
  trend = "OLS of the observed 8 on the code (slope per step, and per step in the 8's SD), Spearman rho, and departure from linearity as the F test of the level factor against the linear term. Monotonic = the level means are non-decreasing; each adjacent step with its Welch 95% CI.",
  latent = "07's MIMIC: the H1 model with S regressed on the covariate(s), S's residual variance 1 (so coefficients are in within-level latent SDs), the agreement factor uncorrelated with every covariate (RI ~ 0*x); observed covariates exogenous with conditional.x = FALSE; cases complete on the covariates. The latent level means are dummy coefficients against level 0.",
  practitioners_only = "The trend is also reported among regular practitioners only (practice_regular = yes; freq5 0-4), so 'no regular practice' does not carry it.",
  spirituality = "Spirituality and religiosity (bg_spiritual, bg_religious, 1-5) as linear covariates, as 07 enters them. Rank correlations: Spearman with the observed 8; latent r with each as a single ordinal indicator (04's single-item treatment).",
  practice_types = "Read from the extraction (01's first complete sf-background submission per participant; option ids only). Types reported by at least 30 regular practitioners. Multi-select, so types overlap; each type is compared with regular practitioners who did not report it: OLS of the 8 on the type indicator plus frequency as a factor (freq5 levels among regular practitioners), then also plus spirituality and religiosity. Descriptive: overlapping, self-selected categories.",
  groups = "Dedicated (01: regular, 3+ a week, a year or more, teacher-led training), current non-dedicated (01's 'current', 1-2 a week or more, not dedicated), never (no regular practice). 07's sf2b_kg for each pair (Student's t, pooled-SD d); latent: MIMIC with two dummies against never.",
  net_of_general_and_observing = "G_obs exactly as 13 builds it: psych::principal, one component, of the 28 other observed scale scores (13's battery without the 8), standardized, complete cases, oriented with the 8; checked against 13's adjusted H6. The dose trend with G_obs, with observed FFMQ Observing, and with both, as covariates, on the same cases; latent: S regressed on frequency and the latent Observing factor (Observing also regressed on frequency), and S on frequency and G_obs (observed).",
  pause = "Take-up and engaged-pause rates among the 710 offered, by freq5, freq6, group and type (n >= 30 offered), Wilson 95% CIs; a chi-square test for trend in proportions across freq5 (prop.trend.test), descriptive only given 73 and 44 events.",
  improper = "10's improper_of() on every fitted model.")

## Practice types from the extraction ---------
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
TY <- local({
  rows <- fromJSON(raw_file, simplifyVector = FALSE)$rows
  bg <- Filter(function(r) identical(r$source, "questionnaire") && identical(r$slug, "sf-background") && is.null(r$resubmission_of) &&
                 isTRUE(r$participant %in% M$participant) && "bg_practice_regular" %in% names(r$responses), rows)
  rm(rows)
  rec <- rbindlist(lapply(bg, function(r) {
    v <- r$responses$bg_practice_types
    ids <- if (is.list(v)) unique(vapply(v, function(o) if (is.list(o) && !is.null(o$option_id)) o$option_id else NA_character_, "")) else character()
    reg <- r$responses$bg_practice_regular
    data.table(participant = r$participant, received = r$received_at %||% NA_character_,
               types = paste(sort(na.omit(ids)), collapse = " "), types_status = if (is.list(v)) "answered" else if (is.null(v)) "absent" else as.character(v),
               regular_check = if (is.list(reg) && !is.null(reg$option_id)) reg$option_id else NA_character_)
  }))
  rec[, received := ts_parse(as.list(received))]
  rec <- rec[order(participant, received)][!duplicated(participant)]
  rec[match(M$participant, participant)]
})
stopifnot(identical(TY$participant, M$participant))
types_check <- list(mismatched_regular_answer = sum(TY$regular_check != M$practice_regular, na.rm = TRUE),
                    status = as.list(table(TY$types_status, useNA = "ifany")),
                    status_by_regular = as.data.table(table(regular = M$practice_regular, status = TY$types_status, useNA = "ifany")))
stopifnot(types_check$mismatched_regular_answer == 0)
TYPES <- names(TYPE_LABEL)
TYM <- sapply(TYPES, function(k) as.integer(vapply(strsplit(TY$types, " "), function(x) k %in% x, TRUE)))
TYM[!(TY$types_status %in% "answered") | !(M$practice_regular %in% "yes"), ] <- NA

## The analysis frame ---------
PF <- data.table(sf8 = S$sf8, freq5 = M$freq_code,
                freq6 = fifelse(M$practice_regular %in% "no", 0, fifelse(M$practice_regular %in% "yes",
                        unname(c(none = 1, lt_weekly = 2, `1_2_week` = 3, `3_5_week` = 4, daily = 5)[M$practice_freq]), NA_real_)),
                regular = fifelse(is.na(M$practice_regular), NA_integer_, as.integer(M$practice_regular == "yes")),
                spiritual = M$spiritual, religious = M$religious, observing = S$ffmq_observing,
                group3 = fifelse(M$dedicated, "dedicated", fifelse(M$practice_group %in% "current", "current, not dedicated",
                                 fifelse(M$practice_group %in% "never", "never", NA_character_))),
                took = M$pause_taken, engaged = M$engaged_pause, offered = M$has_pause)
for (k in TYPES) PF[[paste0("type_", k)]] <- TYM[, k]
stopifnot(PF[group3 == "dedicated", .N] == 94, PF[group3 == "never", .N] == 164, PF[group3 == "current, not dedicated", .N] == 349)

# G_obs, 13's construction.
J13 <- fromJSON(sf2b_out("13_exploratory_general_factor.json"))
BAT <- J13$matrix$battery; B28 <- setdiff(BAT, "sf8")
OBSD <- as.data.frame(S[, ..BAT])
pc28 <- principal(as.data.frame(scale(OBSD[, B28])), nfactors = 1, scores = TRUE, missing = FALSE)
Gobs <- as.numeric(pc28$scores[, 1])
if (cor(Gobs, S$sf8, use = "complete.obs") < 0) Gobs <- -Gobs
PF$G_obs <- Gobs
g_cur <- fifelse(M$practice_group == "current", 1L, fifelse(M$practice_group == "never", 0L, NA_integer_))
chk13 <- list(here = sf2b_kg(S$sf8, g_cur, data.table(G = Gobs))[, .(n_group, n_reference, d)],
              in_13 = tryCatch(fread(sf2b_out("13_observed_known_groups.csv"))[analysis == "observed 8, adjusted for G_obs", .(n_group, n_reference, d)],
                               error = function(e) NULL))
stopifnot(is.null(chk13$in_13) || abs(chk13$here$d - chk13$in_13$d) < .001)

## Helpers ---------
num <- function(DT, k = 3) { cols <- names(DT)[vapply(DT, is.numeric, TRUE)]; DT[, (cols) := lapply(.SD, function(x) round(x, k)), .SDcols = cols]; DT }
rnd <- function(x) if (is.data.frame(x)) x else if (is.list(x)) lapply(x, rnd) else if (is.numeric(x)) signif(x, 4) else x
ci_mean <- function(x) { x <- x[!is.na(x)]; n <- length(x); m <- mean(x); s <- sd(x); h <- qt(.975, n - 1) * s / sqrt(n); c(n = n, mean = m, sd = s, lower = m - h, upper = m + h) }
wilson <- function(k, n) if (n > 0) prop.test(k, n, correct = FALSE)$conf.int else c(NA, NA)
level_table <- function(y, lev, labels, what) {
  ok <- !is.na(y) & !is.na(lev)
  rbindlist(lapply(sort(unique(lev[ok])), function(l) { x <- ci_mean(y[ok & lev == l])
    data.table(outcome = what, level = l, label = labels[l + 1], n = x[["n"]], mean = x[["mean"]], sd = x[["sd"]], lower = x[["lower"]], upper = x[["upper"]]) }))
}
steps <- function(y, lev) {
  ok <- !is.na(y) & !is.na(lev); L <- sort(unique(lev[ok]))
  rbindlist(lapply(seq_len(length(L) - 1), function(i) { a <- y[ok & lev == L[i]]; b <- y[ok & lev == L[i + 1]]
    t <- t.test(b, a); data.table(from = L[i], to = L[i + 1], difference = mean(b) - mean(a), lower = t$conf.int[1], upper = t$conf.int[2], p = t$p.value) }))
}
trend <- function(y, x, covars = NULL, label, sample = "all", keep = rep(TRUE, length(y))) {
  d <- data.table(y = y, x = x); if (!is.null(covars)) d <- cbind(d, covars)
  d <- d[keep & complete.cases(d)]
  cv <- setdiff(names(d), c("y", "x"))
  m1 <- lm(reformulate(c("x", cv), "y"), data = d)
  mf <- lm(reformulate(c("factor(x)", cv), "y"), data = d)
  m0 <- lm(reformulate(c("x", cv), "y"), data = d)
  co <- summary(m1)$coefficients["x", ]; ci <- confint(m1)["x", ]
  nl <- anova(m0, mf)
  data.table(analysis = label, sample = sample, n = nrow(d), covariates = paste(cv, collapse = " + "),
             slope = co[["Estimate"]], lower = ci[[1]], upper = ci[[2]], p = co[["Pr(>|t|)"]],
             slope_per_sd_of_8 = co[["Estimate"]] / sd(d$y), spearman_rho = if (!length(cv)) cor(d$y, d$x, method = "spearman") else NA_real_,
             nonlinearity_F = nl$F[2], nonlinearity_df = nl$Df[2], nonlinearity_p = nl$`Pr(>F)`[2])
}
fit_mimic <- function(covs, rhs, extra = character(), with_obs = FALSE, keep = rep(TRUE, nrow(DAT)), label) {
  d <- DAT; for (v in names(covs)) d[[v]] <- covs[[v]]
  d <- d[keep & complete.cases(as.data.frame(covs)), , drop = FALSE]
  # Covariates mean-centred in the analysed cases. Without it the fit with the freq5 dummies
  # plus spirituality and religiosity (1-5) converged to a degenerate solution (S's
  # loadings near 0, coefficients near 100): S's intercept is fixed at 0, so an uncentred
  # covariate's mean has to be carried by the thresholds. Centring leaves the slopes of a
  # proper solution unchanged.
  for (v in names(covs)) d[[v]] <- d[[v]] - mean(d[[v]])
  ord <- c(I8, if (with_obs) IO)
  syn <- c(sf2b_scale_syntax(I8), paste0("S ~ ", paste(rhs, collapse = " + ")), paste0("RI ~ ", paste0("0*", names(covs), collapse = " + ")), extra)
  if (with_obs) { x <- sf2b_validity_syntax(OBS, D); syn <- c(syn, x[!grepl("^S ~~", x)]) }
  f <- sf2b_try(sf2b_sem(syn, d, ordered = ord))
  if (sf2b_failed(f)) return(data.table(model = label, n = nrow(d), error = sf2b_msg(f)))
  pe <- parameterEstimates(f); ss <- standardizedSolution(f)
  keepl <- pe$label[nzchar(pe$label) & pe$op %in% c("~", ":=")]
  x <- pe[pe$label %in% keepl, c("label", "lhs", "op", "rhs", "est", "se", "ci.lower", "ci.upper", "pvalue")]
  x$std_all <- ss$est.std[match(x$label, ss$label)]
  bad <- if (with_obs) sf2b_orientation(f, I8, OBS, D) else sf2b_orientation(f, I8, character(), D)
  fs <- sf2b_fit_stats(f)
  # A degenerate solution shows as S almost wholly predicted by the covariates.
  r2 <- tryCatch(lavInspect(f, "r2")[["S"]], error = function(e) NA_real_)
  lam <- lavInspect(f, "std")$lambda[I8, "S"] * ifelse(I8 %in% D$reversed, -1, 1)
  cbind(data.table(model = label, n = nrow(d), error = if (length(bad)) paste("not oriented:", paste(bad, collapse = ", ")) else "",
                   improper = paste(c(improper_list(f), if (isTRUE(r2 > .95)) sprintf("degenerate? R2 of S = %.3f", r2)), collapse = "; "),
                   r2_S = r2, min_keyed_std_loading = min(lam), max_keyed_std_loading = max(lam),
                   cfi = fs$cfi.scaled, rmsea = fs$rmsea.scaled), as.data.table(x))
}

## 1-2. Dose, observed ---------
sf2b_log("16: dose")
LEV <- rbind(cbind(code = "freq5", level_table(PF$sf8, PF$freq5, c("0: none", "1: less than weekly", "2: 1-2 a week", "3: 3-5 a week", "4: daily"), "the 8 (observed)")),
             cbind(code = "freq6", level_table(PF$sf8, PF$freq6, FREQ6, "the 8 (observed)")),
             cbind(code = "freq5", level_table(PF$observing, PF$freq5, c("0: none", "1: less than weekly", "2: 1-2 a week", "3: 3-5 a week", "4: daily"), "FFMQ Observing (observed)")),
             cbind(code = "freq5", level_table(PF$G_obs, PF$freq5, c("0: none", "1: less than weekly", "2: 1-2 a week", "3: 3-5 a week", "4: daily"), "G_obs (13)")))
# Spirituality-adjusted level means (OLS with the level factor, at the covariate means).
adj_means <- function(code) {
  d <- PF[!is.na(sf8) & !is.na(get(code)) & !is.na(spiritual) & !is.na(religious), .(y = sf8, lev = factor(get(code)), spiritual, religious)]
  m <- lm(y ~ lev + spiritual + religious, data = d)
  nd <- data.table(lev = factor(levels(d$lev), levels = levels(d$lev)), spiritual = mean(d$spiritual), religious = mean(d$religious))
  p <- predict(m, nd, interval = "confidence")
  data.table(code = code, level = as.numeric(levels(d$lev)), adjusted_mean = p[, "fit"], adjusted_lower = p[, "lwr"], adjusted_upper = p[, "upr"])
}
AM <- rbind(adj_means("freq5"), adj_means("freq6"))
LEV <- merge(LEV, AM, by = c("code", "level"), all.x = TRUE, sort = FALSE)
LEV[outcome != "the 8 (observed)", (c("adjusted_mean", "adjusted_lower", "adjusted_upper")) := NA_real_]
MONO <- rbindlist(lapply(c("freq5", "freq6"), function(cd) {
  m <- LEV[code == cd & outcome == "the 8 (observed)"][order(level)]
  st <- steps(PF$sf8, PF[[cd]])
  cbind(code = cd, st, monotonic_means = all(diff(m$mean) >= 0), monotonic_adjusted_means = all(diff(m$adjusted_mean) >= 0))
}))
REG <- PF$regular %in% 1
TR <- rbind(
  trend(PF$sf8, PF$freq5, label = "the 8 on freq5"),
  trend(PF$sf8, PF$freq6, label = "the 8 on freq6"),
  trend(PF$sf8, PF$freq5, label = "the 8 on freq5", sample = "regular practitioners", keep = REG),
  trend(PF$sf8, PF$freq5, PF[, .(spiritual, religious)], label = "the 8 on freq5, spirituality and religiosity held constant"),
  trend(PF$sf8, PF$freq6, PF[, .(spiritual, religious)], label = "the 8 on freq6, spirituality and religiosity held constant"),
  trend(PF$sf8, PF$freq5, PF[, .(spiritual, religious)], label = "the 8 on freq5, spirituality and religiosity held constant", sample = "regular practitioners", keep = REG),
  # 5. net of G_obs and of Observing, on the cases complete on all of them
  { k <- complete.cases(PF[, .(sf8, freq5, G_obs, observing, spiritual, religious)])
    rbind(trend(PF$sf8, PF$freq5, label = "the 8 on freq5 (same cases as below)", sample = "complete on G_obs, Observing, spirituality", keep = k),
          trend(PF$sf8, PF$freq5, PF[, .(G_obs)], label = "the 8 on freq5, net of G_obs", sample = "complete on G_obs, Observing, spirituality", keep = k),
          trend(PF$sf8, PF$freq5, PF[, .(observing)], label = "the 8 on freq5, net of observed FFMQ Observing", sample = "complete on G_obs, Observing, spirituality", keep = k),
          trend(PF$sf8, PF$freq5, PF[, .(G_obs, observing)], label = "the 8 on freq5, net of G_obs and Observing", sample = "complete on G_obs, Observing, spirituality", keep = k),
          trend(PF$sf8, PF$freq5, PF[, .(G_obs, observing, spiritual, religious)], label = "the 8 on freq5, net of G_obs, Observing, spirituality and religiosity",
                sample = "complete on G_obs, Observing, spirituality", keep = k),
          trend(PF$observing, PF$freq5, label = "observed FFMQ Observing on freq5 (context)", sample = "complete on G_obs, Observing, spirituality", keep = k),
          trend(PF$G_obs, PF$freq5, label = "G_obs on freq5 (context)", sample = "complete on G_obs, Observing, spirituality", keep = k)) })

## 1-2, 4-5. Latent (MIMIC) ---------
sf2b_log("16: latent models")
dum <- function(x, levs) setNames(lapply(levs, function(l) as.integer(x == l)), paste0("lev", levs))
D5 <- as.data.table(dum(PF$freq5, 1:4)); D6 <- as.data.table(dum(PF$freq6, 1:5))
LAT <- rbindlist(list(
  fit_mimic(list(freq5 = PF$freq5), "b_freq*freq5", label = "S on freq5 (linear)"),
  fit_mimic(as.list(D5), paste0("b_", names(D5), "*", names(D5)), label = "S on freq5 levels (dummies vs 0)"),
  fit_mimic(as.list(D6), paste0("b_", names(D6), "*", names(D6)), label = "S on freq6 levels (dummies vs no regular practice)"),
  fit_mimic(list(freq5 = PF$freq5), "b_freq*freq5", keep = REG, label = "S on freq5 (linear), regular practitioners"),
  fit_mimic(list(freq5 = PF$freq5, spiritual = PF$spiritual, religious = PF$religious), c("b_freq*freq5", "b_spiritual*spiritual", "b_religious*religious"),
            label = "S on freq5, spirituality and religiosity"),
  fit_mimic(c(as.list(D5), list(spiritual = PF$spiritual, religious = PF$religious)), c(paste0("b_", names(D5), "*", names(D5)), "b_spiritual*spiritual", "b_religious*religious"),
            label = "S on freq5 levels, spirituality and religiosity"),
  fit_mimic(list(freq5 = PF$freq5, G_obs = PF$G_obs), c("b_freq*freq5", "b_G*G_obs"), keep = !is.na(PF$G_obs), label = "S on freq5 and G_obs"),
  fit_mimic(list(freq5 = PF$freq5), c("b_freq*freq5", paste0("b_obs*", OBS)), extra = paste0(OBS, " ~ a_freq*freq5"), with_obs = TRUE,
            label = "S on freq5 and latent FFMQ Observing (Observing on freq5)"),
  fit_mimic(list(freq5 = PF$freq5, G_obs = PF$G_obs), c("b_freq*freq5", "b_G*G_obs", paste0("b_obs*", OBS)), extra = paste0(OBS, " ~ a_freq*freq5 + a_G*G_obs"),
            with_obs = TRUE, keep = !is.na(PF$G_obs), label = "S on freq5, G_obs and latent FFMQ Observing"),
  # 4. groups
  fit_mimic(list(ded = as.integer(PF$group3 == "dedicated"), cnd = as.integer(PF$group3 == "current, not dedicated")),
            c("b_dedicated*ded", "b_current_not_dedicated*cnd"), extra = "dif_ded_cnd := b_dedicated - b_current_not_dedicated",
            keep = !is.na(PF$group3), label = "S on group (dedicated, current not dedicated; vs never)")), fill = TRUE)

## 2. Spirituality and religiosity alone ---------
SPIR <- rbindlist(lapply(c("spiritual", "religious"), function(v) {
  r <- sf2b_r_ci(PF$sf8, PF[[v]], "spearman")
  d <- DAT; d[[v]] <- PF[[v]]
  syn <- c(sf2b_scale_syntax(I8), paste0("L_", v, " =~ 1*", v), paste0(v, " ~~ 0*", v), paste0("L_", v, " ~~ 1*L_", v),
           paste0("S ~~ r_", v, "*L_", v), paste0("RI ~~ 0*L_", v))
  f <- sf2b_try(sf2b_sem(syn, d, ordered = c(I8, v)))
  lr <- if (sf2b_failed(f)) list(est = NA, ci.lower = NA, ci.upper = NA) else { pe <- parameterEstimates(f); pe[pe$label == paste0("r_", v), ] }
  data.table(variable = v, n = r$n, spearman_rho_with_8 = r$r, lower = r$lower, upper = r$upper, p = r$p,
             spearman_rho_with_freq5 = cor(PF[[v]], PF$freq5, method = "spearman", use = "complete.obs"),
             latent_r_with_8 = lr$est, latent_lower = lr$ci.lower, latent_upper = lr$ci.upper,
             latent_error = if (sf2b_failed(f)) sf2b_msg(f) else "", latent_improper = if (sf2b_failed(f)) "" else paste(improper_list(f), collapse = "; "))
}))
SPIR <- rbind(SPIR, data.table(variable = "spiritual with religious", n = sum(complete.cases(PF$spiritual, PF$religious)),
                               spearman_rho_with_8 = cor(PF$spiritual, PF$religious, method = "spearman", use = "complete.obs")), fill = TRUE)

## 3. Practice types ---------
sf2b_log("16: practice types")
RP <- PF[regular %in% 1]
type_n <- vapply(TYPES, function(k) sum(RP[[paste0("type_", k)]] %in% 1), 0)
TYP <- rbindlist(lapply(TYPES, function(k) {
  v <- paste0("type_", k); has <- RP[[v]]
  base <- data.table(type = k, label = TYPE_LABEL[[k]], n_reporting = sum(has %in% 1), n_only_this_type = sum(has %in% 1 & rowSums(RP[, paste0("type_", TYPES), with = FALSE], na.rm = TRUE) == 1),
                     reported = sum(has %in% 1) >= TYPE_MIN_N)
  if (!base$reported) return(base)
  x <- ci_mean(RP$sf8[has %in% 1]); x0 <- ci_mean(RP$sf8[has %in% 0])
  d1 <- RP[!is.na(has) & !is.na(sf8) & !is.na(freq5), .(y = sf8, h = has, f = factor(freq5), spiritual, religious)]
  m1 <- lm(y ~ h + f, data = d1); c1 <- summary(m1)$coefficients["h", ]; ci1 <- confint(m1)["h", ]
  d2 <- d1[complete.cases(d1)]
  m2 <- lm(y ~ h + f + spiritual + religious, data = d2); c2 <- summary(m2)$coefficients["h", ]; ci2 <- confint(m2)["h", ]
  sdp <- sd(RP$sf8, na.rm = TRUE)
  cbind(base, data.table(mean_8 = x[["mean"]], sd_8 = x[["sd"]], mean_8_others = x0[["mean"]], n_others = x0[["n"]],
                         mean_freq5 = mean(RP$freq5[has %in% 1], na.rm = TRUE), mean_freq5_others = mean(RP$freq5[has %in% 0], na.rm = TRUE),
                         mean_spiritual = mean(RP$spiritual[has %in% 1], na.rm = TRUE), mean_religious = mean(RP$religious[has %in% 1], na.rm = TRUE),
                         diff_adj_freq = c1[["Estimate"]], diff_adj_freq_lower = ci1[[1]], diff_adj_freq_upper = ci1[[2]], p_adj_freq = c1[["Pr(>|t|)"]],
                         d_adj_freq = c1[["Estimate"]] / sdp, n_adj_freq = nrow(d1),
                         diff_adj_freq_spirit = c2[["Estimate"]], diff_adj_freq_spirit_lower = ci2[[1]], diff_adj_freq_spirit_upper = ci2[[2]],
                         p_adj_freq_spirit = c2[["Pr(>|t|)"]], d_adj_freq_spirit = c2[["Estimate"]] / sdp, n_adj_freq_spirit = nrow(d2)))
}), fill = TRUE)

## 4. Groups, observed ---------
GRPM <- PF[!is.na(group3) & !is.na(sf8), .(n = .N, mean = mean(sf8), sd = sd(sf8), spiritual_mean = mean(spiritual, na.rm = TRUE),
                                         religious_mean = mean(religious, na.rm = TRUE)), by = group3][order(match(group3, c("never", "current, not dedicated", "dedicated")))]
pair <- function(a, b, covars = NULL) { g <- fifelse(PF$group3 == a, 1L, fifelse(PF$group3 == b, 0L, NA_integer_))
  cbind(contrast = paste(a, "vs", b), sf2b_kg(PF$sf8, g, covars)) }
GRP <- rbind(pair("dedicated", "never"), pair("current, not dedicated", "never"), pair("dedicated", "current, not dedicated"),
             cbind(pair("dedicated", "never", PF[, .(spiritual, religious)])), cbind(pair("current, not dedicated", "never", PF[, .(spiritual, religious)])),
             cbind(pair("dedicated", "current, not dedicated", PF[, .(spiritual, religious)])), fill = TRUE)

## 6. The pause ---------
O <- PF[offered %in% TRUE]
prate <- function(v, label, levels_label = NULL) {
  d <- O[!is.na(get(v))]
  d[, .(n = .N, took = sum(took), engaged = sum(engaged)), by = .(level = get(v))][order(level)][, `:=`(
    by = label, took_pct = 100 * took / n, took_lower = 100 * mapply(function(k, n) wilson(k, n)[1], took, n),
    took_upper = 100 * mapply(function(k, n) wilson(k, n)[2], took, n), engaged_pct = 100 * engaged / n,
    engaged_lower = 100 * mapply(function(k, n) wilson(k, n)[1], engaged, n), engaged_upper = 100 * mapply(function(k, n) wilson(k, n)[2], engaged, n))][]
}
PZ <- rbind(prate("freq5", "freq5"), prate("freq6", "freq6"), prate("group3", "group"), fill = TRUE)
PZ[, level := as.character(level)]
PZT <- rbindlist(lapply(TYPES, function(k) { v <- paste0("type_", k); d <- O[get(v) %in% 1]
  if (nrow(d) < TYPE_MIN_N) return(data.table(by = "type (regular practitioners reporting it)", level = k, n = nrow(d), note = "fewer than 30 offered: not reported"))
  data.table(by = "type (regular practitioners reporting it)", level = k, n = nrow(d), took = sum(d$took), engaged = sum(d$engaged),
             took_pct = 100 * mean(d$took), took_lower = 100 * wilson(sum(d$took), nrow(d))[1], took_upper = 100 * wilson(sum(d$took), nrow(d))[2],
             engaged_pct = 100 * mean(d$engaged), engaged_lower = 100 * wilson(sum(d$engaged), nrow(d))[1], engaged_upper = 100 * wilson(sum(d$engaged), nrow(d))[2]) }), fill = TRUE)
PZ <- rbind(PZ, PZT, fill = TRUE)
setcolorder(PZ, c("by", "level"))
ptrend <- lapply(c("took", "engaged"), function(o) { t5 <- O[!is.na(freq5), .(k = sum(get(o)), n = .N), by = freq5][order(freq5)]
  x <- suppressWarnings(prop.trend.test(t5$k, t5$n, score = t5$freq5)); list(outcome = o, chisq = unname(x$statistic), p = x$p.value) })

## Write ---------
fwrite(num(copy(LEV)), sf2b_out("16_dose_levels.csv"))
fwrite(num(copy(MONO)), sf2b_out("16_dose_steps.csv"))
fwrite(num(copy(TR), 4), sf2b_out("16_trend.csv"))
fwrite(num(copy(LAT), 4), sf2b_out("16_latent.csv"))
fwrite(num(copy(SPIR)), sf2b_out("16_spirituality.csv"))
fwrite(num(copy(TYP)), sf2b_out("16_types.csv"))
fwrite(num(copy(GRPM)), sf2b_out("16_groups_means.csv"))
fwrite(num(copy(GRP), 4), sf2b_out("16_groups.csv"))
fwrite(num(copy(PZ), 1), sf2b_out("16_pause.csv"))
out <- list(
  analysis = "16: practice as a dose, with spirituality held constant. EXPLORATORY, NOT REGISTERED (osf.io/d4a5u does not include it); specified after H6, 11 and 13 had been seen.",
  mode = SF2B$mode, tag = SF2B$tag, written = format(Sys.time(), "%Y-%m-%d %H:%M"),
  software = list(R = R.version.string, lavaan = as.character(packageVersion("lavaan")), psych = as.character(packageVersion("psych"))),
  n = nrow(PF), choices = CHOICES,
  practice_types_read = list(status = types_check$status, status_by_regular = types_check$status_by_regular,
                             regular_answer_mismatches_with_01 = types_check$mismatched_regular_answer,
                             reporting_by_type_regular_practitioners = as.list(type_n)),
  g_obs_check_against_13 = chk13,
  dose_levels = LEV, dose_steps = MONO, trend = TR, latent = LAT, spirituality = SPIR, types = TYP,
  groups = list(means = GRPM, contrasts = GRP), pause = list(rates = PZ, trend_test_freq5 = ptrend),
  improper = LAT[nzchar(improper %||% ""), unique(paste(model, improper, sep = ": "))],
  run_seconds = round(as.numeric(difftime(Sys.time(), T0, units = "secs"))))
write_json(rnd(out), sf2b_out("16_exploratory_practice.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
sf2b_log("16 done in ", out$run_seconds, " s; outputs ", sf2b_out("16_*"))
print(num(copy(LEV))[outcome == "the 8 (observed)"]); print(num(copy(MONO)))
print(num(copy(TR), 3)[, .(analysis, sample, n, slope, lower, upper, p, slope_per_sd_of_8, spearman_rho, nonlinearity_p)])
print(num(copy(LAT), 3)[, .(model, n, label, est, ci.lower, ci.upper, pvalue, std_all, improper, error)])
print(num(copy(SPIR))); print(num(copy(TYP))); print(num(copy(GRPM))); print(num(copy(GRP), 3)[, .(contrast, covariates, n_group, n_reference, d, d_lower, d_upper, p)])
print(num(copy(PZ), 1)); print(ptrend)
