# 22: how the demographics relate to the 8 (and the 10), to the agreement factor and to
# FFMQ Observing. EXPLORATORY, NOT REGISTERED.
#
# Not part of the filed registration (osf.io/d4a5u) and not in its list of exploratory
# analyses. Specified on 16 November 2026, after the registered analyses and 10-21 had
# been seen, for the Psychological Assessment manuscript (sample description and
# Constraints on Generality). Every number here is hypothesis-generating: no correction
# for multiplicity is applied to the decisions (there are none), but a Holm-adjusted p is
# given within each family (below). No registered cell, test or output is changed.
#
# What is already reported elsewhere, and is cited rather than repeated:
#   08_invariance.csv        gender (men / women) invariance of the 8: thresholds and
#                            loadings both held (H8's procedure).
#   08_exploratory.json      country_rule: country invariance not run (the six countries
#                            173, elsewhere 545; the registered minimum is 200 each).
#   06_h5a.csv, 06_h5b.csv   H5 enters age and gender (man / woman / other) as step-1
#                            covariates; their own coefficients are not reported there.
#   16_spirituality.csv      Spearman and latent r of spirituality and religiosity with the 8.
#   07_known_groups.csv      H6 adjusted for spirituality and religiosity.
#
# Variables (from 01's object; race/ethnicity from the extraction, as 23 reads it):
#   age            z-scored in the model's cases (per SD; the sample SD is about 10 years)
#   gender         woman vs man, other (non-binary or self-described) vs man (H5's coding);
#                  "prefer not to say" leaves the model's cases
#   race/ethnicity the mutually exclusive coding of 23: White only as the reference.
#                  (a) not White (every other answer, incl. more than one selected and
#                  self-described); (b) each category with n >= 30 against White only,
#                  the smaller categories pooled in one dummy that is not interpreted.
#                  Prefer not to say leaves the model's cases. A social category, not biology.
#   education      bachelor's degree or higher vs not (the only SES indicator); also the
#                  five levels as a linear 0-4 code (z), alone
#   country        lives outside the six planned countries (US, UK, CA, AU, IE, NZ) vs in
#                  them; and each named country with n >= 30 against everywhere else
#                  (the other countries were typed as free text and cannot be split)
#   religiosity, spirituality  1-5, z-scored
#   practice       01's practice group: current and past-or-occasional against never
#
# Model sets:
#   alone           each variable (block of dummies) on its own; practice group too, as
#                   the baseline for the two practice sets below
#   joint           age, gender, not White, degree, outside the six, religiosity and
#                   spirituality together
#   joint_practice  joint plus practice group: whether a demographic relation survives
#                   holding practice constant, and practice holding the demographics
#   practice_race   practice group and not White only: the overlap of the two (current
#                   practitioners are mostly not White; see the crosstabs)
#
# Outcomes and estimation:
#   observed   OLS of the observed 8, the 10 and FFMQ Observing, each z-scored in the
#              model's cases: coefficients in outcome SDs (for a 0/1 term a d-like
#              difference over the total SD; for a z-scored term a standardized beta),
#              95% CIs from the t distribution. The 8 minus Observing: the difference of
#              the two coefficients on the same cases, 2,000-resample percentile bootstrap
#              (seed 2026111622).
#   latent     MIMIC in the registered style (07, 16): ordinal items, WLSMV, theta,
#              missing = "pairwise", conditional.x = FALSE; the H1 model (S, loadings
#              free; the agreement factor RI, loadings 1 on the unrecoded statements,
#              variance free) with S regressed on the covariates and S's residual
#              variance 1, so a coefficient is in S's residual SDs (for a 0/1 term a
#              within-group latent d, as 07). Covariates mean-centred in the model's
#              cases (16's finding), cases complete on the covariates.
#                S8           the 8, RI ~ 0*x (registered: RI uncorrelated with covariates)
#                S10          the 10, the same
#                S8_RI_free   the 8 with RI also regressed on the covariates: does the
#                             agreement factor relate to demographics? RI's coefficients are
#                             also given over RI's residual SD (ri_std, delta-method CI)
#                S8_Observing the 8 and FFMQ Observing (04's latent factor of its three
#                             items, residual variance 1) both regressed on the covariates,
#                             their residuals correlated, RI ~ 0 with both; the difference
#                             of the two coefficients is a defined parameter (delta method)
# Holm families: each (outcome or factor, method, model set). In "alone" the family is
# every term of every variable, for that outcome.
# Practice overlap: practice group by each demographic, counts and row %, chi-square
# (Monte Carlo p, 10,000 tables, when an expected count is below 5) and Cramer's V; age,
# religiosity and spirituality by practice group (means, F, eta squared). Race by country
# group, the same way.
#
# Data. The pipeline's saved real-mode objects (04's top-level code, as 10-16 do). One
# fact the derived object does not hold, race/ethnicity, is read from the same extraction
# 01 read (named in the object), from the same submission 01 uses, option ids only: free
# text is never read, printed or written. Nothing from the extraction is written except
# aggregates.
#
# Writes output/real/22_*.csv and 22_exploratory_demographics.json (aggregates only).
# Run from the repo root:  SF2B_MODE=real Rscript study2b/analysis/22_exploratory_demographics.R

# Set Up ---------
## Load libraries ---------
packages <- c("jsonlite", "data.table", "lavaan")
new_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(new_packages)) install.packages(new_packages)
for (thispack in packages) {
  library(thispack, character.only = TRUE, quietly = TRUE, verbose = FALSE)
}

## Settings ---------
source("study2b/analysis/config.R")
source("study2b/analysis/sf2b_models.R")
stopifnot(SF2B$mode == "real", SF2B$first_n == 0)
T0 <- Sys.time()
ex04 <- parse("study2b/analysis/04_validity.R", keep.source = FALSE)
first_fit <- which(vapply(ex04, function(e) grepl("^sf2b_log\\(\"fitting", deparse(e)[1]), TRUE))
stopifnot(length(first_fit) == 1)
for (i in seq_len(first_fit - 1)) eval(ex04[[i]], envir = globalenv())
stopifnot(exists("DAT"), exists("score_items"), !exists("RES"))
ex10 <- parse("study2b/analysis/10_deviation_maia.R", keep.source = FALSE)
is_def <- function(e, nm) is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name(nm))
eval(Filter(function(e) is_def(e, "improper_of"), as.list(ex10))[[1]], envir = globalenv())
# lavaan's phantom factor for an observed covariate (x =~ x) has a standardized loading
# of 1 by construction: not improper (as 16).
improper_list <- function(f) { x <- improper_of(f)$list
  x[!grepl("^std loading ([[:alnum:]_]+) =~ \\1 = ", x)] }
stopifnot(identical(M$participant, S$participant), nrow(DAT) == nrow(M))

I8 <- D$sf8; I10 <- D$sf10; OBS <- "ffmq_observing"; IO <- score_items[[OBS]]
MIN_N <- 30
B_BOOT <- 2000; SEED_BOOT <- 2026111622; SEED_MC <- 2026111623
num <- function(DT, k = 4) { cols <- names(DT)[vapply(DT, is.numeric, TRUE)]; DT[, (cols) := lapply(.SD, function(x) round(x, k)), .SDcols = cols]; DT }
rnd <- function(x) if (is.data.frame(x)) x else if (is.list(x)) lapply(x, rnd) else if (is.numeric(x)) signif(x, 4) else x

## Race/ethnicity from the extraction (option ids only; as 23_demographics.R) ---------
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
RC <- local({
  rows <- fromJSON(raw_file, simplifyVector = FALSE)$rows
  bg <- Filter(function(r) identical(r$source, "questionnaire") && identical(r$slug, "sf-background") && is.null(r$resubmission_of) &&
                 isTRUE(r$participant %in% M$participant) && "bg_practice_regular" %in% names(r$responses), rows)
  rm(rows)
  rec <- rbindlist(lapply(bg, function(r) {
    rc <- r$responses$bg_race; multi <- is.list(rc) && is.null(names(rc))
    ids <- if (multi) unique(na.omit(vapply(rc, function(o) if (is.list(o) && !is.null(o$option_id)) as.character(o$option_id) else NA_character_, ""))) else character()
    g <- r$responses$bg_gender
    data.table(participant = r$participant, received = r$received_at %||% NA_character_, race_ids = paste(sort(ids), collapse = " "),
               race_status = if (multi && length(ids)) "answered" else if (is.character(rc) && identical(rc, "pna")) "prefer_not" else "missing",
               gender_check = if (is.list(g) && !is.null(g$option_id)) as.character(g$option_id) else NA_character_)
  }))
  rec[, received_t := ts_parse(as.list(received))]
  rec <- rec[order(participant, received_t)][!duplicated(participant)]
  rec[match(M$participant, participant)]
})
stopifnot(identical(RC$participant, M$participant), all((is.na(RC$gender_check) & is.na(M$gender)) | RC$gender_check == M$gender))
race_me_of <- function(s, status) {
  if (status != "answered") return(status)
  sub <- setdiff(s, "prefer_not")
  if (!length(sub)) return("prefer_not")
  if (length(sub) >= 2 || "mixed" %in% sub) return("more_than_one")
  if (identical(sub, "other")) return("self_described")
  sub
}
race_me <- mapply(race_me_of, strsplit(RC$race_ids, " ", fixed = TRUE), RC$race_status)
RACE_LABEL <- c(white = "White", black = "Black", mena = "Middle Eastern or North African", south_asian = "South Asian",
                hispanic_latin = "Hispanic or Latino/a/x", east_asian = "East Asian", southeast_asian = "Southeast Asian",
                indigenous = "Indigenous", pacific = "Pacific Islander", more_than_one = "More than one selected",
                self_described = "Self-described only")

## The analysis frame ---------
race_n <- table(race_me[!race_me %in% c("prefer_not", "missing")])
RACE_BIG <- setdiff(names(race_n)[race_n >= MIN_N], "white")
RACE_BIG <- RACE_BIG[order(-race_n[RACE_BIG])]
COUNTRY_BIG <- names(which(table(M$country[M$country %in% SF2B_SIX_COUNTRIES]) >= MIN_N))
COUNTRY_BIG <- COUNTRY_BIG[order(-table(M$country)[COUNTRY_BIG])]
ok_race <- !race_me %in% c("prefer_not", "missing")
dm <- function(cond, ok) fifelse(ok, as.integer(cond), NA_integer_)
EDUC_LEVEL <- c(lt_hs = 0, hs = 1, some_college = 2, bachelors = 3, graduate = 4)
PF <- data.table(
  age = M$age,
  woman = dm(M$gender3 %in% "woman", !is.na(M$gender3)), gender_other = dm(M$gender3 %in% "other", !is.na(M$gender3)),
  not_white = dm(race_me != "white", ok_race),
  degree = dm(M$education %in% c("bachelors", "graduate"), M$education %in% names(EDUC_LEVEL)),
  educ_level = unname(EDUC_LEVEL[M$education]),
  outside_six = dm(M$country %in% "other", M$country %in% c(SF2B_SIX_COUNTRIES, "other")),
  religious = M$religious, spiritual = M$spiritual,
  practice_current = dm(M$practice_group %in% "current", !is.na(M$practice_group)),
  practice_past = dm(M$practice_group %in% "past_or_occasional", !is.na(M$practice_group)))
for (k in RACE_BIG) PF[[paste0("race_", k)]] <- dm(race_me == k, ok_race)
PF[["race_smaller_pooled"]] <- dm(!(race_me %in% c("white", RACE_BIG)), ok_race)
for (k in COUNTRY_BIG) PF[[paste0("country_", k)]] <- dm(M$country %in% k, M$country %in% c(SF2B_SIX_COUNTRIES, "other"))
CONT <- c("age", "educ_level", "religious", "spiritual")

LABEL <- c(age = "age (per SD)", woman = "woman vs man", gender_other = "other gender (non-binary) vs man",
           not_white = "not White vs White only", degree = "bachelor's degree or higher vs not", educ_level = "education level (0-4, per SD)",
           outside_six = "outside the six countries vs in them", religious = "religiosity (per SD)", spiritual = "spirituality (per SD)",
           practice_current = "current practice vs never", practice_past = "past or occasional practice vs never",
           setNames(paste(RACE_LABEL[RACE_BIG], "vs White only"), paste0("race_", RACE_BIG)),
           race_smaller_pooled = "smaller race/ethnicity categories pooled vs White only (not interpreted)",
           setNames(paste(toupper(COUNTRY_BIG), "vs everywhere else"), paste0("country_", COUNTRY_BIG)))
BLOCKS <- list(age = "age", gender = c("woman", "gender_other"), race_white_vs_not = "not_white",
               race_categories = c(paste0("race_", RACE_BIG), "race_smaller_pooled"),
               education_degree = "degree", education_level = "educ_level", country_six = "outside_six",
               country_named = paste0("country_", COUNTRY_BIG), religiosity = "religious", spirituality = "spiritual",
               practice = c("practice_current", "practice_past"))
JOINT <- c("age", "woman", "gender_other", "not_white", "degree", "outside_six", "religious", "spiritual")
SETS <- c(lapply(setNames(names(BLOCKS), names(BLOCKS)), function(b) list(set = "alone", block = b, x = BLOCKS[[b]])),
          list(joint = list(set = "joint", block = "joint", x = JOINT),
               joint_practice = list(set = "joint_practice", block = "joint_practice", x = c(JOINT, "practice_current", "practice_past")),
               practice_race = list(set = "practice_race", block = "practice_race", x = c("practice_current", "practice_past", "not_white"))))

# Covariates for one model: cases complete on them; continuous ones z-scored, every one
# mean-centred, in those cases.
covs_of <- function(x) {
  keep <- complete.cases(PF[, x, with = FALSE])
  X <- as.data.frame(PF[keep, x, with = FALSE])
  for (v in x) X[[v]] <- if (v %in% CONT) (X[[v]] - mean(X[[v]])) / sd(X[[v]]) else X[[v]] - mean(X[[v]])
  list(keep = keep, X = X, n1 = vapply(x, function(v) if (v %in% CONT) NA_integer_ else sum(PF[[v]][keep] == 1L), 0L))
}

## Observed ---------
sf2b_log("22: observed models")
Y <- data.table(sf8 = S$sf8, sf10 = S$sf10, ffmq_observing = S$ffmq_observing)
z <- function(v) (v - mean(v)) / sd(v)
ols <- function(y, X) {
  d <- cbind(y = z(y), X); m <- lm(y ~ ., data = d)
  co <- summary(m)$coefficients; ci <- confint(m)
  list(coef = co[-1, , drop = FALSE], ci = ci[-1, , drop = FALSE], r2 = summary(m)$r.squared)
}
OBSV <- rbindlist(lapply(SETS, function(st) rbindlist(lapply(names(Y), function(o) {
  cv <- covs_of(st$x); ok <- cv$keep & !is.na(Y[[o]])
  X <- cv$X[ok[cv$keep], , drop = FALSE]
  for (v in st$x) X[[v]] <- if (v %in% CONT) z(X[[v]]) else X[[v]] - mean(X[[v]])
  r <- ols(Y[[o]][ok], X)
  data.table(set = st$set, block = st$block, outcome = o, term = st$x, label = unname(LABEL[st$x]), n = sum(ok),
             n_term = vapply(st$x, function(v) if (v %in% CONT) NA_integer_ else sum(PF[[v]][ok] == 1L), 0L),
             est = r$coef[st$x, "Estimate"], lower = r$ci[st$x, 1], upper = r$ci[st$x, 2], p = r$coef[st$x, "Pr(>|t|)"], r2_model = r$r2)
}))))
OBSV[, p_holm := p.adjust(p, "holm"), by = .(set, outcome)]

# The 8 minus Observing, same cases, percentile bootstrap.
sf2b_log("22: observed 8 minus Observing (bootstrap)")
DIFF <- rbindlist(lapply(SETS, function(st) {
  cv <- covs_of(st$x); ok <- cv$keep & !is.na(Y$sf8) & !is.na(Y$ffmq_observing)
  X0 <- as.data.frame(PF[ok, st$x, with = FALSE]); y8 <- Y$sf8[ok]; yo <- Y$ffmq_observing[ok]
  est_of <- function(ix) {
    X <- X0[ix, , drop = FALSE]
    for (v in st$x) X[[v]] <- if (v %in% CONT) z(X[[v]]) else X[[v]] - mean(X[[v]])
    b8 <- coef(lm(y ~ ., data = cbind(y = z(y8[ix]), X)))[st$x]; bo <- coef(lm(y ~ ., data = cbind(y = z(yo[ix]), X)))[st$x]
    rbind(b8, bo, b8 - bo)
  }
  e <- est_of(seq_len(sum(ok)))
  set.seed(SEED_BOOT)
  bt <- replicate(B_BOOT, est_of(sample.int(sum(ok), replace = TRUE))[3, ], simplify = "array")
  bt <- matrix(bt, nrow = length(st$x))
  data.table(set = st$set, block = st$block, term = st$x, label = unname(LABEL[st$x]), n = sum(ok),
             est_8 = e[1, ], est_observing = e[2, ], diff = e[3, ],
             lower = apply(bt, 1, quantile, .025, na.rm = TRUE), upper = apply(bt, 1, quantile, .975, na.rm = TRUE),
             p_boot = apply(bt, 1, function(b) min(1, 2 * min(mean(b <= 0, na.rm = TRUE), mean(b >= 0, na.rm = TRUE)))),
             failed_resamples = apply(bt, 1, function(b) sum(is.na(b))))
}))
DIFF[, p_holm := p.adjust(p_boot, "holm"), by = set]

## Latent (MIMIC) ---------
sf2b_log("22: latent models")
fit_latent <- function(type, st) {
  items <- if (type == "S10") I10 else I8
  cv <- covs_of(st$x); x <- st$x
  d <- cbind(DAT[cv$keep, , drop = FALSE], cv$X)
  syn <- sf2b_scale_syntax(items)
  syn <- c(syn, paste0("S ~ ", paste0("b_", x, "*", x, collapse = " + ")))
  if (type == "S8_RI_free") {
    syn[syn == "RI ~~ RI"] <- "RI ~~ v_ri*RI"
    syn <- c(syn, paste0("RI ~ ", paste0("a_", x, "*", x, collapse = " + ")), paste0("ri_std_", x, " := a_", x, " / sqrt(v_ri)"))
  } else syn <- c(syn, paste0("RI ~ ", paste0("0*", x, collapse = " + ")))
  ord <- items
  if (type == "S8_Observing") {
    syn <- c(syn, sf2b_validity_syntax(OBS, D), paste0(OBS, " ~ ", paste0("o_", x, "*", x, collapse = " + ")),
             paste0("dif_", x, " := b_", x, " - o_", x))
    ord <- c(items, IO)
  }
  f <- sf2b_try(sf2b_sem(syn, d, ordered = ord))
  base <- data.table(model = type, set = st$set, block = st$block, n = nrow(d))
  if (sf2b_failed(f)) return(cbind(base, error = sf2b_msg(f)))
  bad <- sf2b_orientation(f, items, if (type == "S8_Observing") OBS else character(), D)
  pe <- parameterEstimates(f); ss <- standardizedSolution(f)
  pfx <- c(b = "S", a = "RI", ri_std = "RI (over its residual SD)", o = "FFMQ Observing", dif = "S minus Observing")
  rows <- rbindlist(lapply(names(pfx), function(p) {
    lab <- paste0(p, "_", x); r <- pe[match(lab, pe$label), ]
    if (all(is.na(r$est))) return(NULL)
    data.table(factor = pfx[[p]], term = x, label = unname(LABEL[x]), n_term = cv$n1,
               est = r$est, se = r$se, lower = r$ci.lower, upper = r$ci.upper, p = r$pvalue,
               std_all = if (p %in% c("b", "a", "o")) ss$est.std[match(lab, ss$label)] else NA_real_)
  }))
  fs <- sf2b_fit_stats(f)
  cbind(base, rows, error = if (length(bad)) paste("not oriented:", paste(bad, collapse = ", ")) else "",
        improper = paste(improper_list(f), collapse = "; "), cfi = fs$cfi.scaled, rmsea = fs$rmsea.scaled)
}
LAT <- rbindlist(unlist(lapply(c("S8", "S10", "S8_RI_free", "S8_Observing"), function(type)
  lapply(SETS, function(st) { sf2b_log("  ", type, " / ", st$block); fit_latent(type, st) })), recursive = FALSE), fill = TRUE)
LAT[!is.na(p), p_holm := p.adjust(p, "holm"), by = .(model, factor, set)]

## Practice overlap ---------
sf2b_log("22: practice overlap")
age_band <- as.character(cut(M$age, c(18, 25, 35, 45, Inf), c("18-24", "25-34", "35-44", "45+"), right = FALSE))
race_cat <- fifelse(race_me %in% c("white", RACE_BIG), race_me, fifelse(ok_race, "smaller_pooled", NA_character_))
CATS <- list(gender = M$gender3, race_white_vs_not = fifelse(ok_race, fifelse(race_me == "white", "white_only", "not_white"), NA_character_),
             race_categories = race_cat, education_degree = fifelse(is.na(PF$degree), NA_character_, fifelse(PF$degree == 1, "degree", "no_degree")),
             country_six = fifelse(is.na(PF$outside_six), NA_character_, fifelse(PF$outside_six == 1, "outside_six", "six_countries")),
             age_band = age_band)
xtab <- function(a, b, name_a, name_b) {
  ok <- !is.na(a) & !is.na(b); tb <- table(a[ok], b[ok])
  ex <- suppressWarnings(chisq.test(tb, correct = FALSE))
  small <- any(ex$expected < 5)
  ct <- if (small) { set.seed(SEED_MC); suppressWarnings(chisq.test(tb, correct = FALSE, simulate.p.value = TRUE, B = 10000)) } else ex
  v <- sqrt(unname(ex$statistic) / (sum(tb) * (min(dim(tb)) - 1)))
  cells <- as.data.table(tb); setnames(cells, c("level", "group", "n"))
  cells[, `:=`(variable = name_a, against = name_b, row_pct = 100 * n / sum(n)), by = level]
  cells[, col_pct := 100 * n / sum(n), by = group]
  list(cells = cells[, .(variable, level, against, group, n, row_pct, col_pct)],
       test = data.table(variable = name_a, against = name_b, n = sum(tb), chisq = unname(ex$statistic), df = unname(ex$parameter),
                         p = ct$p.value, p_method = if (small) "Monte Carlo (10,000), an expected count < 5" else "asymptotic",
                         cramers_v = v))
}
XT <- lapply(names(CATS), function(v) xtab(CATS[[v]], M$practice_group, v, "practice_group"))
XT <- c(XT, list(xtab(CATS$race_categories, CATS$country_six, "race_categories", "country_six"),
                 xtab(CATS$race_white_vs_not, CATS$country_six, "race_white_vs_not", "country_six")))
XCELLS <- rbindlist(lapply(XT, `[[`, "cells")); XTEST <- rbindlist(lapply(XT, `[[`, "test"))
XTEST[, p_holm := p.adjust(p, "holm")]
CONTP <- rbindlist(lapply(c("age", "religious", "spiritual", "educ_level"), function(v) {
  d <- data.table(y = PF[[v]], g = M$practice_group)[!is.na(y) & !is.na(g)]
  a <- anova(lm(y ~ g, data = d))
  cbind(d[, .(n = .N, mean = mean(y), sd = sd(y)), by = g][order(match(g, c("current", "past_or_occasional", "never")))][, variable := v],
        F = a$`F value`[1], df1 = a$Df[1], df2 = a$Df[2], p = a$`Pr(>F)`[1], eta_sq = a$`Sum Sq`[1] / sum(a$`Sum Sq`))
}))
setcolorder(CONTP, "variable")

## Write ---------
CHOICES <- list(
  race = sprintf("Mutually exclusive coding (23): White only is the reference; categories with n >= %d (%s) each against it; the smaller ones (%s) pooled in one dummy, not interpreted. Prefer not to say (and none answered) leave the model's cases.",
                 MIN_N, paste(RACE_LABEL[RACE_BIG], collapse = ", "), paste(setdiff(names(race_n), c("white", RACE_BIG)), collapse = ", ")),
  gender = "H5's coding (01's gender3): man (reference), woman, other (non-binary or self-described). The 'other' group is small (n in n_term); its estimates are reported but are not interpretable.",
  country = sprintf("Outside the six planned countries vs in them; and the named countries with n >= %d (%s) each against everywhere else. Other countries were typed as free text and are not split.",
                    MIN_N, paste(toupper(COUNTRY_BIG), collapse = ", ")),
  education = "Bachelor's degree or higher vs not (prefer not to say leaves the cases); the five levels as a linear code (z), alone only.",
  scaling = "Observed: outcome z-scored in the model's cases; continuous covariates z-scored in them; 0/1 terms raw (a d-like difference over the outcome's total SD). Latent: S's (and Observing's) residual variance 1, so 0/1 terms are within-group latent ds; continuous covariates z-scored. Observed and latent units differ slightly (total vs residual SD).",
  agreement = "S8_RI_free: RI regressed on the covariates, RI's residual variance labelled v_ri; ri_std = a / sqrt(v_ri), delta-method CI. Not the registered treatment, in which RI is uncorrelated with every covariate.",
  observing = "S8_Observing: Observing is 04's latent factor of ffmq_01, ffmq_06, ffmq_11 (residual variance 1), correlated with S's residual, uncorrelated with RI; dif = b - o (delta method). Observed: the difference of coefficients on z-scored outcomes, same cases, percentile bootstrap.",
  holm = "Holm within each (outcome or factor, method, model set); in 'alone' the family is every term of every variable for that outcome. Exploratory: no decision rests on any p.",
  practice = "Practice group as 01 (current, past or occasional, never); one participant without a regular-practice answer is left out of the crosstabs.")
cited <- list(
  gender_invariance = "08_invariance.csv: gender (man / woman), thresholds then loadings, both held (d_cfi >= -.010, d_rmsea <= .015).",
  country_invariance = sprintf("08_exploratory.json country_rule: not analysed (six countries %d, elsewhere %d; 200 each required).",
                               A$country_rule$groups$six_countries %||% NA, A$country_rule$groups$other %||% NA),
  h5 = "06_h5a.csv, 06_h5b.csv: age and gender are step-1 covariates in H5; their coefficients are not reported there.",
  spirituality = "16_spirituality.csv: Spearman and latent r of spirituality and religiosity with the 8.",
  h6_adjusted = "07_known_groups.csv: H6 adjusted for spirituality and religiosity.")
n_terms <- list(race = as.list(race_n), country_named = as.list(table(M$country)[COUNTRY_BIG]),
                gender = as.list(table(M$gender3, useNA = "ifany")))
num_cols <- function(DT) num(copy(DT))
fwrite(num_cols(OBSV), sf2b_out("22_observed.csv"))
fwrite(num_cols(DIFF), sf2b_out("22_observed_8_minus_observing.csv"))
fwrite(num_cols(LAT), sf2b_out("22_latent.csv"))
fwrite(num_cols(XCELLS), sf2b_out("22_practice_crosstabs.csv"))
fwrite(num_cols(XTEST), sf2b_out("22_practice_tests.csv"))
fwrite(num_cols(CONTP), sf2b_out("22_practice_continuous.csv"))
out <- list(
  analysis = "22: demographics and the 8, the 10, the agreement factor and FFMQ Observing. EXPLORATORY, NOT REGISTERED (osf.io/d4a5u does not include it); specified 16 November 2026, after the registered analyses had been seen.",
  mode = SF2B$mode, tag = SF2B$tag, extraction = A$extraction, written = format(Sys.time(), "%Y-%m-%d %H:%M"),
  software = list(R = R.version.string, lavaan = as.character(packageVersion("lavaan"))),
  n = nrow(M), choices = CHOICES, cited_not_repeated = cited, group_sizes = n_terms,
  seeds = list(bootstrap = SEED_BOOT, monte_carlo = SEED_MC), bootstrap_resamples = B_BOOT,
  observed = OBSV, observed_8_minus_observing = DIFF, latent = LAT,
  practice_overlap = list(tests = XTEST, cells = XCELLS, continuous = CONTP),
  latent_problems = LAT[nzchar(paste0(error, improper)), unique(paste(model, block, error, improper, sep = ": "))],
  run_seconds = round(as.numeric(difftime(Sys.time(), T0, units = "secs"))))
write_json(rnd(out), sf2b_out("22_exploratory_demographics.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
sf2b_log("22 done in ", out$run_seconds, " s; outputs ", sf2b_out("22_*"))
print(num(copy(OBSV), 3)[, .(set, outcome, term, n, n_term, est, lower, upper, p, p_holm)], nrows = 500)
print(num(copy(DIFF), 3)[, .(set, term, n, est_8, est_observing, diff, lower, upper, p_boot, p_holm)], nrows = 500)
print(num(copy(LAT), 3)[, .(model, set, factor, term, n, est, lower, upper, p, p_holm, error, improper)], nrows = 500)
print(num(copy(XTEST), 3)); print(num(copy(CONTP), 3))
