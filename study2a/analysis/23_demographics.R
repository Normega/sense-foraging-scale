# 23: demographics of the analysed sample (Table 1). DESCRIPTIVE, NOT PREREGISTERED.
#
# Written 16 November 2026 for the Psychological Assessment manuscript, which requires
# the sex/gender and race/ethnicity of every sample (and SES where available). Nothing
# here changes a registered analysis.
#
# Sample: the analysed 500 (01's `meta`: everyone who submitted the scale and did not
# fail two or more of the three attention checks), the same sample as every registered
# analysis. Study 2a used the same background form as Study 1 and 2b (sf-background,
# the same items and options); 01 keeps every variable except race/ethnicity.
#
# Data. 01's real-mode object, plus the one fact it does not keep (race/ethnicity),
# read from the same extraction 01 read (named in
# the object) and from the same submission 01 uses (each participant's first complete
# sf-background by received_at, flagged resubmissions dropped). Option ids only: the
# free text of "Other (please specify)" and "Prefer to self-describe" is never read into
# the analysis, printed or written; self-described answers are counted. The variables 01
# does keep (age, gender, education, country, regular practice, religiosity,
# spirituality, tradition) are re-derived the same way and checked equal to 01's, so
# this is the same submission.
#
# Race/ethnicity (bg_race, select all that apply) is reported two ways:
#   any      each category as selected by anyone (percentages can sum above 100);
#   mutually exclusive
#            one category per person: a single named category; "more than one" (two
#            or more selections, or "Mixed or multiple"); self-described only; prefer
#            not to say; missing. "Prefer not to say" ticked with a category counts as
#            that category (the count is reported).
# Education is the only SES indicator the background form asks. Country of residence:
# the six planned countries by name, everywhere else as "Other" (typed as free text,
# so the specific countries are not tabulated); the posting was not country-filtered.
# Practice group as 01 defines it.
#
# Writes output/real/study2a_demographics.json and study2a_demographics.csv (aggregates only).
# Run from the repo root:  SF2A_MODE=real Rscript study2a/analysis/23_demographics.R

# Set Up ---------
## Load libraries ---------
packages <- c("jsonlite", "data.table")
new_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(new_packages)) install.packages(new_packages)
for (thispack in packages) {
  library(thispack, character.only = TRUE, quietly = TRUE, verbose = FALSE)
}

## Settings ---------
source("study2a/analysis/config.R")
stopifnot(SF2A$mode == "real")
STUDY <- "study2a"
A <- readRDS(sf2a_derived("sf2a"))
M <- A$meta
N <- nrow(M)
raw_file <- sf2a_path(A$extraction)
stopifnot(file.exists(raw_file), !anyDuplicated(M$participant))

## Codes and labels (sf-background.json) ---------
RACE <- c(east_asian = "East Asian", south_asian = "South Asian", southeast_asian = "Southeast Asian", black = "Black",
          hispanic_latin = "Hispanic or Latino/a/x", indigenous = "Indigenous", mena = "Middle Eastern or North African",
          pacific = "Pacific Islander", white = "White", mixed = "Mixed or multiple", other = "Self-described (Other)",
          prefer_not = "Prefer not to say")
RACE_ME <- c(RACE[setdiff(names(RACE), c("mixed", "other", "prefer_not"))],
             more_than_one = "More than one selected (incl. Mixed or multiple)", self_described = "Self-described only",
             prefer_not = "Prefer not to say")
GENDER <- c(man = "Man", woman = "Woman", non_binary = "Non-binary", self_describe = "Self-described", prefer_not = "Prefer not to say")
EDUC <- c(lt_hs = "Less than high school", hs = "High school or equivalent", some_college = "Some college",
          bachelors = "Bachelor's degree", graduate = "Graduate or professional degree", prefer_not = "Prefer not to say")
COUNTRY <- c(us = "United States", uk = "United Kingdom", ca = "Canada", au = "Australia", ie = "Ireland", nz = "New Zealand",
             other = "Other country (free text, not tabulated)", prefer_not = "Prefer not to answer")
TRADITION <- c(none = "No religious/spiritual affiliation", christian = "Christian", muslim = "Muslim", jewish = "Jewish",
               buddhist = "Buddhist", hindu = "Hindu", sikh = "Sikh", sbnr = "Spiritual but not religious", other = "Other",
               prefer_not = "Prefer not to say")
PRACTICE <- c(current = "Current (regular, 1-2 a week or more)", past_or_occasional = "Past or occasional (regular, less than weekly or none recently)",
              never = "Never practised regularly")
SIX <- c("us", "uk", "ca", "au", "ie", "nz")
AGE_BANDS <- c(18, 25, 35, 45, 55, 65, Inf)
AGE_LABELS <- c("18-24", "25-34", "35-44", "45-54", "55-64", "65+")

## Background records from the extraction (option ids only) ---------
# The token a declined answer carries; any other string is not read.
oid <- function(v) {
  if (is.list(v) && !is.null(v$option_id)) return(as.character(v$option_id))
  if (is.character(v) && length(v) == 1 && v %in% c("pna", "prefer_not")) return("prefer_not")
  NA_character_
}
likert <- function(v) if (is.numeric(v) && length(v) == 1) as.numeric(v) else if (identical(oid(v), "prefer_not")) -1 else NA_real_
bg_record <- function(r) {
  b <- r$responses; rc <- b$bg_race
  multi <- is.list(rc) && is.null(names(rc))
  ids <- if (multi) unique(na.omit(vapply(rc, function(o) if (is.list(o) && !is.null(o$option_id)) as.character(o$option_id) else NA_character_, ""))) else character()
  data.table(participant = r$participant, received = r$received_at %||% NA_character_,
             age = if (is.list(b$bg_age) && identical(oid(b$bg_age), "age")) suppressWarnings(as.numeric(b$bg_age$value)) else NA_real_,
             age_status = if (is.null(b$bg_age)) NA_character_ else oid(b$bg_age),
             gender = oid(b$bg_gender), education = oid(b$bg_education), country = oid(b$bg_country),
             tradition = oid(b$bg_tradition), practice_regular = oid(b$bg_practice_regular),
             religious = likert(b$bg_religious), spiritual = likert(b$bg_spiritual),
             race_ids = paste(sort(ids), collapse = " "),
             race_status = if (multi && length(ids)) "answered" else if (identical(oid(rc), "prefer_not")) "prefer_not" else "missing")
}
# 01's pick(): slug sf-background, not a flagged resubmission, bg_practice_regular
# present, the earliest received_at (as a string, as 01 orders it).
BG <- local({
  rows <- fromJSON(raw_file, simplifyVector = FALSE)$rows
  bg <- Filter(function(r) !is.null(r$slug) && identical(r$slug, "sf-background") && is.null(r$resubmission_of) &&
                 isTRUE(r$participant %in% M$participant) && "bg_practice_regular" %in% names(r$responses), rows)
  rm(rows)
  rec <- rbindlist(lapply(bg, bg_record))
  rec <- rec[order(participant, received)][!duplicated(participant)]
  rec[match(M$participant, participant)]
})
stopifnot(identical(BG$participant, M$participant))

## Same submission as 01 ---------
same <- function(a, b) all((is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b))
pna <- function(x) fifelse(x %in% "pna", "prefer_not", x)   # 01 keeps a declined answer's token as is
stopifnot(same(BG$gender, pna(M$gender)), same(BG$education, pna(M$education)), same(BG$country, pna(M$country)),
          same(BG$practice_regular, pna(M$practice_regular)), same(BG$age, M$age), same(BG$tradition, pna(M$tradition)),
          same(fifelse(BG$religious < 0, NA_real_, BG$religious), M$religious),
          same(fifelse(BG$spiritual < 0, NA_real_, BG$spiritual), M$spiritual))

## Tables ---------
pct <- function(n) round(100 * n / N, 1)
cat_table <- function(x, labels, variable) {
  x[is.na(x)] <- "missing"
  bad <- setdiff(unique(x), c(names(labels), "missing"))
  if (length(bad)) stop(variable, ": unexpected codes ", paste(bad, collapse = ", "))
  labels <- c(labels, missing = "Missing (no answer)")
  n <- vapply(names(labels), function(k) sum(x == k), 0L)
  data.table(variable = variable, code = names(labels), category = unname(labels), n = n, pct = pct(n), value = NA_real_)
}
cont_table <- function(x, variable, declined) {
  ok <- x[!is.na(x) & x >= 0]
  data.table(variable = variable, code = c("n_answered", "mean", "sd", "median", "min", "max", "prefer_not", "missing"),
             category = c("Answered", "Mean", "SD", "Median", "Minimum", "Maximum", "Prefer not to answer", "Missing (no answer)"),
             n = c(length(ok), rep(NA_integer_, 5), sum(declined), N - length(ok) - sum(declined)),
             pct = c(pct(length(ok)), rep(NA_real_, 5), pct(sum(declined)), pct(N - length(ok) - sum(declined))),
             value = c(NA, round(mean(ok), 2), round(sd(ok), 2), median(ok), min(ok), max(ok), NA, NA))
}

race_sets <- strsplit(BG$race_ids, " ", fixed = TRUE)
race_any <- rbindlist(lapply(names(RACE), function(k) {
  n <- sum(vapply(race_sets, function(s) k %in% s, TRUE))
  data.table(variable = "race_any", code = k, category = RACE[[k]], n = n, pct = pct(n), value = NA_real_)
}))
race_any <- rbind(race_any, data.table(variable = "race_any", code = "missing", category = "Missing (no answer)",
                                       n = sum(BG$race_status == "missing"), pct = pct(sum(BG$race_status == "missing")), value = NA_real_))
race_me_of <- function(s, status) {
  if (status == "prefer_not") return("prefer_not")
  if (status == "missing") return("missing")
  sub <- setdiff(s, "prefer_not")
  if (!length(sub)) return("prefer_not")
  if (length(sub) >= 2 || "mixed" %in% sub) return("more_than_one")
  if (identical(sub, "other")) return("self_described")
  sub
}
BG[, race_me := mapply(race_me_of, race_sets, race_status)]
BG[, n_race_selected := vapply(race_sets, function(s) length(setdiff(s, c("prefer_not", ""))), 0L)]
race_me <- cat_table(fifelse(BG$race_me == "missing", NA_character_, BG$race_me), RACE_ME, "race_exclusive")

BG[, age_declined := age_status %in% "prefer_not" & is.na(age)]
age_band <- as.character(cut(BG$age, AGE_BANDS, AGE_LABELS, right = FALSE))
age_band[BG$age_declined] <- "prefer_not"
practice <- M$practice_group
degree <- fifelse(BG$education %in% c("bachelors", "graduate"), "degree",
          fifelse(BG$education %in% c("lt_hs", "hs", "some_college"), "no_degree", BG$education))
country_six <- fifelse(BG$country %in% SIX, "six", BG$country)

TAB <- rbind(
  cont_table(BG$age, "age", BG$age_declined),
  cat_table(age_band, setNames(c(AGE_LABELS, "Prefer not to answer"), c(AGE_LABELS, "prefer_not")), "age_band"),
  cat_table(BG$gender, GENDER, "gender"),
  race_any, race_me,
  cat_table(BG$education, EDUC, "education"),
  cat_table(degree, c(degree = "Bachelor's degree or higher", no_degree = "No degree", prefer_not = "Prefer not to say"), "education_degree"),
  cat_table(BG$country, COUNTRY, "country"),
  cat_table(country_six, c(six = "The six planned countries (US, UK, CA, AU, IE, NZ)", other = "Elsewhere", prefer_not = "Prefer not to answer"), "country_six"),
  cat_table(practice, PRACTICE, "practice_group"),
  cont_table(BG$religious, "religiosity_1to5", BG$religious %in% -1),
  cont_table(BG$spiritual, "spirituality_1to5", BG$spiritual %in% -1),
  cat_table(BG$tradition, TRADITION, "tradition"))
TAB[, study := STUDY]
setcolorder(TAB, "study")

## Checks on the tables ---------
for (v in setdiff(unique(TAB$variable), c("race_any", "age", "religiosity_1to5", "spirituality_1to5")))
  stopifnot(sum(TAB[variable == v, n]) == N)
qa <- list(
  race_prefer_not_with_a_category = sum(vapply(race_sets, function(s) "prefer_not" %in% s && length(setdiff(s, "prefer_not")) > 0, TRUE)),
  race_selections_per_person = as.list(table(factor(BG$n_race_selected, levels = 0:max(3, max(BG$n_race_selected))))),
  self_described_counts = list(gender = sum(BG$gender %in% "self_describe"), race = sum(vapply(race_sets, function(s) "other" %in% s, TRUE)),
                               country_other = sum(BG$country %in% "other")))

## Write ---------
as_list <- function(v) { x <- TAB[variable == v]
  if (all(is.na(x$value))) setNames(lapply(seq_len(nrow(x)), function(i) list(category = x$category[i], n = x$n[i], pct = x$pct[i])), x$code)
  else setNames(lapply(seq_len(nrow(x)), function(i) if (is.na(x$value[i])) list(n = x$n[i], pct = x$pct[i]) else x$value[i]), x$code) }
out <- list(
  study = STUDY, description = "Demographics of the analysed sample (Table 1). Descriptive, not preregistered; aggregates only.",
  script = "study2a/analysis/23_demographics.R", extraction = A$extraction, analysed_n = N,
  sample = "01's analysed sample: submitted the scale and did not fail two or more of the three attention checks (500).",
  notes = list(
    percentages = "Of the analysed sample, including missing and prefer-not rows.",
    race = "race_any: each category as selected by anyone (select all that apply; sums above 100%). race_exclusive: one category per person; 'more than one' combines two or more selections with 'Mixed or multiple'; 'Prefer not to say' ticked with a category counts as that category.",
    ses = "Education is the only SES indicator the background form asks.",
    country = "Country of residence: the six planned countries by name; other countries were typed as free text and are counted, not tabulated. Study 2a's posting was not country-filtered.",
    practice = "01's practice_group: current = regular practice and 1-2 a week or more over the past 3 months; past_or_occasional = regular practice, less often; never = no regular practice.",
    free_text = "Self-described answers (gender, race/ethnicity 'Other', country 'Other') are counted only; their text is never read."),
  qa = qa)
for (v in unique(TAB$variable)) out[[v]] <- as_list(v)
write_json(out, sf2a_out("study2a_demographics.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
fwrite(TAB, sf2a_out("study2a_demographics.csv"))
sf2a_log("23: demographics of ", N, " written (aggregates only)")
print(TAB[, .(variable, code, n, pct, value)], nrows = 500)
