# Sense Foraging Study 2a: settings shared by every script.
#
# The analysis registered at osf.io/vt7gn (study2a/preregistration_draft.md is the
# text it was filed from; OSF governs): the factor structure of the revised
# 58-statement set, a check that the two modes were understood, and the rules
# that choose the statements Study 2b will test. Section numbers in comments
# refer to that registration. PREREG: marks a reading of the text that it does
# not settle; each is listed in the README.
#
# MODE. "permuted" (the default) shuffles every item's answers independently
# across participants (statements, vignettes, IUS-12, the Study 1 wording items):
# each keeps its distribution, correlations are destroyed. "simulated" replaces
# the statements with data from a KNOWN model (the 5 hypothesized factors,
# reverse-keyed statements loading negatively, factor correlations .3, a
# random-intercept agreement factor with variance .10), cut to the real
# statements' marginal proportions, and the vignettes with ratings from a known
# separation; 99_check_recovery.R compares what is recovered with the truth.
# "real" uses the answers.
#
# DATA. Never in the repo: SF2A_DATA_DIR (default F:/radlab_data/sf_study2a),
# local and unsynced. The extraction is the newest responses_<date>.json there
# (extract.sql via the Supabase CLI); derived/ and output/<mode>/ are written by
# these scripts.
#
# Run from the repo root:  Rscript study2a/analysis/run_all.R

suppressPackageStartupMessages({
  library(jsonlite); library(data.table); library(psych); library(lavaan)
  library(GPArotation); library(EFAtools)
})

SF2A <- list(
  mode     = Sys.getenv("SF2A_MODE", "permuted"),
  data_dir = Sys.getenv("SF2A_DATA_DIR", "F:/radlab_data/sf_study2a"),
  seed     = 20261030,
  # The live questionnaire definitions (radlab scripts/sense_foraging/study2a/ at
  # 1e23b565, the build that ran): statement wording, hypothesized factors,
  # keying, look-alike pairs, ordering rules.
  scale     = "study2a/materials/sf2a-scale.json",
  vignettes = "study2a/materials/sf2a-vignettes.json",
  items     = "study2a/materials/items.json",
  prolific_study = c("6ae3b35e1cebf4b7733f555e" = "main"),   # the one general posting (pilot + main)
  # Registered thresholds (Section 5, "Data exclusion" and "Decision rules").
  exclude_checks_failed = 2,        # two or more of the three checks
  speeder_secs_per_item = 2,        # on the 61 scale positions
  longstring_flag       = 15,
  pna_flag              = 11,       # more than 11 of the 58 statements declined
  separation_cut        = 1,        # a person tells the modes apart at >= 1 point
  vignette_share        = .75,      # ... and the distinction came across at >= 75%
  min_vignettes_each    = 3,
  load_min = .40, load_pref = .50, cross_max = .30, gap_min = .20,
  pna_max_pct = 5, agree_lo = 15, agree_hi = 85,
  per_factor = c(5, 6), min_reversed = 2,
  # Number of content factors in the primary model: the decision K, made from
  # 04's output by the registered criteria and recorded here with its reason.
  # DECIDED 2026-10-30 (Norm): K = 4. Criterion 1 (every factor >= 3 statements
  # with primary loadings >= .40) is met only at k = 4, in the geomin EFAs
  # (k = 4-7; parallel analysis said 7) and in the ESEM with the random
  # intercept (at k = 5 the hypothesized Valuing factor has 2). Criterion 3, fit,
  # is about equal (k = 4 CFI .953, RMSEA .033; k = 5 .958, .032). On criterion
  # 2: a broad factor of 26 forward statements from all five hypothesized
  # factors, and three smaller factors defined mostly by reverse-keyed
  # statements. The k = 5 ESEM is kept as an exploratory comparison.
  k_primary = as.integer(Sys.getenv("SF2A_K", "4")),
  rstarts = 30,
  pa_datasets = 500
)
stopifnot(SF2A$mode %in% c("permuted", "simulated", "real"))

sf2a_path <- function(...) file.path(SF2A$data_dir, ...)
sf2a_out  <- function(...) {
  d <- sf2a_path("output", SF2A$mode)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  file.path(d, ...)
}
sf2a_log <- function(...) cat(sprintf("[%s] ", SF2A$mode), ..., "\n", sep = "")
sf2a_derived <- function(name) {
  dir.create(sf2a_path("derived"), showWarnings = FALSE)
  sf2a_path("derived", sprintf("%s_%s.rds", name, SF2A$mode))
}

# The definitions: statements, factors, keying, look-alike pairs, the statements
# that name a mode, vignette kinds, the trial's Access + Valuing form.
sf2a_definition <- function() {
  s <- fromJSON(SF2A$scale, simplifyVector = FALSE)
  comps <- unlist(lapply(s$pages, `[[`, "components"), recursive = FALSE)
  comps <- Filter(function(c) c$type != "information", comps)
  text <- setNames(vapply(comps, `[[`, "", "question"), vapply(comps, `[[`, "", "id"))
  meta <- rbindlist(lapply(s$items, function(i) data.table(item = i$id, source = i$source, factor = i$factor,
                                                           reversed = isTRUE(i$reverse_keyed),
                                                           study1_item = i$study1_item %||% NA_character_)))
  items <- meta$item
  pairs <- t(vapply(s$interleave$apart_pages, function(p) c(a = p[[1]], b = p[[2]]), c("", "")))
  # The statements that name a mode: the ordering rule's group with min_gap 3.
  names_mode <- unlist(Filter(function(g) isTRUE(g$min_gap == 3), s$interleave$groups)[[1]]$items)
  v <- fromJSON(SF2A$items, simplifyVector = FALSE)
  vkind <- setNames(ifelse(grepl("^Sensing", vapply(v$vignettes, `[[`, "", "kind")), "sensing", "doing"),
                    sprintf("vig_%s", sub("^V", "", vapply(v$vignettes, `[[`, "", "id"))))
  list(items = items, text = text[items], factor = setNames(meta$factor, items),
       reversed = meta[reversed == TRUE, item], source = setNames(meta$source, items),
       study1_item = setNames(meta$study1_item, items), pairs = pairs, names_mode = names_mode,
       factors = names(s$factors), checks = unlist(s$attention_checks),
       positions = c(items, names(s$attention_checks)),
       vignettes = c(names(vkind), "vig_self"), vkind = vkind,
       ius = sprintf("ius12_%02d", 1:12),
       end = c("sf6_7", "sf6_8", "sf6_12", "sf6_18", "sf6_27", "sf6_11"),
       # The trial's 10-statement Access + Valuing form in Study 1 wording (Section 4,
       # "Indices"): four are in the scale in that wording, six come at the end.
       av10 = c(sf6_9 = "sf2a_01", sf6_36 = "sf2a_02", sf6_35 = "sf2a_03", sf6_44 = "sf2a_27",
                sf6_7 = "sf6_7", sf6_8 = "sf6_8", sf6_12 = "sf6_12", sf6_18 = "sf6_18", sf6_27 = "sf6_27", sf6_11 = "sf6_11"),
       av10_reversed = c("sf6_35", "sf6_36", "sf6_44"))
}
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
