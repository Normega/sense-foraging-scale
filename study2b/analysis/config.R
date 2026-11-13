# Sense Foraging Study 2b: settings shared by every script.
#
# The analysis the Study 2b preregistration specifies (study2/preregistration_2b.md;
# once filed, OSF governs). Section numbers in comments refer to it. The registration
# says this code is the specification wherever its text is less precise; every place
# where the text allows more than one reading carries a comment starting "PREREG:"
# that names the reading taken, and README.md lists them.
#
# MODE (SF2B_MODE).
#   "permuted" (the default) reads the real extraction and shuffles every item's
#   answers independently across the participants who answered that questionnaire,
#   before anything structural is computed: each item keeps its distribution, every
#   correlation between items and between measures is destroyed. Exclusions,
#   attention checks, timing, background answers and the pause stay real.
#   "simulated" reads raw/simulated_extract.json, written by 00_simulate.R from a
#   KNOWN model in exactly the shape of the real extraction, so the code path from 01
#   on is identical. 99_check_recovery.R compares what is recovered with the truth.
#   "real" reads the real extraction. Only after the registration is filed.
#
# DATA. Participant data never live in the repo. They sit in SF2B_DATA_DIR (default
# F:/radlab_data/sf_study2b, local and unsynced):
#   raw/       extract_<date>.json (extract.sql via the Supabase CLI);
#              simulated_extract.json (00_simulate.R)
#   derived/   R objects written by the scripts, per mode
#   output/<mode>/   tables and the aggregates-only report JSON
#
# Run from the repo root:  SF2B_MODE=simulated Rscript study2b/analysis/run_all.R

suppressPackageStartupMessages({
  library(jsonlite); library(data.table); library(psych); library(lavaan); library(MASS); library(parallel)
})

SF2B <- list(
  mode     = Sys.getenv("SF2B_MODE", "permuted"),
  data_dir = Sys.getenv("SF2B_DATA_DIR", "F:/radlab_data/sf_study2b"),
  study_id = "1864d663-d01c-5d82-91f0-3047462470c3",   # radlab scripts/sense_foraging/study2b/ids.json "study"
  # The platform definitions as built (radlab scripts/sense_foraging/study2b/ at c016b10d;
  # the shared validity measures from study2/materials, identical to the platform's copies).
  materials = "study2b/materials",

  # ---- seeds (recorded; cited in the registration) ---------------------------------
  seed      = 20261110,   # general (permutation)
  seed_boot = 202611103,  # H3d bootstrap
  seed_sim  = as.integer(Sys.getenv("SF2B_SEED_SIM", "20261199")),   # simulated mode only

  # ---- registered thresholds (Section 5) ---------------------------------------------
  checks = c(sf2b_attn_agree = 5, sf2b_attn_disagree = 2, bfi2_attn_check = 4),
  exclude_checks_failed = 2,       # two or more of the three
  speeder_secs_per_item = 2,       # flag: < 2 s per statement on the scale's 14 positions
  longstring_flag       = 10,      # flag: a run of 10+ identical answers on the 14 positions
  pna_flag              = 2,       # flag: "Prefer not to answer" on MORE than 2 of the 13 statements
  min_answered = c(sf8 = 7, sf10 = 9),   # observed scores (Section 4)
  min_prop_answered = 0.8,         # PREREG: validity scales, Study 2's 80% rule (the text says "per their manuals")
  separation_cut = 1, min_vignettes_each = 3,
  pause_min_secs = 20, pause_max_away = .5,   # engaged pause (decision 22)
  # H1 (decisions 17, 33)
  h1_cfi = .95, h1_resid = .10, h1_d_omega = .03, h1_d_lambda = .05, h1_closefit_alpha = .05,
  # H2-H3
  h2_lower = .30, h3a_upper = .80, h3c_bound = .30, h3d_upper = .80, h3e_bsss = .40, h3e_bound = .30,
  h3b_needed = 4,                  # of the five contrasts
  h5_delta_r2 = .02, h6_d = .50,
  inv_d_cfi = -.010, inv_d_rmsea = .015, min_category_n = 5, partial_max_free = 4,
  country_min = 200,

  # ---- computation ------------------------------------------------------------------
  boot_h3d = as.integer(Sys.getenv("SF2B_BOOT", "2000")),   # registered: 2,000
  cores    = as.integer(Sys.getenv("SF2B_CORES", as.character(max(1L, min(20L, parallel::detectCores() - 2L)))))
)
stopifnot(SF2B$mode %in% c("permuted", "simulated", "real"))

sf2b_path <- function(...) file.path(SF2B$data_dir, ...)
sf2b_out  <- function(...) {
  d <- sf2b_path("output", SF2B$mode)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  file.path(d, ...)
}
sf2b_derived <- function(name) {
  dir.create(sf2b_path("derived"), recursive = TRUE, showWarnings = FALSE)
  sf2b_path("derived", sprintf("%s_%s.rds", name, SF2B$mode))
}
sf2b_log <- function(...) cat(sprintf("[%s %s] ", SF2B$mode, format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
r3 <- function(x) round(x, 3)

sf2b_raw_file <- function() {
  if (SF2B$mode == "simulated") return(sf2b_path("raw", "simulated_extract.json"))
  f <- tail(sort(list.files(sf2b_path("raw"), "^extract_.*\\.json$", full.names = TRUE)), 1)
  if (!length(f)) stop("no raw/extract_*.json in ", sf2b_path("raw"), ": run extract.sql first (README)")
  f
}
sf2b_mat <- function(f) fromJSON(file.path(SF2B$materials, f), simplifyVector = FALSE)
sf2b_ids <- function(d) {
  comps <- unlist(lapply(d$pages, `[[`, "components"), recursive = FALSE)
  vapply(Filter(function(c) c$type != "information", comps), `[[`, "", "id")
}

# The session (radlab build_study2b.py): step 0 background, 1 vignettes, 2 the scale,
# 3 the trial statements, 4-12 the first half of the validity block (shuffled), 13 the
# pause offer, 14-19 the second half (shuffled), 20 debrief.
SF2B_STEPS <- list(half_a = 4:12, pause = 13L, half_b = 14:19)

# ---- the scale ---------------------------------------------------------------------------
sf2b_definitions <- function() {
  s <- sf2b_mat("sf2b-scale.json")
  meta <- rbindlist(lapply(s$items, function(i) data.table(item = i$id, source = i$source, role = i$role,
                                                           reversed = isTRUE(i$reverse_keyed))))
  positions <- sf2b_ids(s)                                   # the 14 positions, presented (fixed) order
  items <- meta$item
  sf8 <- unlist(s$scoring$sf8); sf10 <- unlist(s$scoring$sf10)
  # Section 4 / Appendix A, checked against the definition.
  stopifnot(length(items) == 13, length(positions) == 14,
            identical(sf8, sprintf("sf2b_%02d", c(1, 2, 3, 4, 6, 7, 9, 10))),
            identical(sf10, sprintf("sf2b_%02d", c(1, 2, 3, 4, 6, 7, 8, 9, 10, 12))),
            setequal(meta[reversed == TRUE, item], sprintf("sf2b_%02d", c(2, 5, 7, 9, 11, 13))),
            setequal(meta[role == "exploratory", item], sprintf("sf2b_%02d", c(5, 11, 13))),
            identical(names(s$attention_checks), "sf2b_attn_agree"), s$attention_checks[[1]] == SF2B$checks[["sf2b_attn_agree"]])
  e <- sf2b_mat("sf2b-end.json")
  end_items <- setdiff(sf2b_ids(e), "sf2b_attn_disagree")
  stopifnot(setequal(end_items, c("sf6_9", "sf6_35", "sf6_7", "sf6_8", "sf6_12", "sf6_18", "sf6_27", "sf6_11")),
            identical(names(e$attention_checks), "sf2b_attn_disagree"), e$attention_checks[[1]] == SF2B$checks[["sf2b_attn_disagree"]])
  v <- sf2b_mat("sf2b-vignettes.json")
  vk <- setNames(ifelse(vapply(v$vignettes, function(x) isTRUE(x$sensing), TRUE), "sensing", "doing"),
                 vapply(v$vignettes, `[[`, "", "id"))
  list(items = items, positions = positions, source = setNames(meta$source, items), role = setNames(meta$role, items),
       reversed = meta[reversed == TRUE, item], sf8 = sf8, sf10 = sf10, exploratory = meta[role == "exploratory", item],
       end_items = end_items, vignettes = c(names(vk), "vig_self"), vkind = vk,
       # The trial's 10-statement Access + Valuing score in its own wording (Section 4):
       # two in the scale (sf6_36 = sf2b_02, sf6_44 = sf2b_11), eight at the end.
       av10 = c(sf6_9 = "sf6_9", sf6_36 = "sf2b_02", sf6_35 = "sf6_35", sf6_7 = "sf6_7", sf6_8 = "sf6_8",
                sf6_12 = "sf6_12", sf6_18 = "sf6_18", sf6_27 = "sf6_27", sf6_11 = "sf6_11", sf6_44 = "sf2b_11"),
       av10_reversed = c("sf6_36", "sf6_35", "sf6_44"),
       validity = sf2b_validity_scales(), past_week = sf2b_past_week_map())
}

# ---- validity measures ----------------------------------------------------------------------
# Item ids, keying and scoring per the manuals. Keys the platform definition states are
# checked against it. method "sum" = manual sum (prorated as mean x k when items are
# missing), "mean" = manual mean.
sf2b_validity_scales <- function() {
  ids <- function(fmt, v) sprintf(fmt, v)
  maia <- function(v) ids("maia_%02d", v)
  s <- list(
    ffmq15 = list(slug = "sf2-ffmq15", range = c(1, 5), method = "sum", items = ids("ffmq_%02d", 1:15),
                  reverse = ids("ffmq_%02d", c(3, 4, 7, 8, 9, 13, 14)),
                  scores = list(ffmq_observing = ids("ffmq_%02d", c(1, 6, 11)), ffmq_describing = ids("ffmq_%02d", c(2, 7, 12)),
                                ffmq_acting = ids("ffmq_%02d", c(3, 8, 13)), ffmq_nonjudging = ids("ffmq_%02d", c(4, 9, 14)),
                                ffmq_nonreactivity = ids("ffmq_%02d", c(5, 10, 15)))),
    maia = list(slug = "sf2b-maia2-brief", range = c(0, 5), method = "mean",
                items = maia(c(2, 3, 4, 6, 9, 10, 11, 12, 15, 17, 19, 20, 25, 26, 27, 29:37)),
                reverse = maia(c(6, 9, 10, 11, 12, 15)),
                scores = list(maia_noticing = maia(2:4), maia_notdistracting = maia(c(6, 9, 10)),
                              maia_notworrying = maia(c(11, 12, 15)), maia_attnregulation = maia(c(17, 19, 20)),
                              maia_emoawareness = maia(25:27), maia_selfregulation = maia(29:31),
                              maia_bodylistening = maia(32:34), maia_trusting = maia(35:37))),
    bidr16 = list(slug = "sf2-bidr16", range = c(1, 7), method = "mean", items = ids("bidr_%02d", 1:16),
                  reverse = ids("bidr_%02d", c(1, 3, 5, 8, 9, 11, 12, 13)),
                  scores = list(bidr_sde = ids("bidr_%02d", 1:8), bidr_im = ids("bidr_%02d", 9:16))),
    phq4 = list(slug = "sf2-phq4", range = c(0, 3), method = "sum", items = ids("phq4_%d", 1:4), reverse = character(),
                scores = list(phq4_total = ids("phq4_%d", 1:4))),
    sbi = list(slug = "sf2-sbi-moment", range = c(1, 7), method = "mean", items = ids("sbi_%02d", c(2, 5, 8, 11, 14, 17, 20, 23)),
               reverse = ids("sbi_%02d", c(2, 8, 14, 20)), scores = list(sbi_moment = ids("sbi_%02d", c(2, 5, 8, 11, 14, 17, 20, 23)))),
    bfi2 = list(slug = "sf2-bfi2-om", range = c(1, 5), method = "mean", items = ids("bfi2_%02d", seq(5, 60, 5)),
                reverse = ids("bfi2_%02d", c(5, 25, 30, 45, 50, 55)), check = "bfi2_attn_check",
                scores = list(bfi_aesthetic = ids("bfi2_%02d", c(5, 20, 35, 50)), bfi_intellectual = ids("bfi2_%02d", c(10, 25, 40, 55)),
                              bfi_creative = ids("bfi2_%02d", c(15, 30, 45, 60)), bfi_openmindedness = ids("bfi2_%02d", seq(5, 60, 5)))),
    eq = list(slug = "sf2-eq-decentering", range = c(1, 5), method = "sum", items = ids("eq_%02d", c(3, 15, 9, 14, 10, 16, 6, 17, 18, 12, 20)),
              reverse = character(), scores = list(eq_decentering = ids("eq_%02d", c(3, 15, 9, 14, 10, 16, 6, 17, 18, 12, 20)))),
    pss4 = list(slug = "sf2-pss4", range = c(0, 4), method = "sum", items = ids("pss4_%d", 1:4), reverse = ids("pss4_%d", 2:3),
                scores = list(pss4 = ids("pss4_%d", 1:4))),
    bsss = list(slug = "sf2b-bsss8", range = c(1, 5), method = "mean", items = ids("bsss_%d", 1:8), reverse = character(),
                scores = list(bsss8 = ids("bsss_%d", 1:8))),
    ncs = list(slug = "sf2b-ncs6", range = c(1, 5), method = "mean", items = ids("ncs_%02d", c(1, 2, 3, 4, 11, 15)),
               reverse = ids("ncs_%02d", c(3, 4)), scores = list(ncs6 = ids("ncs_%02d", c(1, 2, 3, 4, 11, 15)))),
    spsq = list(slug = "sf2b-spsq-sd", range = c(1, 7), method = "mean", items = ids("spsq_%02d", c(33, 7, 42, 25)),
                reverse = character(), scores = list(spsq_sd = ids("spsq_%02d", c(33, 7, 42, 25)))),
    asi = list(slug = "sf2b-asi3-pc", range = c(0, 4), method = "sum", items = ids("asi_%02d", c(3, 4, 7, 8, 12, 15)),
               reverse = character(), scores = list(asi_pc = ids("asi_%02d", c(3, 4, 7, 8, 12, 15)))),
    samd = list(slug = "sf2-samd-md", range = c(1, 7), method = "mean", items = ids("samd_md_%02d", 1:13),
                reverse = character(), scores = list(samd_md = ids("samd_md_%02d", 1:13)))
  )
  for (k in names(s)) {
    d <- sf2b_mat(paste0(s[[k]]$slug, ".json"))
    present <- sf2b_ids(d)
    stopifnot(all(s[[k]]$items %in% present), is.null(s[[k]]$check) || s[[k]]$check %in% present)
    if (!is.null(d$scoring$reverse_items)) stopifnot(setequal(unlist(d$scoring$reverse_items), s[[k]]$reverse))
    subs <- d$scoring$subscales
    if (is.list(subs) && length(subs) && is.null(names(subs)))     # [{name, item_ids}]
      for (sub in subs) stopifnot(any(vapply(s[[k]]$scores, function(v) setequal(v, unlist(sub$item_ids)), TRUE)))
    if (k == "maia") for (sub in subs) stopifnot(any(vapply(s[[k]]$scores, function(v) setequal(v, unlist(sub)), TRUE)))
  }
  s
}

# Past-week items are multiple choice: each option's number comes from the definition;
# "nothing_unexpected" has none and becomes missing (Section 4).
sf2b_past_week_map <- function() {
  d <- sf2b_mat("sf2-past-week.json")
  comps <- unlist(lapply(d$pages, `[[`, "components"), recursive = FALSE)
  setNames(lapply(comps, function(c) {
    v <- vapply(c$options, function(o) if (is.numeric(o$value)) as.numeric(o$value) else NA_real_, 0)
    setNames(v, vapply(c$options, `[[`, "", "id"))
  }), vapply(comps, `[[`, "", "id"))
}

# ---- latent validity factors (Section 5, H2-H3) -------------------------------------------------
# Each validity measure enters as latent factors of its items: a measure with subscales
# as correlated factors of that measure. Factor names are the score names. The BFI-2
# Open-Mindedness total (H3b contrast 3) is a single factor of all 12 items, in a model
# of its own (PREREG: the text says "BFI-2 Open-Mindedness, total" and does not say how
# a total enters a latent model; one factor of its items is the reading taken).
SF2B_MEASURES <- list(
  ffmq = c("ffmq_observing", "ffmq_describing", "ffmq_acting", "ffmq_nonjudging", "ffmq_nonreactivity"),
  maia = c("maia_noticing", "maia_notdistracting", "maia_notworrying", "maia_attnregulation", "maia_emoawareness",
           "maia_selfregulation", "maia_bodylistening", "maia_trusting"),
  bfi = c("bfi_aesthetic", "bfi_intellectual", "bfi_creative"),
  bfi_total = "bfi_openmindedness",
  bidr = c("bidr_sde", "bidr_im"),
  sbi = "sbi_moment", eq = "eq_decentering", samd = "samd_md", bsss = "bsss8", ncs = "ncs6",
  spsq = "spsq_sd", asi = "asi_pc", phq = "phq4_total", pss = "pss4")

# The predicted validity matrix (Section 5): latent r with the 8. Band = (lower, upper)
# of |r| implied by the label, for description only; the decisions use the cells below.
SF2B_PREDICTED <- data.table::rbindlist(list(
  list("ffmq_observing", "FFMQ Observing", "+ moderate (.30-.60)"),
  list("ffmq_describing", "FFMQ Describing", "+ small"),
  list("ffmq_acting", "FFMQ Acting with awareness", "+ small"),
  list("ffmq_nonjudging", "FFMQ Nonjudging", "+ small"),
  list("ffmq_nonreactivity", "FFMQ Nonreactivity", "+ small to moderate"),
  list("sbi_moment", "SBI Savoring the Moment", "+ moderate"),
  list("bfi_aesthetic", "BFI-2 Aesthetic Sensitivity", "+ moderate"),
  list("bfi_intellectual", "BFI-2 Intellectual Curiosity", "+ small"),
  list("bfi_openmindedness", "BFI-2 Open-Mindedness, total", "+ small to moderate, below SBI"),
  list("bfi_creative", "BFI-2 Creative Imagination", "+ small"),
  list("maia_noticing", "Brief MAIA-2 Noticing", "+ moderate"),
  list("maia_notworrying", "Brief MAIA-2 Not-Worrying", "near zero"),
  list("maia_attnregulation", "Brief MAIA-2 Attention Regulation", "+ small to moderate"),
  list("maia_bodylistening", "Brief MAIA-2 Body Listening", "+ small to moderate"),
  list("maia_emoawareness", "Brief MAIA-2 Emotional Awareness", "+ small to moderate"),
  list("maia_selfregulation", "Brief MAIA-2 Self-Regulation", "+ small to moderate"),
  list("maia_notdistracting", "Brief MAIA-2 Not-Distracting", "+ small"),
  list("maia_trusting", "Brief MAIA-2 Trusting", "+ small"),
  list("eq_decentering", "EQ Decentering", "+ small to moderate"),
  list("samd_md", "SAMD Meaning-Detecting", "+ small to moderate, below FFMQ Observing"),
  list("bsss8", "BSSS-8", "+ small"),
  list("ncs6", "NCS-6", "near zero"),
  list("spsq_sd", "SPSQ-SF Sensory Discomfort", "near zero"),
  list("asi_pc", "ASI-3 Physical Concerns", "near zero"),
  list("bidr_sde", "BIDR Self-Deceptive Enhancement", "near zero to + small"),
  list("bidr_im", "BIDR Impression Management", "near zero to + small"),
  list("phq4_total", "PHQ-4", "near zero to - small"),
  list("pss4", "PSS-4", "near zero to - small"),
  list("lifesat", "Life satisfaction", "+ small"),
  list("crit_pause_days", "Past week: pause days", "+ small to moderate"),
  list("engaged_pause", "Engaged pause", "+ small"),
  list("crit_notice_first", "Past week: notice first", "+ small"),
  list("crit_outdoors", "Past week: time outdoors", "+ small")))
data.table::setnames(SF2B_PREDICTED, c("score", "label", "prediction"))

# Confirmatory cells.
SF2B_H2 <- c("ffmq_observing", "sbi_moment", "bfi_aesthetic", "maia_noticing")
# H3b: first > second; `model` names the model holding both (a measure's own model, or a pair).
SF2B_H3B <- data.table(contrast = 1:5,
  first  = c("bfi_aesthetic", "ffmq_observing", "sbi_moment", "maia_noticing", "ffmq_observing"),
  second = c("bfi_intellectual", "ffmq_describing", "bfi_openmindedness", "maia_notworrying", "samd_md"),
  model  = c("bfi", "ffmq", "pair", "maia", "pair"))
SF2B_H3C <- c("bidr_sde", "bidr_im")
SF2B_H3E <- data.table(score = c("bsss8", "ncs6", "spsq_sd", "asi_pc"),
                       lower = c(-Inf, -.30, -.30, -.30), upper = c(.40, .30, .30, .30))

# Practice frequency (H5b): 0 = no regular practice or none in the past three months,
# 1 = less than weekly, 2 = 1-2 a week, 3 = 3-5 a week, 4 = daily or near-daily.
SF2B_FREQ_CODE <- c(none = 0, lt_weekly = 1, `1_2_week` = 2, `3_5_week` = 3, daily = 4)
SF2B_SIX_COUNTRIES <- c("us", "uk", "ca", "au", "ie", "nz")
