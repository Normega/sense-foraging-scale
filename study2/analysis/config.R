# Sense Foraging Study 2: settings shared by every script.
#
# The preregistration (study2/preregistration_draft.md as of commit 51c41d5; the
# file is now study2/preregistration_2b.md, rewritten for 2b) says this code is the
# specification wherever its text is less precise. Every place where the text
# allowed more than one reading carries a comment starting "PREREG:" that names
# the reading taken; README.md lists them.
#
# MODE (SF2_MODE).
#   "permuted" (the default) reads the real extraction and shuffles every item's
#   answers independently across the participants who answered it, before anything
#   structural is computed. Each item keeps its distribution; every correlation
#   between items, between measures and between occasions is destroyed. Exclusions,
#   attention checks, timing, posting and practice group stay real.
#   "simulated" reads raw/simulated_extract.json, written by 00_simulate.R from a
#   KNOWN model in exactly the shape of the real extraction, so the code path from
#   01 on is identical. 99_check_recovery.R compares what the pipeline recovers
#   with the generating values.
#   "real" reads the real extraction. Only after the registration is filed.
#
# DATA. Participant data never live in the repo. They sit in SF2_DATA_DIR
# (default F:/radlab_data/sf_study2, local and unsynced):
#   raw/       extract_<date>.json (extract.sql via the Supabase CLI),
#              retest_invited_<date>.csv (draw_retest_sample.R),
#              simulated_extract.json and simulated_retest_invited.csv (00_simulate.R)
#   derived/   R objects written by the scripts, per mode
#   output/<mode>/   tables and the aggregates-only report JSON
#
# Run from the repo root:  SF2_MODE=simulated Rscript study2/analysis/run_all.R

# Packages: lavaan (every model), psych, irr (ICC), MASS (polr for H5b; mvrnorm in the
# simulation), data.table, jsonlite, parallel; logistf (H5b's Firth sensitivity
# analysis, called as logistf::logistf) and GPArotation (simulation truth only).
# semTools, which Section 5 lists, is loaded but not needed: the Wu & Estabrook
# invariance models are written out explicitly (sf2_inv_syntax), because the random
# intercept's loadings stay fixed at 1 and its variance free in every group.
suppressPackageStartupMessages({
  library(jsonlite); library(data.table); library(psych); library(lavaan)
  library(semTools); library(irr); library(MASS); library(parallel)
})

SF2 <- list(
  mode      = Sys.getenv("SF2_MODE", "permuted"),
  data_dir  = Sys.getenv("SF2_DATA_DIR", "F:/radlab_data/sf_study2"),
  study_id  = "528674d9-f691-5cff-a81c-1f4f0996b08b",   # study2 ids.json "study"
  retest_id = "65ffadea-be05-56c7-922a-20b78258c6ab",   # study2 ids.json "retest"
  # Materials: the platform definitions, exactly as administered.
  short_form = "study1/materials/sf-short-20.json",
  background = "study1/materials/sf-background.json",
  materials  = "study2/materials",
  # Prolific study id -> posting. The two Study 2 postings do not exist yet: fill
  # these in when they are created (01 stops in real/permuted mode while empty).
  postings_real = c(),
  postings_sim  = c("sim-general-posting" = "general", "sim-contemplative-posting" = "contemplative"),

  # ---- seeds (recorded; cited in the preregistration) ------------------------
  seed            = 20261012,   # general
  seed_split      = 202610127,  # step 7: selection / holdout halves
  seed_retest     = 202610128,  # follow-up invitations (draw_retest_sample.R)
  seed_boot_h3    = 202610123,  # H3b and H3d bootstraps
  seed_boot_step7 = 202610124,  # step 7 bootstrap stability
  seed_sim        = as.integer(Sys.getenv("SF2_SEED_SIM", "20261099")),   # simulated mode only (98 varies it)

  # ---- registered thresholds (Section 5) ---------------------------------------
  speeder_secs_per_item = 2,    # flag: fewer than 2 s per short-form item
  longstring_flag       = 10,   # flag: a run of 10 or more identical answers in the short form
  pna_flag              = 4,    # flag: "Prefer not to answer" on MORE than 4 of the 20 items
  min_prop_answered     = 0.8,  # observed scores need at least 80% of a scale's items
  retest_window         = c(14, 21),  # days; outside it a pair is kept and flagged
  min_category_n        = 5,    # sparse-category rule (invariance and step 7)
  checks = c(maia_attn_check = 2, bfi2_attn_check = 4),   # instructed answer: MAIA "2", BFI-2 "Agree a little" = 4

  # ---- computation ----------------------------------------------------------------
  rstarts    = as.integer(Sys.getenv("SF2_RSTARTS", "30")),      # registered: 30 random starts
  boot_h3    = as.integer(Sys.getenv("SF2_BOOT_H3", "2000")),    # registered: 2,000 (H3b, H3d)
  boot_step7 = as.integer(Sys.getenv("SF2_BOOT_STEP7", "200")),  # registered: 200 (step 7)
  cores      = as.integer(Sys.getenv("SF2_CORES", as.character(max(1L, min(12L, parallel::detectCores() - 2L))))),
  partial_max_free = as.integer(Sys.getenv("SF2_PARTIAL_MAX", "6"))  # partial invariance: most items released per step
)
stopifnot(SF2$mode %in% c("permuted", "simulated", "real"))
SF2$postings <- if (SF2$mode == "simulated") SF2$postings_sim else SF2$postings_real

sf2_path <- function(...) file.path(SF2$data_dir, ...)
sf2_out  <- function(...) {
  d <- sf2_path("output", SF2$mode)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  file.path(d, ...)
}
sf2_derived <- function(name) {
  dir.create(sf2_path("derived"), recursive = TRUE, showWarnings = FALSE)
  sf2_path("derived", sprintf("%s_%s.rds", name, SF2$mode))
}
sf2_log <- function(...) cat(sprintf("[%s %s] ", SF2$mode, format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# The raw extraction this mode reads.
sf2_raw_file <- function() {
  if (SF2$mode == "simulated") return(sf2_path("raw", "simulated_extract.json"))
  f <- tail(sort(list.files(sf2_path("raw"), "^extract_.*\\.json$", full.names = TRUE)), 1)
  if (!length(f)) stop("no raw/extract_*.json in ", sf2_path("raw"), ": run extract.sql first (README)")
  f
}
sf2_invited_file <- function() {
  if (SF2$mode == "simulated") return(sf2_path("raw", "simulated_retest_invited.csv"))
  f <- tail(sort(list.files(sf2_path("raw"), "^retest_invited_.*\\.csv$", full.names = TRUE)), 1)
  if (length(f)) f else NA_character_
}

# ---- the instruments ---------------------------------------------------------------

# Content-factor names used everywhere in the code, in the order of Section 4.
SF2_FACTORS <- c(Access = "Access", `Understanding the modes` = "Understanding",
                 `Valuing and savouring` = "Valuing", `Letting be` = "LettingBe", Forage = "Forage")

sf2_definitions <- function() {
  sf <- fromJSON(SF2$short_form, simplifyVector = FALSE)
  comps <- unlist(lapply(sf$pages, `[[`, "components"), recursive = FALSE)
  sf_items <- vapply(comps, `[[`, "", "id")                          # the fixed presented order
  sf_text  <- setNames(vapply(comps, `[[`, "", "question"), sf_items)
  sf_factor <- unlist(lapply(names(sf$factors), function(f)
    setNames(rep(unname(SF2_FACTORS[f]), length(sf$factors[[f]])), unlist(sf$factors[[f]]))))[sf_items]
  sf_rev <- unlist(sf$reverse_keyed)
  # Section 4's table, checked against the definition.
  stopifnot(setequal(names(sf_factor)[sf_factor == "Access"], sprintf("sf6_%d", c(35, 36, 9, 8, 7))),
            setequal(names(sf_factor)[sf_factor == "Understanding"], sprintf("sf6_%d", c(22, 3, 15, 31, 37))),
            setequal(names(sf_factor)[sf_factor == "Valuing"], sprintf("sf6_%d", c(12, 18, 11, 27, 44))),
            setequal(names(sf_factor)[sf_factor == "LettingBe"], sprintf("sf6_%d", c(40, 46, 43, 47, 21))),
            setequal(sf_rev, sprintf("sf6_%d", c(35, 36, 37, 44, 40, 46, 43, 47))))

  cd <- fromJSON(file.path(SF2$materials, "sf2-sfs-candidates.json"), simplifyVector = FALSE)
  cand <- rbindlist(lapply(cd$candidates, function(c) data.table(
    id = c$id, factor = unname(SF2_FACTORS[c$factor]), reversed = isTRUE(c$reverse_keyed),
    study1_item = c$study1_item %||% NA_character_)))
  ccomps <- unlist(lapply(cd$pages, `[[`, "components"), recursive = FALSE)
  cand[, text := setNames(vapply(ccomps, `[[`, "", "question"), vapply(ccomps, `[[`, "", "id"))[id]]
  stopifnot(nrow(cand) == 24, !anyNA(cand$factor),
            setequal(cand[factor == "Understanding", id], sprintf("sfc_%d", c(2, 3, 7, 8, 9))),
            setequal(cand[factor == "LettingBe", id], sprintf("sfc_%d", c(5, 6, 10:14))),
            setequal(cand[factor == "Valuing" & is.na(study1_item), id], sprintf("sfc_%d", c(4, 15))),
            setequal(cand[factor == "Forage", id], sprintf("sfc_%d", 16:20)),
            setequal(cand[reversed == TRUE, id], sprintf("sfc_%d", c(2, 3, 9, 12, 13, 14, 4, 15, 19, 20))))

  list(sf_items = sf_items, sf_text = sf_text, sf_factor = sf_factor, sf_rev = sf_rev,
       cand = cand, cand_items = cand$id,
       factors4 = c("Access", "Understanding", "Valuing", "LettingBe"),
       validity = sf2_validity_scales(), past_week = sf2_past_week_map())
}

# Validity scales: item ids, keying and scoring per the manuals. Keys that the
# platform definition carries (BFI-2, SBI, PSS-4, IUS-12) are read from it and
# checked; FFMQ-15, MAIA-2 and BIDR-16 keys are from the REB-approved appendix
# (Validity_Scales_Appendix.docx, "Scoring" paragraphs), which the definitions
# were built from. method "sum" = manual sum (prorated as mean x k when items are
# missing, see sf2_score); "mean" = manual mean.
sf2_validity_scales <- function() {
  ids <- function(fmt, v) sprintf(fmt, v)
  s <- list(
    ffmq15 = list(slug = "sf2-ffmq15", range = c(1, 5), method = "sum", items = ids("ffmq_%02d", 1:15),
                  reverse = ids("ffmq_%02d", c(3, 4, 7, 8, 9, 13, 14)),
                  scores = list(ffmq_observing = ids("ffmq_%02d", c(1, 6, 11)), ffmq_describing = ids("ffmq_%02d", c(2, 7, 12)),
                                ffmq_acting = ids("ffmq_%02d", c(3, 8, 13)), ffmq_nonjudging = ids("ffmq_%02d", c(4, 9, 14)),
                                ffmq_nonreactivity = ids("ffmq_%02d", c(5, 10, 15)))),
    maia2 = list(slug = "sf2-maia2", range = c(0, 5), method = "mean", items = ids("maia_%02d", 1:37),
                 reverse = ids("maia_%02d", c(5:12, 15)), check = "maia_attn_check",
                 scores = list(maia_noticing = ids("maia_%02d", 1:4), maia_notdistracting = ids("maia_%02d", 5:10),
                               maia_notworrying = ids("maia_%02d", 11:15), maia_attnregulation = ids("maia_%02d", 16:22),
                               maia_emoawareness = ids("maia_%02d", 23:27), maia_selfregulation = ids("maia_%02d", 28:31),
                               maia_bodylistening = ids("maia_%02d", 32:34), maia_trusting = ids("maia_%02d", 35:37))),
    bidr16 = list(slug = "sf2-bidr16", range = c(1, 7), method = "mean", items = ids("bidr_%02d", 1:16),
                  reverse = ids("bidr_%02d", c(1, 3, 5, 8, 9, 11, 12, 13)),
                  scores = list(bidr_sde = ids("bidr_%02d", 1:8), bidr_im = ids("bidr_%02d", 9:16))),
    phq4 = list(slug = "sf2-phq4", range = c(0, 3), method = "sum", items = ids("phq4_%d", 1:4), reverse = character(),
                scores = list(phq4_total = ids("phq4_%d", 1:4), phq4_anxiety = ids("phq4_%d", 1:2), phq4_depression = ids("phq4_%d", 3:4))),
    ius12 = list(slug = "sf2-ius12", range = c(1, 5), method = "sum", items = ids("ius12_%02d", 1:12), reverse = character(),
                 scores = list(ius12_total = ids("ius12_%02d", 1:12), ius12_prospective = ids("ius12_%02d", c(1, 2, 4, 5, 8, 9, 11)),
                               ius12_inhibitory = ids("ius12_%02d", c(3, 6, 7, 10, 12)))),
    sbi = list(slug = "sf2-sbi-moment", range = c(1, 7), method = "mean", items = ids("sbi_%02d", c(2, 5, 8, 11, 14, 17, 20, 23)),
               reverse = ids("sbi_%02d", c(2, 8, 14, 20)), scores = list(sbi_moment = ids("sbi_%02d", c(2, 5, 8, 11, 14, 17, 20, 23)))),
    bfi2 = list(slug = "sf2-bfi2-om", range = c(1, 5), method = "mean", items = ids("bfi2_%02d", seq(5, 60, 5)),
                reverse = ids("bfi2_%02d", c(5, 25, 30, 45, 50, 55)), check = "bfi2_attn_check",
                scores = list(bfi_openmindedness = ids("bfi2_%02d", seq(5, 60, 5)), bfi_aesthetic = ids("bfi2_%02d", c(5, 20, 35, 50)),
                              bfi_intellectual = ids("bfi2_%02d", c(10, 25, 40, 55)), bfi_creative = ids("bfi2_%02d", c(15, 30, 45, 60)))),
    eq = list(slug = "sf2-eq-decentering", range = c(1, 5), method = "sum", items = ids("eq_%02d", c(3, 15, 9, 14, 10, 16, 6, 17, 18, 12, 20)),
              reverse = character(), scores = list(eq_decentering = ids("eq_%02d", c(3, 15, 9, 14, 10, 16, 6, 17, 18, 12, 20)))),
    pss4 = list(slug = "sf2-pss4", range = c(0, 4), method = "sum", items = ids("pss4_%d", 1:4), reverse = ids("pss4_%d", 2:3),
                scores = list(pss4 = ids("pss4_%d", 1:4)))
  )
  # Check against the platform definitions: every id exists, and keys the definition states agree.
  for (k in names(s)) {
    d <- fromJSON(file.path(SF2$materials, paste0(s[[k]]$slug, ".json")), simplifyVector = FALSE)
    present <- vapply(unlist(lapply(d$pages, `[[`, "components"), recursive = FALSE), `[[`, "", "id")
    stopifnot(all(s[[k]]$items %in% present), is.null(s[[k]]$check) || s[[k]]$check %in% present)
    if (!is.null(d$scoring$reverse_items)) stopifnot(setequal(unlist(d$scoring$reverse_items), s[[k]]$reverse))
    if (!is.null(d$scoring$subscales)) for (sub in d$scoring$subscales)
      stopifnot(any(vapply(s[[k]]$scores, function(v) setequal(v, unlist(sub$item_ids)), TRUE)))
  }
  s
}

# The past-week items are multiple choice: each option's numeric value comes from
# the definition; "nothing_unexpected" has no number and becomes missing (Section 4).
sf2_past_week_map <- function() {
  d <- fromJSON(file.path(SF2$materials, "sf2-past-week.json"), simplifyVector = FALSE)
  comps <- unlist(lapply(d$pages, `[[`, "components"), recursive = FALSE)
  setNames(lapply(comps, function(c) {
    v <- vapply(c$options, function(o) if (is.numeric(o$value)) as.numeric(o$value) else NA_real_, 0)
    setNames(v, vapply(c$options, `[[`, "", "id"))
  }), vapply(comps, `[[`, "", "id"))
}

# Every validity score, its label, and the predicted cell of the validity matrix
# (Section 5, "Predicted validity matrix"; Access, Understanding, Valuing, Letting be).
# Codes: "+mod" moderate, "+small", "+sm-mod" small to moderate, "+mod-lg" moderate
# to large, "0" near zero (|r| < .20), "-small", "-mod", "-sm-mod". NA = not in the
# predicted matrix (reported, no prediction).
SF2_VALIDITY_SCORES <- data.table::rbindlist(list(
  list("ffmq_observing", "FFMQ Observing", "+mod", "+small", "+mod", "+small"),
  list("ffmq_describing", "FFMQ Describing", "+small", "+small", "0", "0"),
  list("ffmq_acting", "FFMQ Acting with awareness", "+small", "0", "0", "+small"),
  list("ffmq_nonjudging", "FFMQ Nonjudging", "+small", "0", "0", "+sm-mod"),
  list("ffmq_nonreactivity", "FFMQ Nonreactivity", "+mod", "+small", "+small", "+mod"),
  list("maia_noticing", "MAIA-2 Noticing", "+mod", "+small", "+small", "0"),
  list("maia_notdistracting", "MAIA-2 Not-Distracting", "0", "0", "0", "+small"),
  list("maia_notworrying", "MAIA-2 Not-Worrying", "0", "0", "0", "+small"),
  list("maia_attnregulation", "MAIA-2 Attention Regulation", "+mod", "+small", "+small", "+small"),
  list("maia_emoawareness", "MAIA-2 Emotional Awareness", NA, NA, NA, NA),
  list("maia_selfregulation", "MAIA-2 Self-Regulation", NA, NA, NA, NA),
  list("maia_bodylistening", "MAIA-2 Body Listening", "+sm-mod", "+small", "+mod", "0"),
  list("maia_trusting", "MAIA-2 Trusting", "+sm-mod", "+small", "+mod", "0"),
  list("eq_decentering", "EQ Decentering", "+mod", "+sm-mod", "+small", "+mod"),
  list("ius12_total", "IUS-12", "-small", "0", "0", "-mod"),
  list("ius12_prospective", "IUS-12 Prospective", NA, NA, NA, NA),
  list("ius12_inhibitory", "IUS-12 Inhibitory", NA, NA, NA, NA),
  list("sbi_moment", "SBI Savoring the Moment", "+mod", "+small", "+mod-lg", "+small"),
  list("bfi_openmindedness", "BFI-2 Open-Mindedness", "+small", "+small", "+small", "+small"),
  list("bfi_aesthetic", "BFI-2 Aesthetic Sensitivity", "+small", "+small", "+mod", "+small"),
  list("bfi_intellectual", "BFI-2 Intellectual Curiosity", NA, NA, NA, NA),
  list("bfi_creative", "BFI-2 Creative Imagination", NA, NA, NA, NA),
  list("bidr_sde", "BIDR SDE", "+small", "0", "0", "+small"),
  list("bidr_im", "BIDR IM", "0", "0", "0", "0"),
  list("phq4_total", "PHQ-4", "-small", "0", "0", "-small"),
  list("phq4_anxiety", "PHQ-4 anxiety", NA, NA, NA, NA),
  list("phq4_depression", "PHQ-4 depression", NA, NA, NA, NA),
  list("pss4", "PSS-4", "-small", "0", "0", "-sm-mod"),
  list("lifesat", "Life satisfaction", "+small", "0", "+small", "+small"),
  list("crit_pause_days", "Past week: pause days", "+sm-mod", "0", "+sm-mod", "0"),
  list("crit_notice_first", "Past week: notice first", "+small", "+small", "0", "+sm-mod"),
  list("crit_outdoors", "Past week: time outdoors", "+small", "0", "+small", "0")
))
data.table::setnames(SF2_VALIDITY_SCORES, c("score", "label", "Access", "Understanding", "Valuing", "LettingBe"))
# Band each code allows (|r| lower, upper, sign); descriptive only.
SF2_BANDS <- list(`+small` = c(.10, .30, 1), `+mod` = c(.30, .50, 1), `+sm-mod` = c(.10, .50, 1), `+mod-lg` = c(.30, 1, 1),
                  `0` = c(0, .20, 0), `-small` = c(.10, .30, -1), `-mod` = c(.30, .50, -1), `-sm-mod` = c(.10, .50, -1))

# Confirmatory cells. H2: factor, measure. H3d: neighbour sets. H4: factor, criterion.
SF2_H2 <- data.table(factor = c("Access", "Access", "LettingBe", "Valuing"),
                     score = c("ffmq_observing", "maia_attnregulation", "eq_decentering", "sbi_moment"))
SF2_H3D <- list(Access = c("ffmq_observing", "maia_noticing", "maia_attnregulation", "eq_decentering"),
                Valuing = c("sbi_moment", "bfi_aesthetic", "maia_bodylistening", "ffmq_observing"))
SF2_H4 <- data.table(factor = c("Access", "Valuing", "LettingBe", "Access"),
                     score = c("crit_pause_days", "crit_pause_days", "crit_notice_first", "crit_notice_first"))

# Agreement index (Section 4, "Indices"): the reversed-partner pairs available among
# the short form and the candidates. PREREG: the pairs are the DESIGNATED ones only:
# Study 1's reversal partners (sf-pool6 interleave.apart_pages) whose both members
# were administered, and the opposite pairs the candidate block was built with
# (sf2-sfs-candidates interleave.apart_pages; its other apart pairs are same-direction
# look-alikes, not reversals). 9 reversed + 9 forward items. Two further opposites are
# obvious but were never designated and are NOT used: sf6_47 / sfc_5 and sf6_40 / sfc_6.
SF2_AGREE_PAIRS <- data.table(
  reversed = c("sf6_36", "sf6_37", "sf6_43", "sf6_44", "sf6_46", "sfc_9", "sfc_12", "sfc_20", "sfc_19"),
  forward  = c("sf6_9",  "sf6_3",  "sf6_11", "sf6_12", "sfc_25", "sfc_8", "sfc_11", "sfc_17", "sfc_16"),
  source   = c(rep("Study 1 partners", 5), rep("candidate opposites", 4)))

# Practice frequency (Section 4, "Indices"; H5b, the H6 trend, step 9c): five ordered
# levels, 0 = no regular practice or none in the past three months, 1 = less than
# weekly, 2 = 1-2 times a week, 3 = 3-5 times a week, 4 = daily or near-daily.
# (01 codes "no regular practice" as 0 too.)
SF2_FREQ_CODE <- c(none = 0, lt_weekly = 1, `1_2_week` = 2, `3_5_week` = 3, daily = 4)

# Country group (Section 4): the six countries the REB application named vs everywhere
# else; if either group has fewer than 200 participants, country is not analysed.
SF2_SIX_COUNTRIES <- c("us", "uk", "ca", "au", "ie", "nz")

# Study 1 values for the five re-run items (step 10): primary keyed loading in Study
# 1's 4-factor ESEM with the random intercept (output/real/04_primary_loadings_k04.csv,
# F1 = Access, F3 = Valuing) and the item mean (02_item_descriptives.csv).
SF2_STUDY1_RERUN <- data.table(
  candidate = c("sfc_21", "sfc_22", "sfc_23", "sfc_24", "sfc_25"),
  study1_item = c("sf6_6", "sf6_2", "sf6_14", "sf6_25", "sf6_28"),
  factor = c("Access", "Access", "Access", "Valuing", "Valuing"),
  study1_loading = c(0.647, 0.573, 0.579, 0.462, 0.421),
  study1_mean = c(3.60, 3.98, 4.30, 4.12, 3.94))

# The 14-item Access + Valuing form (Section 4; step 10), in Study 2's ids.
SF2_AV14 <- list(Access = c("sf6_9", "sf6_7", "sfc_21", "sfc_22", "sf6_8", "sfc_23", "sf6_36", "sf6_35"),
                 Valuing = c("sf6_18", "sf6_12", "sfc_24", "sf6_27", "sfc_25", "sf6_44"))
