# 01: from the raw extraction to one analysis row per participant, applying the
# registered rules (Section 5, "Data exclusion" and "Missing data").
#
#   excluded  staff test enrollments; sessions without a submitted short form;
#             participants who answered BOTH instructed-response checks (MAIA-2, BFI-2)
#             and failed both. A check never reached is not a failure.
#   first     each questionnaire's first complete submission is used; rows the platform
#             flagged as resubmissions (resubmission_of set) are dropped first.
#   flagged   failing one check, or the only check reached; < 2 s per short-form item;
#             a longstring of 10 or more within the short form (fixed presented order;
#             "Prefer not to answer" breaks a run); "Prefer not to answer" on more than
#             4 of the 20 items. Flags are for sensitivity analyses, never exclusions.
#   samples   structure (H1, invariance, steps 7-10): everyone not excluded;
#             validity (H2-H5): those who reached the end of the validity block;
#             known groups (H6): background answers and the short form.
#   follow-up retest rows are paired through external_identity_links (the extract does
#             the join); pairs outside 14-21 days are kept and flagged.
# In permuted mode every item's answers are shuffled across the participants who
# answered that questionnaire BEFORE any item-derived index (longstring, scores,
# agreement index) is computed. The "Prefer not to answer" count, exclusions, checks,
# timing, posting and background answers stay real, as in Study 1.

source("study2/analysis/config.R")
source("study2/analysis/sf2_models.R")
D <- sf2_definitions()
if (SF2$mode != "simulated" && !length(SF2$postings))
  stop("config.R: set SF2$postings_real (Prolific study id -> posting) once the two postings exist")

raw_file <- sf2_raw_file()
ext <- fromJSON(raw_file, simplifyVector = FALSE)
rows <- ext$rows
sf2_log("extraction: ", basename(raw_file), ", ", length(rows), " rows")

# The session's structure (build_study2.py): step 0 background, 1 short form, 2 candidates,
# 3-13 the validity block (11 steps in a per-participant order), 14 debrief.
VAL_FIRST <- 3L; VAL_LAST <- 13L
VAL_SLUGS <- c(vapply(D$validity, `[[`, "", "slug"), "sf2-past-week")
PW_ITEMS <- names(D$past_week)

ts_parse <- function(x) {
  x <- unlist(lapply(x, function(v) if (is.null(v)) NA_character_ else v))
  off <- ifelse(grepl("[+-][0-9]{2}(:?[0-9]{2})?$", x), sub(".*([+-][0-9]{2})(:?([0-9]{2}))?$", "\\1\\3", x), "+00")
  base <- sub("[+-][0-9]{2}(:?[0-9]{2})?$", "", sub("Z$", "", x))
  t <- as.POSIXct(base, tz = "UTC", format = "%Y-%m-%d %H:%M:%OS")
  h <- as.numeric(substr(off, 1, 3)); m <- suppressWarnings(as.numeric(substr(off, 4, 5))); m[is.na(m)] <- 0
  t - sign(h) * (abs(h) * 3600 + m * 60)
}
num <- function(v) if (is.null(v) || is.list(v) || is.character(v)) NA_real_ else as.numeric(v)
opt <- function(v) if (is.list(v) && !is.null(v$option_id)) v$option_id else NA_character_
fld <- function(r, k) { v <- r[[k]]; if (is.null(v)) NA else v }

src <- vapply(rows, function(r) fld(r, "source"), ""); stu <- vapply(rows, function(r) fld(r, "study"), "")

# ---- main study ------------------------------------------------------------------------------
enr <- rows[src == "enrollment" & stu == "main"]
meta_all <- rbindlist(lapply(enr, function(r) data.table(
  participant = r$participant, is_test = isTRUE(r$is_test), prolific_study_id = fld(r, "prolific_study_id"),
  repository_consent = isTRUE(r$repository_consent), consented = !is.null(r$consent_date),
  schedule_id = fld(r, "schedule_id"), session_status = fld(r, "session_status"),
  session_completed = !is.null(r$session_completed_at), step_crashes = as.numeric(fld(r, "step_crashes")))))
meta_all <- meta_all[!duplicated(participant)]
meta_all[, posting := unname(SF2$postings[prolific_study_id])]

main_rows <- rows[src %in% c("questionnaire", "vas", "step") & stu == "main"]
by_person <- split(main_rows, vapply(main_rows, function(r) r$participant, ""))

first_complete <- function(rs, slug, need) {
  cand <- Filter(function(r) identical(r$source, "questionnaire") && identical(r$slug, slug) &&
                   is.null(r$resubmission_of) && all(need %in% names(r$responses)), rs)
  if (!length(cand)) return(NULL)
  cand[[order(ts_parse(lapply(cand, `[[`, "received_at")))[1]]]
}
vas_slugs <- unique(unlist(lapply(main_rows, function(r) if (identical(r$source, "vas")) r$slug)))
if (length(vas_slugs) > 1) stop("more than one VAS scale in Study 2 sessions (", paste(vas_slugs, collapse = ", "),
                                "): the session has one, life satisfaction; check the extraction")

subs_log <- list()
recs <- lapply(meta_all$participant, function(p) {
  rs <- by_person[[p]] %||% list()
  qs <- Filter(function(r) identical(r$source, "questionnaire"), rs)
  slugs <- vapply(qs, `[[`, "", "slug")
  subs_log[[p]] <<- data.table(participant = p, slug = slugs, flagged_resubmission = vapply(qs, function(r) !is.null(r$resubmission_of), TRUE))
  out <- list(participant = p)
  bg <- first_complete(rs, "sf-background", "bg_practice_regular")
  sf <- first_complete(rs, "sf-short-20", D$sf_items)
  cd <- first_complete(rs, "sf2-sfs-candidates", D$cand_items)
  out$has_background <- !is.null(bg); out$has_sf <- !is.null(sf); out$has_candidates <- !is.null(cd)
  if (!is.null(bg)) {
    b <- bg$responses
    out$age <- if (identical(opt(b$bg_age), "age")) num(b$bg_age$value) else NA_real_
    out$gender <- opt(b$bg_gender); out$country <- opt(b$bg_country); out$education <- opt(b$bg_education)
    out$practice_regular <- opt(b$bg_practice_regular); out$practice_freq <- opt(b$bg_practice_freq)
    out$practice_years <- if (identical(opt(b$bg_practice_years), "years")) num(b$bg_practice_years$value) else NA_real_
    out$practice_training <- opt(b$bg_practice_training)
    types <- if (is.list(b$bg_practice_types)) vapply(b$bg_practice_types, function(t) t$option_id %||% NA_character_, "") else character()
    out$practice_types <- paste(types, collapse = " ")
    out$spiritual <- num(b$bg_spiritual); out$religious <- num(b$bg_religious)
    out$tradition <- opt(b$bg_tradition)
  }
  if (!is.null(sf)) {
    s <- sf$responses
    out$sf <- vapply(D$sf_items, function(i) num(s[[i]]), 0)
    out$sf_pna <- sum(vapply(D$sf_items, function(i) identical(s[[i]], "pna"), TRUE))
    out$sf_secs <- as.numeric(s[["_timing"]]$total_ms %||% NA) / 1000
    out$sf_received <- sf$received_at
    out$presentation_format <- sf$presentation_format %||% NA_character_
  }
  if (!is.null(cd)) {
    s <- cd$responses
    out$cand <- vapply(D$cand_items, function(i) num(s[[i]]), 0)
    out$cand_pos <- vapply(D$cand_items, function(i) as.numeric(s[["_order"]]$positions[[i]] %||% NA), 0)
    out$cand_secs <- as.numeric(s[["_timing"]]$total_ms %||% NA) / 1000
  }
  out$val <- list(); out$val_step <- c(); out$val_secs <- c()
  for (k in names(D$validity)) {
    v <- D$validity[[k]]
    q <- first_complete(rs, v$slug, c(v$items, v$check))
    if (is.null(q)) next
    out$val[[k]] <- vapply(v$items, function(i) num(q$responses[[i]]), 0)
    if (!is.null(v$check)) out$val[[paste0(k, "_check")]] <- num(q$responses[[v$check]])
    out$val_step[v$slug] <- as.numeric(q$step_index %||% NA)
    out$val_secs[v$slug] <- as.numeric(q$responses[["_timing"]]$total_ms %||% NA) / 1000
  }
  pw <- first_complete(rs, "sf2-past-week", PW_ITEMS)
  if (!is.null(pw)) {
    out$val$past_week <- vapply(PW_ITEMS, function(i) { o <- opt(pw$responses[[i]]); if (is.na(o)) NA_real_ else unname(D$past_week[[i]][o]) }, 0)
    out$val_step["sf2-past-week"] <- as.numeric(pw$step_index %||% NA)
    out$val_secs["sf2-past-week"] <- as.numeric(pw$responses[["_timing"]]$total_ms %||% NA) / 1000
  }
  vas <- Filter(function(r) identical(r$source, "vas") && is.null(r$resubmission_of), rs)
  if (length(vas)) {
    v1 <- vas[[order(ts_parse(lapply(vas, `[[`, "received_at")))[1]]]
    out$lifesat <- num(v1$value); out$val_step["lifesat"] <- as.numeric(v1$step_index %||% NA)
  }
  steps <- Filter(function(r) identical(r$source, "step"), rs)
  out$max_step_exited <- if (length(steps)) max(vapply(steps, function(r) as.numeric(r$step_index %||% NA), 0), na.rm = TRUE) else NA_real_
  vs <- Filter(function(r) grepl("^vas_", r$slug %||% ""), steps)
  if (length(vs)) out$val_secs["lifesat"] <- as.numeric(difftime(ts_parse(vs[[1]]$exited_at), ts_parse(vs[[1]]$entered_at), units = "secs"))
  out
})
names(recs) <- meta_all$participant
subs_log <- rbindlist(subs_log)

scalar <- c("has_background", "has_sf", "has_candidates", "age", "gender", "country", "education", "practice_regular",
            "practice_freq", "practice_years", "practice_training", "practice_types", "spiritual", "religious", "tradition",
            "sf_pna", "sf_secs", "sf_received", "presentation_format", "cand_secs", "lifesat", "max_step_exited")
M <- rbindlist(lapply(recs, function(r) { x <- r[intersect(scalar, names(r))]; x$participant <- r$participant; as.data.table(x) }), fill = TRUE)
M <- merge(meta_all, M, by = "participant", sort = FALSE)
setcolorder(M, "participant")

# ---- exclusions -------------------------------------------------------------------------------------
chk <- function(r, k) { v <- r$val[[paste0(k, "_check")]]; if (is.null(v)) NA_real_ else v }
M[, maia_check := vapply(recs[participant], chk, 0, k = "maia2")]
M[, bfi_check := vapply(recs[participant], chk, 0, k = "bfi2")]
M[, checks_reached := (!is.na(maia_check)) + (!is.na(bfi_check))]
M[, checks_failed := (!is.na(maia_check) & maia_check != SF2$checks[["maia_attn_check"]]) +
                     (!is.na(bfi_check) & bfi_check != SF2$checks[["bfi2_attn_check"]])]
M[, `:=`(excl_test = is_test, excl_no_sf = !is_test & !has_sf,
         excl_attention = !is_test & has_sf & checks_reached == 2 & checks_failed == 2)]
M[, analysable := !excl_test & !excl_no_sf & !excl_attention]
M[, reached_end_validity := vapply(recs[participant], function(r) {
  any(r$val_step >= VAL_LAST, na.rm = TRUE) || isTRUE(r$max_step_exited >= VAL_LAST) }, TRUE) | session_completed]

flow <- rbindlist(list(
  list("enrollments in the main study", nrow(M)),
  list("excluded: staff test enrollment", sum(M$excl_test)),
  list("excluded: no submitted short form", sum(M$excl_no_sf)),
  list("excluded: answered both checks and failed both", sum(M$excl_attention)),
  list("analysable: structure, invariance, steps 7-10 (H1)", sum(M$analysable)),
  list("  of whom stopped before the end of the validity block", sum(M$analysable & !M$reached_end_validity)),
  list("analysable: validity (H2-H5), reached the end of the validity block", sum(M$analysable & M$reached_end_validity))))
setnames(flow, c("step", "n"))

# ---- item matrices ------------------------------------------------------------------------------------
mat <- function(field, cols, who) {
  X <- matrix(NA_real_, length(who), length(cols), dimnames = list(who, cols))
  for (p in who) { v <- recs[[p]][[field]]; if (!is.null(v)) X[p, ] <- v[cols] }
  X
}
vmat <- function(k, cols, who) {
  X <- matrix(NA_real_, length(who), length(cols), dimnames = list(who, cols))
  for (p in who) { v <- recs[[p]]$val[[k]]; if (!is.null(v)) X[p, ] <- v[cols] }
  X
}
who <- M[has_sf == TRUE & !is_test, participant]
SFX <- mat("sf", D$sf_items, who); CX <- mat("cand", D$cand_items, who); CPOS <- mat("cand_pos", D$cand_items, who)
VX <- lapply(setNames(names(D$validity), names(D$validity)), function(k) vmat(k, D$validity[[k]]$items, who))
VX$past_week <- vmat("past_week", PW_ITEMS, who)
LS <- vapply(who, function(p) recs[[p]]$lifesat %||% NA_real_, 0)

if (SF2$mode == "permuted") {
  set.seed(SF2$seed)
  perm <- function(X, has) { for (j in seq_len(ncol(X))) X[has, j] <- X[has, j][sample.int(sum(has))]; X }
  SFX <- perm(SFX, rep(TRUE, nrow(SFX)))
  CX <- perm(CX, vapply(who, function(p) !is.null(recs[[p]]$cand), TRUE))
  for (k in names(VX)) VX[[k]] <- perm(VX[[k]], vapply(who, function(p) !is.null(recs[[p]]$val[[k]]), TRUE))
  hl <- !is.na(LS); LS[hl] <- LS[hl][sample.int(sum(hl))]
  sf2_log("every item's answers permuted across the participants who answered it (structure destroyed by design)")
}

# Item-derived flags, on whichever matrix the mode uses.
ls_run <- apply(SFX, 1, function(v) { r <- rle(ifelse(is.na(v), -1, v)); max(c(0, r$lengths[r$values != -1])) })
M[, longstring := NA_real_][match(who, participant), longstring := ls_run]
M[, `:=`(flag_one_check = checks_failed == 1,
         flag_speeder = !is.na(sf_secs) & sf_secs / length(D$sf_items) < SF2$speeder_secs_per_item,
         flag_longstring = !is.na(longstring) & longstring >= SF2$longstring_flag,
         flag_pna = !is.na(sf_pna) & sf_pna > SF2$pna_flag)]
M[, flagged := flag_one_check | flag_speeder | flag_longstring | flag_pna]

# ---- background indices -----------------------------------------------------------------------------------
# PREREG: someone who practises regularly but whose years, frequency or training is
# unknown cannot be shown to be dedicated and is "other" (everyone else); a declined
# or missing regular-practice answer leaves the group missing.
M[, practice_group := fifelse(practice_regular %in% "no", "never",
                       fifelse(practice_regular %in% "yes",
                               fifelse(practice_freq %in% c("3_5_week", "daily") & !is.na(practice_years) & practice_years >= 1 &
                                         practice_training %in% "yes", "dedicated", "other"), NA_character_))]
M[, meditation_or_yoga := grepl("\\b(meditation|yoga)\\b", practice_types)]
M[, freq_code := fifelse(practice_regular %in% "no", 0, unname(SF2_FREQ_CODE[practice_freq]))]
M[, gender2 := fifelse(gender %in% c("man", "woman"), gender, NA_character_)]
M[, gender3 := fifelse(gender %in% c("man", "woman"), gender, fifelse(gender %in% c("non_binary", "self_describe"), "other", NA_character_))]
M[, country2 := fifelse(country %in% SF2_SIX_COUNTRIES, "six_countries", fifelse(country %in% "other", "other", NA_character_))]
M[, posting_contemplative := as.integer(posting == "contemplative")]

A <- M[analysable == TRUE]
# Country group (Section 4): the six countries the REB application named vs everywhere
# else ("other"); if either group has fewer than 200 participants, country is not
# analysed. Counted in the analysable sample; a declined country answer is missing.
cn2 <- table(A$country2)
country_rule <- list(groups = as.list(cn2), analysed = length(cn2) == 2 && all(cn2 >= 200))
A[, country_group := if (country_rule$analysed) country2 else NA_character_]

ids <- A$participant
SFX <- SFX[ids, , drop = FALSE]; CX <- CX[ids, , drop = FALSE]; CPOS <- CPOS[ids, , drop = FALSE]
VX <- lapply(VX, function(X) X[ids, , drop = FALSE]); LS <- LS[ids]

# ---- observed scores ---------------------------------------------------------------------------------------
fac_items <- function(f) D$sf_items[D$sf_factor == f]
S <- data.table(participant = ids,
  Access = sf2_score(SFX, fac_items("Access"), D$sf_rev), Understanding = sf2_score(SFX, fac_items("Understanding"), D$sf_rev),
  Valuing = sf2_score(SFX, fac_items("Valuing"), D$sf_rev), LettingBe = sf2_score(SFX, fac_items("LettingBe"), D$sf_rev),
  AccessValuing = sf2_score(SFX, c(fac_items("Access"), fac_items("Valuing")), D$sf_rev))
for (k in names(D$validity)) { v <- D$validity[[k]]
  for (s in names(v$scores)) S[[s]] <- sf2_score(VX[[k]], v$scores[[s]], v$reverse, v$range, v$method) }
# BIDR-16 dichotomous (Paulhus): after reversal, 6 or 7 counts 1; the count per subscale (0-8), 80% rule.
for (s in c("bidr_sde", "bidr_im")) {
  it <- D$validity$bidr16$scores[[s]]; Xb <- VX$bidr16[, it]; rv <- intersect(it, D$validity$bidr16$reverse)
  Xb[, rv] <- 8 - Xb[, rv]; Xb <- (Xb >= 6) * 1
  S[[paste0(s, "_dichotomous")]] <- ifelse(rowSums(!is.na(Xb)) >= ceiling(.8 * length(it)), rowMeans(Xb, na.rm = TRUE) * length(it), NA)
}
for (i in PW_ITEMS) S[[i]] <- VX$past_week[, i]
S$lifesat <- LS
# Agreement index: mean raw answer over the designated reversed-partner pairs where
# both members were answered (SF2_AGREE_PAIRS); variants for H7 (short-form pairs
# only, the follow-up has no candidates) and step 9a (without the Forage pairs).
RAW <- cbind(SFX, CX)
agree <- function(pairs) {
  both <- sapply(seq_len(nrow(pairs)), function(k) !is.na(RAW[, pairs$reversed[k]]) & !is.na(RAW[, pairs$forward[k]]))
  num <- sapply(seq_len(nrow(pairs)), function(k) RAW[, pairs$reversed[k]] + RAW[, pairs$forward[k]])
  num[!both] <- NA
  ifelse(rowSums(both) > 0, rowSums(num, na.rm = TRUE) / (2 * rowSums(both)), NA)
}
sf_pairs <- SF2_AGREE_PAIRS[grepl("^sf6_", reversed) & grepl("^sf6_", forward)]
S$agreement <- agree(SF2_AGREE_PAIRS)
S$agreement_sf <- agree(sf_pairs)
S$agreement_noforage <- agree(SF2_AGREE_PAIRS[!(forward %in% D$cand[factor == "Forage", id] | reversed %in% D$cand[factor == "Forage", id])])

# Presented position of each validity step (1-11) and seconds per questionnaire.
VPOS <- rbindlist(lapply(ids, function(p) { s <- recs[[p]]$val_step
  as.data.table(as.list(c(participant = p, setNames(as.character(s - VAL_FIRST + 1), names(s))))) }), fill = TRUE)
for (j in setdiff(names(VPOS), "participant")) set(VPOS, j = j, value = as.numeric(VPOS[[j]]))
VSECS <- rbindlist(lapply(ids, function(p) as.data.table(as.list(c(participant = p, recs[[p]]$val_secs)))), fill = TRUE)
for (j in setdiff(names(VSECS), "participant")) set(VSECS, j = j, value = as.numeric(VSECS[[j]]))

# ---- the follow-up -----------------------------------------------------------------------------------------
re <- rows[stu == "retest"]
renr <- rbindlist(lapply(Filter(function(r) identical(r$source, "enrollment"), re), function(r) data.table(
  retest_code = r$retest_code, participant = fld(r, "participant"), n_links = as.numeric(fld(r, "n_links")),
  is_test = isTRUE(r$is_test))))
renr <- renr[!duplicated(retest_code)]
rq <- split(Filter(function(r) identical(r$source, "questionnaire"), re), vapply(Filter(function(r) identical(r$source, "questionnaire"), re), `[[`, "", "retest_code"))
renr[, has_sf := retest_code %in% names(Filter(function(rs) !is.null(first_complete(rs, "sf-short-20", D$sf_items)), rq))]
RT <- t(vapply(renr$retest_code, function(rc) {
  s <- first_complete(rq[[rc]] %||% list(), "sf-short-20", D$sf_items)
  if (is.null(s)) rep(NA_real_, length(D$sf_items)) else vapply(D$sf_items, function(i) num(s$responses[[i]]), 0) }, numeric(length(D$sf_items))))
colnames(RT) <- D$sf_items; rownames(RT) <- renr$retest_code
renr[, received := vapply(retest_code, function(rc) { s <- first_complete(rq[[rc]] %||% list(), "sf-short-20", D$sf_items)
  if (is.null(s)) NA_character_ else s$received_at }, "")]
renr[, `:=`(excl_test = is_test, excl_unlinked = !is_test & (is.na(participant) | n_links != 1),
            excl_no_sf = !is_test & !has_sf)]
renr[, excl_baseline := !excl_test & !excl_unlinked & !(participant %in% ids)]
renr[, paired := !excl_test & !excl_unlinked & !excl_no_sf & !excl_baseline]
P2 <- renr[paired == TRUE]
P2 <- P2[!duplicated(participant)]            # a person linked twice keeps the first retest code
P2[, days := as.numeric(difftime(ts_parse(received), ts_parse(A$sf_received[match(participant, A$participant)]), units = "days"))]
# PREREG: "14-21 days" read as 14 <= days < 22 (whole days 14 to 21 after the main short form).
P2[, in_window := days >= SF2$retest_window[1] & days < SF2$retest_window[2] + 1]
T2 <- RT[P2$retest_code, , drop = FALSE]; rownames(T2) <- P2$participant
if (SF2$mode == "permuted") { set.seed(SF2$seed + 1); for (j in seq_len(ncol(T2))) T2[, j] <- T2[sample.int(nrow(T2)), j] }
T1 <- SFX[P2$participant, , drop = FALSE]
# The follow-up also repeats the three past-week items (after the short form): first
# complete submission, option values from the definition, as at baseline.
PW2 <- t(vapply(P2$retest_code, function(rc) {
  q <- first_complete(rq[[rc]] %||% list(), "sf2-past-week", PW_ITEMS)
  if (is.null(q)) return(rep(NA_real_, length(PW_ITEMS)))
  vapply(PW_ITEMS, function(i) { o <- opt(q$responses[[i]]); if (is.na(o)) NA_real_ else unname(D$past_week[[i]][o]) }, 0) },
  numeric(length(PW_ITEMS))))
PW2 <- matrix(PW2, nrow = nrow(P2), dimnames = list(P2$participant, PW_ITEMS))
if (SF2$mode == "permuted") { set.seed(SF2$seed + 2); for (j in seq_len(ncol(PW2))) PW2[, j] <- PW2[sample.int(nrow(PW2)), j] }
PW1 <- VX$past_week[P2$participant, , drop = FALSE]
retest_flow <- rbindlist(list(
  list("follow-up enrollments", nrow(renr)), list("excluded: staff test", sum(renr$excl_test)),
  list("excluded: not linked to exactly one Study 2 code", sum(renr$excl_unlinked)),
  list("excluded: no submitted short form", sum(renr$excl_no_sf & !renr$excl_test & !renr$excl_unlinked)),
  list("excluded: main-session record excluded or absent", sum(renr$excl_baseline & !renr$excl_no_sf)),
  list("paired", nrow(P2)), list("  of whom outside 14-21 days (kept, flagged)", sum(!P2$in_window))))
setnames(retest_flow, c("step", "n"))
invited_file <- sf2_invited_file()
invited <- if (!is.na(invited_file) && file.exists(invited_file)) fread(invited_file)$participant else NULL

analysis <- list(
  mode = SF2$mode, extraction = basename(raw_file), defs = D,
  meta = A, meta_all = M, sf_raw = SFX, cand_raw = CX, cand_positions = CPOS, val_items = VX, scores = S,
  validity_positions = VPOS, validity_secs = VSECS, country_rule = country_rule,
  retest = list(pairs = P2, t1 = T1, t2 = T2, past_week_t1 = PW1, past_week_t2 = PW2, enrollments = renr, invited = invited),
  flow = flow, retest_flow = retest_flow, submissions = subs_log)
saveRDS(analysis, sf2_derived("sf2"))

fwrite(flow, sf2_out("01_flow.csv")); fwrite(retest_flow, sf2_out("01_retest_flow.csv"))
flags <- M[, .(n = .N, one_check = sum(flag_one_check, na.rm = TRUE), speeder = sum(flag_speeder, na.rm = TRUE),
               longstring = sum(flag_longstring, na.rm = TRUE), pna_over_4 = sum(flag_pna, na.rm = TRUE),
               any_flag = sum(flagged, na.rm = TRUE)), by = .(analysable)]
fwrite(flags, sf2_out("01_flags.csv"))
resub <- subs_log[, .(submissions = .N, flagged_resubmissions = sum(flagged_resubmission),
                      participants_with_repeat = sum(table(participant[!flagged_resubmission]) > 1)), by = slug]
fwrite(resub, sf2_out("01_submissions.csv"))
print(flow); print(retest_flow); print(flags)
sf2_log("country groups: ", paste(names(cn2), cn2, sep = " = ", collapse = ", "), if (country_rule$analysed) "" else " (country not analysed)")
