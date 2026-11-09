# 01: from the raw extraction to one analysis row per participant, applying the
# registered rules (Section 5, "Data exclusion" and "Missing data"; Section 4, "Indices").
#
#   excluded  staff test enrollments; sessions without a submitted scale; failing two
#             or more of the three attention checks (sf2b_attn_agree, sf2b_attn_disagree,
#             bfi2_attn_check). A check not reached is not a failure.
#   first     each questionnaire's first complete submission is used; rows the platform
#             flagged as resubmissions (resubmission_of set) are dropped first.
#   flagged   failing exactly one check; < 2 s per position on the scale (14 positions);
#             a longstring of 10 or more on the 14 positions in presented order ("Prefer
#             not to answer" breaks a run); "Prefer not to answer" on more than 2 of the
#             13 statements. Flags are for sensitivity analyses, never exclusions.
#   partial   a participant who submitted the scale is kept for H1; each later analysis
#             uses those who reached its measures; the pause analyses those who reached
#             the offer.
# In permuted mode every item's answers are shuffled across the participants who
# answered that questionnaire BEFORE any answer-derived index (longstring, scores) is
# computed. Exclusions, checks, timing, background answers and the pause stay real.

source("study2b/analysis/config.R")
source("study2b/analysis/sf2b_models.R")
D <- sf2b_definitions()

raw_file <- sf2b_raw_file()
rows <- fromJSON(raw_file, simplifyVector = FALSE)$rows
sf2b_log("extraction: ", basename(raw_file), ", ", length(rows), " rows")

ts_parse <- function(x) {
  x <- unlist(lapply(x, function(v) if (is.null(v)) NA_character_ else v))
  off <- ifelse(grepl("[+-][0-9]{2}(:?[0-9]{2})?$", x), sub(".*([+-][0-9]{2})(:?([0-9]{2}))?$", "\\1\\3", x), "+00")
  base <- sub("[+-][0-9]{2}(:?[0-9]{2})?$", "", sub("Z$", "", sub("T", " ", x)))
  t <- as.POSIXct(base, tz = "UTC", format = "%Y-%m-%d %H:%M:%OS")
  h <- as.numeric(substr(off, 1, 3)); m <- suppressWarnings(as.numeric(substr(off, 4, 5))); m[is.na(m)] <- 0
  t - sign(h) * (abs(h) * 3600 + m * 60)
}
num <- function(v) if (is.null(v) || is.list(v) || is.character(v)) NA_real_ else as.numeric(v)
opt <- function(v) if (is.list(v) && !is.null(v$option_id)) v$option_id else NA_character_
fld <- function(r, k) { v <- r[[k]]; if (is.null(v)) NA else v }
chr <- function(v) if (is.null(v) || is.list(v)) NA_character_ else as.character(v)

src <- vapply(rows, function(r) fld(r, "source"), "")
enr <- rows[src == "enrollment"]
meta_all <- rbindlist(lapply(enr, function(r) data.table(
  participant = r$participant, is_test = isTRUE(r$is_test), repository_consent = isTRUE(r$repository_consent),
  consented = !is.null(r$consent_date), schedule_id = fld(r, "schedule_id"), session_status = fld(r, "session_status"),
  session_completed = !is.null(r$session_completed_at), step_crashes = as.numeric(fld(r, "step_crashes")))))
meta_all <- meta_all[!duplicated(participant)]

data_rows <- rows[src %in% c("questionnaire", "vas", "step")]
by_person <- split(data_rows, vapply(data_rows, function(r) r$participant, ""))

first_complete <- function(rs, slug, need) {
  cand <- Filter(function(r) identical(r$source, "questionnaire") && identical(r$slug, slug) &&
                   is.null(r$resubmission_of) && all(need %in% names(r$responses)), rs)
  if (!length(cand)) return(NULL)
  cand[[order(ts_parse(lapply(cand, `[[`, "received_at")))[1]]]
}
vas_slugs <- unique(unlist(lapply(data_rows, function(r) if (identical(r$source, "vas")) r$slug)))
if (length(vas_slugs) > 1) stop("more than one VAS scale in Study 2b sessions (", paste(vas_slugs, collapse = ", "),
                                "): the session has one, life satisfaction; check the extraction")
secs_of <- function(q) as.numeric(q$responses[["_timing"]]$total_ms %||% NA) / 1000

subs_log <- list()
recs <- lapply(meta_all$participant, function(p) {
  rs <- by_person[[p]] %||% list()
  qs <- Filter(function(r) identical(r$source, "questionnaire"), rs)
  subs_log[[p]] <<- data.table(participant = p, slug = vapply(qs, `[[`, "", "slug"),
                               flagged_resubmission = vapply(qs, function(r) !is.null(r$resubmission_of), TRUE))
  out <- list(participant = p)
  bg <- first_complete(rs, "sf-background", "bg_practice_regular")
  sc <- first_complete(rs, "sf2b-scale", D$positions)
  vg <- first_complete(rs, "sf2b-vignettes", D$vignettes)
  en <- first_complete(rs, "sf2b-end", c(D$end_items, "sf2b_attn_disagree"))
  pz <- first_complete(rs, "sf2b-pause", "pause_1")
  out$has_background <- !is.null(bg); out$has_scale <- !is.null(sc); out$has_vignettes <- !is.null(vg)
  out$has_end <- !is.null(en); out$has_pause <- !is.null(pz)
  if (!is.null(bg)) {
    b <- bg$responses
    out$age <- if (identical(opt(b$bg_age), "age")) num(b$bg_age$value) else NA_real_
    out$gender <- opt(b$bg_gender); out$country <- opt(b$bg_country); out$education <- opt(b$bg_education)
    out$practice_regular <- opt(b$bg_practice_regular); out$practice_freq <- opt(b$bg_practice_freq)
    out$practice_years <- if (identical(opt(b$bg_practice_years), "years")) num(b$bg_practice_years$value) else NA_real_
    out$practice_training <- opt(b$bg_practice_training)
    out$spiritual <- num(b$bg_spiritual); out$religious <- num(b$bg_religious)
  }
  if (!is.null(sc)) {
    s <- sc$responses
    out$scale <- vapply(D$items, function(i) num(s[[i]]), 0)
    out$scale_pna <- sum(vapply(D$items, function(i) identical(s[[i]], "pna"), TRUE))
    out$check_agree <- num(s$sf2b_attn_agree)
    out$scale_secs <- secs_of(sc)
    out$presentation_format <- sc$presentation_format %||% NA_character_
  }
  if (!is.null(vg)) out$vig <- vapply(D$vignettes, function(i) num(vg$responses[[i]]), 0)
  if (!is.null(en)) { out$end <- vapply(D$end_items, function(i) num(en$responses[[i]]), 0)
    out$check_disagree <- num(en$responses$sf2b_attn_disagree) }
  if (!is.null(pz)) {
    rec <- pz$responses$pause_1; pa <- rec$pause
    out$pause_choice <- chr(rec$choice); out$pause_decision_s <- num(rec$decision_ms) / 1000
    out$pause_decline_reason <- chr(rec$decline_reason)
    out$pause_total_s <- if (is.list(pa)) num(pa$total_ms) / 1000 else NA_real_
    out$pause_on_page_s <- if (is.list(pa)) num(pa$on_page_ms) / 1000 else NA_real_
    out$pause_unfocused_s <- if (is.list(pa)) num(pa$unfocused_ms) / 1000 else NA_real_
    out$pause_hidden_s <- if (is.list(pa)) num(pa$hidden_ms) / 1000 else NA_real_
    out$pause_ended_by <- if (is.list(pa)) chr(pa$ended_by) else NA_character_
    out$pause_nudged <- if (is.list(pa)) isTRUE(pa$nudged) else NA
    out$pause_has_note <- !is.null(rec$note) && nzchar(chr(rec$note))
    out$pause_setting <- chr(rec$setting); out$pause_tired <- num(rec$tired)
    out$pause_device <- chr(rec$device); out$pause_button_side <- chr(rec$pause_button_side)
    out$pause_note <- chr(rec$note)
  }
  out$val <- list(); out$val_step <- c(); out$val_secs <- c()
  for (k in names(D$validity)) {
    v <- D$validity[[k]]
    q <- first_complete(rs, v$slug, c(v$items, v$check))
    if (is.null(q)) next
    out$val[[k]] <- vapply(v$items, function(i) num(q$responses[[i]]), 0)
    if (!is.null(v$check)) out$check_bfi <- num(q$responses[[v$check]])
    out$val_step[v$slug] <- as.numeric(q$step_index %||% NA)
    out$val_secs[v$slug] <- secs_of(q)
  }
  pw <- first_complete(rs, "sf2-past-week", names(D$past_week))
  if (!is.null(pw)) {
    out$val$past_week <- vapply(names(D$past_week), function(i) { o <- opt(pw$responses[[i]]); if (is.na(o)) NA_real_ else unname(D$past_week[[i]][o]) }, 0)
    out$val_step["sf2-past-week"] <- as.numeric(pw$step_index %||% NA)
    out$val_secs["sf2-past-week"] <- secs_of(pw)
  }
  vas <- Filter(function(r) identical(r$source, "vas") && is.null(r$resubmission_of), rs)
  if (length(vas)) {
    v1 <- vas[[order(ts_parse(lapply(vas, `[[`, "received_at")))[1]]]
    out$lifesat <- num(v1$value); out$val_step["lifesat"] <- as.numeric(v1$step_index %||% NA)
  }
  steps <- Filter(function(r) identical(r$source, "step"), rs)
  out$max_step_exited <- if (length(steps)) max(vapply(steps, function(r) as.numeric(r$step_index %||% NA), 0), na.rm = TRUE) else NA_real_
  out
})
names(recs) <- meta_all$participant
subs_log <- rbindlist(subs_log)

scalar <- c("has_background", "has_scale", "has_vignettes", "has_end", "has_pause", "age", "gender", "country", "education",
            "practice_regular", "practice_freq", "practice_years", "practice_training", "spiritual", "religious",
            "scale_pna", "check_agree", "check_disagree", "check_bfi", "scale_secs", "presentation_format",
            grep("^pause_", unique(unlist(lapply(recs, names))), value = TRUE), "lifesat", "max_step_exited")
scalar <- setdiff(scalar, "pause_note")
M <- rbindlist(lapply(recs, function(r) { x <- r[intersect(scalar, names(r))]; x$participant <- r$participant; as.data.table(x) }), fill = TRUE)
M <- merge(meta_all, M, by = "participant", sort = FALSE)
setcolorder(M, "participant")

# ---- exclusions ------------------------------------------------------------------------------------
fails <- function(x, correct) !is.na(x) & x != correct
M[, checks_reached := (!is.na(check_agree)) + (!is.na(check_disagree)) + (!is.na(check_bfi))]
M[, checks_failed := fails(check_agree, SF2B$checks[["sf2b_attn_agree"]]) + fails(check_disagree, SF2B$checks[["sf2b_attn_disagree"]]) +
                     fails(check_bfi, SF2B$checks[["bfi2_attn_check"]])]
M[, `:=`(excl_test = is_test, excl_no_scale = !is_test & !has_scale,
         excl_attention = !is_test & has_scale & checks_failed >= SF2B$exclude_checks_failed)]
M[, analysable := !excl_test & !excl_no_scale & !excl_attention]

flow <- rbindlist(list(
  list("enrollments", nrow(M)),
  list("excluded: staff test enrollment", sum(M$excl_test)),
  list("excluded: no submitted scale", sum(M$excl_no_scale)),
  list("excluded: failed two or more attention checks", sum(M$excl_attention)),
  list("analysable (H1)", sum(M$analysable)),
  list("  reached the first half of the validity block", sum(M$analysable & vapply(recs[M$participant], function(r) any(r$val_step %in% SF2B_STEPS$half_a), TRUE))),
  list("  reached the pause offer (H4b)", sum(M$analysable & M$has_pause)),
  list("  reached the end of the session", sum(M$analysable & (M$session_completed | M$max_step_exited >= 19)))))
setnames(flow, c("step", "n"))

# ---- item matrices -----------------------------------------------------------------------------------
who <- M[has_scale == TRUE & !is_test, participant]
mat <- function(field, cols) {
  X <- matrix(NA_real_, length(who), length(cols), dimnames = list(who, cols))
  for (p in who) { v <- recs[[p]][[field]]; if (!is.null(v)) X[p, ] <- v[cols] }
  X
}
vmat <- function(k, cols) {
  X <- matrix(NA_real_, length(who), length(cols), dimnames = list(who, cols))
  for (p in who) { v <- recs[[p]]$val[[k]]; if (!is.null(v)) X[p, ] <- v[cols] }
  X
}
SCX <- mat("scale", D$items); VIG <- mat("vig", D$vignettes); END <- mat("end", D$end_items)
VX <- lapply(setNames(names(D$validity), names(D$validity)), function(k) vmat(k, D$validity[[k]]$items))
VX$past_week <- vmat("past_week", names(D$past_week))
LS <- vapply(who, function(p) recs[[p]]$lifesat %||% NA_real_, 0)

if (SF2B$mode == "permuted") {
  set.seed(SF2B$seed)
  answered <- function(X) rowSums(!is.na(X)) > 0
  perm <- function(X, has = answered(X)) { for (j in seq_len(ncol(X))) X[has, j] <- X[has, j][sample.int(sum(has))]; X }
  SCX <- perm(SCX); VIG <- perm(VIG); END <- perm(END)
  for (k in names(VX)) VX[[k]] <- perm(VX[[k]])
  hl <- !is.na(LS); LS[hl] <- LS[hl][sample.int(sum(hl))]
  sf2b_log("every item's answers permuted across the participants who answered it (structure destroyed by design)")
}

# Answer-derived flags, on whichever matrix the mode uses. The 14 positions in presented
# (fixed) order: the 13 statements and the check.
POS <- cbind(SCX, sf2b_attn_agree = M$check_agree[match(who, M$participant)])[, D$positions]
ls_run <- apply(POS, 1, function(v) { r <- rle(ifelse(is.na(v), -1, v)); max(c(0, r$lengths[r$values != -1])) })
M[, longstring := NA_real_][match(who, participant), longstring := ls_run]
M[, `:=`(flag_one_check = checks_failed == 1,
         flag_speeder = !is.na(scale_secs) & scale_secs / length(D$positions) < SF2B$speeder_secs_per_item,
         flag_longstring = !is.na(longstring) & longstring >= SF2B$longstring_flag,
         flag_pna = !is.na(scale_pna) & scale_pna > SF2B$pna_flag)]
M[, flagged := flag_one_check | flag_speeder | flag_longstring | flag_pna]

# ---- background indices (Section 4) -----------------------------------------------------------------------
# Current = regular practice "yes" and current frequency 1-2 a week or more; never =
# regular practice "no" (2a's definitions). Dedicated = 3+ a week, a year or more,
# teacher-led training. PREREG: someone who practises regularly but whose frequency is
# unknown is in neither H6 group; a missing regular-practice answer leaves every group missing.
M[, practice_group := fifelse(practice_regular %in% "no", "never",
                       fifelse(practice_regular %in% "yes",
                               fifelse(practice_freq %in% c("1_2_week", "3_5_week", "daily"), "current", "past_or_occasional"), NA_character_))]
M[, dedicated := practice_regular %in% "yes" & practice_freq %in% c("3_5_week", "daily") & !is.na(practice_years) &
                 practice_years >= 1 & practice_training %in% "yes"]
M[, freq_code := fifelse(practice_regular %in% "no", 0, unname(SF2B_FREQ_CODE[practice_freq]))]
# PREREG: gender in H5 is man / woman / other (non-binary or self-described); "prefer
# not to say" leaves the complete-case sample.
M[, gender3 := fifelse(gender %in% c("man", "woman"), gender, fifelse(gender %in% c("non_binary", "self_describe"), "other", NA_character_))]
M[, gender2 := fifelse(gender %in% c("man", "woman"), gender, NA_character_)]
M[, country2 := fifelse(country %in% SF2B_SIX_COUNTRIES, "six_countries", fifelse(country %in% "other", "other", NA_character_))]

# Engaged pause (decision 22): chose the pause, at least 20 s on the page, and away
# (unfocused + hidden) under half the pause's total time; 0 otherwise; missing when the
# offer was not reached. PREREG: a record whose choice is "pause" but which carries no
# pause timing (should not occur) counts as not engaged and is counted in 01_pause_check.csv.
M[, pause_taken := fifelse(has_pause, as.integer(pause_choice %in% "pause"), NA_integer_)]
M[, engaged_pause := fifelse(!has_pause, NA_integer_,
       as.integer(pause_choice %in% "pause" & !is.na(pause_on_page_s) & pause_on_page_s >= SF2B$pause_min_secs &
                  !is.na(pause_total_s) & (pause_unfocused_s + pause_hidden_s) < SF2B$pause_max_away * pause_total_s))]
M[, pause_device := factor(pause_device, levels = c("desktop", "phone", "tablet"))]

A <- M[analysable == TRUE]
cn2 <- table(A$country2)
country_rule <- list(groups = as.list(cn2), analysed = length(cn2) == 2 && all(cn2 >= SF2B$country_min))
A[, country_group := if (country_rule$analysed) country2 else NA_character_]

ids <- A$participant
SCX <- SCX[ids, , drop = FALSE]; VIG <- VIG[ids, , drop = FALSE]; END <- END[ids, , drop = FALSE]
VX <- lapply(VX, function(X) X[ids, , drop = FALSE]); LS <- LS[ids]

# ---- observed scores -------------------------------------------------------------------------------------------
S <- data.table(participant = ids,
  sf8 = sf2b_score(SCX, D$sf8, D$reversed, need = SF2B$min_answered[["sf8"]]),
  sf10 = sf2b_score(SCX, D$sf10, D$reversed, need = SF2B$min_answered[["sf10"]]))
for (k in names(D$validity)) { v <- D$validity[[k]]
  for (s in names(v$scores)) S[[s]] <- sf2b_score(VX[[k]], v$scores[[s]], v$reverse, v$range, v$method) }
# BIDR-16 dichotomous (Paulhus; exploratory, H3c): after reversal, 6 or 7 counts 1; count per subscale (0-8), 80% rule.
for (s in c("bidr_sde", "bidr_im")) {
  it <- D$validity$bidr16$scores[[s]]; Xb <- VX$bidr16[, it]; rv <- intersect(it, D$validity$bidr16$reverse)
  Xb[, rv] <- 8 - Xb[, rv]; Xb <- (Xb >= 6) * 1
  S[[paste0(s, "_dichotomous")]] <- ifelse(rowSums(!is.na(Xb)) >= ceiling(.8 * length(it)), rowMeans(Xb, na.rm = TRUE) * length(it), NA)
}
for (i in names(D$past_week)) S[[i]] <- VX$past_week[, i]
S$lifesat <- LS
S$engaged_pause <- A$engaged_pause
# The trial's Access + Valuing score in its own wording: its 10 statements, mean of those answered.
AV <- cbind(END, sf2b_02 = SCX[, "sf2b_02"], sf2b_11 = SCX[, "sf2b_11"])[, unname(D$av10)]
colnames(AV) <- names(D$av10)
S$trial_av10 <- sf2b_score(AV, names(D$av10), D$av10_reversed, need = 1)

# Vignettes (Section 4): separation = mean sensing - mean doing, at least three of each.
sens <- names(D$vkind)[D$vkind == "sensing"]; doing <- names(D$vkind)[D$vkind == "doing"]
sep <- rowMeans(VIG[, sens, drop = FALSE], na.rm = TRUE) - rowMeans(VIG[, doing, drop = FALSE], na.rm = TRUE)
sep[rowSums(!is.na(VIG[, sens, drop = FALSE])) < SF2B$min_vignettes_each | rowSums(!is.na(VIG[, doing, drop = FALSE])) < SF2B$min_vignettes_each] <- NA
A[, `:=`(separation = sep, tells_apart = fifelse(is.na(sep), NA_character_, fifelse(sep >= SF2B$separation_cut, "yes", "no")),
         vig_self = VIG[, "vig_self"])]

# Presented position of each validity measure within its half (1-9, 1-6), and seconds per questionnaire.
pos_in_half <- function(st) ifelse(st %in% SF2B_STEPS$half_a, st - min(SF2B_STEPS$half_a) + 1,
                                   ifelse(st %in% SF2B_STEPS$half_b, st - min(SF2B_STEPS$half_b) + 1, NA))
VPOS <- rbindlist(lapply(ids, function(p) { s <- recs[[p]]$val_step
  as.data.table(as.list(c(participant = p, setNames(as.character(pos_in_half(s)), names(s))))) }), fill = TRUE)
for (j in setdiff(names(VPOS), "participant")) set(VPOS, j = j, value = as.numeric(VPOS[[j]]))
VSECS <- rbindlist(lapply(ids, function(p) as.data.table(as.list(c(participant = p, recs[[p]]$val_secs)))), fill = TRUE)
for (j in setdiff(names(VSECS), "participant")) set(VSECS, j = j, value = as.numeric(VSECS[[j]]))
# Pause notes: free text, kept only in derived/ (never in output/), for blind coding.
NOTES <- data.table(participant = ids, note = vapply(ids, function(p) recs[[p]]$pause_note %||% NA_character_, ""))[!is.na(note)]

analysis <- list(mode = SF2B$mode, extraction = basename(raw_file), defs = D, meta = A, meta_all = M,
                 scale_raw = SCX, end_raw = END, vignettes = VIG, val_items = VX, scores = S,
                 validity_positions = VPOS, validity_secs = VSECS, country_rule = country_rule, flow = flow,
                 submissions = subs_log)
saveRDS(analysis, sf2b_derived("sf2b"))
saveRDS(NOTES, sf2b_derived("pause_notes"))

fwrite(flow, sf2b_out("01_flow.csv"))
flags <- M[, .(n = .N, one_check = sum(flag_one_check, na.rm = TRUE), speeder = sum(flag_speeder, na.rm = TRUE),
               longstring = sum(flag_longstring, na.rm = TRUE), pna_over_2 = sum(flag_pna, na.rm = TRUE),
               any_flag = sum(flagged, na.rm = TRUE), step_crash = sum(step_crashes > 0, na.rm = TRUE)), by = .(analysable)]
fwrite(flags, sf2b_out("01_flags.csv"))
resub <- subs_log[, .(submissions = .N, flagged_resubmissions = sum(flagged_resubmission),
                      participants_with_repeat = sum(table(participant[!flagged_resubmission]) > 1)), by = slug]
fwrite(resub, sf2b_out("01_submissions.csv"))
fwrite(A[has_pause == TRUE, .(reached_offer = .N, took_pause = sum(pause_taken), engaged = sum(engaged_pause),
                              pause_without_timing = sum(pause_choice %in% "pause" & is.na(pause_total_s)),
                              ended_by_cap = sum(pause_ended_by %in% "cap"), nudged = sum(pause_nudged %in% TRUE),
                              wrote_a_note = sum(pause_has_note))], sf2b_out("01_pause_check.csv"))
print(flow); print(flags)
sf2b_log("country groups: ", paste(names(cn2), cn2, sep = " = ", collapse = ", "), if (country_rule$analysed) "" else " (country not analysed)")
