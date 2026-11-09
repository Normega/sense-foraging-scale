# 00 (simulated mode only): a raw extraction from a KNOWN model, in exactly the shape
# extract.sql returns, so everything from 01 on runs the code path the real data will.
#
# The truth:
#   the scale     Study 2a's fitted one-factor model of the 13 statements with the
#                 agreement factor (sim_population.json): standardized content and
#                 agreement loadings and thresholds. The model of H1 is true.
#   validity      each score is a latent factor F = rho * S + sqrt(1 - rho^2) * (g * W + sqrt(1 - g^2) * U),
#                 W shared within a battery (FFMQ facets, Brief MAIA-2 subscales, BFI-2
#                 facets, BIDR, wellbeing), so rho is F's correlation with the scale's
#                 factor S and the battery correlations are known. Items load .60-.80 on
#                 their factor (negatively when reverse-keyed) and are cut into their
#                 response categories. No agreement factor in the validity items.
#   criteria      the past-week items and the pause from latents with known rho.
#   practice      a latent practice propensity correlating .35 with S; current
#                 practitioners are its top 55%, never-practised its bottom 27% (2a's shares).
#   vignettes     each person's separation ~ N(1.5, 1), as in 2a's simulation.
#   invariance    for women, V10 (sf2b_03) loads .20 instead of its population loading,
#                 for the exploratory invariance analyses to find.
# Planted mess: staff tests; enrollments that never reach the scale; two and one failed
# attention checks; sessions stopping after the scale, before the pause and in the
# second half; speeders, longstrings, many "Prefer not to answer"; a flagged
# resubmission and a later re-entry with different answers.
# Writes raw/simulated_extract.json and derived/simulation_truth_simulated.rds.

source("study2b/analysis/config.R")
source("study2b/analysis/sf2b_models.R")
stopifnot(SF2B$mode == "simulated")
D <- sf2b_definitions()
set.seed(SF2B$seed_sim)
pop <- fromJSON("study2b/analysis/sim_population.json")
stopifnot(identical(pop$items, D$items))

N_MAIN <- 700           # complete, analysable sessions before the planted mess below
plant <- list(test = 8, no_scale = 14, two_checks = 12, one_check = 15, stop_after_end = 18, stop_before_pause = 12,
              stop_in_b = 10, speeder = 6, longstring = 6, many_pna = 6, resubmission = 3, reentry = 2)
n <- N_MAIN + plant$test + plant$no_scale + plant$two_checks + plant$one_check
codes <- sprintf("P-%06d", sample(1e5:9e5, n))
role <- c(rep("main", N_MAIN), rep("test", plant$test), rep("no_scale", plant$no_scale),
          rep("two_checks", plant$two_checks), rep("one_check", plant$one_check))
role <- sample(role)
codes[role == "test"] <- sprintf("test-%02d", seq_len(plant$test))

# ---- latent scores ---------------------------------------------------------------------------
S <- rnorm(n); AGR <- rnorm(n)
BATTERY <- list(ffmq = .50, maia = .60, bfi = .60, bidr = .40)
RHO <- c(ffmq_observing = .50, ffmq_describing = .25, ffmq_acting = .15, ffmq_nonjudging = .10, ffmq_nonreactivity = .30,
         maia_noticing = .45, maia_notdistracting = .10, maia_notworrying = .05, maia_attnregulation = .35,
         maia_emoawareness = .30, maia_selfregulation = .30, maia_bodylistening = .35, maia_trusting = .20,
         bfi_aesthetic = .45, bfi_intellectual = .25, bfi_creative = .20,
         bidr_sde = .15, bidr_im = .05, sbi_moment = .50, eq_decentering = .30, samd_md = .25, bsss8 = .20,
         ncs6 = .10, spsq_sd = .05, asi_pc = .00, phq4_total = -.15, pss4 = -.15, lifesat = .20,
         crit_pause_days = .30, crit_notice_first = .20, crit_outdoors = .15, pause = .25)
battery_of <- function(sc) {
  b <- sub("_.*", "", sc)
  if (b %in% names(BATTERY)) b else if (sc %in% c("phq4_total", "pss4", "lifesat")) "wellbeing" else NA_character_
}
W <- list(ffmq = rnorm(n), maia = rnorm(n), bfi = rnorm(n), bidr = rnorm(n), wellbeing = rnorm(n))
G_WB <- c(phq4_total = .70, pss4 = .70, lifesat = -.60)
gamma_of <- function(sc) { b <- battery_of(sc); if (is.na(b)) 0 else if (b == "wellbeing") G_WB[[sc]] else BATTERY[[b]] }
F <- sapply(names(RHO), function(sc) {
  r <- RHO[[sc]]; g <- gamma_of(sc); b <- battery_of(sc)
  r * S + sqrt(1 - r^2) * ((if (is.na(b)) 0 else g * W[[b]]) + sqrt(1 - g^2) * rnorm(n))
})
# The population correlation matrix of S and every factor (the estimands).
fac <- names(RHO)
Rtrue <- diag(length(fac) + 1); dimnames(Rtrue) <- list(c("S", fac), c("S", fac))
for (a in fac) { Rtrue["S", a] <- Rtrue[a, "S"] <- RHO[[a]]
  for (b in fac) if (a != b) { ba <- battery_of(a); bb <- battery_of(b)
    Rtrue[a, b] <- RHO[[a]] * RHO[[b]] + if (!is.na(ba) && identical(ba, bb)) sqrt((1 - RHO[[a]]^2) * (1 - RHO[[b]]^2)) * gamma_of(a) * gamma_of(b) else 0 } }

# ---- items -------------------------------------------------------------------------------------------
cut_cat <- function(z, th) findInterval(z, th) + 1L
# The scale: 2a's population, except one planted non-invariance for the exploratory
# invariance analyses to find: for women, V10 (sf2b_03, noticing the body) loads .20
# instead of its population loading (residual variance raised to keep its variance 1).
gender <- sample(c("man", "woman", "non_binary", "self_describe", "prefer_not"), n, TRUE, c(.47, .48, .03, .01, .01))
NONINV <- list(item = "sf2b_03", group = "woman", loading = .20)
thr <- lapply(pop$items, function(i) unlist(pop$thresholds[[i]]))
SC <- sapply(seq_along(pop$items), function(j) {
  l <- rep(pop$loading_S[j], n)
  if (pop$items[j] == NONINV$item) l[gender == NONINV$group] <- NONINV$loading
  e <- sqrt(1 - l^2 - pop$loading_RI[j]^2)
  cut_cat(l * S + pop$loading_RI[j] * AGR + e * rnorm(n), thr[[j]]) })
colnames(SC) <- pop$items
# The trial's statements in their own wording: the same S and agreement, loadings .65 (forward) / -.55 (reversed).
th6 <- qnorm(c(.04, .12, .30, .62, .88))
EN <- sapply(D$end_items, function(i) { l <- if (i %in% D$av10_reversed) -.55 else .65
  cut_cat(l * S + .25 * AGR + sqrt(1 - l^2 - .0625) * rnorm(n), if (l < 0) -rev(th6) else th6) })
# Validity items: each on its score's factor; reverse-keyed items load negatively.
score_items <- sf2b_score_items(D)
VI <- list()
for (k in names(D$validity)) {
  v <- D$validity[[k]]; ncat <- diff(v$range) + 1
  base <- qnorm(seq_len(ncat - 1) / ncat)
  M <- sapply(v$items, function(i) {
    sc <- names(score_items)[vapply(score_items, function(s) i %in% s, TRUE)]
    sc <- setdiff(sc, "bfi_openmindedness")[1]
    lam <- runif(1, .60, .80) * if (i %in% v$reverse) -1 else 1
    th <- base + rnorm(1, 0, .3)
    cut_cat(lam * F[, sc] + sqrt(1 - lam^2) * rnorm(n), th) + v$range[1] - 1L })
  VI[[k]] <- M
}
LIFESAT <- cut_cat(F[, "lifesat"] + rnorm(n, 0, .5), qnorm(c(.05, .15, .35, .65, .88), 0, sqrt(1.25)))
PW <- cbind(crit_pause_days = cut_cat(F[, "crit_pause_days"], qnorm(c(.20, .35, .50, .62, .72, .80, .88))) - 1L,
            crit_notice_first = cut_cat(F[, "crit_notice_first"], qnorm(c(.08, .25, .60, .88))),
            crit_outdoors = cut_cat(F[, "crit_outdoors"], qnorm(c(.10, .35, .70, .92))) - 1L)
nothing_unexp <- runif(n) < .10

# ---- the pause -----------------------------------------------------------------------------------------
take <- F[, "pause"] > qnorm(.70)
on_page <- ifelse(take, pmin(600, exp(rnorm(n, log(60), .9))), NA)
away_frac <- ifelse(take, ifelse(runif(n) < .80, rbeta(n, 1, 12), rbeta(n, 5, 3)), NA)
away <- on_page * away_frac / (1 - away_frac)
total <- pmin(600, on_page + away)
capped <- take & (on_page + away) >= 600
scale_away <- ifelse(take, total / (on_page + away), NA)        # trim proportionally at the 10-minute cap
on_page <- on_page * scale_away; away <- away * scale_away
unf <- away * runif(n); hid <- away - unf
engaged_true <- take & on_page >= SF2B$pause_min_secs & (unf + hid) < SF2B$pause_max_away * total
device <- sample(c("phone", "tablet", "desktop"), n, TRUE, c(.45, .05, .50))

# ---- background -------------------------------------------------------------------------------------------
Pz <- .35 * S + sqrt(1 - .35^2) * rnorm(n)
q <- quantile(Pz, c(.27, .45))
prac <- ifelse(Pz <= q[1], "never", ifelse(Pz <= q[2], "occasional", "current"))
freq <- ifelse(prac == "never", NA, ifelse(prac == "occasional", sample(c("none", "lt_weekly"), n, TRUE),
               sample(c("1_2_week", "3_5_week", "daily"), n, TRUE, c(.45, .30, .25))))
years <- ifelse(prac == "never", NA, round(rexp(n, 1 / 4), 1))
training <- ifelse(prac == "never", NA, sample(c("yes", "no"), n, TRUE, c(.4, .6)))
spiritual <- cut_cat(.3 * Pz + rnorm(n), qnorm(c(.25, .45, .70, .88))); religious <- cut_cat(.15 * Pz + rnorm(n), qnorm(c(.45, .65, .80, .92)))
age <- round(pmin(80, 18 + rexp(n, 1 / 18)))
country <- sample(c("uk", "us", "ca", "au", "ie", "nz", "other"), n, TRUE, c(.40, .30, .06, .04, .02, .01, .17))
sep <- rnorm(n, 1.5, 1); centre <- rnorm(n, 3.5, .5)

# ---- who does what ---------------------------------------------------------------------------------------
main_ix <- which(role == "main")
pick <- function(k, from) { s <- sample(from, k); s }
avail <- main_ix
take_from <- function(k) { s <- sample(avail, k); avail <<- setdiff(avail, s); s }
stop_end <- take_from(plant$stop_after_end); stop_pause <- take_from(plant$stop_before_pause); stop_b <- take_from(plant$stop_in_b)
speeders <- take_from(plant$speeder); longs <- take_from(plant$longstring); pnas <- take_from(plant$many_pna)
resub <- take_from(plant$resubmission); reent <- take_from(plant$reentry)
no_scale_bg <- which(role == "no_scale")[seq_len(plant$no_scale) %% 2 == 0]   # half of them answered the background

# ---- rows ------------------------------------------------------------------------------------------------
t0 <- as.POSIXct("2026-11-20 09:00:00", tz = "UTC") + sort(runif(n, 0, 6 * 86400))
iso <- function(t) format(t, "%Y-%m-%d %H:%M:%OS3+00", tz = "UTC")
uuid <- function() paste(sprintf("%08x", sample.int(2^31 - 1, 1)), sprintf("%04x", sample.int(65535, 1)), "4000",
                         sprintf("%04x", sample.int(65535, 1)), sprintf("%012.0f", runif(1, 1e10, 9e11)), sep = "-")
rows <- vector("list", 80000L); n_rows <- 0L
add <- function(r) { n_rows <<- n_rows + 1L; rows[[n_rows]] <<- r }
base_row <- function(i, sched, status, completed) list(participant = codes[i], is_test = role[i] == "test",
  prolific_study_id = "sim-posting", repository_consent = runif(1) < .8, enrollment_status = "active",
  consent_date = iso(t0[i]), enrolled_at = iso(t0[i]), schedule_id = sched, session_status = status,
  session_completed_at = if (completed) iso(t0[i] + 1800) else NA, step_crashes = 0L)
likert <- function(x, pna_p = .01) lapply(x, function(v) if (runif(1) < pna_p) "pna" else as.integer(v))
mc <- function(id, value = NULL) if (is.null(value)) list(option_id = id) else list(option_id = id, value = value)

for (i in seq_len(n)) {
  sched <- uuid(); r0 <- role[i]
  last_step <- if (r0 == "no_scale") (if (i %in% no_scale_bg) 1L else -1L) else
    if (i %in% stop_end) 3L else if (i %in% stop_pause) 12L else if (i %in% stop_b) 16L else 20L
  completed <- last_step == 20L
  add(c(list(source = "enrollment"), base_row(i, sched, if (completed) "completed" else "pending", completed),
        list(response_id = NA, slug = NA, step_index = NA, received_at = NA, resubmission_of = NA,
             presentation_format = NA, responses = NA, value = NA, entered_at = NA, exited_at = NA)))
  if (last_step < 0) next
  tt <- t0[i] + 60
  q_row <- function(slug, step, resp, at = tt, resub_of = NA) c(list(source = "questionnaire"), base_row(i, sched, if (completed) "completed" else "pending", completed),
    list(response_id = uuid(), slug = slug, step_index = step, received_at = iso(at), resubmission_of = resub_of,
         presentation_format = "stacked", responses = resp, value = NA, entered_at = NA, exited_at = NA))
  step_row <- function(slug, step, secs) c(list(source = "step"), base_row(i, sched, if (completed) "completed" else "pending", completed),
    list(response_id = NA, slug = slug, step_index = step, received_at = NA, resubmission_of = NA,
         presentation_format = NA, responses = NA, value = NA, entered_at = iso(tt - secs), exited_at = iso(tt)))
  timing <- function(secs) list(total_ms = round(secs * 1000), v = 1L)
  # 0 background
  b <- list(bg_age = mc("age", age[i]), bg_gender = mc(gender[i]), bg_country = mc(country[i]), bg_education = mc("bachelors"),
            bg_race = mc("white"), bg_practice_regular = mc(if (prac[i] == "never") "no" else "yes"),
            bg_spiritual = as.integer(spiritual[i]), bg_religious = as.integer(religious[i]), bg_tradition = mc("none"),
            `_timing` = timing(90))
  if (prac[i] != "never") { b$bg_practice_freq <- mc(freq[i]); b$bg_practice_years <- mc("years", years[i])
    b$bg_practice_training <- mc(training[i]); b$bg_practice_types <- list(mc("meditation")) }
  tt <- tt + 90; add(q_row("sf-background", 0L, b)); add(step_row("questionnaire", 0L, 90))
  if (last_step < 2) next
  # 1 vignettes
  vv <- vapply(names(D$vkind), function(v) pmin(6, pmax(1, round(centre[i] + (if (D$vkind[[v]] == "sensing") .5 else -.5) * sep[i] + rnorm(1, 0, .8)))), 0)
  vr <- likert(c(vv, vig_self = pmin(6, pmax(1, round(3.5 + .6 * S[i] + rnorm(1))))))
  vr$`_order` <- list(method = "interleave", version = 1L, positions = as.list(setNames(sample(9), names(D$vkind))))
  vr$`_timing` <- timing(150)
  tt <- tt + 150; add(q_row("sf2b-vignettes", 1L, vr)); add(step_row("questionnaire", 1L, 150))
  # 2 the scale (14 positions, fixed order)
  x <- SC[i, ]
  if (i %in% longs) x[] <- 5L
  chk_agree <- if (r0 == "two_checks" || (r0 == "one_check" && i %% 2 == 0)) 3L else 5L
  sr <- likert(x, if (i %in% pnas) .35 else .01); sr$sf2b_attn_agree <- chk_agree
  if (i %in% pnas) { k <- sample(names(x), 4); for (j in k) sr[[j]] <- "pna" }
  secs <- if (i %in% speeders) 14 * 1.4 else 14 * runif(1, 3, 7)
  sr$`_timing` <- timing(secs)
  tt <- tt + secs; add(q_row("sf2b-scale", 2L, sr))
  if (i %in% resub) { first_id <- rows[[n_rows]]$response_id
    add(q_row("sf2b-scale", 2L, sr, at = tt + 0.6, resub_of = first_id)) }
  if (i %in% reent) { sr2 <- sr; for (j in names(x)) sr2[[j]] <- as.integer(sample(1:6, 1)); add(q_row("sf2b-scale", 2L, sr2, at = tt + 900)) }
  add(step_row("questionnaire", 2L, secs))
  # 3 the trial statements, with the second check
  er <- likert(EN[i, ]); er$sf2b_attn_disagree <- if (r0 == "two_checks" || (r0 == "one_check" && i %% 2 == 1)) 4L else 2L
  er$`_timing` <- timing(40)
  tt <- tt + 40; add(q_row("sf2b-end", 3L, er)); add(step_row("questionnaire", 3L, 40))
  if (last_step < 4) next
  # 4-12 first half, 13 pause, 14-19 second half (per-participant orders)
  half_a <- sample(c("sf2-ffmq15", "sf2-bidr16", "sf2-sbi-moment", "sf2b-bsss8", "sf2b-ncs6", "sf2-phq4", "sf2b-spsq-sd", "sf2-pss4", "lifesat"))
  half_b <- sample(c("sf2b-maia2-brief", "sf2-bfi2-om", "sf2-eq-decentering", "sf2-samd-md", "sf2b-asi3-pc", "sf2-past-week"))
  slug_key <- setNames(names(D$validity), vapply(D$validity, `[[`, "", "slug"))
  give <- function(slug, step) {
    secs <- runif(1, 30, 120)
    tt <<- tt + secs
    if (slug == "lifesat") {
      add(c(list(source = "vas"), base_row(i, sched, if (completed) "completed" else "pending", completed),
            list(response_id = uuid(), slug = "vas_life_satisfaction", step_index = step, received_at = iso(tt), resubmission_of = NA,
                 presentation_format = NA, responses = NA, value = as.integer(LIFESAT[i]), entered_at = NA, exited_at = NA)))
      add(step_row("vas_life_satisfaction", step, secs)); return(invisible())
    }
    if (slug == "sf2-past-week") {
      ids <- names(D$past_week)
      resp <- lapply(ids, function(id) { val <- PW[i, id]
        if (id == "crit_notice_first" && nothing_unexp[i]) return(mc("nothing_unexpected"))
        mc(names(D$past_week[[id]])[match(val, D$past_week[[id]])]) })
      names(resp) <- ids
    } else {
      k <- slug_key[[slug]]; v <- D$validity[[k]]
      resp <- likert(VI[[k]][i, ])
      if (!is.null(v$check)) resp[[v$check]] <- if (r0 == "two_checks") 2L else 4L
    }
    resp$`_timing` <- timing(secs)
    add(q_row(slug, step, resp)); add(step_row("questionnaire", step, secs))
  }
  for (s in seq_along(half_a)) give(half_a[s], 3L + s)
  if (last_step < 13) next
  tt <- tt + 20
  rec <- list(v = 1L, slot = 1L, pause_button_side = sample(c("left", "right"), 1), device = device[i],
              viewport_w = if (device[i] == "phone") 390L else 1440L, local_hour = sample(8:22, 1),
              choice = if (take[i]) "pause" else "continue", decision_ms = round(runif(1, 2000, 12000)),
              decline_reason = if (take[i]) NA else sample(c("no_time", "not_good_place", "not_interested", "other", "skipped"), 1),
              pause = if (take[i]) list(total_ms = round(total[i] * 1000), on_page_ms = round(on_page[i] * 1000),
                                        unfocused_ms = round(unf[i] * 1000), hidden_ms = round(hid[i] * 1000),
                                        hidden_episodes = as.integer(hid[i] > 0), unfocused_episodes = as.integer(unf[i] > 0),
                                        taps = 1L, keys = 0L, scrolls = 0L, ended_by = if (capped[i]) "cap" else "button",
                                        nudged = on_page[i] >= 180) else NA,
              note = NA, setting = sample(c("home", "work", "outdoors", "transit", "other", "pna"), 1, prob = c(.7, .12, .04, .04, .06, .04)),
              tired = if (runif(1) < .03) "pna" else sample(1:6, 1), done = TRUE)
  tt <- tt + if (take[i]) total[i] else 5
  add(q_row("sf2b-pause", 13L, list(pause_1 = rec, `_timing` = timing(if (take[i]) total[i] + 20 else 15))))
  add(step_row("questionnaire", 13L, 20))
  for (s in seq_along(half_b)) { if (13L + s > last_step) break; give(half_b[s], 13L + s) }
  if (completed) add(step_row("activity", 20L, 30))
}

rows <- rows[seq_len(n_rows)]
ext <- list(rows = rows)
dir.create(sf2b_path("raw"), recursive = TRUE, showWarnings = FALSE)
write_json(ext, sf2b_path("raw", "simulated_extract.json"), auto_unbox = TRUE, null = "null", na = "null", digits = NA)

# What 01 should find, for 99_check_recovery.R.
planted <- list(
  enrolled = n, tests = plant$test, no_scale = plant$no_scale,
  two_checks = plant$two_checks, one_check = plant$one_check,
  analysable = n - plant$test - plant$no_scale - plant$two_checks,
  speeders = codes[speeders], longstrings = codes[longs], many_pna = codes[pnas],
  resubmission = codes[resub], reentry = codes[reent], reached_pause = sum(role %in% c("main", "one_check")) - length(c(stop_end, stop_pause)))
an <- role %in% c("main", "one_check")
truth <- list(population = pop, rho = RHO, noninvariance = NONINV,
              sample_r = setNames(vapply(fac, function(a) cor(S[an], F[an, a]), 0), fac), R = Rtrue, battery = BATTERY, gamma_wellbeing = G_WB,
              engaged_rate = mean(engaged_true[role %in% c("main", "one_check")]), take_rate = mean(take),
              practice_d_latent = NA, planted = planted,
              agreement_variance_std = mean(pop$loading_RI[D$sf8]^2))
saveRDS(truth, sf2b_derived("simulation_truth"))
sf2b_log("simulated extraction: ", n, " enrollments, ", length(rows), " rows -> ", sf2b_path("raw", "simulated_extract.json"))
