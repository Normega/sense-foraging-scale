# 24: measurement invariance of the 8 by race/ethnicity (and by practice group), and the
# latent mean differences it licenses. EXPLORATORY, NOT REGISTERED.
#
# Not part of the filed registration (osf.io/d4a5u) and not in its list of exploratory
# analyses. Specified on 17 November 2026, after the registered analyses and 10-23 had
# been seen, because 22 found the 8 higher among participants who are not White (latent
# d .70 alone, .33 adjusted) and practice confounded with race: a group difference is
# only interpretable if the 8 measures the same thing in both groups.
#
# Procedure: H8's, unchanged (sf2b_run_mg in sf2b_models.R): the H1 model (S + the
# agreement factor RI, loadings 1, variance free per group, orthogonal to S), WLSMV,
# theta; configural, then thresholds, then loadings (Wu & Estabrook 2016); a step is
# non-invariant when delta scaled CFI < -.010 AND delta scaled RMSEA > .015; partial
# invariance by releasing up to SF2B$partial_max_free items, best CFI first; categories
# with fewer than 5 answers in any group merged toward the centre.
#
# Added here (H8 stops at loadings): a fourth step, intercepts. Every latent-response
# intercept is fixed at 0 in every group (except items released at this step), and the
# means of S and of RI are freed in groups 2..G; the same delta rule judges it against
# the loadings step. Two versions of the means model:
#   means_ri_free    S and RI means both free (primary: 22 found agreement differs)
#   means_ri_fixed   RI mean fixed at 0 in every group (all shift goes to S)
# Latent d of group k vs the reference = mean_k / sqrt((n_1 * var_1 + n_k * var_k) /
# (n_1 + n_k)), S's reference variance 1; a defined parameter, delta-method 95% CI.
#
# Comparisons (reference group first; labels sort alphabetically, so the reference is 1_):
#   race2      White only vs not White (every other answered category, incl. more than
#              one selected and self-described; prefer not / missing excluded)
#   race_bw    White only vs Black only (the only other category with n >= 200,
#              SF2B$country_min, the registered minimum for a group comparison)
#   race3      White only / Black only / every other answered category pooled
#   practice   never practised vs current practice (H6's groups). Never has n = 164,
#              below the 200 minimum: reported, flagged as under-powered.
#
# Race/ethnicity is read from the extraction exactly as 22 and 23 do (option ids only;
# free text is never read, printed or written). Writes aggregates only:
#   output/real/24_invariance.csv             every step of every comparison
#   output/real/24_invariance_merges.csv      categories merged per comparison
#   output/real/24_latent_means.csv           mean, variance and d of S and RI per group
#   output/real/24_exploratory_invariance.json
# Run from the repo root:  SF2B_MODE=real Rscript study2b/analysis/24_exploratory_invariance_race.R

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
A <- readRDS(sf2b_derived("sf2b")); D <- A$defs; M <- A$meta; X <- A$scale_raw
I8 <- D$sf8
stopifnot(nrow(X) == nrow(M))

## Race/ethnicity from the extraction (option ids only; as 22 and 23) ---------
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
    data.table(participant = r$participant, received = r$received_at %||% NA_character_, race_ids = paste(sort(ids), collapse = " "),
               race_status = if (multi && length(ids)) "answered" else if (is.character(rc) && identical(rc, "pna")) "prefer_not" else "missing")
  }))
  rec[, received_t := ts_parse(as.list(received))]
  rec <- rec[order(participant, received_t)][!duplicated(participant)]
  rec[match(M$participant, participant)]
})
stopifnot(identical(RC$participant, M$participant))
race_me_of <- function(s, status) {
  if (status != "answered") return(status)
  sub <- setdiff(s, "prefer_not")
  if (!length(sub)) return("prefer_not")
  if (length(sub) >= 2 || "mixed" %in% sub) return("more_than_one")
  if (identical(sub, "other")) return("self_described")
  sub
}
race_me <- mapply(race_me_of, strsplit(RC$race_ids, " ", fixed = TRUE), RC$race_status)
ok_race <- !race_me %in% c("prefer_not", "missing")

## Groupings ---------
GROUPS <- list(
  race2    = list(label = "race/ethnicity: White only vs not White",
                  g = fifelse(!ok_race, NA_character_, fifelse(race_me == "white", "1_white", "2_not_white"))),
  race_bw  = list(label = "race/ethnicity: White only vs Black only",
                  g = fifelse(race_me == "white", "1_white", fifelse(race_me == "black", "2_black", NA_character_))),
  race3    = list(label = "race/ethnicity: White only / Black only / other not White",
                  g = fifelse(!ok_race, NA_character_, fifelse(race_me == "white", "1_white", fifelse(race_me == "black", "2_black", "3_other_not_white")))),
  practice = list(label = "practice: never vs current (H6's groups)",
                  g = fifelse(M$practice_group == "never", "1_never", fifelse(M$practice_group == "current", "2_current", NA_character_)))
)
for (k in names(GROUPS)) { n <- table(GROUPS[[k]]$g); GROUPS[[k]]$n <- as.list(n)
  GROUPS[[k]]$below_min <- names(n)[n < SF2B$country_min]
  sf2b_log(k, ": ", paste(names(n), n, sep = " = ", collapse = ", ")) }

## The intercepts step and the latent means (added to H8's procedure) ---------
# Takes the loadings-step syntax (with the partial releases sf2b_run_mg found) and edits
# it: intercepts fixed at 0 in every group, S and RI means and variances labelled.
means_syntax <- function(items, ncat, G, free, free_int = character(), ri_mean_free = TRUE) {
  syn <- strsplit(sf2b_inv_syntax(items, ncat, G, "loadings", free), "\n", fixed = TRUE)[[1]]
  g0 <- function(x) paste0("c(", paste(x, collapse = ", "), ")")
  for (i in setdiff(items, free_int)) {
    w <- which(startsWith(syn, paste0(i, " ~ ")))
    stopifnot(length(w) == 1)
    syn[w] <- paste0(i, " ~ ", g0(rep("0", G)), "*1")
  }
  rep_line <- function(pattern, new) { w <- grep(pattern, syn); stopifnot(length(w) == 1); syn[w] <<- new }
  gk <- seq(2, G)
  # lavaan: a c() modifier that mixes numbers and labels makes every element a label
  # (so c(1, vS_g2) would leave group 1 free). Fixed values and labels therefore go in
  # two modifiers on the same parameter: c(1, NA)*S + c(vS_g1, vS_g2)*S.
  lab <- function(p) g0(paste0(p, "_g", seq_len(G)))
  rep_line("^S ~ c[(].*[*]1$", paste0("S ~ ", g0(c("0", rep("NA", G - 1))), "*1 + ", lab("mS"), "*1"))
  rep_line("^S ~~ c[(].*[*]S$", paste0("S ~~ ", g0(c("1", rep("NA", G - 1))), "*S + ", lab("vS"), "*S"))
  rep_line("^RI ~ c[(].*[*]1$", if (ri_mean_free) paste0("RI ~ ", g0(c("0", rep("NA", G - 1))), "*1 + ", lab("mR"), "*1")
                                 else paste0("RI ~ ", g0(rep("0", G)), "*1"))
  rep_line("^RI ~~ c[(].*[*]RI$", paste0("RI ~~ ", lab("vR"), "*RI"))
  syn
}
defined_d <- function(n, ri_mean_free) {
  G <- length(n); out <- character()
  for (k in seq(2, G)) {
    out <- c(out, sprintf("dS_g%d := mS_g%d / sqrt((%d + %d * vS_g%d) / %d)", k, k, n[1], n[k], k, n[1] + n[k]))
    if (ri_mean_free) out <- c(out, sprintf("dR_g%d := mR_g%d / sqrt((%d * vR_g1 + %d * vR_g%d) / %d)", k, k, n[1], n[k], k, n[1] + n[k]))
  }
  out
}

run_means <- function(key, inv) {
  gr <- GROUPS[[key]]; keep <- !is.na(gr$g)
  mg <- sf2b_merge_sparse(X[keep, I8, drop = FALSE], gr$g[keep])
  ncat <- apply(mg$X, 2, function(x) length(unique(x[!is.na(x)])))
  d <- as.data.frame(mg$X); d$.group <- as.character(gr$g[keep]); d <- d[order(d$.group), ]
  G <- length(unique(d$.group)); n <- as.integer(table(d$.group)); labs <- names(table(d$.group))
  last <- tail(inv$table, 1)
  ref <- list(cfi = last$cfi, rmsea = last$rmsea)
  # Two starts per model (the loadings-step solution and lavaan's default), a check that
  # the means model is identified: both should reach the same fmin AND the same S mean.
  # Kept: the converged fit with the lower fmin (on a tie, the one whose parameter
  # covariance matrix is closer to positive definite). Recorded per fit: its start, both
  # fmins, the smallest eigenvalue of its covariance matrix relative to the largest
  # (equality labels make it singular, so ~0 is expected; clearly negative is flagged),
  # and the gap in the S mean between starts (a gap at equal fmin = not identified).
  fl <- sf2b_try(sf2b_sem(sf2b_inv_syntax(I8, ncat, G, "loadings", inv$free), d, I8, group = ".group",
                          auto.fix.first = FALSE, control = list(iter.max = 3000)))
  eig_rel <- function(f) { ev <- tryCatch(eigen(lavInspect(f, "vcov"), symmetric = TRUE, only.values = TRUE)$values, error = function(e) NA_real_)
    min(ev) / max(abs(ev)) }
  fmin_of <- function(f) if (sf2b_failed(f)) NA_real_ else fitMeasures(f, "fmin")[[1]]
  mS_of <- function(f) if (sf2b_failed(f)) NA_real_ else { pe <- parameterEstimates(f); pe$est[pe$label == "mS_g2"][1] }
  fit_m <- function(free_int, ri_free) {
    syn <- c(means_syntax(I8, ncat, G, inv$free, free_int, ri_free), defined_d(n, ri_free))
    one <- function(start) { a <- list(syn, d, I8, group = ".group", auto.fix.first = FALSE, control = list(iter.max = 3000))
      if (inherits(start, "lavaan")) a$start <- start
      sf2b_try(do.call(sf2b_sem, a)) }
    fs <- list(loadings = if (!sf2b_failed(fl)) one(fl) else simpleError("no loadings fit"), default = one(NULL))
    fm <- vapply(fs, fmin_of, 0); er <- vapply(fs, function(f) if (sf2b_failed(f)) NA_real_ else eig_rel(f), 0)
    ok <- !is.na(fm)
    pick <- if (!any(ok)) "default" else if (sum(ok) == 1) names(fs)[ok] else
      if (abs(diff(fm)) < 1e-7) names(fs)[which.max(er)] else names(fs)[which.min(fm)]
    f <- fs[[pick]]
    attr(f, "check") <- list(start = pick, fmin_loadings_start = fm[["loadings"]], fmin_default_start = fm[["default"]],
                             vcov_min_eig_rel = er[[pick]], vcov_not_pd = isTRUE(er[[pick]] < -1e-6),
                             mS_gap_between_starts = abs(mS_of(fs$loadings) - mS_of(fs$default)))
    f
  }
  row <- function(f, st, released) {
    fs <- sf2b_fit_stats(f)
    r <- data.table(analysis = gr$label, step = st, released = paste(released, collapse = " "), converged = !sf2b_failed(f),
                    error = sf2b_msg(f), chisq = fs$chisq.scaled, df = fs$df, cfi = fs$cfi.scaled, rmsea = fs$rmsea.scaled,
                    srmr = fs$srmr, d_cfi = round(fs$cfi.scaled - ref$cfi, 4), d_rmsea = round(fs$rmsea.scaled - ref$rmsea, 4))
    r[, noninvariant := !converged | (!is.na(d_cfi) & d_cfi < SF2B$inv_d_cfi & d_rmsea > SF2B$inv_d_rmsea)]
    ck <- attr(f, "check"); if (!is.null(ck)) r <- cbind(r, as.data.table(ck))
    r
  }
  rows <- list(); means <- list()
  for (ri_free in c(TRUE, FALSE)) {
    st <- if (ri_free) "intercepts, means (RI mean free)" else "intercepts, means (RI mean fixed 0)"
    free_int <- character(); f <- fit_m(free_int, ri_free); r <- row(f, st, free_int); rows[[length(rows) + 1]] <- r
    while (isTRUE(r$noninvariant) && length(free_int) < SF2B$partial_max_free) {
      cand <- setdiff(I8, free_int)
      cf <- vapply(cand, function(i) { f2 <- fit_m(c(free_int, i), ri_free); if (sf2b_failed(f2)) NA_real_ else fitMeasures(f2, "cfi.scaled")[[1]] }, 0)
      if (all(is.na(cf))) break
      free_int <- c(free_int, cand[which.max(cf)])
      f <- fit_m(free_int, ri_free); r <- row(f, paste0(st, " (partial)"), free_int); rows[[length(rows) + 1]] <- r
    }
    if (!sf2b_failed(f)) {
      pe <- as.data.table(parameterEstimates(f))
      get <- function(lab) { x <- pe[label == lab]; if (!nrow(x)) return(c(NA, NA, NA, NA)); c(x$est, x$se, x$ci.lower, x$ci.upper) }
      for (k in seq_len(G)) {
        if (k == 1) { mS <- c(0, NA, NA, NA); vS <- c(1, NA, NA, NA); mR <- c(0, NA, NA, NA); dS <- dR <- c(NA, NA, NA, NA) }
        else { mS <- get(paste0("mS_g", k)); vS <- get(paste0("vS_g", k)); mR <- if (ri_free) get(paste0("mR_g", k)) else c(0, NA, NA, NA)
          dS <- get(paste0("dS_g", k)); dR <- if (ri_free) get(paste0("dR_g", k)) else c(NA, NA, NA, NA) }
        vR <- get(paste0("vR_g", k))
        means[[length(means) + 1]] <- data.table(comparison = key, analysis = gr$label, model = st, intercepts_released = paste(free_int, collapse = " "),
          model_noninvariant = r$noninvariant, start = r$start, vcov_not_pd = r$vcov_not_pd, group = labs[k], n = n[k], reference = labs[1],
          S_mean = mS[1], S_mean_se = mS[2], S_var = vS[1], S_d = dS[1], S_d_lo = dS[3], S_d_hi = dS[4],
          RI_mean = mR[1], RI_mean_se = mR[2], RI_var = vR[1], RI_d = dR[1], RI_d_lo = dR[3], RI_d_hi = dR[4])
      }
    }
  }
  list(table = rbindlist(rows, fill = TRUE), means = rbindlist(means, fill = TRUE))
}

## Run ---------
INV <- list(); STEPS <- list(); MEANS <- list(); MERGES <- list()
for (k in names(GROUPS)) {
  inv <- sf2b_run_mg(X, I8, GROUPS[[k]]$g, GROUPS[[k]]$label)
  if (!"noninvariant" %in% names(inv$table)) inv$table[, noninvariant := NA]
  INV[[k]] <- inv; STEPS[[k]] <- cbind(comparison = k, inv$table)
  if (!is.null(inv$merges)) MERGES[[k]] <- cbind(comparison = k, inv$merges)
  last <- tail(inv$table, 1)
  # The means step needs the loadings step reached and held (fully or partially).
  if (grepl("^loadings", last$step) && !isTRUE(last$noninvariant)) {
    mm <- run_means(k, inv); STEPS[[paste0(k, "_means")]] <- cbind(comparison = k, mm$table); MEANS[[k]] <- mm$means
  } else sf2b_log(k, ": loadings step not held; means not compared")
}
steps <- rbindlist(STEPS, fill = TRUE); means <- rbindlist(MEANS, fill = TRUE)
num <- function(DT, k = 4) { cols <- names(DT)[vapply(DT, is.numeric, TRUE)]; DT[, (cols) := lapply(.SD, function(x) round(x, k)), .SDcols = cols]; DT }
fwrite(num(steps), sf2b_out("24_invariance.csv"))
if (length(MERGES)) fwrite(rbindlist(MERGES, fill = TRUE), sf2b_out("24_invariance_merges.csv"))
fwrite(num(means), sf2b_out("24_latent_means.csv"))

out <- list(
  status = "EXPLORATORY, NOT REGISTERED. Specified 17 November 2026 after the registered analyses and 10-23 were seen.",
  procedure = "H8's (sf2b_run_mg): configural, thresholds, loadings; non-invariant when delta scaled CFI < -.010 and delta scaled RMSEA > .015; partial invariance up to 4 items. Added: intercepts with S and RI means free (and RI mean fixed 0), judged against the loadings step by the same rule.",
  groups = lapply(GROUPS, function(g) list(label = g$label, n = g$n, below_registered_minimum = g$below_min)),
  results = lapply(names(INV), function(k) {
    t <- steps[comparison == k]
    list(comparison = k, released = INV[[k]]$free,
         last_registered_step = tail(INV[[k]]$table$step, 1),
         noninvariant_steps = t[noninvariant == TRUE, step],
         steps = t[, .(step, released, converged, cfi, rmsea, d_cfi, d_rmsea, noninvariant)])
  }),
  latent_means = means,
  runtime_min = round(as.numeric(difftime(Sys.time(), T0, units = "mins")), 1),
  session = list(R = R.version.string, lavaan = as.character(packageVersion("lavaan")))
)
names(out$results) <- names(INV)
write_json(out, sf2b_out("24_exploratory_invariance.json"), auto_unbox = TRUE, pretty = TRUE, digits = 4, na = "null")
print(steps[, .(comparison, step, released, cfi, rmsea, d_cfi, d_rmsea, noninvariant)])
print(means[, .(comparison, model, group, n, S_d, S_d_lo, S_d_hi, RI_d, RI_d_lo, RI_d_hi, S_var, RI_var)])
sf2b_log("24 done in ", out$runtime_min, " min")
