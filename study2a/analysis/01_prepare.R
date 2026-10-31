# 01: from the raw extraction to one analysis row per participant.
#
# Registered rules applied here (Section 5, "Data exclusion"):
#   excluded  - staff test enrollments (already out of the extraction);
#               sessions without a submitted scale;
#               failing two or more of the three attention checks
#   first     - each questionnaire's first complete submission is used
#   flagged   - failing exactly one check; < 2 s per statement on the scale;
#               a longstring >= 15 on the 61 scale positions; "Prefer not to
#               answer" on more than 11 of the 58 statements
# Also recorded, not a registered flag: a step crash in the session (Study 1's
# translation-crash group; the browser was machine-translating the page).
# In permuted and simulated modes the answers are replaced BEFORE any
# answer-derived index (longstring, separation) is computed.

source("study2a/analysis/config.R")
def <- sf2a_definition()

raw_file <- tail(sort(list.files(SF2A$data_dir, "^responses_.*\\.json$", full.names = TRUE)), 1)
stopifnot(length(raw_file) == 1)
rows <- fromJSON(raw_file, simplifyVector = FALSE)$rows
sf2a_log("extraction: ", basename(raw_file), ", ", length(rows), " rows")

num <- function(v) if (is.null(v) || is.list(v) || is.character(v)) NA_real_ else as.numeric(v)
opt <- function(v) if (is.list(v) && !is.null(v$option_id)) v$option_id else if (is.character(v)) v else NA_character_
is_pna <- function(v) identical(v, "pna")

by_person <- split(rows, vapply(rows, `[[`, "", "participant"))
recs <- lapply(by_person, function(rs) {
  p <- rs[[1]]
  rs <- Filter(function(r) !is.null(r$slug), rs)      # a session row with no submissions
  pick <- function(slug, need) {
    cand <- Filter(function(r) identical(r$slug, slug) && is.null(r$resubmission_of) &&
                     all(need %in% names(r$responses)), rs)
    if (!length(cand)) return(NULL)
    cand[[order(vapply(cand, `[[`, "", "received_at"))[1]]]
  }
  n_subs <- function(slug) sum(vapply(rs, function(r) identical(r$slug, slug), TRUE))
  scale <- pick("sf2a-scale", def$positions)
  vig <- pick("sf2a-vignettes", def$vignettes)
  bg <- pick("sf-background", "bg_practice_regular")
  ius <- pick("sf2-ius12", def$ius)
  end <- pick("sf2a-end", def$end)
  crashes <- if (is.null(scale)) 0 else sum(vapply(Filter(function(r) identical(r$schedule_id, scale$schedule_id), rs),
                                                   function(r) as.numeric(r$step_crashes %||% 0), 0) > 0)
  out <- list(participant = p$participant, posting = unname(SF2A$prolific_study[p$prolific_study_id %||% ""]) %||% NA_character_,
              repository_consent = isTRUE(p$repository_consent), has_scale = !is.null(scale),
              n_scale_submissions = n_subs("sf2a-scale"), n_vignette_submissions = n_subs("sf2a-vignettes"),
              n_sessions = length(unique(vapply(rs, `[[`, "", "schedule_id"))), step_crash = crashes > 0,
              has_vignettes = !is.null(vig), has_ius = !is.null(ius), has_end = !is.null(end))
  if (!is.null(bg)) {
    b <- bg$responses
    out$practice_regular <- opt(b$bg_practice_regular); out$practice_freq <- opt(b$bg_practice_freq)
    out$age <- if (is.list(b$bg_age)) num(b$bg_age$value) else NA_real_
    out$gender <- opt(b$bg_gender); out$education <- opt(b$bg_education); out$country <- opt(b$bg_country)
    out$religious <- num(b$bg_religious); out$spiritual <- num(b$bg_spiritual); out$tradition <- opt(b$bg_tradition)
  }
  if (!is.null(scale)) {
    s <- scale$responses
    out$items <- vapply(def$items, function(i) num(s[[i]]), 0)
    out$pna <- sum(vapply(def$items, function(i) is_pna(s[[i]]), TRUE))
    out$checks <- vapply(names(def$checks), function(i) num(s[[i]]), 0)
    out$positions <- vapply(def$positions, function(i) as.numeric(s[["_order"]]$positions[[i]] %||% NA), 0)
    out$order_method <- s[["_order"]]$method %||% NA_character_
    out$scale_secs <- as.numeric(s[["_timing"]]$total_ms %||% NA) / 1000
    out$format <- scale$presentation_format %||% NA_character_
  }
  out$vig <- if (!is.null(vig)) vapply(def$vignettes, function(i) num(vig$responses[[i]]), 0) else setNames(rep(NA_real_, 10), def$vignettes)
  out$vig_pna <- if (!is.null(vig)) sum(vapply(def$vignettes, function(i) is_pna(vig$responses[[i]]), TRUE)) else NA_integer_
  out$ius <- if (!is.null(ius)) vapply(def$ius, function(i) num(ius$responses[[i]]), 0) else setNames(rep(NA_real_, 12), def$ius)
  out$end <- if (!is.null(end)) vapply(def$end, function(i) num(end$responses[[i]]), 0) else setNames(rep(NA_real_, 6), def$end)
  out
})
mats <- c("items", "checks", "positions", "vig", "ius", "end")
meta <- rbindlist(lapply(recs, function(r) r[setdiff(names(r), mats)]), fill = TRUE)
flow <- data.table(step = "enrolled (not staff tests)", n = nrow(meta))
flow <- rbind(flow, data.table(step = "submitted the scale", n = sum(meta$has_scale)))

with_scale <- recs[vapply(recs, `[[`, TRUE, "has_scale")]
M <- meta[has_scale == TRUE]
grab <- function(nm) { Z <- do.call(rbind, lapply(with_scale, `[[`, nm)); rownames(Z) <- M$participant; Z }
X <- grab("items"); CH <- grab("checks"); POS <- grab("positions")
V <- grab("vig"); IUS <- grab("ius"); END <- grab("end")

M[, checks_failed := rowSums(sweep(CH, 2, def$checks[colnames(CH)], `!=`) | is.na(CH))]
M[, excluded_attention := checks_failed >= SF2A$exclude_checks_failed]
flow <- rbind(flow, data.table(step = "excluded: failed two or more attention checks", n = sum(M$excluded_attention)))

set.seed(SF2A$seed)
if (SF2A$mode == "simulated") {
  n <- nrow(X); k <- length(def$factors); p <- length(def$items)
  Phi <- matrix(.3, k, k); diag(Phi) <- 1
  Lambda <- matrix(0, p, k, dimnames = list(def$items, def$factors))
  for (j in seq_len(p)) Lambda[j, def$factor[def$items[j]]] <- runif(1, .55, .75) * ifelse(def$items[j] %in% def$reversed, -1, 1)
  eta <- MASS::mvrnorm(n, rep(0, k), Phi)
  acq <- rnorm(n, 0, sqrt(.10))
  uniq <- sqrt(pmax(.05, 1 - rowSums((Lambda %*% Phi) * Lambda)))
  Y <- eta %*% t(Lambda) + acq + sweep(matrix(rnorm(n * p), n), 2, uniq, `*`)
  # Cut each statement at its real marginal proportions; real gaps stay missing.
  for (j in seq_len(p)) {
    miss <- is.na(X[, j])
    pr <- cumsum(tabulate(X[, j], nbins = 6)) / sum(!miss)
    X[, j] <- findInterval(Y[, j], quantile(Y[, j], pr[1:5])) + 1
    X[miss, j] <- NA
  }
  # Vignettes: each person's true separation ~ N(1.5, 1), centred anywhere on the
  # scale; ratings = centre +/- separation / 2 + noise, rounded into 1-6.
  sep <- rnorm(n, 1.5, 1); centre <- rnorm(n, 3.5, .5)
  sgn <- ifelse(def$vkind[def$vignettes[1:9]] == "sensing", .5, -.5)
  for (j in 1:9) { miss <- is.na(V[, j]); V[, j] <- pmin(6, pmax(1, round(centre + sgn[j] * sep + rnorm(n, 0, .8)))); V[miss, j] <- NA }
  for (Z in c("IUS", "END")) { A0 <- get(Z); A0[] <- apply(A0, 2, sample); assign(Z, A0) }
  saveRDS(list(Lambda = Lambda, Phi = Phi, acq_var = .10, separation = sep), sf2a_derived("simulation_truth"))
  sf2a_log("statements simulated from the known 5-factor + agreement model; vignettes from a known separation")
}
if (SF2A$mode == "permuted") {
  for (Z in c("X", "V", "IUS", "END")) { A0 <- get(Z); A0[] <- apply(A0, 2, sample); assign(Z, A0) }
  sf2a_log("answers permuted within every item (structure destroyed by design)")
}

# Answer-derived indices, on whichever matrices are in use.
ALL <- cbind(X, CH)[, def$positions]
longstring <- vapply(seq_len(nrow(ALL)), function(i) {
  v <- ALL[i, order(POS[i, ])]; r <- rle(ifelse(is.na(v), -1, v)); max(r$lengths[r$values != -1], 0)
}, 0)
sens <- names(def$vkind)[def$vkind == "sensing"]; doing <- names(def$vkind)[def$vkind == "doing"]
n_s <- rowSums(!is.na(V[, sens])); n_d <- rowSums(!is.na(V[, doing]))
separation <- rowMeans(V[, sens], na.rm = TRUE) - rowMeans(V[, doing], na.rm = TRUE)
separation[n_s < SF2A$min_vignettes_each | n_d < SF2A$min_vignettes_each] <- NA
M[, `:=`(longstring = longstring, secs_per_item = scale_secs / length(def$positions),
         sensing_mean = rowMeans(V[, sens], na.rm = TRUE), doing_mean = rowMeans(V[, doing], na.rm = TRUE),
         separation = separation, tells_apart = separation >= SF2A$separation_cut, vig_self = V[, "vig_self"])]
M[, `:=`(flag_one_check = checks_failed == 1,
         flag_speeder = !is.na(secs_per_item) & secs_per_item < SF2A$speeder_secs_per_item,
         flag_longstring = longstring >= SF2A$longstring_flag,
         flag_pna = pna > SF2A$pna_flag)]
M[, flagged := flag_one_check | flag_speeder | flag_longstring | flag_pna]
M[, practice_group := fifelse(practice_regular == "no", "never",
                       fifelse(practice_regular == "yes" & practice_freq %in% c("1_2_week", "3_5_week", "daily"), "current",
                       fifelse(practice_regular == "yes", "past_or_occasional", NA_character_)))]
M[, ius_total := rowSums(IUS)]

keep <- !M$excluded_attention
flow <- rbind(flow, data.table(step = "analysable", n = sum(keep)))
key <- function(Z, rv) { Z[, rv] <- 7 - Z[, rv]; Z }
analysis <- list(
  mode = SF2A$mode, extraction = basename(raw_file), definition = def,
  items_raw = X[keep, , drop = FALSE],                       # unrecoded: for models with the random intercept
  items_keyed = key(X[keep, , drop = FALSE], def$reversed),  # reverse-keyed recoded 7 - x
  positions = POS[keep, , drop = FALSE], vignettes = V[keep, , drop = FALSE],
  ius = IUS[keep, , drop = FALSE], end = END[keep, , drop = FALSE],
  meta = M[keep], meta_all = M, flow = flow)
saveRDS(analysis, sf2a_derived("sf2a"))
fwrite(flow, sf2a_out("01_flow.csv"))
flags <- M[, .(n = .N, one_check = sum(flag_one_check), speeder = sum(flag_speeder), longstring = sum(flag_longstring),
               pna_over_11 = sum(flag_pna), any_flag = sum(flagged), step_crash = sum(step_crash),
               repeat_scale = sum(n_scale_submissions > 1), two_sessions = sum(n_sessions > 1),
               any_statement_missing = sum(rowSums(is.na(X[M$participant, , drop = FALSE])) > 0)),
           by = .(excluded_attention)]
fwrite(flags, sf2a_out("01_flags.csv"))
print(flow); print(flags)
