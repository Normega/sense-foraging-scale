# 11: follow-up analyses after the registered analysis. EXPLORATORY and NOT PREREGISTERED: run on
# 11-12 October 2026, with Norm, to decide how to score the short form and what to put in Study 2's
# candidate block. Writes 11_followup_study1.json for senseforaging.com/research/study1.
#
# Real mode only: it reads background answers (practice years, religiosity, spirituality, tradition)
# and per-answer timing that the registered pipeline does not use. Uses 04_primary_loadings_k04 (the
# primary model) and, when present, the k05/k06 ESEMs (run 04 with SF1_K=5 and 6 first) and the
# plain EFAs from 03.

source("study1/analysis/config.R")
suppressMessages({ library(psych); library(data.table); library(jsonlite) })
stopifnot(SF1$mode == "real")
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
set.seed(20261012)

A <- readRDS(sf1_path("derived", "sf1_real.rds"))
m <- as.data.table(A$meta); keep <- !m$excluded_attention %in% TRUE
k <- A$items_keyed[keep, ]; raw <- A$items_raw[keep, ]; m <- m[keep]
stopifnot(identical(as.character(m$participant), rownames(k)))
txt <- A$definition$text; rev <- A$definition$reversed; pairs <- A$definition$pairs
out_real <- function(f) fread(sf1_path("output", "real", f))

# ---- background fields and answer timing, from the raw extraction -------------------------------
raw_file <- tail(sort(list.files(sf1_path("raw"), "^extract_.*\\.json$", full.names = TRUE)), 1)
rows <- fromJSON(raw_file, simplifyVector = FALSE)$rows
first_sub <- function(slug) {
  r <- Filter(function(x) identical(x$slug, slug) && is.null(x$resubmission_of), rows)
  r <- r[order(vapply(r, `[[`, "", "received_at"))]
  r[!duplicated(vapply(r, `[[`, "", "participant"))]
}
num1 <- function(v) if (is.numeric(v)) v else NA_real_
bg <- rbindlist(lapply(first_sub("sf-background"), function(x) {
  b <- x$responses
  y <- b$bg_practice_years
  data.table(participant = x$participant,
             years = if (is.list(y) && identical(y$option_id, "years")) num1(y$value) else NA_real_,
             spiritual = num1(b$bg_spiritual), religious = num1(b$bg_religious),
             tradition = if (is.list(b$bg_tradition)) b$bg_tradition$option_id %||% NA_character_ else NA_character_)
}))
m <- merge(m, bg, by = "participant", all.x = TRUE, sort = FALSE)
stopifnot(identical(as.character(m$participant), rownames(k)))

# Per-item latency: the gap since the previous first answer on the same page (the first answer on a
# page also carries the page's reading time, so it is dropped); changed = answered more than once.
timing <- rbindlist(lapply(first_sub("sf-pool6"), function(x) {
  t <- x$responses[["_timing"]]; if (is.null(t)) return(NULL)
  ans <- t$answers; pages <- t$pages
  rbindlist(lapply(pages, function(pg) {
    st <- pg[[2]]; en <- pg[[3]]
    first <- vapply(ans, function(v) as.numeric(v[[1]]), 0); n <- vapply(ans, length, 0L)
    on <- which(first >= st & first < en); if (!length(on)) return(NULL)
    o <- on[order(first[on])]
    data.table(participant = x$participant, item = names(ans)[o], rank = seq_along(o) - 1L,
               latency = first[o] - c(st, head(first[o], -1)), clicks = n[o])
  }))
}))
timing <- timing[participant %in% m$participant & grepl("^sf6_[0-9]+$", item)]

# ---- shared definitions ---------------------------------------------------------------------------
SF <- list(Access = c(35, 36, 9, 8, 7), Understanding = c(22, 3, 15, 31, 37),
           Valuing = c(12, 18, 11, 27, 44), LettingBe = c(40, 46, 43, 47, 21))
ids <- function(v) paste0("sf6_", v)
score <- function(v) rowMeans(k[, ids(v), drop = FALSE], na.rm = TRUE)
S <- as.data.table(lapply(SF, score)); S[, Total20 := score(unlist(SF))]; S[, Total15 := score(unlist(SF[1:3]))]
SC <- c("Access", "Understanding", "Valuing", "LettingBe")
alpha_of <- function(X) { C <- cor(X, use = "pair"); p <- ncol(C); (p / (p - 1)) * (1 - p / sum(C)) }
mean_r <- function(X) { C <- cor(X, use = "pair"); mean(C[upper.tri(C)]) }
r_p <- function(x, y) { t <- suppressWarnings(cor.test(x, y)); list(r = unname(round(t$estimate, 3)), p = signif(t$p.value, 3)) }
d_p <- function(x, g) {
  a <- x[g %in% TRUE]; b <- x[g %in% FALSE]; a <- a[!is.na(a)]; b <- b[!is.na(b)]
  sp <- sqrt(((length(a) - 1) * var(a) + (length(b) - 1) * var(b)) / (length(a) + length(b) - 2))
  list(d = round((mean(a) - mean(b)) / sp, 3), p = signif(t.test(a, b)$p.value, 3), n1 = length(a), n0 = length(b))
}
freq_map <- c(none = 0, lt_weekly = 0.5, `1_2_week` = 1.5, `3_5_week` = 4, daily = 6.5)
m[, freq := ifelse(practice_regular == "no", 0, freq_map[as.character(practice_freq)])]
m[, log_years := ifelse(practice_regular == "no", 0, log1p(pmin(years, 40)))]
m[, training := as.numeric(practice_training == "yes")][practice_regular == "no", training := 0]
m[, weekly := ifelse(practice_regular == "no", 0, ifelse(practice_freq %in% c("1_2_week", "3_5_week", "daily"), 1, NA))]
m[, ever := as.numeric(practice_regular == "yes")]
m[, dedicated := ifelse(practice_regular == "no", 0,
                        ifelse(practice_freq %in% c("3_5_week", "daily") & !is.na(years) & years >= 1 & practice_training == "yes", 1, NA))]
anglo <- m$country %in% c("uk", "us", "ca", "au", "ie", "nz")
acq <- rowMeans(raw[, c(pairs[, 1], pairs[, 2])], na.rm = TRUE)   # agreeing with a statement AND its opposite
out <- list(generated = format(Sys.time(), "%Y-%m-%d"), n = nrow(k),
            note = "Exploratory follow-up analyses, not preregistered (11-12 October 2026).")

# ---- 1. reliability of the decided short form -----------------------------------------------------
rel <- rbindlist(lapply(c(names(SF), "Total20", "Total15"), function(n) {
  v <- switch(n, Total20 = unlist(SF), Total15 = unlist(SF[1:3]), SF[[n]]); X <- k[, ids(v)]
  pc <- suppressWarnings(polychoric(X)$rho)
  data.table(scale = n, n_items = length(v), alpha = round(alpha_of(X), 3),
             ordinal_alpha = round((ncol(pc) / (ncol(pc) - 1)) * (1 - ncol(pc) / sum(pc)), 3), mean_r = round(mean_r(X), 3),
             alpha_anglo = round(alpha_of(X[anglo, ]), 3), alpha_other = round(alpha_of(X[!anglo, ]), 3))
}))
om <- suppressWarnings(omega(polychoric(k[, ids(unlist(SF))])$rho, nfactors = 4, n.obs = nrow(k), plot = FALSE))
drop1 <- lapply(SF[c("Understanding", "LettingBe")], function(v) {
  a <- suppressWarnings(alpha(k[, ids(v)], warnings = FALSE))$alpha.drop
  data.table(item = as.integer(sub("sf6_", "", rownames(a))), alpha_if_dropped = round(a$raw_alpha, 3))
})
sc_cor <- round(cor(S[, ..SC], use = "pair"), 3)
out$reliability <- list(scales = rel, omega_total = round(om$omega.tot, 3), omega_h = round(om$omega_h, 3),
                        alpha_if_dropped = drop1, subscale_correlations = as.data.table(sc_cor, keep.rownames = "scale"))

# ---- 2. practice ----------------------------------------------------------------------------------
pv <- list(freq = m$freq, log_years = m$log_years, training = m$training, spiritual = m$spiritual,
           contemplative_posting = as.numeric(m$posting == "contemplative"))
out$practice_r <- rbindlist(lapply(names(pv), function(v)
  c(list(variable = v), unlist(lapply(setNames(c(SC, "Total15"), c(SC, "Total15")), function(s) r_p(pv[[v]], S[[s]])), recursive = FALSE))))
out$dedicated <- rbindlist(lapply(c(SC, "Total15"), function(s) c(list(scale = s), d_p(S[[s]], m$dedicated == 1))))
auc <- function(p, y) { r <- rank(p); n1 <- sum(y == 1); n0 <- sum(y == 0); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
cv <- function(fml, data, binary, reps = 20, K = 10) {
  data <- model.frame(fml, data); y <- model.response(data)
  mean(replicate(reps, { fold <- sample(rep(1:K, length.out = nrow(data))); pred <- numeric(nrow(data))
    for (j in 1:K) { fit <- glm(fml, data = data[fold != j, ], family = if (binary) binomial() else gaussian())
                     pred[fold == j] <- predict(fit, data[fold == j, ], type = "response") }
    if (binary) auc(pred, y) else cor(pred, y)^2 }))
}
PD <- cbind(m[, .(ever, weekly, dedicated, freq, log_years, spiritual)], S)
models <- c(four_subscales = "Access + Understanding + Valuing + LettingBe", access_valuing = "Access + Valuing",
            total15 = "Total15", spiritual = "spiritual", spiritual_plus_four = "spiritual + Access + Understanding + Valuing + LettingBe")
pred <- function(o, binary) {
  dd <- PD[!is.na(PD[[o]]) & !is.na(spiritual)]
  c(list(outcome = o, n = nrow(dd), cases = if (binary) sum(dd[[o]]) else NA),
    lapply(models, function(f) round(cv(as.formula(paste(o, "~", f)), dd, binary), 3)))
}
out$prediction <- list(binary = rbindlist(lapply(c("ever", "weekly", "dedicated"), pred, binary = TRUE)),
                       continuous = rbindlist(lapply(c("freq", "log_years"), pred, binary = FALSE)))
dd <- PD[!is.na(dedicated)]; p <- fitted(glm(dedicated ~ Access + Valuing, data = dd, family = binomial()))
ths <- sort(unique(p)); J <- sapply(ths, function(t) mean(p[dd$dedicated == 1] >= t) + mean(p[dd$dedicated == 0] < t) - 1)
t <- ths[which.max(J)]
out$prediction$dedicated_cut <- list(sensitivity = round(mean(p[dd$dedicated == 1] >= t), 2), specificity = round(mean(p[dd$dedicated == 0] < t), 2))

# ---- 3. demographics ------------------------------------------------------------------------------
dem <- list(
  list("women_vs_men", "d", function(s) d_p(s, ifelse(m$gender == "woman", TRUE, ifelse(m$gender == "man", FALSE, NA)))),
  list("age", "r", function(s) r_p(as.numeric(m$age), s)),
  list("degree_vs_none", "d", function(s) d_p(s, ifelse(m$education %in% c("bachelors", "graduate"), TRUE,
                                                        ifelse(m$education %in% c("lt_hs", "hs", "some_college"), FALSE, NA)))),
  list("anglo_vs_other", "d", function(s) d_p(s, anglo)),
  list("religious", "r", function(s) r_p(m$religious, s)),
  list("spiritual", "r", function(s) r_p(m$spiritual, s)))
out$demographics <- rbindlist(lapply(dem, function(x) c(list(contrast = x[[1]], stat = x[[2]]),
  unlist(lapply(setNames(SC, SC), function(s) { v <- x[[3]](S[[s]]); list(est = v$d %||% v$r, p = v$p) }), recursive = FALSE))))
TR <- cbind(m[, .(tradition)], S[, ..SC]); TR[, AccVal := (S$Access + S$Valuing) / 2]
out$tradition <- TR[, c(.(n = .N), lapply(.SD, function(v) round(mean(v, na.rm = TRUE), 2))), by = tradition][n >= 15][order(-AccVal)]

# ---- 4. general agreement in the averaged scores --------------------------------------------------
ctrl <- function(s, x) { r <- resid(lm(s ~ acq, na.action = na.exclude)); round(cor(r, x, use = "pair"), 3) }
out$agreement <- list(
  index_mean = round(mean(acq), 2), index_sd = round(sd(acq), 2),
  with_scales = lapply(setNames(SC, SC), function(s) round(cor(S[[s]], acq, use = "pair"), 3)),
  with = lapply(list(religious = m$religious, spiritual = m$spiritual, freq = m$freq, age = as.numeric(m$age)),
                function(x) round(cor(acq, x, use = "pair"), 3)),
  anglo_d = round((mean(acq[anglo]) - mean(acq[!anglo])) / sd(acq), 3),
  controlled = rbindlist(lapply(c("Access", "Valuing", "LettingBe"), function(s) data.table(
    scale = s, r_acq = round(cor(S[[s]], acq, use = "pair"), 3),
    spiritual_raw = round(cor(S[[s]], m$spiritual, use = "pair"), 3), spiritual_ctrl = ctrl(S[[s]], m$spiritual),
    religious_raw = round(cor(S[[s]], m$religious, use = "pair"), 3), religious_ctrl = ctrl(S[[s]], m$religious),
    freq_raw = round(cor(S[[s]], m$freq, use = "pair"), 3), freq_ctrl = ctrl(S[[s]], m$freq),
    age_raw = round(cor(S[[s]], as.numeric(m$age), use = "pair"), 3), age_ctrl = ctrl(S[[s]], as.numeric(m$age))))))

# ---- 5. items: wording, agreement rates, timing ---------------------------------------------------
L4 <- out_real("04_primary_loadings_k04.csv"); lk <- as.matrix(L4[, .(F1, F2, F3, F4)]); rownames(lk) <- L4$item
primary <- apply(abs(lk), 1, max)
TT <- timing[rank > 0 & latency > 300 & latency < 60000][, ll := log(latency)][, llc := ll - mean(ll), by = participant]
TI <- TT[, .(latency_s = round(median(latency) / 1000, 1), latency_resid = round(mean(llc), 3)), by = item]
CH <- timing[, .(changed_pct = round(100 * mean(clicks > 1), 1)), by = item]
sfi <- unlist(lapply(names(SF), function(n) setNames(rep(n, 5), ids(SF[[n]]))))
itr <- function(i, g = rep(TRUE, nrow(k))) { if (!i %in% names(sfi)) return(NA_real_)
  s <- setdiff(names(sfi)[sfi == sfi[i]], i); round(cor(k[g, i], rowMeans(k[g, s, drop = FALSE], na.rm = TRUE), use = "pair"), 3) }
IT <- data.table(item = A$definition$items)
IT[, `:=`(num = as.integer(sub("sf6_", "", item)), subscale = unname(sfi[item]), reversed = item %in% rev,
          words = lengths(strsplit(txt[item], " ")), says_mode = grepl("mode", txt[item]),
          mean = round(colMeans(raw[, item], na.rm = TRUE), 2),
          agree_pct = round(100 * colMeans(raw[, item] >= 4, na.rm = TRUE)),
          loading = round(primary[item], 3))]
IT[, `:=`(item_total_r = sapply(item, itr), item_total_r_anglo = sapply(item, itr, g = anglo), item_total_r_other = sapply(item, itr, g = !anglo),
          r_practice = round(sapply(item, function(i) cor(k[, i], m$weekly, use = "pair")), 3),
          r_spiritual_raw = round(sapply(item, function(i) cor(raw[, i], m$spiritual, use = "pair")), 3))]
IT <- merge(merge(IT, TI, by = "item", all.x = TRUE), CH, by = "item", all.x = TRUE); setorder(IT, num)
coefs <- function(f) { s <- summary(lm(f, IT))$coef; as.list(setNames(round(s[, 1], 3), rownames(s))) |>
  c(list(p = as.list(setNames(signif(s[, 4], 3), rownames(s))))) }
out$items <- list(
  table = IT,
  mode_vs_not = IT[, .(n = .N, latency_resid = round(mean(latency_resid), 3), changed_pct = round(mean(changed_pct), 1),
                       loading = round(mean(loading), 3), r_practice_abs = round(mean(abs(r_practice)), 3)), by = says_mode],
  latency_model = coefs(latency_resid ~ words + says_mode + reversed),
  changed_model = coefs(changed_pct ~ words + says_mode + reversed),
  loading_model = coefs(loading ~ words + says_mode + reversed))

# ---- 6. the most reliable item sets ---------------------------------------------------------------
orient <- c(F1 = 1, F2 = 1, F3 = 1, F4 = -1); NM <- c(F1 = "Access", F2 = "Understanding", F3 = "Valuing", F4 = "LettingBe")
lkk <- lk; lkk[rev, ] <- -lkk[rev, ]; lkk <- sweep(lkk, 2, orient, "*")
prim <- colnames(lkk)[apply(abs(lkk), 1, which.max)]; names(prim) <- rownames(lkk)
pval <- apply(lkk, 1, function(r) r[which.max(abs(r))])
greedy <- function(X, kk, need_rev = 0) {
  C <- cor(X, use = "pair"); diag(C) <- NA; ij <- which(C == max(C, na.rm = TRUE), arr.ind = TRUE)[1, ]; sel <- colnames(X)[ij]
  while (length(sel) < kk) {
    left <- setdiff(colnames(X), sel)
    if (need_rev > sum(sel %in% rev) && kk - length(sel) <= need_rev - sum(sel %in% rev)) left <- intersect(left, rev)
    if (!length(left)) break
    sel <- c(sel, names(which.max(sapply(left, function(i) alpha_of(X[, c(sel, i)])))))
  }
  sel
}
cvalpha <- function(X, kk, need_rev, reps = 20) mean(replicate(reps, {
  h <- sample(nrow(X), nrow(X) %/% 2); alpha_of(X[-h, greedy(X[h, ], kk, need_rev)]) }))
setrow <- function(label, s, X = NULL, kk = NULL, nr = 0) { sc <- rowMeans(k[, s], na.rm = TRUE)
  data.table(set = label, k = length(s), items = paste(sub("sf6_", "", s), collapse = " "), n_reversed = sum(s %in% rev),
             alpha = round(alpha_of(k[, s]), 3), cv_alpha = if (is.null(X)) NA_real_ else round(cvalpha(X, kk, nr), 3),
             r_agreement = round(cor(sc, acq, use = "pair"), 3), r_practice = round(cor(sc, m$weekly, use = "pair"), 3)) }
mr <- list()
for (f in names(NM)) {
  pool <- names(prim)[prim == f & pval >= .30]; X <- k[, pool]
  for (kk in c(5, 8)) if (kk <= length(pool)) for (nr in unique(c(0, min(2, sum(pool %in% rev))))) {
    s <- greedy(X, kk, nr); mr[[length(mr) + 1]] <- cbind(data.table(factor = NM[f], pool = length(pool)),
      setrow(if (nr == 0) "most reliable" else paste0("with ", nr, " reversed"), s, X, kk, nr)) }
}
out$max_reliability <- list(
  by_factor = rbindlist(mr),
  current = rbindlist(lapply(names(SF), function(n) cbind(data.table(factor = n, pool = NA_integer_), setrow("current short form", ids(SF[[n]]))))),
  proposed = rbind(cbind(data.table(factor = "Access", pool = NA_integer_), setrow("proposed: 8, 2 reversed, no near-duplicates", ids(c(9, 7, 6, 2, 8, 14, 36, 35)))),
                   cbind(data.table(factor = "Valuing", pool = NA_integer_), setrow("proposed: 6, 1 reversed", ids(c(18, 12, 25, 27, 28, 44))))),
  general = rbindlist(lapply(c(8, 12), function(kk) cbind(data.table(factor = "all 48", pool = 48L), setrow("most reliable", greedy(k, kk), k, kk, 0)))))

# ---- 7. foraging items ----------------------------------------------------------------------------
FI <- ids(c(24, 25, 26, 27, 28, 30, 31, 32)); fs <- rowMeans(k[, FI], na.rm = TRUE)
efa <- function(kk) { f <- sf1_path("output", "real", sprintf("03_efa_k%02d_loadings.csv", kk)); if (!file.exists(f)) return(NULL)
  L <- fread(f)[item %in% FI]; fcols <- grep("^f[0-9]+$", names(L), value = TRUE)
  best <- names(which.max(colSums(abs(as.matrix(L[item %in% ids(c(25, 26, 27, 30, 31)), ..fcols])))))
  data.table(k = kk, item = as.integer(sub("sf6_", "", L$item)), loading_on_discovery_factor = round(L[[best]], 3)) }
esem <- function(kk) { f <- sf1_path("output", "real", sprintf("04_primary_loadings_k%02d.csv", kk)); if (!file.exists(f)) return(NULL)
  L <- fread(f); fc <- grep("^F[0-9]+$", names(L), value = TRUE); M <- as.matrix(L[, ..fc]); rownames(M) <- L$item
  pr <- fc[apply(abs(M), 1, which.max)]; names(pr) <- L$item
  members <- function(f) { it <- names(pr)[pr == f & apply(abs(M), 1, max) >= .30]; paste(sort(as.integer(sub("sf6_", "", it))), collapse = " ") }
  data.table(k = kk, item = as.integer(sub("sf6_", "", FI)), factor = unname(pr[FI]), loading = round(apply(abs(M[FI, ]), 1, max), 3),
             factor_members = vapply(unname(pr[FI]), members, "")) }
out$forage <- list(
  items = FI, alpha = round(alpha_of(k[, FI]), 3), mean_r = round(mean_r(k[, FI]), 3),
  r = list(access = round(cor(fs, S$Access, use = "pair"), 3), valuing = round(cor(fs, S$Valuing, use = "pair"), 3),
           agreement = round(cor(fs, acq, use = "pair"), 3), practice = round(cor(fs, m$weekly, use = "pair"), 3)),
  agree_pct = as.list(IT[item %in% FI, setNames(agree_pct, num)]),
  plain_efa = rbindlist(lapply(c(5, 6, 7), efa)),
  esem_with_agreement = rbindlist(lapply(c(4, 5, 6), esem)),
  fit_k6 = if (file.exists(sf1_path("output", "real", "04_primary_fit_k06.csv"))) out_real("04_primary_fit_k06.csv")[1] else NULL)

f <- sf1_out("11_followup_study1.json")
write_json(out, f, auto_unbox = TRUE, digits = NA, pretty = TRUE, na = "null", dataframe = "rows")
sf1_log("follow-up: ", f)
