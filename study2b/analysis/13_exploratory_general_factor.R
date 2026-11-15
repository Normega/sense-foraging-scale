# 13: a general self-description factor across the validity battery. EXPLORATORY, NOT REGISTERED.
#
# Not part of the filed registration (osf.io/d4a5u) and not in its list of exploratory
# analyses. Run after the confirmatory real-mode run, with its answers seen; nothing here
# changes a registered cell or output. Every number it writes is exploratory.
#
# Why: the 8 correlates .45-.87 with most of the battery, including NCS-6 (.45) and BIDR
# Self-Deceptive Enhancement (.48). This asks how much of that is one general factor
# running through everyone's self-descriptions, and what of the 8 is left net of it.
#
# Method.
#   1. The latent correlation matrix of the battery: the 8 (with its agreement factor, as
#      registered), the five FFMQ facets, SBI, the three BFI-2 facets, EQ Decentering,
#      SAMD, NCS-6, BIDR SDE and IM, BSSS-8, SPSQ-SF, ASI-3, the eight Brief MAIA-2
#      subscales, PSS-4, PHQ-4 and life satisfaction; plus two extension variables that do
#      not define the factor: past-week pause days and practice group (current vs never).
#      Each cell comes from a CFA of the two instruments it joins (registered settings;
#      subscales of one instrument as correlated factors; within-instrument cells from the
#      instrument's own model). A single model of all ~170 ordinal items is not attempted:
#      WLSMV's weight matrix would be ~15,000 square.
#   2. The general factor: an exploratory Schmid-Leiman bifactor of that matrix (minres,
#      oblimin), the number of group factors from parallel analysis. Chosen because it
#      works on a latent matrix assembled cell by cell (no raw-data model needed), gives
#      every measure a general loading net of the clusters the battery plainly has
#      (FFMQ/MAIA, openness, distress, social desirability), and is the standard route
#      to omega-hierarchical. Checked against: a one-factor model of the same matrix;
#      a Jennrich-Bentler bifactor rotation; Schmid-Leiman with one group factor fewer
#      and one more; the general factor defined without the 8 (the 8 as an extension
#      variable); and the first principal component of the observed scale scores.
#   3. Net of the general factor: the 8's correlation with each measure with G partialled
#      out of both, (r - l8*lx) / sqrt((1 - l8^2)(1 - lx^2)); the extension variables'
#      loadings on G by Dwyer's (1937) extension of the oblique first-order solution; the
#      practice group difference net of G as a d (the 8 regressed on G and group; the
#      group coefficient over the 8's residual SD, 07's metric); pause days as a partial r.
#      These latent quantities are point estimates (the matrix is assembled from 140
#      models, so it has no joint SE). The observed-score check gives CIs: G as the first
#      principal component of the other 28 observed scores (the 8 left out of it),
#      complete cases.
#
# Settings as registered (Section 5): ordinal items, WLSMV, theta, missing = "pairwise",
# conditional.x = FALSE, the 8's agreement factor as registered. Brief MAIA-2 enters with
# decision 46's bound on maia_12 (residual variance at least .10 of its total), the repair
# that gave the confirmatory Noticing cell; without it every Brief MAIA-2 model fails. A
# variable whose cells cannot be estimated is dropped and listed. CHOICES below are
# written to the JSON.
#
# Writes output/<tag>/13_*.csv and 13_exploratory_general_factor.json (aggregates only).
# Run from the repo root:  SF2B_MODE=real Rscript study2b/analysis/13_exploratory_general_factor.R

# Set Up ---------
## Load libraries ---------
packages <- c("jsonlite", "data.table", "lavaan", "psych", "GPArotation")
new_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(new_packages)) install.packages(new_packages)
for (thispack in packages) {
  library(thispack, character.only = TRUE, quietly = TRUE, verbose = FALSE)
}

## Settings ---------
source("study2b/analysis/config.R")
source("study2b/analysis/sf2b_models.R")
T0 <- Sys.time()

# The pipeline's own data: 04_validity.R's top-level expressions up to its first fit (as
# 10_deviation_maia.R does): A, D, M, S, DAT, score_items; nothing is fitted.
ex04 <- parse("study2b/analysis/04_validity.R", keep.source = FALSE)
first_fit <- which(vapply(ex04, function(e) grepl("^sf2b_log\\(\"fitting", deparse(e)[1]), TRUE))
stopifnot(length(first_fit) == 1)
for (i in seq_len(first_fit - 1)) eval(ex04[[i]], envir = globalenv())
stopifnot(exists("DAT"), exists("score_items"), !exists("RES"))
ex10 <- parse("study2b/analysis/10_deviation_maia.R", keep.source = FALSE)
is_def <- function(e, nm) is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name(nm))
eval(Filter(function(e) is_def(e, "improper_of"), as.list(ex10))[[1]], envir = globalenv())
improper_list <- function(f) { x <- improper_of(f)$list; x[!grepl("^std loading (L_[[:alnum:]_]+|grp) =~ ", x)] }

# Practice group as 07 codes it, numeric (an observed covariate, as in 07's MIMIC).
DAT$grp <- fifelse(M$practice_group == "current", 1, fifelse(M$practice_group == "never", 0, NA_real_))
P_CUR <- mean(DAT$grp, na.rm = TRUE)
I8 <- D$sf8

BLOCKS <- list(sf8 = "sf8", ffmq = SF2B_MEASURES$ffmq, sbi = "sbi_moment", bfi = c("bfi_aesthetic", "bfi_intellectual", "bfi_creative"),
               eq = "eq_decentering", samd = "samd_md", ncs = "ncs6", bidr = c("bidr_sde", "bidr_im"), bsss = "bsss8",
               spsq = "spsq_sd", asi = "asi_pc", maia = SF2B_MEASURES$maia, pss = "pss4", phq = "phq4_total", lifesat = "lifesat",
               pause = "crit_pause_days", grp = "grp")
EXTENSION <- c("crit_pause_days", "grp")
SINGLE_ORD <- c("lifesat", "crit_pause_days")
BATTERY <- setdiff(unlist(BLOCKS, use.names = FALSE), EXTENSION)
FOCUS <- c("ncs6", "bidr_sde", "eq_decentering", "ffmq_observing")
MAIA_BOUND <- 9   # decision 46 under theta: l^2 <= 9, std loading <= .949, residual >= .10 of total

CHOICES <- list(
  matrix = "Latent correlations, each from a CFA of the two instruments it joins (136 pair models), within-instrument cells from the instrument's own model (FFMQ, BFI-2, BIDR, Brief MAIA-2). Registered settings throughout. A single model of all items was not attempted (size).",
  the_8 = "The 8's content factor with its agreement factor (registered); in every pair model the agreement factor is uncorrelated with the other factors, so the 8's cells are 04's model for that measure.",
  bfi = "The three BFI-2 facets, not the Open-Mindedness total (it is their sum and would duplicate them).",
  maia = "All eight Brief MAIA-2 subscales, with decision 46's bound on maia_12 in every model holding it (without it the registered Brief MAIA-2 models fail). A subscale whose cells fail is dropped.",
  single_items = "Life satisfaction and past-week pause days as single-indicator ordinal factors with no residual variance (04's treatment). Practice group (current 1, never 0; 07's coding) as an observed numeric variable through a phantom factor, so its cells are latent point-biserial correlations; the 8's reproduces 07's latent d.",
  extension = "Pause days and practice group do not define the general factor; their loadings come from Dwyer's (1937) least-squares extension of the oblique first-order solution, carried to G through the second-order loadings.",
  non_positive_definite = "If the assembled matrix is not positive definite it is smoothed (psych::cor.smooth) and the smallest eigenvalue reported.",
  number_of_group_factors = "Parallel analysis (psych::fa.parallel, fa = 'fa', minres, n.obs = 720) on the latent matrix; Schmid-Leiman with that number, and with one fewer and one more as a check.",
  primary_method = "Schmid-Leiman (psych::schmid, minres, oblimin, n.obs = 720). The general factor is oriented so the 8 loads positively.",
  checks = "One-factor minres of the latent matrix; Jennrich-Bentler bifactor rotation (psych::fa, rotate = 'bifactor', group factors + 1; the general factor taken as the column with the largest sum of squared loadings); SL with the 8 left out (the 8 extended); first principal component of the observed scores (pairwise Pearson).",
  net_of_G = "Partial correlation with G out of both: (r - l8*lx)/sqrt((1 - l8^2)(1 - lx^2)); share of r accounted for = l8*lx / r. Practice group net of G: the 8 regressed on G and group (standardized, from the correlations), the group coefficient per unit of group over the 8's residual SD: 07's metric for d. Point estimates only.",
  observed_check = "G_obs = first principal component (psych::principal) of the 28 other observed scale scores, standardized, complete cases on the 28 (and on the variables tested). G_obs contains each focus measure, so its observed partial correlations remove part of that measure itself and can turn negative: read them as a bound, not an estimate. The 8's group difference adjusted for G_obs (07's sf2b_kg, d over the 8's pooled SD), the partial Spearman correlation with pause days (ranks residualized on G_obs's rank), and partial Pearson correlations with the focus measures; CIs by Fisher's z with n - 4.")

## Syntax ---------
has_maia12 <- function(scores) "maia_notworrying" %in% scores
block_syntax <- function(b, has8) {
  sc <- BLOCKS[[b]]
  if (b == "sf8") return(list(syn = sf2b_scale_syntax(I8), ord = I8, fac = c(sf8 = "S")))
  if (b %in% c("lifesat", "pause")) {
    v <- sc
    return(list(syn = c(paste0("L_", v, " =~ 1*", v), paste0(v, " ~~ 0*", v), paste0("L_", v, " ~~ 1*L_", v), if (has8) paste0("RI ~~ 0*L_", v)),
                ord = v, fac = setNames(paste0("L_", v), v)))
  }
  if (b == "grp") return(list(syn = c("L_grp =~ 1*grp", "grp ~~ 0*grp", if (has8) "RI ~~ 0*L_grp"), ord = character(), fac = c(grp = "L_grp")))
  x <- sf2b_validity_syntax(sc, D, ri = has8)
  if (!has8) x <- x[!grepl("^S ~~", x)]
  if (has_maia12(sc)) x <- c(x, "maia_notworrying =~ l_maia12*maia_12", sprintf("l_maia12^2 < %s", MAIA_BOUND))
  list(syn = x, ord = unlist(score_items[sc]), fac = setNames(sc, sc))
}
model_syntax <- function(bs) {
  has8 <- "sf8" %in% bs
  parts <- lapply(bs, block_syntax, has8 = has8)
  list(syn = unlist(lapply(parts, `[[`, "syn")), ord = unique(unlist(lapply(parts, `[[`, "ord"))), fac = unlist(lapply(parts, `[[`, "fac")),
       scores = unlist(lapply(bs, function(b) if (b %in% c("sf8", "lifesat", "pause", "grp")) character() else BLOCKS[[b]])), has8 = has8)
}
# Orientation as keyed (sf2b_orientation when the 8 is in the model; its loop otherwise).
orient_bad <- function(f, m) {
  if (m$has8) return(sf2b_orientation(f, I8, m$scores, D))
  L <- lavInspect(f, "std")$lambda; its <- sf2b_score_items(D); rv <- unlist(lapply(D$validity, `[[`, "reverse"))
  Filter(function(sc) sum(L[its[[sc]], sc] * ifelse(its[[sc]] %in% rv, -1, 1)) <= 0, m$scores)
}
run_model <- function(bs) {
  m <- model_syntax(bs)
  t0 <- Sys.time()
  f <- sf2b_try(sf2b_sem(m$syn, DAT, ordered = m$ord))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  res <- list(blocks = bs, ok = FALSE, secs = secs)
  if (sf2b_failed(f)) return(c(res, error = sf2b_msg(f)))
  bad <- orient_bad(f, m)
  if (length(bad)) return(c(res, error = paste("not oriented as keyed:", paste(bad, collapse = ", "))))
  C <- lavInspect(f, "cor.lv")[m$fac, m$fac, drop = FALSE]; dimnames(C) <- list(names(m$fac), names(m$fac))
  out <- c(res[setdiff(names(res), "ok")], list(ok = TRUE, error = "", cor = C, fit = sf2b_fit_stats(f), improper = paste(improper_list(f), collapse = "; ")))
  if (has_maia12(m$scores)) out$maia12_std_residual <- diag(lavInspect(f, "std")$theta)[["maia_12"]]
  out
}

## 1. Latent matrix ---------
multi <- names(BLOCKS)[lengths(BLOCKS) > 1]
JOBS <- c(lapply(multi, identity), combn(names(BLOCKS), 2, simplify = FALSE))
sf2b_log("13: ", length(JOBS), " CFA models on ", SF2B$cores, " workers")
t1 <- Sys.time()
RES <- sf2b_par_lapply(JOBS, run_model, data = list(DAT = DAT, D = D, score_items = score_items, I8 = I8, BLOCKS = BLOCKS, MAIA_BOUND = MAIA_BOUND,
                                                     has_maia12 = has_maia12, block_syntax = block_syntax, model_syntax = model_syntax,
                                                     orient_bad = orient_bad, improper_of = improper_of, improper_list = improper_list))
secs_models <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
for (r in Filter(function(r) !r$ok, RES)) sf2b_log("FAILED: ", paste(r$blocks, collapse = " + "), ": ", r$error)
fits <- rbindlist(lapply(RES, function(r) c(list(model = paste(r$blocks, collapse = " + "), ok = r$ok, error = r$error, improper = r$improper %||% "",
                                                   maia12_std_residual = r$maia12_std_residual %||% NA_real_, seconds = round(r$secs, 1)),
                                              if (r$ok) r$fit)), fill = TRUE)
fwrite(fits, sf2b_out("13_model_fit.csv"))

VARS <- unlist(BLOCKS, use.names = FALSE)
blk <- setNames(rep(names(BLOCKS), lengths(BLOCKS)), VARS)
RL <- matrix(NA_real_, length(VARS), length(VARS), dimnames = list(VARS, VARS)); diag(RL) <- 1
within_all <- list()
for (r in Filter(function(r) r$ok, RES)) {
  C <- r$cor; v <- rownames(C)
  for (i in v) for (j in v) if (i != j) {
    same <- blk[[i]] == blk[[j]]
    if (same && length(r$blocks) == 2) { within_all[[length(within_all) + 1]] <- data.table(a = i, b = j, r = C[i, j]); next }
    RL[i, j] <- C[i, j]
  }
}
# Stability of within-instrument cells across the pair models that also estimate them.
within_dt <- rbindlist(within_all)
stab <- if (nrow(within_dt)) within_dt[match(a, VARS) < match(b, VARS), .(n_models = .N, min = min(r), max = max(r)), by = .(a, b)] else NULL
if (!is.null(stab)) { stab[, own_model := RL[cbind(a, b)]]; stab[, max_abs_dev := pmax(abs(max - own_model), abs(min - own_model))] }
dropped <- VARS[apply(RL, 1, function(x) any(is.na(x)))]
dropped_battery <- intersect(dropped, BATTERY)
if (length(dropped)) sf2b_log("dropped (cells missing): ", paste(dropped, collapse = ", "))
BAT <- setdiff(BATTERY, dropped)
EXT <- setdiff(EXTENSION, dropped)
fwrite(data.table(var = VARS, round(RL, 4)), sf2b_out("13_latent_matrix.csv"))
if (!is.null(stab)) fwrite(stab[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)], sf2b_out("13_within_cell_stability.csv"))

# The 8's row against 04's registered cells (same models, should agree).
LAT04 <- fread(sf2b_out("04_latent_validity.csv"))[form == "sf8" & sample == "all" & ri_free == FALSE]
chk04 <- merge(data.table(score = setdiff(VARS, "sf8"), r_13 = RL["sf8", setdiff(VARS, "sf8")]), LAT04[, .(score, r_04 = est)], by = "score")
chk04[, diff := round(r_13 - r_04, 3)]
dmaia <- tryCatch(fromJSON(sf2b_out("10_deviation_maia.json"))$sf8$h2_noticing$r, error = function(e) NA_real_)

Rm <- RL[BAT, BAT]
ev <- eigen(Rm, symmetric = TRUE, only.values = TRUE)$values
smoothed <- min(ev) <= 0
if (smoothed) Rm <- cor.smooth(Rm)

## 2. The general factor ---------
N <- nrow(DAT)
set.seed(SF2B$seed)
pa <- fa.parallel(Rm, n.obs = N, fm = "minres", fa = "fa", plot = FALSE)
K <- max(3L, as.integer(pa$nfact))
orient <- function(x, ref = "sf8") if (isTRUE(x[ref] < 0)) -x else x
sl_fit <- function(R, k) {
  s <- schmid(R, nfactors = k, fm = "minres", rotate = "oblimin", n.obs = N)
  g <- setNames(s$sl[, "g"], rownames(s$sl))
  flip <- if (isTRUE(g["sf8"] < 0) || (!"sf8" %in% names(g) && sum(g) < 0)) -1 else 1
  list(s = s, g = g * flip, flip = flip, k = k)
}
SL <- sl_fit(Rm, K)
SLm <- sl_fit(Rm, K - 1L); SLp <- sl_fit(Rm, K + 1L)
one <- fa(Rm, nfactors = 1, fm = "minres", n.obs = N)
g_one <- orient(setNames(as.numeric(one$loadings[, 1]), rownames(Rm)))
jb <- tryCatch(fa(Rm, nfactors = K + 1, fm = "minres", rotate = "bifactor", n.obs = N), error = function(e) e)
g_jb <- if (inherits(jb, "error")) setNames(rep(NA_real_, nrow(Rm)), rownames(Rm)) else {
  Lj <- unclass(jb$loadings); gc <- which.max(colSums(Lj^2)); orient(setNames(Lj[, gc], rownames(Rm))) }

# Dwyer extension of an SL solution to variables `ext` (correlations Roe with the original variables).
sl_extend <- function(sl, Roe) {
  P <- unclass(sl$s$oblique); Phi <- sl$s$phi; Sg <- P %*% Phi
  pe <- Roe[, rownames(P), drop = FALSE] %*% Sg %*% solve(crossprod(Sg))
  gam <- as.numeric(unclass(sl$s$gloading))
  drop(pe %*% gam) * sl$flip
}
# The general factor without the 8; the 8 extended.
B28 <- setdiff(BAT, "sf8")
SL28 <- sl_fit(Rm[B28, B28], K)
if (cor(SL28$g, SL$g[B28]) < 0) { SL28$g <- -SL28$g; SL28$flip <- -SL28$flip }
g8_ext <- unname(sl_extend(SL28, Rm["sf8", B28, drop = FALSE]))[1]

# Heywood check: the primary SL without BFI-2 Intellectual Curiosity (communality > 1 in
# every solution: its latent r with Creative Imagination is about .91).
heywood <- function(sl) { h <- sl$s$sl[, "h2"]; names(h)[h > 1] }
B_noint <- setdiff(BAT, "bfi_intellectual")
SLni <- sl_fit(Rm[B_noint, B_noint], K)
# The number of group factors, swept.
KSWEEP <- rbindlist(lapply(3:9, function(k) { x <- suppressWarnings(sl_fit(Rm, k))
  data.table(group_factors = k, max_communality = max(x$s$sl[, "h2"]), heywood = paste(heywood(x), collapse = " "),
             rms = x$s$rms, g_share_of_common = sum(x$g^2) / sum(x$s$sl[, "h2"]),
             g_sf8 = x$g[["sf8"]], g_ffmq_observing = x$g[["ffmq_observing"]], g_eq_decentering = x$g[["eq_decentering"]],
             g_ncs6 = x$g[["ncs6"]], g_bidr_sde = x$g[["bidr_sde"]], g_bidr_im = x$g[["bidr_im"]]) }))
fwrite(KSWEEP[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)], sf2b_out("13_k_sensitivity.csv"))

# Observed scores: first principal component (pairwise Pearson).
OBSD <- as.data.frame(S[, ..BAT])
Ro <- cor(OBSD, use = "pairwise.complete.obs")
pc <- principal(Ro, nfactors = 1, n.obs = N)
g_pc <- orient(setNames(as.numeric(pc$loadings[, 1]), BAT))

labels <- c(setNames(SF2B_PREDICTED$label, SF2B_PREDICTED$score), sf8 = "The 8", grp = "Practice group (current vs never)")
LOAD <- data.table(measure = BAT, label = labels[BAT], sl_g = SL$g[BAT], sl_g_one_fewer = SLm$g[BAT], sl_g_one_more = SLp$g[BAT],
                   one_factor = g_one[BAT], jennrich_bentler_g = g_jb[BAT],
                   sl_g_without_8 = c(SL28$g, sf8 = g8_ext)[BAT], sl_g_without_bfi_intellectual = SLni$g[BAT], observed_pc1 = g_pc[BAT],
                   sl_h2 = SL$s$sl[BAT, "h2"], sl_g_share_of_common = SL$s$sl[BAT, "p2"])
LOAD <- LOAD[order(-sl_g)]
fwrite(LOAD[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)], sf2b_out("13_general_loadings.csv"))
fwrite(data.table(measure = rownames(SL$s$sl), round(unclass(SL$s$sl), 3)), sf2b_out("13_schmid_leiman.csv"))
agree <- cor(LOAD[, .(sl_g, sl_g_one_fewer, sl_g_one_more, one_factor, jennrich_bentler_g, sl_g_without_8, sl_g_without_bfi_intellectual, observed_pc1)], use = "pairwise")

## 3. Net of G ---------
# Three definitions of G: the primary SL; SL without the 8 (the 8 and the extension
# variables extended onto it, so the 8 does not help define what is partialled from it);
# and the one-factor model (G absorbs the cluster variance too: an upper bound on what
# "general" can remove).
ext1 <- function(Roe, l) drop(Roe[, names(l), drop = FALSE] %*% l) / sum(l^2)
GDEF <- list(
  sl = c(SL$g, if (length(EXT)) sl_extend(SL, RL[EXT, BAT, drop = FALSE])),
  sl_without_8 = c(SL28$g, sf8 = g8_ext, if (length(EXT)) sl_extend(SL28, RL[EXT, B28, drop = FALSE])),
  one_factor = c(g_one, if (length(EXT)) ext1(RL[EXT, BAT, drop = FALSE], g_one)))
lam <- SL$g; ext_g <- GDEF$sl[EXT]
net <- rbindlist(lapply(names(GDEF), function(gd) { la <- GDEF[[gd]]
  rbindlist(lapply(setdiff(c(BAT, EXT), "sf8"), function(x) {
    r <- RL["sf8", x]; l8 <- la[["sf8"]]; lx <- la[[x]]
    data.table(g_definition = gd, measure = x, label = labels[[x]], r_with_8 = r, g_loading_8 = l8, g_loading = lx, accounted_by_G = l8 * lx,
               share_accounted = (l8 * lx) / r, partial_r_net_of_G = (r - l8 * lx) / sqrt((1 - l8^2) * (1 - lx^2)), extension = x %in% EXT)
  })) }))
net[, focus := measure %in% c(FOCUS, EXT)]
fwrite(net[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)][order(g_definition, -focus, -abs(r_with_8))], sf2b_out("13_net_of_general.csv"))
# Practice group, 07's metric, from the correlations.
d_from <- function(r8g, r8G = 0, rGg = 0) {
  b <- solve(matrix(c(1, rGg, rGg, 1), 2), c(r8G, r8g))          # standardized: G, grp
  r2 <- sum(b * c(r8G, r8g))
  coef <- b[2] / sqrt(P_CUR * (1 - P_CUR))
  # in SDs of the 8's residual given G and group (07's metric), and in the 8's within-group SD (sf2b_kg's)
  c(d_residual_sd = coef / sqrt(1 - r2), d_8_within_sd = coef / sqrt(1 - r8g^2))
}
grp_net <- if ("grp" %in% EXT) list(
  d_raw = d_from(RL["sf8", "grp"])[["d_residual_sd"]], r_8_grp = RL["sf8", "grp"], p_current = P_CUR,
  d_07_latent = tryCatch(fread(sf2b_out("07_latent_difference.csv"))[form == "sf8", latent_d], error = function(e) NA_real_),
  net = rbindlist(lapply(names(GDEF), function(gd) data.table(g_definition = gd, g_loading_8 = GDEF[[gd]][["sf8"]], g_loading_grp = GDEF[[gd]][["grp"]],
                                                              d_net_of_G = d_from(RL["sf8", "grp"], GDEF[[gd]][["sf8"]], GDEF[[gd]][["grp"]])[["d_residual_sd"]],
                                                              d_net_of_G_in_8_within_sd = d_from(RL["sf8", "grp"], GDEF[[gd]][["sf8"]], GDEF[[gd]][["grp"]])[["d_8_within_sd"]])))) else NULL
if (!is.null(grp_net)) fwrite(grp_net$net[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)], sf2b_out("13_practice_group_net.csv"))

## Observed-score check, with CIs ---------
pc28 <- principal(as.data.frame(scale(OBSD[, B28])), nfactors = 1, scores = TRUE, missing = FALSE)
Gobs <- as.numeric(pc28$scores[, 1])
if (cor(Gobs, S$sf8, use = "complete.obs") < 0) Gobs <- -Gobs
n_g <- sum(!is.na(Gobs))
fisher <- function(r, n, k = 1) { se <- 1 / sqrt(n - 3 - k); tanh(atanh(r) + c(-1, 1) * qnorm(.975) * se) }
pcor <- function(x, y, z, rank = FALSE) {
  ok <- complete.cases(x, y, z); x <- x[ok]; y <- y[ok]; z <- z[ok]
  if (rank) { x <- rank(x); y <- rank(y); z <- rank(z) }
  r <- cor(resid(lm(x ~ z)), resid(lm(y ~ z))); ci <- fisher(r, sum(ok))
  data.table(n = sum(ok), r = r, lower = ci[1], upper = ci[2], p = 2 * pt(-abs(r * sqrt((sum(ok) - 3) / (1 - r^2))), sum(ok) - 3))
}
g_cur <- DAT$grp
keep <- !is.na(Gobs)
obs_check <- list(
  n_G = n_g, g_pc1_variance = pc28$values[1] / length(B28),
  known_groups = rbind(cbind(analysis = "observed 8, same cases, unadjusted", sf2b_kg(fifelse(keep, S$sf8, NA_real_), g_cur)),
                       cbind(analysis = "observed 8, adjusted for G_obs", sf2b_kg(S$sf8, g_cur, data.table(G = Gobs)))),
  pause = rbind(cbind(analysis = "Spearman, same cases", as.data.table(sf2b_r_ci(fifelse(keep, S$sf8, NA_real_), S$crit_pause_days, "spearman"))[, .(n, r, lower, upper, p)]),
                cbind(analysis = "partial Spearman, net of G_obs", pcor(S$sf8, S$crit_pause_days, Gobs, rank = TRUE))),
  focus = rbindlist(lapply(FOCUS, function(x) cbind(measure = x, rbind(
    cbind(analysis = "Pearson, same cases", as.data.table(sf2b_r_ci(fifelse(keep, S$sf8, NA_real_), S[[x]]))[, .(n, r, lower, upper, p)]),
    cbind(analysis = "partial Pearson, net of G_obs", pcor(S$sf8, S[[x]], Gobs)))))))
fwrite(obs_check$known_groups[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 4) else x)], sf2b_out("13_observed_known_groups.csv"))
fwrite(rbind(cbind(measure = "crit_pause_days", obs_check$pause), obs_check$focus)[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 4) else x)],
       sf2b_out("13_observed_partial.csv"))

## Report ---------
out <- list(
  analysis = "13: a general self-description factor across the battery. EXPLORATORY, NOT REGISTERED (osf.io/d4a5u does not include it); run after the confirmatory run, with its results seen.",
  mode = SF2B$mode, tag = SF2B$tag, written = format(Sys.time(), "%Y-%m-%d %H:%M"),
  software = list(R = R.version.string, lavaan = as.character(packageVersion("lavaan")), psych = as.character(packageVersion("psych")),
                  GPArotation = as.character(packageVersion("GPArotation"))),
  n = N, choices = CHOICES,
  matrix = list(models = length(JOBS), failed = lapply(Filter(function(r) !r$ok, RES), function(r) list(model = paste(r$blocks, collapse = " + "), error = r$error)),
                improper = fits[nzchar(improper), .(model, improper)], dropped = dropped, battery = BAT,
                maia12_bound_active_in = fits[!is.na(maia12_std_residual) & maia12_std_residual < .1001, model],
                maia12_std_residual_range = range(fits$maia12_std_residual, na.rm = TRUE),
                smallest_eigenvalue = min(ev), smoothed = smoothed, check_against_04 = chk04, noticing_10_deviation = dmaia,
                within_cell_max_abs_dev = if (!is.null(stab)) max(stab$max_abs_dev) else NA_real_),
  parallel_analysis = list(nfact = pa$nfact, group_factors_used = K, fa_eigenvalues = round(pa$fa.values, 3)),
  schmid_leiman = list(k = K, omega_h = NULL, g_variance_share_of_common = sum(SL$g^2) / sum(SL$s$sl[, "h2"]),
                       g_variance_share_of_total = mean(SL$g^2), second_order_loadings = as.numeric(unclass(SL$s$gloading)),
                       fit_rms = SL$s$rms),
  one_factor = list(variance_share = mean(g_one^2), rms = one$rms,
                    note = "psych's chi-square-based indices (TLI, RMSEA) are not reported: they are not valid for a latent matrix assembled from separate models"),
  heywood = list(primary = heywood(SL), without_bfi_intellectual = heywood(SLni),
                 note = "communality above 1 (ultra-Heywood) for BFI-2 Intellectual Curiosity in every SL solution, k = 3-9; its latent r with Creative Imagination is about .91. The SL without it is reported as a check."),
  k_sensitivity = KSWEEP,
  jennrich_bentler = if (inherits(jb, "error")) list(error = conditionMessage(jb)) else list(factors = K + 1, rms = jb$rms),
  agreement_between_methods = round(agree, 3),
  loadings = LOAD, the_8 = list(sl_g = lam[["sf8"]], g_without_8_extended = g8_ext, bidr_sde = lam[["bidr_sde"]], bidr_im = lam[["bidr_im"]]),
  net_of_G = net, practice_group_net = grp_net, observed_check = obs_check,
  run_seconds = list(models = round(secs_models), total = round(as.numeric(difftime(Sys.time(), T0, units = "secs")))), cores = SF2B$cores)
# omega-hierarchical of the battery as a whole (descriptive only: the battery is not a scale).
om <- tryCatch(omega(Rm, nfactors = K, fm = "minres", n.obs = N, flip = FALSE, plot = FALSE), error = function(e) NULL)
out$schmid_leiman$omega_h <- if (!is.null(om)) om$omega_h else NA_real_
write_json(out, sf2b_out("13_exploratory_general_factor.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
sf2b_log("13 done in ", out$run_seconds$total, " s; outputs ", sf2b_out("13_*"))
print(LOAD[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])
print(net[focus == TRUE & g_definition == "sl", lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
print(grp_net); print(KSWEEP); print(obs_check$known_groups[, .(analysis, n_group, n_reference, d = round(d, 3), d_lower = round(d_lower, 3), d_upper = round(d_upper, 3))])
print(obs_check$pause); print(chk04)
