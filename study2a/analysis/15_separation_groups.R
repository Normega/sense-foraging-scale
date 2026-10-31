# 15: exploratory, not preregistered (Norm's question, 2026-10-31): do the people who did
# best on the vignettes (the two modes "quiz") answer the statements differently: higher
# internal consistency, a stronger or cleaner factor structure, less agreement? Real mode
# only, after 12 and 14.
#
# Groups: the registered cut (separation >= 1, "tells the modes apart", 346 vs 154) and the
# top vs bottom third of separation (equal n, so eigenvalues and congruence compare fairly).
#   1. alpha of the main sets per group, with a bootstrap CI of the difference
#   2. the keyed polychoric matrix per group: eigenvalues, and how the forward and the
#      reverse-keyed statements correlate
#   3. plain EFA, k = 4 (minres, oblimin) per group; Tucker congruence between the groups,
#      against a baseline of random splits of the same people (what sampling alone gives)
#   4. 12's G + RI + W model per group: agreement variance, the general factor's loadings,
#      content omega of 12's best 8 and the trial's 10 in the scale

source("study2a/analysis/config.R")
source("study2a/analysis/sf2a_models.R")
stopifnot(SF2A$mode == "real")
A <- readRDS(sf2a_derived("sf2a"))
def <- A$definition; X <- A$items_raw; K <- A$items_keyed; M <- A$meta; items <- def$items; rv <- def$reversed
out15 <- function(f) sf2a_out(paste0("15_", f))
set.seed(SF2A$seed + 15)

sep <- M$separation
terc <- cut(sep, quantile(sep, c(0, 1/3, 2/3, 1), na.rm = TRUE), include.lowest = TRUE, labels = c("bottom", "middle", "top"))
contrasts <- list(registered = list(hi = which(M$tells_apart %in% TRUE), lo = which(M$tells_apart %in% FALSE)),
                  thirds = list(hi = which(terc == "top"), lo = which(terc == "bottom")))

# ---- the sets ----------------------------------------------------------------------------------
PL <- fread(sf2a_out("05_primary_loadings_keyed_k04.csv"))
b8 <- fread(sf2a_out("12_strongest_scale.csv"))[n == 8 & min_reversed == 0, strsplit(ids, " ")[[1]]]
und <- strsplit(fread(sf2a_out("14_chosen_summary.csv"))[factor == "Understanding", statements], " ")[[1]]
und <- items[match(und, def$source[items])]
sets <- c(list(`all 58` = items, `best 8 (12)` = b8,
               `trial's 10 in the scale` = c("sf2a_01", "sf2a_02", "sf2a_03", "sf2a_04", "sf2a_05", "sf2a_23", "sf2a_24", "sf2a_25", "sf2a_26", "sf2a_27")),
          setNames(lapply(paste0("F", 1:4), function(f) PL[primary_factor == f & primary >= .40, item]),
                   paste0(c("F1 pull of the doing mode", "F2 modes the same", "F3 broad", "F4 need to explain or act"), " (05)")),
          list(`Understanding 3 (14)` = und, `statements naming a mode` = intersect(items, def$names_mode)))
KE <- cbind(K, A$end)                                         # the trial's 10 in Study 1 wording
sets_ext <- list(`trial's 10, Study 1 wording` = unname(def$av10))
alpha_fast <- function(Z) { S <- cov(Z, use = "pairwise.complete.obs"); k <- ncol(S); k / (k - 1) * (1 - sum(diag(S)) / sum(S)) }

# ---- 1. alpha -----------------------------------------------------------------------------------
alpha_tab <- rbindlist(lapply(names(contrasts), function(cn) { g <- contrasts[[cn]]
  rbindlist(lapply(c(names(sets), names(sets_ext)), function(s) {
    Z <- if (s %in% names(sets)) K[, sets[[s]], drop = FALSE] else KE[, sets_ext[[s]], drop = FALSE]
    ah <- alpha_fast(Z[g$hi, ]); al <- alpha_fast(Z[g$lo, ])
    bd <- replicate(1000, alpha_fast(Z[sample(g$hi, replace = TRUE), ]) - alpha_fast(Z[sample(g$lo, replace = TRUE), ]))
    sc <- rowMeans(Z, na.rm = TRUE)
    data.table(contrast = cn, set = s, n_items = ncol(Z), n_hi = length(g$hi), n_lo = length(g$lo),
               alpha_hi = ah, alpha_lo = al, diff = ah - al, diff_lo = quantile(bd, .025), diff_hi = quantile(bd, .975),
               mean_hi = mean(sc[g$hi], na.rm = TRUE), mean_lo = mean(sc[g$lo], na.rm = TRUE),
               sd_hi = sd(sc[g$hi], na.rm = TRUE), sd_lo = sd(sc[g$lo], na.rm = TRUE))
  }))
}))
num3 <- function(T) T[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 3) else v)]
fwrite(num3(alpha_tab), out15("alpha.csv"))

# ---- 2. eigenvalues, forward vs reversed ----------------------------------------------------------
poly <- function(rows) suppressWarnings(psych::polychoric(K[rows, items], correct = 0)$rho)
fw <- setdiff(items, rv)
struct_tab <- rbindlist(lapply(names(contrasts), function(cn) rbindlist(lapply(c("hi", "lo"), function(side) {
  rows <- contrasts[[cn]][[side]]; R <- poly(rows); ev <- eigen(R, only.values = TRUE)$values
  mr <- function(a, b) { r <- R[a, b]; if (identical(a, b)) mean(r[upper.tri(r)]) else mean(r) }
  data.table(contrast = cn, group = side, n = length(rows), ev1 = ev[1], ev1_pct = 100 * ev[1] / length(items), ev2 = ev[2],
             ev3 = ev[3], ev4 = ev[4], ev1_over_ev2 = ev[1] / ev[2],
             r_fwd_fwd = mr(fw, fw), r_rev_rev = mr(rv, rv), r_fwd_rev = mr(fw, rv),
             r_composites = cor(rowMeans(K[rows, fw], na.rm = TRUE), rowMeans(K[rows, rv], na.rm = TRUE), use = "pair"))
}))))
fwrite(num3(struct_tab), out15("structure.csv"))

# ---- 3. EFA congruence between groups vs random splits ---------------------------------------------
efa4 <- function(rows) unclass(suppressWarnings(psych::fa(poly(rows), nfactors = 4, n.obs = length(rows), fm = "minres",
                                                          rotate = "oblimin"))$loadings)
perms <- function(k) if (k == 1) list(1) else do.call(c, lapply(seq_len(k), function(i) lapply(perms(k - 1), function(p) c(i, setdiff(seq_len(k), i)[p]))))
match_cong <- function(L1, L2) {
  C <- abs(psych::factor.congruence(L1, L2))
  best <- perms(ncol(C))[[which.max(vapply(perms(ncol(C)), function(p) mean(C[cbind(seq_len(ncol(C)), p)]), 0))]]
  C[cbind(seq_len(ncol(C)), best)]
}
library(parallel)
cl <- makeCluster(min(12, detectCores() - 2))
clusterExport(cl, c("K", "items", "poly", "efa4", "perms", "match_cong"), envir = environment())
cong_tab <- rbindlist(lapply(names(contrasts), function(cn) {
  g <- contrasts[[cn]]
  obs <- match_cong(efa4(g$hi), efa4(g$lo))
  pool <- c(g$hi, g$lo); nlo <- length(g$lo)
  base <- parSapply(cl, seq_len(20), function(i, pool, nlo) { s <- sample(pool, nlo); mean(match_cong(efa4(setdiff(pool, s)), efa4(s))) },
                    pool = pool, nlo = nlo)
  data.table(contrast = cn, congruence_by_factor = paste(sprintf("%.2f", sort(obs, decreasing = TRUE)), collapse = " "),
             mean_congruence_groups = mean(obs), mean_congruence_random_splits = mean(base),
             random_split_lo = quantile(base, .05), random_split_hi = quantile(base, .95))
}))
stopCluster(cl)
fwrite(num3(cong_tab), out15("efa_congruence.csv"))

# ---- 4. G + RI + W per group ---------------------------------------------------------------------
g_syntax <- function(items, wording = TRUE) {
  rw <- intersect(items, rv)
  c(paste0("G =~ NA*", items[1], " + ", paste(items[-1], collapse = " + ")), "G ~~ 1*G",
    paste0("RI =~ ", paste0("1*", items, collapse = " + ")), "RI ~~ RI", "G ~~ 0*RI",
    if (wording) c(paste0("W =~ NA*", rw[1], " + ", paste(rw[-1], collapse = " + ")), "W ~~ 1*W", "W ~~ 0*G + 0*RI"))
}
jobs <- unlist(lapply(names(contrasts), function(cn) lapply(c("hi", "lo"), function(s) list(cn = cn, side = s, rows = contrasts[[cn]][[s]]))),
               recursive = FALSE)
cl <- makeCluster(length(jobs))
clusterCall(cl, setwd, getwd())
clusterEvalQ(cl, { Sys.setenv(SF2A_MODE = "real"); source("study2a/analysis/config.R"); source("study2a/analysis/sf2a_models.R") })
clusterExport(cl, c("g_syntax", "rv", "items", "X"), envir = environment())
gf <- parLapply(cl, jobs, function(j) c(j[c("cn", "side")], list(fit = tryCatch(sf2a_sem(g_syntax(items), X[j$rows, ], ordered = items),
                                                                                error = function(e) e))))
stopCluster(cl)
set_omega <- function(L, set) { s <- ifelse(set %in% rv, -1, 1); l <- L[set, , drop = FALSE]
  g <- sum(s * l[, "G"])^2; a <- sum(s * l[, "RI"])^2; w <- sum(s * l[, "W"])^2; e <- sum(1 - rowSums(l^2))
  c(content = g / (g + a + w + e), agreement = a / (g + a + w + e)) }
g_tab <- rbindlist(lapply(gf, function(z) {
  if (sf2a_failed(z$fit)) return(data.table(contrast = z$cn, group = z$side, note = sf2a_msg(z$fit)))
  L <- lavInspect(z$fit, "std")$lambda[items, c("G", "RI", "W")]; s <- ifelse(items %in% rv, -1, 1)
  if (sum(s * L[, "G"]) < 0) L[, "G"] <- -L[, "G"]
  kg <- s * L[, "G"]; o8 <- set_omega(L, b8); o10 <- set_omega(L, sets$`trial's 10 in the scale`)
  data.table(contrast = z$cn, group = z$side, n = lavInspect(z$fit, "ntotal"), improper = sf2a_improper(z$fit),
             agreement_pct_item_variance = 100 * mean(L[, "RI"]^2), G_mean_keyed_forward = mean(kg[!items %in% rv]),
             G_mean_keyed_reversed = mean(kg[items %in% rv]), G_mean_keyed_names_mode = mean(kg[items %in% def$names_mode]),
             n_G_ge_40 = sum(kg >= .40), W_mean = mean(abs(L[rv, "W"])),
             best8_content_omega = o8[["content"]], best8_agreement = o8[["agreement"]],
             trial10_content_omega = o10[["content"]], trial10_agreement = o10[["agreement"]],
             cfi = fitMeasures(z$fit, "cfi.scaled"), rmsea = fitMeasures(z$fit, "rmsea.scaled"))
}), fill = TRUE)
fwrite(num3(g_tab), out15("G_model.csv"))

print(num3(alpha_tab)[, .(contrast, set, n_items, alpha_hi, alpha_lo, diff, diff_lo, diff_hi, sd_hi, sd_lo, mean_hi, mean_lo)])
print(num3(struct_tab)); print(num3(cong_tab)); print(num3(g_tab))
