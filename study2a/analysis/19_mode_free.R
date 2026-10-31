# 19: exploratory, not preregistered (Norm, 2026-10-31): can the final 8 / 10 lose their
# remaining statements that name a mode (sf6_36, sf6_44)? Real mode only, after 18.
#   - the best sets of 8 and 10 from statements that name no mode, by content omega (12's G
#     model, greedy, one per look-alike pair), with 0, 1, 2 or 3 reverse-keyed statements
#   - hand variants: 17's final 8 with sf6_36 and/or sf6_44 replaced by the strongest
#     mode-free reverse-keyed statements, or by forward ones
#   - each described as in 18

source("study2a/analysis/config.R")
source("study2a/analysis/sf2a_models.R")
stopifnot(SF2A$mode == "real")
A <- readRDS(sf2a_derived("sf2a"))
def <- A$definition; K <- A$items_keyed; M <- A$meta; items <- def$items; rv <- def$reversed
F12 <- readRDS(sf2a_derived("followup12"))
out19 <- function(f) sf2a_out(paste0("19_", f))
num3 <- function(T) T[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 3) else v)]

L <- lavInspect(F12$G, "std")$lambda[items, c("G", "RI", "W")]; s_all <- ifelse(items %in% rv, -1, 1)
if (sum(s_all * L[, "G"]) < 0) L[, "G"] <- -L[, "G"]
set_omega <- function(set) { s <- ifelse(set %in% rv, -1, 1); l <- L[set, , drop = FALSE]
  g <- sum(s * l[, "G"])^2; a <- sum(s * l[, "RI"])^2; w <- sum(s * l[, "W"])^2; e <- sum(1 - rowSums(l^2))
  c(content = g / (g + a + w + e), agreement = a / (g + a + w + e)) }
alpha_fast <- function(Z) { S <- cov(Z, use = "pairwise.complete.obs"); k <- ncol(S); k / (k - 1) * (1 - sum(diag(S)) / sum(S)) }
partner <- function(i) { p <- def$pairs[def$pairs[, "a"] == i | def$pairs[, "b"] == i, , drop = FALSE]; setdiff(c(p), i) }
grp <- M$practice_group
dval <- function(y) { a <- y[grp %in% "current"]; b <- y[grp %in% "never"]
  (mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE)) / sqrt((var(a, na.rm = TRUE) + var(b, na.rm = TRUE)) / 2) }
trial10 <- rowMeans(cbind(K, A$end)[, unname(def$av10)], na.rm = TRUE)
ta <- M$tells_apart %in% TRUE
final8 <- strsplit(fread(sf2a_out("17_swap_variants.csv"))[variant == "sfc_x2 -> U4", ids], " ")[[1]]
final10 <- strsplit(fread(sf2a_out("18_eight_vs_ten.csv"))[set == "10, mode-free pair 1", ids], " ")[[1]]
sc8 <- rowMeans(K[, final8], na.rm = TRUE)
id <- function(src) items[match(src, def$source[items])]

# greedy content omega with an exact number of reverse-keyed statements
greedy <- function(n, n_rev, pool) {
  sel <- character()
  while (length(sel) < n) {
    left <- setdiff(pool, c(sel, unlist(lapply(sel, partner))))
    nr <- sum(sel %in% rv)
    if (nr >= n_rev) left <- setdiff(left, rv) else if (n - length(sel) <= n_rev - nr) left <- intersect(left, rv)
    if (!length(left)) break
    sel <- c(sel, names(which.max(vapply(left, function(i) set_omega(c(sel, i))[["content"]], 0))))
  }
  sel
}
describe <- function(set, label) {
  om <- set_omega(set); sc <- rowMeans(K[, set], na.rm = TRUE); a <- alpha_fast(K[, set]); sem <- sd(sc, na.rm = TRUE) * sqrt(1 - a)
  data.table(set = label, n = length(set), n_reversed = sum(set %in% rv), n_names_mode = sum(set %in% def$names_mode),
             omega_content = om[["content"]], share_agreement = om[["agreement"]],
             alpha = a, alpha_told_apart = alpha_fast(K[ta, set]), alpha_not = alpha_fast(K[!ta, set]),
             mdc95 = 1.96 * sqrt(2) * sem, d_practice = dval(sc), r_final8 = cor(sc, sc8, use = "pair"),
             r_trial10 = cor(sc, trial10, use = "pair"), r_ius = cor(sc, M$ius_total, use = "pair"),
             statements = paste(def$source[set], collapse = " "), ids = paste(set, collapse = " "))
}
free <- setdiff(items, def$names_mode)
res <- rbindlist(c(
  list(describe(final8, "final 8 (17)"), describe(final10, "final 10 (18)")),
  unlist(lapply(c(8, 10), function(n) lapply(0:3, function(r) describe(greedy(n, r, free), sprintf("mode-free best %d, %d reversed", n, r)))),
         recursive = FALSE),
  list(describe(c(setdiff(final8, id("sf6_36")), id("N6")), "final 8, sf6_36 -> N6"),
       describe(c(setdiff(final8, id("sf6_44")), id("N6")), "final 8, sf6_44 -> N6"),
       describe(c(setdiff(final8, id(c("sf6_36", "sf6_44"))), id(c("N6", "sfc_31"))), "final 8, both -> N6 + sfc_31"),
       describe(c(setdiff(final8, id(c("sf6_36", "sf6_44"))), id(c("N6", "sf6_43"))), "final 8, both -> N6 + sf6_43"),
       describe(c(setdiff(final8, id(c("sf6_36", "sf6_44"))), id(c("A9", "sfc_24"))), "final 8, both -> A9 + sfc_24 (forward)"))))
res <- unique(res, by = "ids")
fwrite(num3(res), out19("mode_free.csv"))
print(num3(res)[, !"ids"])
