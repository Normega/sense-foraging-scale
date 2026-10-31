# 18: exploratory, not preregistered (Norm, 2026-10-31): what do two more statements add to
# 17's final 8 (12's best 8 with U4 in place of sfc_x2)? Real mode only, after 17.
#   - every pair of statements added to the 8 (one per look-alike pair), and every single
#     one (a 9), by content omega from 12's G model
#   - for the 8, the best 9 and the best 10s: omega, agreement, alpha (overall and by vignette
#     group), standard error of measurement and the smallest detectable change (95%) of the
#     mean score, practitioner d, correlations with the trial's 10, IUS-12 and the self-rating

source("study2a/analysis/config.R")
source("study2a/analysis/sf2a_models.R")
stopifnot(SF2A$mode == "real")
A <- readRDS(sf2a_derived("sf2a"))
def <- A$definition; K <- A$items_keyed; M <- A$meta; items <- def$items; rv <- def$reversed
F12 <- readRDS(sf2a_derived("followup12"))
out18 <- function(f) sf2a_out(paste0("18_", f))
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
pool <- setdiff(items, c(final8, unlist(lapply(final8, partner))))
ok <- function(set) !any(set %in% unlist(lapply(set, partner)))
singles <- rbindlist(lapply(pool, function(i) data.table(add = i, omega_content = set_omega(c(final8, i))[["content"]])))
pairs <- rbindlist(combn(pool, 2, function(p) if (ok(c(final8, p)))
  data.table(add = paste(p, collapse = " "), omega_content = set_omega(c(final8, p))[["content"]]), simplify = FALSE))

describe <- function(set, label) {
  om <- set_omega(set); sc <- rowMeans(K[, set], na.rm = TRUE); a <- alpha_fast(K[, set]); sem <- sd(sc, na.rm = TRUE) * sqrt(1 - a)
  add <- setdiff(set, final8)
  data.table(set = label, n = length(set), n_reversed = sum(set %in% rv), n_names_mode = sum(set %in% def$names_mode),
             added = paste(def$source[add], collapse = " "), omega_content = om[["content"]], share_agreement = om[["agreement"]],
             alpha = a, alpha_told_apart = alpha_fast(K[ta, set]), alpha_not = alpha_fast(K[!ta, set]),
             sd = sd(sc, na.rm = TRUE), sem = sem, mdc95 = 1.96 * sqrt(2) * sem, d_practice = dval(sc),
             r_trial10 = cor(sc, trial10, use = "pair"), r_final8 = cor(sc, rowMeans(K[, final8], na.rm = TRUE), use = "pair"),
             r_ius = cor(sc, M$ius_total, use = "pair"), r_self = cor(sc, M$vig_self, use = "pair"),
             added_text = paste(def$text[add], collapse = " | "), ids = paste(set, collapse = " "))
}
top_pairs <- pairs[order(-omega_content)][1:5]
top_single <- singles[order(-omega_content)][1]
mode_free <- pairs[!vapply(strsplit(add, " "), function(p) any(p %in% def$names_mode), TRUE)][order(-omega_content)][1:3]
res <- rbindlist(c(list(describe(final8, "final 8")),
                   list(describe(c(final8, top_single$add), "best 9")),
                   lapply(seq_len(nrow(top_pairs)), function(i) describe(c(final8, strsplit(top_pairs$add[i], " ")[[1]]), paste("10, pair", i))),
                   lapply(seq_len(nrow(mode_free)), function(i) describe(c(final8, strsplit(mode_free$add[i], " ")[[1]]), paste("10, mode-free pair", i)))))
res <- unique(res, by = "ids")
fwrite(num3(res), out18("eight_vs_ten.csv"))
print(num3(res)[, !c("added_text", "ids")])
cat("\n"); for (i in seq_len(nrow(res))) cat(res$set[i], ":", res$added_text[i], "\n")
