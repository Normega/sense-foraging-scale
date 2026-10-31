# 17: exploratory, not preregistered (Norm, 2026-10-31): 12's best 8 with a reverse-keyed
# statement that does not name a mode. Which statement should it replace? Real mode only,
# after 12 and 15.
#   - candidates: every reverse-keyed statement that names no mode, in place of each of the
#     8 (one statement per look-alike pair)
#   - per variant: content omega and agreement share (12's G model), alpha overall and in
#     the people who did and did not tell the modes apart, the practitioner difference, and
#     correlations with the best 8, the trial's 10, IUS-12 and the self-rating

source("study2a/analysis/config.R")
source("study2a/analysis/sf2a_models.R")
stopifnot(SF2A$mode == "real")
A <- readRDS(sf2a_derived("sf2a"))
def <- A$definition; K <- A$items_keyed; M <- A$meta; items <- def$items; rv <- def$reversed
F12 <- readRDS(sf2a_derived("followup12"))
out17 <- function(f) sf2a_out(paste0("17_", f))
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

b8 <- fread(sf2a_out("12_strongest_scale.csv"))[n == 8 & min_reversed == 0, strsplit(ids, " ")[[1]]]
sc_b8 <- rowMeans(K[, b8], na.rm = TRUE)
trial10 <- rowMeans(cbind(K, A$end)[, unname(def$av10)], na.rm = TRUE)
ta <- M$tells_apart %in% TRUE
cand <- setdiff(rv[!rv %in% def$names_mode], b8)

describe <- function(set, label, out = NA_character_, inn = NA_character_) {
  om <- set_omega(set); sc <- rowMeans(K[, set], na.rm = TRUE)
  data.table(variant = label, drop = out, add = inn, drop_text = if (is.na(out)) NA else def$text[out],
             add_text = if (is.na(inn)) NA else def$text[inn],
             n_reversed = sum(set %in% rv), n_names_mode = sum(set %in% def$names_mode),
             omega_content = om[["content"]], share_agreement = om[["agreement"]],
             alpha = alpha_fast(K[, set]), alpha_told_apart = alpha_fast(K[ta, set]), alpha_not = alpha_fast(K[!ta, set]),
             d_practice = dval(sc), r_best8 = cor(sc, sc_b8, use = "pair"), r_trial10 = cor(sc, trial10, use = "pair"),
             r_ius = cor(sc, M$ius_total, use = "pair"), r_self = cor(sc, M$vig_self, use = "pair"),
             ids = paste(set, collapse = " "))
}
res <- rbindlist(c(list(describe(b8, "best 8 (12)")),
  unlist(lapply(cand, function(i) lapply(b8, function(o) { set <- c(setdiff(b8, o), i)
    if (any(set %in% unlist(lapply(set, partner)))) return(NULL)
    describe(set, paste(def$source[o], "->", def$source[i]), o, i) })), recursive = FALSE)))
res[, `:=`(drop_source = def$source[drop], add_source = def$source[add])]
fwrite(num3(res[order(-omega_content)]), out17("swap_variants.csv"))

# the candidates themselves: keyed G loading, agreement and wording loadings, item-rest r
cand_tab <- data.table(item = cand, source = def$source[cand], G_keyed = -L[cand, "G"], RI = L[cand, "RI"], W = L[cand, "W"],
                       text = def$text[cand])
fwrite(num3(cand_tab[order(-G_keyed)]), out17("candidates.csv"))
print(num3(cand_tab[order(-G_keyed)])[, !"text"])
print(num3(res[order(-omega_content)])[1:15, .(variant, n_reversed, n_names_mode, omega_content, share_agreement, alpha, alpha_told_apart,
                                                alpha_not, d_practice, r_best8, r_trial10, r_ius)])
