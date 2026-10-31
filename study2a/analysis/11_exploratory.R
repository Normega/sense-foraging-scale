# 11: further exploratory analyses (Section 5, "Exploratory analysis"); none is a
# decision rule.
#   - the self-rating's correlation with each factor score
#   - each reworded trial statement's correlation with its Study 1 wording, and
#     the trial's Access + Valuing score's correlations with the Access and
#     Valuing factor scores
#   - current practitioners vs never-practitioners on each statement and factor,
#     standardized (d with the pooled SD), as a guide to Study 2b
# Position effects are in 10.
# PREREG: the trial's Access + Valuing score is the mean of its 10 statements in
# Study 1 wording (sf6_35, sf6_36, sf6_44 recoded), of those answered; the Access
# and Valuing factors are the primary model's factors 05 matched to them.

source("study2a/analysis/config.R")
source("study2a/analysis/sf2a_models.R")
A <- readRDS(sf2a_derived("sf2a"))
def <- A$definition; k <- SF2A$k_primary; M <- A$meta
P <- readRDS(sf2a_derived(sprintf("primary_k%02d", k)))
FS <- P$scores; K <- A$items_keyed; END <- A$end
fz <- function(r, n) tanh(atanh(r) + c(-1, 1) * 1.96 / sqrt(n - 3))
corr_row <- function(a, b, label) {
  ok <- complete.cases(a, b); n <- sum(ok)
  if (n < 10) return(data.table(what = label, n = n))
  r <- cor(a[ok], b[ok]); ci <- fz(r, n)
  data.table(what = label, n = n, r = round(r, 3), ci_lo = round(ci[1], 3), ci_hi = round(ci[2], 3))
}
lab <- function(f) paste0(f, " (", P$labels[[f]] %||% "", ")")

self <- rbindlist(lapply(P$fnames, function(f) corr_row(M$vig_self, FS[, f], paste("self-rating ~", lab(f)))))

reword <- rbindlist(lapply(def$end, function(s1) {
  i <- def$items[def$source[def$items] == s1]
  if (!length(i)) return(NULL)
  pc <- psych::polychoric(cbind(A$items_raw[, i], END[, s1]), correct = 0)$rho[1, 2]
  cbind(corr_row(A$items_raw[, i], END[, s1], paste(i, "~", s1, "(Study 1 wording)")), polychoric = round(pc, 3))
}), fill = TRUE)

av_items <- cbind(K[, def$av10[c("sf6_9", "sf6_36", "sf6_35", "sf6_44")]], END)
av10 <- rowMeans(av_items, na.rm = TRUE)
# Every factor, not only those matched to Access and Valuing: at k = 4 the two
# hypothesized factors' forward statements define one factor (matched to Forage
# by plurality) and their reverse-keyed statements another.
av <- rbindlist(lapply(P$fnames, function(f) corr_row(av10, FS[, f], paste("trial Access + Valuing ~", lab(f)))))
av_desc <- data.table(what = "trial Access + Valuing score", n = sum(!is.na(av10)), mean = round(mean(av10, na.rm = TRUE), 3),
                      sd = round(sd(av10, na.rm = TRUE), 3), alpha = round(psych::alpha(av_items, warnings = FALSE)$total$raw_alpha, 3))

g <- M$practice_group
dfun <- function(y, label) {
  a <- y[g %in% "current"]; b <- y[g %in% "never"]; a <- a[!is.na(a)]; b <- b[!is.na(b)]
  sp <- sqrt(((length(a) - 1) * var(a) + (length(b) - 1) * var(b)) / (length(a) + length(b) - 2))
  data.table(what = label, n_current = length(a), n_never = length(b), mean_current = round(mean(a), 3),
             mean_never = round(mean(b), 3), d = round((mean(a) - mean(b)) / sp, 3))
}
kg <- rbind(rbindlist(lapply(P$fnames, function(f) dfun(FS[, f], paste("factor", lab(f))))),
            dfun(av10, "trial Access + Valuing score"),
            rbindlist(lapply(def$items, function(i) dfun(K[, i], paste(i, def$source[i], def$factor[i])))))

fwrite(rbind(self, av, fill = TRUE), sf2a_out(sprintf("11_correlations_k%02d.csv", k)))
fwrite(reword, sf2a_out("11_reworded_vs_study1.csv"))
fwrite(av_desc, sf2a_out("11_trial_score.csv"))
fwrite(kg, sf2a_out(sprintf("11_practitioners_k%02d.csv", k)))
fwrite(M[, .N, by = practice_group], sf2a_out("11_practice_groups.csv"))
print(self); print(reword); print(av); print(av_desc); print(head(kg, length(P$fnames) + 1))
