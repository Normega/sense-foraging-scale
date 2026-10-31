# 03: did the two modes come across? (Section 5, step 2; decision V)
#
#   separation   a participant's mean rating of the four sensing vignettes minus
#                their mean of the five doing vignettes, for those who rated at
#                least three of each (computed in 01)
#   tells apart  separation >= 1 point
#   decision V   the distinction came across if >= 75% of participants with a
#                separation score tell the modes apart
#   midpoint     halfway between the sample's mean sensing and mean doing rating;
#                a vignette whose mean falls on the wrong side is flagged
# Nobody is excluded on their vignette ratings.
# PREREG: "the sample's mean sensing rating" is read as the mean, over
# participants, of each participant's mean sensing rating (and likewise doing).

source("study2a/analysis/config.R")
A <- readRDS(sf2a_derived("sf2a"))
def <- A$definition; M <- A$meta; V <- A$vignettes

has <- !is.na(M$separation)
n <- sum(has); k <- sum(M$tells_apart[has])
wilson <- function(k, n, z = 1.96) {
  p <- k / n; d <- 1 + z^2 / n; c0 <- (p + z^2 / (2 * n)) / d
  h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d; c(c0 - h, c0 + h)
}
ci <- wilson(k, n)
sens_mean <- mean(M$sensing_mean[has]); doing_mean <- mean(M$doing_mean[has])
midpoint <- (sens_mean + doing_mean) / 2
tt <- t.test(M$separation[has])
summary_tab <- data.table(
  n_with_separation = n, n_without = sum(!has),
  mean_sensing = round(sens_mean, 3), mean_doing = round(doing_mean, 3), midpoint = round(midpoint, 3),
  mean_separation = round(mean(M$separation[has]), 3), sd_separation = round(sd(M$separation[has]), 3),
  separation_ci_lo = round(tt$conf.int[1], 3), separation_ci_hi = round(tt$conf.int[2], 3),
  n_tell_apart = k, pct_tell_apart = round(100 * k / n, 1),
  pct_ci_lo = round(100 * ci[1], 1), pct_ci_hi = round(100 * ci[2], 1),
  threshold_pct = 100 * SF2A$vignette_share, came_across = k / n >= SF2A$vignette_share)
fwrite(summary_tab, sf2a_out("03_vignettes_decision.csv"))

per <- data.table(vignette = def$vignettes[1:9], kind = def$vkind[def$vignettes[1:9]],
                  n = colSums(!is.na(V[has, 1:9])), mean = round(colMeans(V[has, 1:9], na.rm = TRUE), 3),
                  sd = round(apply(V[has, 1:9], 2, sd, na.rm = TRUE), 3),
                  pna_or_missing = colSums(is.na(V[, 1:9])))
per[, wrong_side := (kind == "sensing" & mean <= midpoint) | (kind == "doing" & mean >= midpoint)]
fwrite(per, sf2a_out("03_vignettes_each.csv"))

dist <- data.table(separation = round(M$separation[has], 3))[, .N, by = separation][order(separation)]
fwrite(dist, sf2a_out("03_separation_distribution.csv"))
self <- data.table(vig_self = c(1:6, NA))[, n := vapply(vig_self, function(v) if (is.na(v)) sum(is.na(M$vig_self)) else sum(M$vig_self == v, na.rm = TRUE), 0L)]
fwrite(self, sf2a_out("03_self_rating.csv"))

sf2a_log("decision V: ", k, " of ", n, " (", summary_tab$pct_tell_apart, "%) tell the modes apart; came across = ",
         summary_tab$came_across, "; vignettes on the wrong side: ", paste(per[wrong_side == TRUE, vignette], collapse = " "))
print(summary_tab); print(per)
