# 07: intolerance of uncertainty (Section 5, step 6; expectation 6).
#
# Correlations of the IUS-12 total (sum of 12 items, 12-60) with the mean of
# Letting be's 7 forward statements and with the mean of its 7 reverse-keyed
# statements (recoded 7 - x), compared with a test for dependent correlations
# (Steiger, 1980: psych::r.test with r12, r13, r23, Williams' T2 as Steiger
# recommends). 95% CIs by Fisher's z.
# PREREG: "more strongly" is read as larger in absolute value (both are expected
# to be negative, higher scores meaning more letting be); participants with an
# IUS total and both means are used.

source("study2a/analysis/config.R")
source("study2a/analysis/sf2a_models.R")
A <- readRDS(sf2a_derived("sf2a"))
def <- A$definition; M <- A$meta; X <- A$items_raw

lb <- def$items[def$factor[def$items] == "Letting be"]
fwd <- setdiff(lb, def$reversed); rv <- intersect(lb, def$reversed)
d <- data.table(ius = M$ius_total, forward = sf2a_score(X, fwd), reversed = sf2a_score(X, rv, rv))
d <- d[complete.cases(d)]
fz <- function(r, n) tanh(atanh(r) + c(-1, 1) * 1.96 / sqrt(n - 3))
r_f <- cor(d$ius, d$forward); r_r <- cor(d$ius, d$reversed); r_fr <- cor(d$forward, d$reversed); n <- nrow(d)
st <- psych::r.test(n = n, r12 = r_f, r13 = r_r, r23 = r_fr)
res <- data.table(n = n, r_ius_forward = round(r_f, 3), ci_forward = paste(round(fz(r_f, n), 3), collapse = " to "),
                  r_ius_reversed = round(r_r, 3), ci_reversed = paste(round(fz(r_r, n), 3), collapse = " to "),
                  r_forward_reversed = round(r_fr, 3), t2 = round(st$t, 3), p = signif(st$p, 3),
                  reversed_stronger = abs(r_r) > abs(r_f),
                  ius_mean = round(mean(d$ius), 2), ius_sd = round(sd(d$ius), 2),
                  ius_alpha = round(psych::alpha(A$ius, warnings = FALSE)$total$raw_alpha, 3))
fwrite(res, sf2a_out("07_ius.csv"))
print(res)
