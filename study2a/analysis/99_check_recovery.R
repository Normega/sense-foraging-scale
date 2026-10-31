# 99: simulated mode only. Does the pipeline recover the model the data came from?
# Tucker congruence of each recovered factor with its true factor (keyed, matched
# by the target), the random-intercept variance against the true .10, and the
# vignette separation against the true separations.

source("study2a/analysis/config.R")
source("study2a/analysis/sf2a_models.R")
stopifnot(SF2A$mode == "simulated", SF2A$k_primary == 5)
truth <- readRDS(sf2a_derived("simulation_truth"))
A <- readRDS(sf2a_derived("sf2a")); def <- A$definition
P <- readRDS(sf2a_derived("primary_k05"))
Lt <- truth$Lambda; Lt[def$reversed, ] <- -Lt[def$reversed, ]
cong <- diag(psych::factor.congruence(P$keyed[def$items, ], Lt[def$items, ]))
ri <- sf2a_ri(P$with_ri)
keep <- !A$meta_all$excluded_attention
sep_r <- cor(A$meta$separation, truth$separation[keep], use = "complete.obs")
# The RI variance is on the theta parameterization's scale (residual variances 1),
# not the generating scale (total variance 1); the comparable figure is its share
# of each statement's latent-response variance, true 10%.
res <- rbind(data.table(what = paste("congruence", P$fnames), value = round(cong, 3)),
             data.table(what = c("RI variance (theta scale)", "RI % of statement variance (true 10)", "separation r with truth"),
                        value = round(c(ri$variance, ri$mean_pct_item_variance, sep_r), 3)))
print(res)
fwrite(res, sf2a_out("99_recovery.csv"))
