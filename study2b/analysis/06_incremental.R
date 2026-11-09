# 06: H5, incremental validity (Section 5, model 4). Observed scores, complete cases.
#   H5a  for each of life satisfaction, PHQ-4 and PSS-4: step 1 age, gender, the five
#        FFMQ-15 facets and the eight Brief MAIA-2 subscales; step 2 the 8. Supported for
#        an outcome when delta R^2 >= .02 and the step-2 F-test is significant after Holm
#        correction across the three; H5a holds when at least one outcome is supported.
#   H5b  proportional-odds regression of practice frequency (0-4): step 1 age, gender,
#        spirituality, religiosity, the FFMQ-15 facets and the Brief MAIA-2 subscales;
#        step 2 the 8. Supported when its coefficient is positive with likelihood-ratio p < .05.
# PREREG: gender is man / woman / other (non-binary or self-described), a factor with
# "man" as reference; both steps of a regression use the cases complete on every
# variable of step 2. The 10 is reported beside the 8.

source("study2b/analysis/config.R")
source("study2b/analysis/sf2b_models.R")
A <- readRDS(sf2b_derived("sf2b")); M <- A$meta; S <- A$scores

FFMQ <- SF2B_MEASURES$ffmq; MAIA <- SF2B_MEASURES$maia
base <- data.table(age = M$age, gender = factor(M$gender3, levels = c("man", "woman", "other")),
                   spiritual = M$spiritual, religious = M$religious, freq_code = M$freq_code)
base <- cbind(base, S[, c(FFMQ, MAIA, "sf8", "sf10", "lifesat", "phq4_total", "pss4"), with = FALSE])

h5a <- function(form) {
  rows <- rbindlist(lapply(c("lifesat", "phq4_total", "pss4"), function(y) {
    d <- base[, c(y, "age", "gender", FFMQ, MAIA, form), with = FALSE]; d <- d[complete.cases(d)]
    d[, gender := droplevels(gender)]
    m1 <- lm(reformulate(c("age", if (nlevels(d$gender) > 1) "gender", FFMQ, MAIA), y), data = d)
    m2 <- update(m1, as.formula(paste(". ~ . +", form)))
    a <- anova(m1, m2)
    co <- summary(m2)$coefficients[form, ]
    data.table(outcome = y, n = nrow(d), r2_step1 = summary(m1)$r.squared, r2_step2 = summary(m2)$r.squared,
               delta_r2 = summary(m2)$r.squared - summary(m1)$r.squared, F = a$F[2], df1 = a$Df[2], df2 = a$Res.Df[2],
               p = a$`Pr(>F)`[2], b = co[["Estimate"]], beta = co[["Estimate"]] * sd(d[[form]]) / sd(d[[y]]))
  }))
  rows[, p_holm := p.adjust(p, "holm")]
  rows[, supported := delta_r2 >= SF2B$h5_delta_r2 & p_holm < .05]
  rows
}
h5b <- function(form) {
  d <- base[, c("freq_code", "age", "gender", "spiritual", "religious", FFMQ, MAIA, form), with = FALSE]; d <- d[complete.cases(d)]
  d[, `:=`(gender = droplevels(gender), freq = factor(freq_code, levels = 0:4, ordered = TRUE))]
  d[, freq := droplevels(freq)]
  z <- function(x) (x - mean(x)) / sd(x)
  d[, x := z(get(form))]
  f1 <- reformulate(c("age", if (nlevels(d$gender) > 1) "gender", "spiritual", "religious", FFMQ, MAIA), "freq")
  m1 <- MASS::polr(f1, data = d, Hess = TRUE); m2 <- MASS::polr(update(f1, . ~ . + x), data = d, Hess = TRUE)
  lr <- 2 * (as.numeric(logLik(m2)) - as.numeric(logLik(m1)))
  co <- coef(summary(m2))["x", ]
  m0 <- MASS::polr(freq ~ 1, data = d)
  mcf <- function(m) 1 - as.numeric(logLik(m)) / as.numeric(logLik(m0))
  data.table(n = nrow(d), coef_per_sd = co[["Value"]], se = co[["Std. Error"]], odds_ratio_per_sd = exp(co[["Value"]]),
             lr_chisq = lr, p = pchisq(lr, 1, lower.tail = FALSE), mcfadden_r2_step1 = mcf(m1), mcfadden_r2_step2 = mcf(m2),
             supported = co[["Value"]] > 0 && pchisq(lr, 1, lower.tail = FALSE) < .05)
}
res_a <- rbindlist(lapply(c("sf8", "sf10"), function(f) cbind(form = f, h5a(f))))
res_b <- rbindlist(lapply(c("sf8", "sf10"), function(f) cbind(form = f, h5b(f))))
num <- function(DT) { cols <- names(DT)[vapply(DT, is.numeric, TRUE)]; DT[, (cols) := lapply(.SD, signif, 4), .SDcols = cols]; DT }
fwrite(num(copy(res_a)), sf2b_out("06_h5a.csv")); fwrite(num(copy(res_b)), sf2b_out("06_h5b.csv"))
dec <- list(H5a = any(res_a[form == "sf8", supported]), H5a_outcomes_supported = res_a[form == "sf8" & supported == TRUE, outcome],
            H5b = isTRUE(res_b[form == "sf8", supported]))
write_json(dec, sf2b_out("06_h5_decision.json"), auto_unbox = TRUE, pretty = TRUE)
print(res_a[, .(form, outcome, n, delta_r2 = r3(delta_r2), p_holm = signif(p_holm, 3), supported)]); print(res_b[, .(form, n, odds_ratio_per_sd, p, supported)])
