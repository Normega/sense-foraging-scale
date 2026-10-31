# 16: exploratory, not preregistered (Norm's question, 2026-10-31): who does well on the
# vignettes (separation = mean sensing rating of the sensing vignettes minus mean of the
# doing ones), and does doing well go with higher scores on 12's best 8 or the trial's 10?
# Real mode only, after 12.
#   1. separation and "tells the modes apart" by each background variable, and both
#      regressed on all of them together, with effort (seconds per statement, checks failed)
#   2. the scores' correlations with separation and with its two halves, means by third of
#      separation, and the scores regressed on separation with the background variables

source("study2a/analysis/config.R")
source("study2a/analysis/sf2a_models.R")
stopifnot(SF2A$mode == "real")
A <- readRDS(sf2a_derived("sf2a"))
def <- A$definition; K <- A$items_keyed; M <- copy(A$meta)
out16 <- function(f) sf2a_out(paste0("16_", f))
num3 <- function(T) T[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 3) else v)]

b8 <- fread(sf2a_out("12_strongest_scale.csv"))[n == 8 & min_reversed == 0, strsplit(ids, " ")[[1]]]
M[, `:=`(best8 = rowMeans(K[, b8], na.rm = TRUE),
         trial10 = rowMeans(cbind(K, A$end)[, unname(def$av10)], na.rm = TRUE),
         log_secs = log(secs_per_item), checks = checks_failed,
         gender3 = factor(ifelse(gender %in% c("man", "woman"), gender, "other or not said"), c("woman", "man", "other or not said")),
         educ = factor(fcase(education %in% c("lt_hs", "hs"), "high school or less", education == "some_college", "some college",
                             education == "bachelors", "bachelor's", education == "graduate", "graduate"),
                       c("high school or less", "some college", "bachelor's", "graduate")),
         country6 = ifelse(country == "other", "elsewhere", "the six planned countries"),
         practice = factor(practice_group, c("never", "past_or_occasional", "current")),
         age_band = cut(age, c(0, 25, 35, 45, 120), c("18-25", "26-35", "36-45", "46+")))]
M[, thirds := cut(separation, quantile(separation, c(0, 1/3, 2/3, 1), na.rm = TRUE), include.lowest = TRUE,
                  labels = c("bottom third", "middle third", "top third"))]

# ---- 1. who does well ----------------------------------------------------------------------------
by_cat <- function(v) M[!is.na(get(v)), .(n = .N, separation = mean(separation, na.rm = TRUE), sd = sd(separation, na.rm = TRUE),
                                          pct_tell_apart = 100 * mean(tells_apart, na.rm = TRUE),
                                          best8 = mean(best8, na.rm = TRUE), trial10 = mean(trial10, na.rm = TRUE)),
                        by = .(level = as.character(get(v)))][, variable := v][order(level)]
cats <- rbindlist(lapply(c("gender3", "age_band", "educ", "country6", "practice", "tradition", "practice_freq"), by_cat), fill = TRUE)
setcolorder(cats, "variable")
fwrite(num3(cats), out16("by_background.csv"))

cr <- function(x, y) { ok <- complete.cases(x, y); r <- cor(x[ok], y[ok]); p <- cor.test(x[ok], y[ok])$p.value
  data.table(n = sum(ok), r = r, p = p) }
conts <- c("age", "religious", "spiritual", "ius_total", "log_secs", "checks", "longstring", "vig_self")
cont_tab <- rbindlist(lapply(conts, function(v) cbind(variable = v, cr(M[[v]], M$separation))))
fwrite(num3(cont_tab), out16("separation_correlations.csv"))

z <- function(x) as.numeric(scale(x))
D <- M[, .(sep_z = z(separation), tells_apart, age_z = z(age), gender3, educ, country6, practice, spiritual_z = z(spiritual),
           ius_z = z(ius_total), log_secs_z = z(log_secs), checks)]
f_sep <- lm(sep_z ~ age_z + gender3 + educ + country6 + practice + spiritual_z + ius_z + log_secs_z + checks, data = D)
f_ta <- glm(tells_apart ~ age_z + gender3 + educ + country6 + practice + spiritual_z + ius_z + log_secs_z + checks, data = D, family = binomial)
tidy_lm <- function(f, what) { s <- summary(f)$coefficients
  data.table(model = what, term = rownames(s), estimate = s[, 1], se = s[, 2], p = s[, 4]) }
reg1 <- rbind(tidy_lm(f_sep, "separation (SD units)"), tidy_lm(f_ta, "tells apart (log odds)"))
fwrite(num3(reg1), out16("separation_regression.csv"))
r2_1 <- data.table(model = c("separation", "separation, effort only", "separation, background only"),
                   r2 = c(summary(f_sep)$r.squared, summary(lm(sep_z ~ log_secs_z + checks, data = D))$r.squared,
                          summary(lm(sep_z ~ age_z + gender3 + educ + country6 + practice + spiritual_z + ius_z, data = D))$r.squared),
                   n = c(nobs(f_sep), NA, NA))

# ---- 2. scores ----------------------------------------------------------------------------------
score_cor <- rbindlist(lapply(c("best8", "trial10"), function(s) rbindlist(lapply(c("separation", "sensing_mean", "doing_mean"), function(v)
  cbind(score = s, with = v, cr(M[[v]], M[[s]]))))))
fwrite(num3(score_cor), out16("score_correlations.csv"))
thirds_tab <- M[!is.na(thirds), .(n = .N, separation = mean(separation), best8 = mean(best8, na.rm = TRUE), best8_sd = sd(best8, na.rm = TRUE),
                                  trial10 = mean(trial10, na.rm = TRUE), trial10_sd = sd(trial10, na.rm = TRUE)), by = thirds][order(thirds)]
fwrite(num3(thirds_tab), out16("scores_by_third.csv"))
D[, `:=`(best8_z = z(M$best8), trial10_z = z(M$trial10), sep2 = z(M$separation)^2)]
reg2 <- rbindlist(lapply(c("best8_z", "trial10_z"), function(s) rbind(
  tidy_lm(lm(as.formula(paste(s, "~ sep_z + sep2")), data = D), paste(s, "~ separation (+ squared)")),
  tidy_lm(lm(as.formula(paste(s, "~ sep_z + age_z + gender3 + educ + country6 + practice + spiritual_z + ius_z + log_secs_z + checks")), data = D),
          paste(s, "~ separation + background")))))
fwrite(num3(reg2), out16("score_regression.csv"))
# does practice predict the score differently in good and poor separators?
inter <- rbindlist(lapply(c("best8_z", "trial10_z"), function(s)
  tidy_lm(lm(as.formula(paste(s, "~ practice * tells_apart")), data = D), paste(s, "~ practice x tells apart"))))
fwrite(num3(inter), out16("practice_by_separation.csv"))

print(num3(cats)); print(num3(cont_tab)); print(num3(reg1)); print(num3(r2_1)); print(num3(score_cor)); print(num3(thirds_tab))
print(num3(reg2)); print(num3(inter))
