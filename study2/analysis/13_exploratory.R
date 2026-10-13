# 13: exploratory analyses named in Section 5 ("Exploratory analysis").
#   Order effects: does a validity measure's position in the shuffled block (1-11)
#     affect its mean (score ~ position) or its correlation with the short form
#     (score ~ Access+Valuing x position)? Validity sample.
#   The bifactor alternative: a general factor plus four specific factors, with the
#     random intercept, fitted as a bifactor ESEM (one orthogonal efa block, target
#     rotation: the general factor free on every item, each specific factor targeted to
#     its items). omega_h = (sum of keyed general loadings)^2 / variance of the keyed
#     sum of latent responses implied by the model (k' Sigma k). Not psych::omega:
#     its flip = TRUE default re-reverses keyed items and inflated Study 1's omega_h.
#   Further outcomes: time outdoors and life satisfaction with each factor (the latent
#     correlations of 04). Latent means across country and posting are in 06.
#   The following week (follow-up sample): each past-week item at the follow-up
#     regressed on baseline Access + Valuing, controlling for the same item at baseline.
#     PREREG: linear regression for all three items (time outdoors, 0-4, as numeric);
#     every linked pair, in or out of the 14-21-day window; a pair missing the item at
#     either occasion ("nothing unexpected happened" included) drops out for that item.
#     Reported: the raw and standardized coefficient of Access + Valuing, its 95% CI.

source("study2/analysis/config.R")
source("study2/analysis/sf2_models.R")
A <- readRDS(sf2_derived("sf2"))
D <- A$defs; M <- A$meta; S <- A$scores; items <- D$sf_items; F4 <- D$factors4
VS <- M$reached_end_validity

# ---- order effects ----------------------------------------------------------------------------------
slug_of <- c(unlist(unname(lapply(D$validity, function(v) setNames(rep(v$slug, length(v$scores)), names(v$scores))))),
             crit_pause_days = "sf2-past-week", crit_notice_first = "sf2-past-week", crit_outdoors = "sf2-past-week", lifesat = "lifesat")
POS <- A$validity_positions
OE <- rbindlist(lapply(SF2_VALIDITY_SCORES$score, function(s) {
  pos <- POS[[slug_of[[s]]]][VS]; y <- S[[s]][VS]; av <- S$AccessValuing[VS]
  d <- data.table(y = y, pos = pos, av = av)[complete.cases(y, pos, av)]
  if (nrow(d) < 30 || length(unique(d$pos)) < 3) return(data.table(score = s, n = nrow(d)))
  m1 <- summary(lm(y ~ pos, d))$coefficients; m2 <- summary(lm(scale(y) ~ scale(av) * pos, d))$coefficients
  data.table(score = s, n = nrow(d), mean_slope_per_position = m1["pos", 1], mean_slope_in_sd = m1["pos", 1] / sd(d$y),
             p_mean = m1["pos", 4], r_slope_per_position = m2["scale(av):pos", 1], p_moderation = m2["scale(av):pos", 4])
}), fill = TRUE)
fwrite(OE[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 4) else x)], sf2_out("13_order_effects.csv"))

# ---- bifactor ESEM ------------------------------------------------------------------------------------
FB <- c("G", F4)
member <- sf2_member(items, D$sf_factor, F4)
memb <- cbind(G = TRUE, member)
tgb <- list(target = matrix(0, length(items), 5, dimnames = list(items, FB)), mask = ifelse(memb, 0, 1))
bf <- sf2_try(sf2_sem(sf2_esem_syntax(items, FB, TRUE), as.data.frame(A$sf_raw), items, target = tgb, orthogonal = TRUE))
bif <- if (sf2_failed(bf)) list(error = sf2_msg(bf)) else {
  L <- sf2_loadings(bf, FB, items); key <- ifelse(items %in% D$sf_rev, -1, 1)
  O <- sf2_orient(L, memb, D$sf_rev)
  Sig <- lavInspect(bf, "cor.ov")[items, items]
  tot <- drop(t(key) %*% Sig %*% key)
  lg <- O$keyed[, "G"]
  spec <- sapply(F4, function(f) { it <- items[member[, f]]; sum(O$keyed[it, f])^2 })
  list(fit = sf2_fit_stats(bf), omega_h = sum(lg)^2 / tot,
       ecv_general = sum(lg^2) / sum(O$keyed^2),
       omega_h_note = "(sum keyed general loadings)^2 / k' Sigma k; Sigma = model-implied latent-response correlations, k = +-1 keying",
       loadings = data.table(item = items, factor = D$sf_factor[items], r3(O$keyed)),
       ri = sf2_ri(bf))
}
# ---- the following week -------------------------------------------------------------------------------
P2 <- A$retest$pairs; PW1 <- A$retest$past_week_t1; PW2 <- A$retest$past_week_t2
av1 <- S$AccessValuing[match(P2$participant, M$participant)]
FW <- rbindlist(lapply(colnames(PW2), function(v) {
  d <- data.table(y2 = PW2[, v], av = av1, y1 = PW1[, v]); d <- d[complete.cases(d)]
  if (nrow(d) < 20) return(data.table(item = v, n = nrow(d)))
  cs <- summary(lm(y2 ~ av + y1, d))$coefficients
  ms <- lm(scale(y2) ~ scale(av) + scale(y1), d); css <- summary(ms)$coefficients; ci <- confint(ms)["scale(av)", ]
  data.table(item = v, n = nrow(d), b = cs["av", 1], se = cs["av", 2], p = cs["av", 4], beta = css["scale(av)", 1],
             beta_lower = ci[[1]], beta_upper = ci[[2]], beta_baseline_item = css["scale(y1)", 1])
}), fill = TRUE)
fwrite(FW[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 4) else x)], sf2_out("13_following_week.csv"))

further <- if (file.exists(sf2_out("04_latent_validity.csv")))
  fread(sf2_out("04_latent_validity.csv"))[score %in% c("crit_outdoors", "lifesat"), .(score, factor, r, lower, upper, predicted)] else NULL

write_json(list(order_effects = OE, bifactor = bif, further_outcomes = further, following_week = FW), sf2_out("13_exploratory.json"),
           auto_unbox = TRUE, pretty = TRUE, digits = 4, na = "null", dataframe = "rows")
print(OE[order(p_mean)][1:5]); print(FW); if (is.null(bif$error)) sf2_log("bifactor ESEM omega_h = ", round(bif$omega_h, 3)) else sf2_log("bifactor: ", bif$error)
