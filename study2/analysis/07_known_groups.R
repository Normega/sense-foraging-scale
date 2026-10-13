# 07: H6, known groups (Section 5, analysis 5), on everyone analysable with background
# answers and the short form.
#
# Dedicated practitioners vs never-practitioners on each observed score, a linear model
# with posting as a covariate; Cohen's d = the group coefficient over the pooled SD,
# with a 95% CI; Holm across Access and Valuing and savouring (two-tailed). H6 holds
# for a factor when the difference favours the dedicated, p_Holm < .05 and d >= 0.50.
# Understanding, Letting be and Access + Valuing are compared with no prediction.
# Pre-specified sensitivity analyses, beside the primary test:
#   (a) + spirituality and religiosity as covariates
#   (b) dedicated restricted to those who name meditation or yoga among their practices
#   (c) latent scores with agreement modelled: the H1 ESEM in the dedicated + never
#       subsample, each content factor regressed on the group and posting (a MIMIC);
#       the coefficient is in factor-residual-SD units, the latent analogue of d.
# Also the primary test without flagged participants (Section 5, "Data exclusion").
# Secondary: a linear trend across the five frequency categories (codes 0-4) among all
# who practise regularly, posting as a covariate; and latent mean differences in the
# invariance models (four-factor CFA, Access + Valuing CFA) across the two practice
# groups, reported only where the intercepts step holds.

source("study2/analysis/config.R")
source("study2/analysis/sf2_models.R")
A <- readRDS(sf2_derived("sf2"))
D <- A$defs; M <- A$meta; S <- A$scores; items <- D$sf_items; F4 <- D$factors4
SC <- c("Access", "Valuing", "AccessValuing", "Understanding", "LettingBe")
grp <- ifelse(M$practice_group %in% "dedicated", 1L, ifelse(M$practice_group %in% "never", 0L, NA_integer_))
grp_b <- ifelse(M$practice_group %in% "dedicated" & M$meditation_or_yoga, 1L, ifelse(M$practice_group %in% "never", 0L, NA_integer_))
post <- data.table(posting = M$posting_contemplative)

kg_set <- function(label, g, covars) {
  t <- rbindlist(lapply(SC, function(s) cbind(analysis = label, score = s, sf2_kg(S[[s]], g, covars))))
  t[score %in% c("Access", "Valuing"), p_holm := p.adjust(p, "holm")]
  t[, h6_rule := score %in% c("Access", "Valuing") & coef > 0 & p_holm < .05 & d >= .50]
  t
}
KG <- rbind(kg_set("primary", grp, post),
            kg_set("(a) + spirituality, religiosity", grp, cbind(post, spiritual = M$spiritual, religious = M$religious)),
            kg_set("(b) dedicated = meditation or yoga", grp_b, post),
            kg_set("primary, without flagged participants", ifelse(M$flagged, NA_integer_, grp), post), fill = TRUE)

# (c) latent, agreement modelled.
sub <- !is.na(grp)
dl <- data.frame(A$sf_raw[sub, ], dedicated = grp[sub], posting = M$posting_contemplative[sub])
member <- sf2_member(items, D$sf_factor, F4)
mimic <- sf2_try(sf2_sem(c(sf2_esem_syntax(items, F4, TRUE), paste0(paste(F4, collapse = " + "), " ~ dedicated + posting")),
                         dl, items, target = sf2_target(member)))
LAT <- if (sf2_failed(mimic)) data.table(error = sf2_msg(mimic)) else {
  o <- sf2_orient(sf2_loadings(mimic, F4, items), member, D$sf_rev)$orient
  pe <- as.data.table(parameterEstimates(mimic))[op == "~" & rhs == "dedicated" & lhs %in% F4]
  pe[, .(analysis = "(c) latent, agreement modelled", score = lhs, n_group = sum(dl$dedicated), n_reference = sum(dl$dedicated == 0),
         coef = est * o[lhs], se, p = pvalue, d = est * o[lhs],
         d_lower = ifelse(o[lhs] > 0, ci.lower, -ci.upper), d_upper = ifelse(o[lhs] > 0, ci.upper, -ci.lower))]
}
if (!"error" %in% names(LAT)) {
  LAT[score %in% c("Access", "Valuing"), p_holm := p.adjust(p, "holm")]
  LAT[, h6_rule := score %in% c("Access", "Valuing") & coef > 0 & p_holm < .05 & d >= .50]
}
KG <- rbind(KG, LAT, fill = TRUE)
fwrite(KG[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 4) else x)], sf2_out("07_known_groups.csv"))

# Secondary: frequency trend among those who practise regularly.
pr <- M$practice_regular %in% "yes" & !is.na(M$freq_code)
TR <- rbindlist(lapply(SC, function(s) {
  d <- data.table(y = S[[s]][pr], f = M$freq_code[pr], posting = M$posting_contemplative[pr])[complete.cases(y, f, posting)]
  m <- lm(y ~ f + posting, d); cf <- summary(m)$coefficients["f", ]
  data.table(score = s, n = nrow(d), slope_per_category = cf[["Estimate"]], se = cf[["Std. Error"]], p = cf[["Pr(>|t|)"]],
             slope_in_sd = cf[["Estimate"]] / sd(d$y))
}))
fwrite(TR[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 4) else x)], sf2_out("07_frequency_trend.csv"))

# Secondary: latent means across practice groups, where scalar invariance holds.
pg <- ifelse(M$practice_group %in% c("dedicated", "never"), M$practice_group, NA)
pg <- ifelse(is.na(pg), NA, ifelse(pg == "never", "1_never", "2_dedicated"))   # never = reference group
sf20 <- items; av10 <- items[D$sf_factor[items] %in% c("Access", "Valuing")]
inv <- list(cfa20 = sf2_run_mg(A$sf_raw, sf20, sf2_pattern_cfa(sf20, D$sf_factor, F4), pg, "cfa20 / practice group"),
            av10 = sf2_run_mg(A$sf_raw, av10, sf2_pattern_cfa(av10, D$sf_factor, c("Access", "Valuing")), pg, "av10 / practice group"))
IT <- rbindlist(lapply(inv, `[[`, "table"), fill = TRUE)
fwrite(IT, sf2_out("07_practice_invariance.csv"))
LMN <- rbindlist(lapply(names(inv), function(k) { r <- inv[[k]]; f <- r$fits$intercepts
  if (!isTRUE(r$completed) || is.null(f) || sf2_failed(f)) return(data.table(model = k, note = "intercepts step does not hold: no latent means"))
  fn <- colnames(r$pattern); lam <- lavInspect(f, "est")[[1]]$lambda[rownames(r$pattern), fn, drop = FALSE]
  o <- sf2_orient(lam, sf2_member(rownames(r$pattern), D$sf_factor, fn), D$sf_rev)$orient
  cbind(model = k, note = "dedicated minus never, in never-group SD units", sf2_latent_means(f, fn, o)) }), fill = TRUE)
fwrite(LMN, sf2_out("07_practice_latent_means.csv"))

dec <- list(h6 = KG[analysis == "primary" & score %in% c("Access", "Valuing"), .(score, n_group, n_reference, d, d_lower, d_upper, p, p_holm, supported = h6_rule)],
            sensitivity = KG[analysis != "primary" & score %in% c("Access", "Valuing"), .(analysis, score, d, p_holm, h6_rule)],
            n_dedicated = sum(grp == 1, na.rm = TRUE), n_never = sum(grp == 0, na.rm = TRUE))
write_json(dec, sf2_out("07_h6_decision.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null", dataframe = "rows")
print(KG[, .(analysis, score, n_group, n_reference, d = round(d, 3), p_holm = signif(p_holm, 3), h6_rule)])
print(TR); print(LMN)
