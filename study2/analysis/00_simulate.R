# 00 (simulated mode only): a raw extraction generated from a KNOWN model, written in
# exactly the shape extract.sql returns through the Supabase CLI
# ({"boundary", "rows", "warning"}, one row per enrollment, questionnaire submission,
# VAS answer and exited step), so that 01 onward run unchanged on real data.
#
# The generating model:
#   - five correlated content factors (Access, Understanding, Valuing, Letting be,
#     Forage) for the 20 short-form items and 24 candidates, in the theta metric
#     (residual variance 1) that the registered models use; loadings on the a-priori
#     factor plus cross-loadings of ~.2-.3 (standardized) on the six items Study 1's
#     ESEM showed them for and ~.1 on six more (so the simple-structure CFA misfits);
#     Understanding has two weak short-form items (so it should NOT be confirmed) and
#     Letting be is confirmed but below omega .70;
#   - a random-intercept agreement factor (variance .22, ~10% of item variance),
#     loading 1 on every unrecoded item;
#   - 6-point answers cut at item-specific thresholds; "Prefer not to answer" at random;
#   - validity constructs with set correlations to the factors (chosen so that some
#     hypotheses pass and some fail: H2 fails for Letting be with EQ Decentering, H3c
#     fails for Letting be with BIDR SDE, H4 fails for Letting be with noticing first,
#     H6 fails for Access on d >= .50, H7 fails for Valuing), their items cut to each
#     instrument's scale;
#   - practice groups at Study 1's rates (about 68 dedicated, 286 never in 1,000) with
#     latent shifts; one item with intercept non-invariance across postings (sf6_15);
#   - a two-week retest of the short form and the past-week items for about 300 of
#     400 invited;
#   - country of residence: 31% in the six listed countries, 69% elsewhere (Study 1's split).
# Planted mess (each case's code is recorded in the truth file, for 99):
#   staff test enrollments, enrollments without a short form, partial sessions that
#   stop before or inside the validity block, both / one / the only reached attention
#   check failed, people who never reached a check, speeders, longstrings, many
#   "Prefer not to answer", one flagged resubmission, one later re-entry with
#   different answers, retests outside 14-21 days, an unlinked retest, a retest test.
# Truth values (population quantities from 60,000 simulated people under the same
# model) go to derived/simulation_truth_simulated.rds: the generating values, and the
# estimand the target-rotated ESEM converges to (see "The ESTIMAND" below), which is
# what 99 checks recovery against.

source("study2/analysis/config.R")
source("study2/analysis/sf2_models.R")
stopifnot(SF2$mode == "simulated")
D <- sf2_definitions()
set.seed(SF2$seed_sim)

FN5 <- c("Access", "Understanding", "Valuing", "LettingBe", "Forage")
Phi <- diag(5); dimnames(Phi) <- list(FN5, FN5)
setr <- function(a, b, r) { Phi[a, b] <<- r; Phi[b, a] <<- r }
setr("Access", "Understanding", .35); setr("Access", "Valuing", .60); setr("Access", "LettingBe", .30); setr("Access", "Forage", .45)
setr("Understanding", "Valuing", .35); setr("Understanding", "LettingBe", .15); setr("Understanding", "Forage", .25)
setr("Valuing", "LettingBe", .25); setr("Valuing", "Forage", .55); setr("LettingBe", "Forage", .20)
RI_VAR <- 0.22

# Keyed loadings, theta metric (standardized ~ value / sqrt(1 + .22 + common variance)).
KEYED <- list(
  sf6_35 = c(Access = 1.20), sf6_36 = c(Access = 1.30, LettingBe = .20), sf6_9 = c(Access = 1.30, Understanding = .20),
  sf6_8 = c(Access = 1.10), sf6_7 = c(Access = 1.00, Valuing = .50),
  sf6_22 = c(Understanding = .90, Valuing = .55), sf6_3 = c(Understanding = .85, Valuing = .20), sf6_15 = c(Understanding = .75, LettingBe = .20),
  sf6_31 = c(Understanding = .35), sf6_37 = c(Understanding = .32),
  sf6_12 = c(Valuing = 1.20, Access = .25), sf6_18 = c(Valuing = 1.00), sf6_11 = c(Valuing = 1.10),
  sf6_27 = c(Valuing = .95, Understanding = .55), sf6_44 = c(Valuing = .90, Understanding = .50),
  sf6_40 = c(LettingBe = .75, Valuing = .55), sf6_46 = c(LettingBe = .85, Access = .45), sf6_43 = c(LettingBe = .70),
  sf6_47 = c(LettingBe = .65), sf6_21 = c(LettingBe = .70, Valuing = .25),
  sfc_2 = c(Understanding = .45), sfc_3 = c(Understanding = .40), sfc_7 = c(Understanding = 1.10),
  sfc_8 = c(Understanding = 1.00), sfc_9 = c(Understanding = .85),
  sfc_5 = c(LettingBe = 1.10), sfc_6 = c(LettingBe = 1.05), sfc_10 = c(LettingBe = 1.00), sfc_11 = c(LettingBe = .90),
  sfc_12 = c(LettingBe = .75), sfc_13 = c(LettingBe = .60), sfc_14 = c(LettingBe = .95),
  sfc_4 = c(Valuing = .80), sfc_15 = c(Valuing = .90),
  sfc_16 = c(Forage = 1.10), sfc_17 = c(Forage = 1.00), sfc_18 = c(Forage = .90, Valuing = .30),
  sfc_19 = c(Forage = .70), sfc_20 = c(Forage = .60),
  sfc_21 = c(Access = 1.10), sfc_22 = c(Access = 1.00), sfc_23 = c(Access = .85, Understanding = .55),
  sfc_24 = c(Valuing = 1.00), sfc_25 = c(Valuing = .80, LettingBe = .30))
ALL_ITEMS <- c(D$sf_items, D$cand_items)
REV <- c(D$sf_rev, D$cand[reversed == TRUE, id])
stopifnot(setequal(names(KEYED), ALL_ITEMS))
LAM <- matrix(0, length(ALL_ITEMS), 5, dimnames = list(ALL_ITEMS, FN5))
for (i in ALL_ITEMS) LAM[i, names(KEYED[[i]])] <- KEYED[[i]] * if (i %in% REV) -1 else 1
DIF <- list(item = "sf6_15", posting = "contemplative", shift = 1.6)   # intercept non-invariance (06 should release it)
PNA_RATE <- setNames(rep(.008, length(ALL_ITEMS)), ALL_ITEMS); PNA_RATE["sfc_10"] <- .07   # sfc_10 fails the <= 5% rule
SHIFT <- rbind(dedicated = c(.45, .10, 1.00, .10, .50), other = c(.15, 0, .30, 0, .15), never = rep(0, 5))
colnames(SHIFT) <- FN5
RHO <- c(Access = .90, Understanding = .85, Valuing = .82, LettingBe = .85, Forage = .85)   # two-week latent stability

# Item thresholds: agree-leaning for forward items, mirrored for reversed (raw scale).
p_fwd <- c(.03, .08, .15, .34, .28, .12); p_rev <- c(.10, .25, .25, .22, .12, .06)
CUTS <- lapply(setNames(ALL_ITEMS, ALL_ITEMS), function(i) {
  p <- (if (i %in% REV) p_rev else p_fwd) * runif(6, .8, 1.2); p <- p / sum(p)
  sd_y <- sqrt(drop(t(LAM[i, ]) %*% Phi %*% LAM[i, ]) + RI_VAR + 1)
  qnorm(cumsum(p)[1:5]) * sd_y
})

# Validity constructs: correlation with (Access, Understanding, Valuing, Letting be,
# Forage), loadings on three auxiliary traits (mindfulness, negative affect, openness)
# and on the standardized RI.
TGT <- rbindlist(list(
  list("ffmq_observing", .50, .20, .45, .15, .45, .30, 0, 0, 0), list("ffmq_describing", .20, .15, .05, .05, .05, .40, 0, 0, 0),
  list("ffmq_acting", .20, .05, .05, .20, 0, .40, -.20, 0, 0), list("ffmq_nonjudging", .15, 0, 0, .25, 0, .30, -.30, 0, 0),
  list("ffmq_nonreactivity", .40, .15, .20, .40, .10, .40, 0, 0, 0),
  list("maia_noticing", .45, .20, .25, .05, .30, .20, 0, 0, 0), list("maia_notdistracting", .05, 0, .05, .15, 0, 0, 0, 0, 0),
  list("maia_notworrying", .05, 0, 0, .20, 0, 0, -.30, 0, 0), list("maia_attnregulation", .50, .20, .20, .20, .20, .30, 0, 0, 0),
  list("maia_emoawareness", .30, .20, .25, .05, .20, .20, 0, 0, 0), list("maia_selfregulation", .40, .15, .25, .25, .10, .30, -.20, 0, 0),
  list("maia_bodylistening", .30, .15, .40, .05, .25, .20, 0, 0, 0), list("maia_trusting", .30, .10, .25, .10, .10, 0, -.20, 0, 0),
  list("bidr_sde", .20, .05, .05, .33, .05, 0, -.20, 0, .20), list("bidr_im", .05, 0, .05, .05, 0, 0, 0, 0, .20),
  list("phq_anx", -.20, 0, 0, -.20, 0, 0, .70, 0, 0), list("phq_dep", -.20, 0, -.05, -.15, 0, 0, .70, 0, 0),
  list("ius", -.15, 0, 0, -.40, 0, 0, .40, 0, 0), list("sbi", .45, .15, .62, .15, .35, 0, -.20, 0, 0),
  list("bfi_as", .20, .15, .45, .15, .35, 0, 0, .50, 0), list("bfi_ic", .15, .15, .20, .15, .30, 0, 0, .50, 0),
  list("bfi_ci", .15, .10, .15, .15, .25, 0, 0, .50, 0), list("eq", .45, .25, .20, .31, .15, .30, 0, 0, 0),
  list("pss", -.20, 0, 0, -.25, 0, 0, .60, 0, 0), list("lifesat", .25, 0, .35, .15, .10, 0, -.30, 0, 0),
  list("pause", .35, .05, .35, .05, .35, 0, 0, 0, 0), list("notice", .25, .15, .05, -.10, .05, .20, 0, 0, 0),
  list("outdoors", .15, 0, .20, 0, .15, 0, 0, 0, 0)))
setnames(TGT, c("T", FN5, "mind", "neg", "open", "ri"))
TB <- t(solve(Phi, t(as.matrix(TGT[, ..FN5]))))                     # weights on the factors
TC <- as.matrix(TGT[, .(mind, neg, open, ri)])
TE <- 1 - rowSums((TB %*% Phi) * TB) - rowSums(TC^2)
stopifnot(all(TE > .1)); rownames(TB) <- rownames(TC) <- names(TE) <- TGT$T

# Each validity scale: its items, the construct each item measures, keying relative
# to the construct, and construct-direction answer proportions.
V <- D$validity
item_T <- list(
  ffmq15 = c(setNames(rep("ffmq_observing", 3), V$ffmq15$scores$ffmq_observing), setNames(rep("ffmq_describing", 3), V$ffmq15$scores$ffmq_describing),
             setNames(rep("ffmq_acting", 3), V$ffmq15$scores$ffmq_acting), setNames(rep("ffmq_nonjudging", 3), V$ffmq15$scores$ffmq_nonjudging),
             setNames(rep("ffmq_nonreactivity", 3), V$ffmq15$scores$ffmq_nonreactivity)),
  maia2 = unlist(lapply(names(V$maia2$scores), function(s) setNames(rep(s, length(V$maia2$scores[[s]])), V$maia2$scores[[s]]))),
  bidr16 = c(setNames(rep("bidr_sde", 8), V$bidr16$scores$bidr_sde), setNames(rep("bidr_im", 8), V$bidr16$scores$bidr_im)),
  phq4 = c(setNames(rep("phq_anx", 2), V$phq4$scores$phq4_anxiety), setNames(rep("phq_dep", 2), V$phq4$scores$phq4_depression)),
  ius12 = setNames(rep("ius", 12), V$ius12$items), sbi = setNames(rep("sbi", 8), V$sbi$items),
  bfi2 = c(setNames(rep("bfi_as", 4), V$bfi2$scores$bfi_aesthetic), setNames(rep("bfi_ic", 4), V$bfi2$scores$bfi_intellectual),
           setNames(rep("bfi_ci", 4), V$bfi2$scores$bfi_creative)),
  eq = setNames(rep("eq", 11), V$eq$items), pss4 = setNames(rep("pss", 4), V$pss4$items))
PROPS <- list(ffmq15 = c(.08, .20, .32, .27, .13), maia2 = c(.05, .12, .22, .28, .22, .11), bidr16 = c(.08, .12, .16, .20, .20, .15, .09),
              phq4 = c(.45, .32, .14, .09), ius12 = c(.25, .30, .22, .15, .08), sbi = c(.04, .07, .12, .18, .25, .22, .12),
              bfi2 = c(.07, .17, .25, .31, .20), eq = c(.06, .18, .34, .30, .12), pss4 = c(.15, .25, .32, .20, .08))
# PSS-4 items 2 and 3 are positively worded: raw answers run against perceived stress.
item_rev <- lapply(names(V), function(k) V[[k]]$reverse); names(item_rev) <- names(V)
VLOAD <- lapply(item_T, function(v) setNames(runif(length(v), .65, .82), names(v)))
SINGLE <- list(pause = c(.20, .12, .14, .14, .10, .08, .06, .16), notice = c(.06, .18, .38, .28, .10),
               outdoors = c(.10, .30, .35, .17, .08), lifesat = c(.04, .10, .22, .34, .22, .08))

# ---- people ----------------------------------------------------------------------------
gen_people <- function(n, posting) {
  p_ded <- ifelse(posting == "general", .03, .125); p_nev <- ifelse(posting == "general", .38, .145)
  u <- runif(n)
  grp <- ifelse(u < p_ded, "dedicated", ifelse(u < p_ded + p_nev, "never", "other"))
  eta <- mvrnorm(n, rep(0, 5), Phi); colnames(eta) <- FN5
  eta <- eta + SHIFT[grp, , drop = FALSE]
  ri <- rnorm(n, 0, sqrt(RI_VAR))
  aux <- cbind(matrix(rnorm(n * 3), n, 3), ri / sqrt(RI_VAR))
  E <- matrix(rnorm(n * nrow(TB)), n); colnames(E) <- rownames(TB)
  Tm <- eta %*% t(TB) + aux %*% t(TC) + sweep(E, 2, sqrt(TE), `*`)
  colnames(Tm) <- rownames(TB)
  list(n = n, posting = posting, group = grp, eta = eta, ri = ri, T = Tm, aux = aux, E = E)
}

gen_sfc <- function(P, eta = P$eta, ri = P$ri) {
  Y <- eta %*% t(LAM) + ri + matrix(rnorm(P$n * length(ALL_ITEMS)), P$n)
  colnames(Y) <- ALL_ITEMS
  Y[P$posting == DIF$posting, DIF$item] <- Y[P$posting == DIF$posting, DIF$item] + DIF$shift
  X <- sapply(ALL_ITEMS, function(i) findInterval(Y[, i], CUTS[[i]]) + 1L)
  list(X = X, Y = Y)
}

gen_validity <- function(P) {
  out <- list()
  for (k in names(item_T)) {
    its <- names(item_T[[k]]); rg <- V[[k]]$range; p <- PROPS[[k]]
    M <- sapply(its, function(i) {
      a <- VLOAD[[k]][[i]]; s <- if (i %in% item_rev[[k]]) -1 else 1
      y <- s * a * P$T[, item_T[[k]][[i]]] + sqrt(1 - a^2) * rnorm(P$n)
      pp <- if (s < 0) rev(p) else p
      rg[1] + findInterval(y, qnorm(cumsum(pp)[-length(pp)]))
    })
    out[[k]] <- matrix(M, P$n, dimnames = list(NULL, its))
  }
  single <- function(t, p, lo) lo + findInterval(sqrt(.65) * P$T[, t] + sqrt(.35) * rnorm(P$n), qnorm(cumsum(p)[-length(p)]))
  out$crit_pause_days <- single("pause", SINGLE$pause, 0)
  out$crit_notice_first <- single("notice", SINGLE$notice, 1)
  out$crit_outdoors <- single("outdoors", SINGLE$outdoors, 0)
  out$lifesat <- single("lifesat", SINGLE$lifesat, 1)
  out
}

# The past week again at the follow-up: the three behaviours' constructs two weeks on,
# from the latent states then (eta2, ri2), auxiliary traits that are stable (.8), and
# unique parts that carry over by half; cut as at baseline.
gen_pastweek2 <- function(P, eta2, ri2) {
  aux2 <- cbind(.8 * P$aux[, 1:3] + .6 * matrix(rnorm(P$n * 3), P$n, 3), ri2 / sqrt(RI_VAR))
  k <- c(crit_pause_days = "pause", crit_notice_first = "notice", crit_outdoors = "outdoors")
  lo <- c(crit_pause_days = 0, crit_notice_first = 1, crit_outdoors = 0)
  sapply(names(k), function(v) { t <- k[[v]]
    T2 <- drop(eta2 %*% TB[t, ]) + drop(aux2 %*% TC[t, ]) + sqrt(TE[[t]]) * (.5 * P$E[, t] + sqrt(.75) * rnorm(P$n))
    lo[[v]] + findInterval(sqrt(.65) * T2 + sqrt(.35) * rnorm(P$n), qnorm(cumsum(SINGLE[[t]])[-length(SINGLE[[t]])])) })
}

# Observed scores exactly as 01 computes them.
score_all <- function(X, VI) {
  S <- data.table(
    Access = sf2_score(X, D$sf_items[D$sf_factor == "Access"], D$sf_rev),
    Understanding = sf2_score(X, D$sf_items[D$sf_factor == "Understanding"], D$sf_rev),
    Valuing = sf2_score(X, D$sf_items[D$sf_factor == "Valuing"], D$sf_rev),
    LettingBe = sf2_score(X, D$sf_items[D$sf_factor == "LettingBe"], D$sf_rev),
    AccessValuing = sf2_score(X, D$sf_items[D$sf_factor %in% c("Access", "Valuing")], D$sf_rev))
  if (!is.null(VI)) {
    for (k in names(item_T)) for (s in names(V[[k]]$scores))
      S[[s]] <- sf2_score(VI[[k]], V[[k]]$scores[[s]], V[[k]]$reverse, V[[k]]$range, V[[k]]$method)
    for (s in c("crit_pause_days", "crit_notice_first", "crit_outdoors", "lifesat")) S[[s]] <- as.numeric(VI[[s]])
  }
  S
}

# ---- the truth: population quantities -------------------------------------------------------
sf2_log("truth: simulating a population of 60,000 under the same model")
NP <- 60000
POP <- gen_people(NP, sample(rep(c("general", "contemplative"), c(.6 * NP, .4 * NP))))
popI <- gen_sfc(POP); popV <- gen_validity(POP); popS <- score_all(popI$X, popV)
sd_eta <- apply(POP$eta, 2, sd); sd_y <- apply(popI$Y, 2, sd)
STD <- sweep(LAM, 2, sd_eta, `*`) / sd_y                       # std.all loadings, raw direction
KSTD <- STD; KSTD[REV, ] <- -KSTD[REV, ]                         # keyed
F4 <- D$factors4
truth <- list(
  Lambda_theta = LAM, Phi_generating = Phi, ri_variance = RI_VAR, keyed_std_generating = KSTD,
  factor_cor_generating = cor(POP$eta), dif = DIF, shift = SHIFT, rho = RHO, pna_rate = PNA_RATE)
vs <- intersect(SF2_VALIDITY_SCORES$score, names(popS))
R_ev <- cor(POP$eta[, F4], as.matrix(popS[, ..vs]), use = "pair")      # generating factors x validity scores

# The ESTIMAND. With cross-loadings in the generating model, the target-rotated ESEM
# does not return the generating factors: rotation is indeterminate, and the target
# criterion picks the rotation with the smallest non-target loadings. What a correct
# implementation converges to in the population is that rotation of the generating
# solution, computed here analytically (GPArotation pst, oblique, started at the
# generating solution; checked against lavaan on 30,000 cases: loadings within .02).
# With orthonormal xi (eta = Lc xi) and A0 = Lambda Lc, the rotated factors are
# f = T' xi, so cor(f, v) = T' Lc^-1 cor(eta, v) and Phi* = T'T.
estimand <- function(its) {
  R4 <- cor(POP$eta[, F4]); Lc <- t(chol(R4))
  mem <- sf2_member(its, setNames(c(D$sf_factor, setNames(D$cand$factor, D$cand$id))[its], its), F4)
  rot <- suppressWarnings(GPArotation::GPFoblq(STD[its, F4] %*% Lc, Tmat = t(Lc), method = "pst", eps = 1e-7, maxit = 20000,
                                               methodArgs = list(W = ifelse(mem, 0, 1), Target = matrix(0, length(its), 4))))
  L <- rot$loadings; dimnames(L) <- list(its, F4); K <- L; K[its %in% REV, ] <- -K[its %in% REV, ]
  Phi_s <- rot$Phi; dimnames(Phi_s) <- list(F4, F4)
  list(raw = L, keyed = K, phi = Phi_s, Tm = rot$Th, Lc = Lc, member = mem)
}
E20 <- estimand(D$sf_items)
R_fv <- t(E20$Tm) %*% solve(E20$Lc) %*% R_ev; rownames(R_fv) <- F4         # rotated factors x validity scores
truth$keyed_std <- E20$keyed; truth$factor_cor <- E20$phi
truth$omega <- sapply(F4, function(f) sf2_omega(E20$keyed[D$sf_items[D$sf_factor == f], f]))
truth$omega_generating <- sapply(F4, function(f) sf2_omega(KSTD[D$sf_items[D$sf_factor == f], f]))
truth$latent_observed <- rbindlist(lapply(vs, function(s) data.table(score = s, factor = c(F4, "RI"),
  r = c(R_fv[, s], cor(POP$ri, popS[[s]], use = "pair")), r_generating = c(R_ev[, s], NA))))
truth$observed_observed <- rbindlist(lapply(vs, function(s) data.table(score = s, factor = c("Access", "Understanding", "Valuing", "LettingBe", "AccessValuing"),
  r = sapply(c("Access", "Understanding", "Valuing", "LettingBe", "AccessValuing"), function(f) cor(popS[[f]], popS[[s]], use = "pair")))))
S_v <- cor(as.matrix(popS[, ..vs]), use = "pair")
mR <- function(r, nb) sqrt(drop(t(r[nb]) %*% solve(S_v[nb, nb], r[nb])))
truth$h3d <- sapply(names(SF2_H3D), function(f) mR(R_fv[f, ], SF2_H3D[[f]]))
truth$h3d_generating <- sapply(names(SF2_H3D), function(f) mR(R_ev[f, ], SF2_H3D[[f]]))
truth$h3b <- R_fv["Valuing", "sbi_moment"] - R_fv["Valuing", "bfi_openmindedness"]
# Step 7's 39-item model: the rotated population loadings, for the expected selections.
truth$keyed_std_39 <- estimand(c(D$sf_items, D$cand[factor != "Forage", id]))$keyed
kg_truth <- function(s) { g <- POP$group; keep <- g %in% c("dedicated", "never")
  sf2_kg(popS[[s]][keep], g[keep] == "dedicated", data.table(posting = as.integer(POP$posting[keep] == "contemplative")))$d }
truth$h6_d <- sapply(c("Access", "Valuing", "AccessValuing", "Understanding", "LettingBe"), kg_truth)
# Retest: latent states two weeks on, fresh item noise; population ICC(2,1), absolute agreement.
eta2 <- sweep(POP$eta, 2, RHO, `*`) + sweep(mvrnorm(NP, rep(0, 5), Phi), 2, sqrt(1 - RHO^2), `*`)
ri2 <- .8 * POP$ri + .6 * rnorm(NP, 0, sqrt(RI_VAR))
popI2 <- gen_sfc(POP, eta2, ri2); popS2 <- score_all(popI2$X, NULL)
icc_pop <- function(a, b) { ok <- complete.cases(a, b); irr::icc(head(cbind(a[ok], b[ok]), 20000), "twoway", "agreement", "single")$value }
truth$icc <- sapply(c("Access", "Valuing", "AccessValuing", "Understanding", "LettingBe"), function(s) icc_pop(popS[[s]], popS2[[s]]))
agree_sf <- SF2_AGREE_PAIRS[grepl("^sf6_", reversed) & grepl("^sf6_", forward)]
truth$icc_agreement <- icc_pop(rowMeans(popI$X[, c(agree_sf$reversed, agree_sf$forward)]), rowMeans(popI2$X[, c(agree_sf$reversed, agree_sf$forward)]))
# The following week (exploratory): each past-week item at the follow-up on baseline
# Access + Valuing, controlling for the same item at baseline (standardized coefficient).
popPW2 <- gen_pastweek2(POP, eta2, ri2)
truth$following_week <- sapply(colnames(popPW2), function(v)
  coef(lm(scale(popPW2[, v]) ~ scale(popS$AccessValuing) + scale(popS[[v]])))[2])
names(truth$following_week) <- colnames(popPW2)
rm(popI, popV, popS, popI2, popS2, eta2, ri2, popPW2); invisible(gc())

# ---- the sample --------------------------------------------------------------------------------
N_COMP <- c(general = 600, contemplative = 400)
plan <- data.table(role = "completer", posting = rep(names(N_COMP), N_COMP))
plan <- rbind(plan,
  data.table(role = "partial_after_sf", posting = sample(names(N_COMP), 8, TRUE)),
  data.table(role = "partial_after_candidates", posting = sample(names(N_COMP), 7, TRUE)),
  data.table(role = "partial_mid_validity", posting = sample(names(N_COMP), 14, TRUE)),
  data.table(role = "no_sf_background_only", posting = sample(names(N_COMP), 6, TRUE)),
  data.table(role = "no_sf_nothing", posting = sample(names(N_COMP), 4, TRUE)),
  data.table(role = "staff_test", posting = sample(names(N_COMP), 3, TRUE)))
plan <- plan[sample(.N)]
N <- nrow(plan)
code <- function(n, prefix) vapply(seq_len(n), function(i) paste0(prefix, paste(sample(c(letters, 0:9), 12, TRUE), collapse = "")), "")
plan[, participant := code(.N, "P-")]
plan[role == "staff_test", participant := sprintf("test-%08d", sample(1e7, .N))]
P <- gen_people(N, plan$posting)
plan[, practice_group := P$group]
SI <- gen_sfc(P); VI <- gen_validity(P)

# Mess on the completers (disjoint sets).
comp <- which(plan$role == "completer")
taken <- integer()
take <- function(k) { x <- sample(setdiff(comp, taken), k); taken <<- c(taken, x); x }
plan[, `:=`(attention = "pass", flag = "", resubmit = FALSE, reentry = FALSE)]
plan[take(6), attention := "fail_both"]
plan[take(15), attention := "fail_one"]
plan[take(8), flag := "speeder"]
plan[take(6), flag := "longstring"]
plan[take(5), flag := "pna_many"]
plan[take(1), resubmit := TRUE]
plan[take(1), reentry := TRUE]

# Validity order (the 11 shuffled steps) and where each person stops.
VSTEPS <- unname(c(vapply(V, `[[`, "", "slug")[c("ffmq15", "maia2", "bidr16", "phq4", "ius12", "sbi", "bfi2", "eq", "pss4")], "sf2-past-week", "VAS"))
ORDER <- t(replicate(N, sample(VSTEPS)))
plan[, last_step := 14L]                                 # 0 background, 1 short form, 2 candidates, 3-13 validity, 14 debrief
plan[role == "partial_after_sf", last_step := 1L]
plan[role == "partial_after_candidates", last_step := 2L]
plan[role == "no_sf_background_only", last_step := 0L]
plan[role == "no_sf_nothing", last_step := -1L]
mid <- which(plan$role == "partial_mid_validity")
plan[mid, last_step := 2L + sample(1:9, length(mid), TRUE)]
# Three partial people reach exactly one check and fail it; two never reach a check.
checks <- c("sf2-maia2", "sf2-bfi2-om")
for (j in seq_along(mid)) {
  i <- mid[j]; k <- plan$last_step[i] - 2L
  if (j <= 3) {                                           # exactly one check reached, failed
    o <- sample(setdiff(VSTEPS, checks)); ch <- sample(checks, 1); other <- setdiff(checks, ch)
    k <- max(k, 2L); o <- c(o[seq_len(k - 1)], ch, other, o[-seq_len(k - 1)])[1:11]
    ORDER[i, ] <- o; plan[i, `:=`(last_step = 2L + k, attention = "fail_only_reached")]
  } else if (j <= 5) {                                    # never reached a check
    o <- c(sample(setdiff(VSTEPS, checks)), checks)[c(1:2, 10:11, 3:9)]
    ORDER[i, ] <- o; plan[i, `:=`(last_step = 4L, attention = "never_reached")]
  }
}
plan[, reached_end := last_step >= 13L]

# Answers with mess applied.
X <- SI$X; storage.mode(X) <- "character"
for (i in ALL_ITEMS) X[runif(N) < PNA_RATE[[i]], i] <- "pna"
for (r in which(plan$flag == "pna_many")) X[r, sample(D$sf_items, sample(5:9, 1))] <- "pna"
for (r in which(plan$flag == "longstring")) X[r, D$sf_items[4:15]] <- "4"
for (k in names(VI)) if (is.matrix(VI[[k]])) {
  M <- VI[[k]]; storage.mode(M) <- "character"; M[runif(length(M)) < .006] <- "pna"; VI[[k]] <- M
}

# ---- rows in the extract's shape ----------------------------------------------------------------
t0 <- as.POSIXct("2026-11-02 14:00:00", tz = "UTC")
ts <- function(t) if (is.null(t)) NULL else paste0(format(t, "%Y-%m-%d %H:%M:%OS6", tz = "UTC"), "+00")
uuid <- function() paste(vapply(c(8, 4, 4, 4, 12), function(k) paste(sample(c(0:9, letters[1:6]), k, TRUE), collapse = ""), ""), collapse = "-")
row <- function(source, study, e, ...) {
  x <- list(source = source, study = study, participant = e$participant, retest_code = e$retest_code, n_links = e$n_links,
            is_test = e$is_test, prolific_study_id = e$prolific_study_id, repository_consent = e$repository_consent,
            enrollment_status = e$enrollment_status, consent_date = e$consent_date, enrolled_at = e$enrolled_at,
            schedule_id = e$schedule_id, session_status = e$session_status, session_completed_at = e$session_completed_at,
            response_id = NULL, slug = NULL, step_index = NULL, received_at = NULL, resubmission_of = NULL,
            presentation_format = NULL, responses = NULL, value = NULL, entered_at = NULL, exited_at = NULL,
            step_crashes = e$step_crashes)
  extra <- list(...); for (nm in names(extra)) x[nm] <- list(extra[[nm]])
  x
}
mc <- function(id, value = NULL) list(option_id = id, value = value)
num_or_pna <- function(v) { v <- unname(v); if (identical(v, "pna")) "pna" else as.numeric(v) }
timing <- function(total_ms, npages) { cuts <- round(sort(c(0, sort(runif(npages - 1)) * total_ms, total_ms)))
  list(v = 1, unit = "ms since the first page was shown", total_ms = total_ms,
       pages = lapply(seq_len(npages), function(p) list(p - 1, cuts[p], cuts[p + 1]))) }
posting_id <- setNames(names(SF2$postings_sim), SF2$postings_sim)

background <- function(i) {
  g <- plan$practice_group[i]
  regular <- if (runif(1) < .01) "pna" else if (g == "never") mc("no") else mc("yes")
  b <- list(bg_age = if (runif(1) < .01) mc("prefer_not") else mc("age", as.integer(pmin(80, 18 + round(rgamma(1, 4, 0.2))))),
            bg_gender = mc(sample(c("man", "woman", "non_binary", "self_describe", "prefer_not"), 1, prob = c(.45, .52, .02, .005, .005))),
            bg_race = list(mc(sample(c("white", "black", "east_asian", "south_asian", "hispanic_latin", "mixed"), 1, prob = c(.7, .08, .07, .06, .05, .04)))),
            bg_country = { cc <- sample(c("uk", "us", "ca", "au", "ie", "nz", "other"), 1, prob = c(.10, .12, .03, .03, .015, .015, .69))
                           mc(cc, if (cc == "other") sample(c("India", "South Africa", "Poland", "Portugal", "Mexico", "Kenya"), 1) else NULL) },
            bg_education = mc(sample(c("hs", "some_college", "bachelors", "graduate"), 1, prob = c(.15, .25, .38, .22))),
            bg_practice_regular = regular)
  if (identical(regular, "pna") || identical(regular$option_id, "no")) {
    b[c("bg_practice_types", "bg_practice_years", "bg_practice_freq", "bg_practice_training")] <- list("not_applicable")
  } else if (g == "dedicated") {
    b$bg_practice_types <- if (runif(1) < .88) list(mc(sample(c("meditation", "yoga"), 1))) else list(mc("taichi_qigong"))
    b$bg_practice_years <- mc("years", sample(1:20, 1)); b$bg_practice_freq <- mc(sample(c("3_5_week", "daily"), 1))
    b$bg_practice_training <- mc("yes")
  } else {
    repeat {
      # how often other practitioners practise rises with Valuing and Access (H5b's signal)
      w <- exp(.5 * (P$eta[i, "Valuing"] + .5 * P$eta[i, "Access"]) * (0:4)); fr <- sample(names(SF2_FREQ_CODE), 1, prob = w)
      yr <- sample(c(0.5, 1:10), 1); tr <- sample(c("yes", "no"), 1, prob = c(.3, .7))
      if (!(fr %in% c("3_5_week", "daily") && yr >= 1 && tr == "yes")) break
    }
    b$bg_practice_types <- list(mc(sample(c("meditation", "yoga", "contemplative_prayer", "other"), 1)))
    b$bg_practice_years <- mc("years", yr); b$bg_practice_freq <- mc(fr); b$bg_practice_training <- mc(tr)
  }
  sp <- 3 + .3 * P$eta[i, "Valuing"] + (g == "dedicated") + rnorm(1)
  b$bg_spiritual <- if (runif(1) < .02) "pna" else as.integer(min(5, max(1, round(sp))))
  b$bg_religious <- if (runif(1) < .02) "pna" else as.integer(min(5, max(1, round(.5 * sp + 1 + rnorm(1)))))
  b$bg_tradition <- mc(sample(c("none", "christian", "buddhist", "sbnr", "other"), 1))
  b$`_timing` <- timing(round(runif(1, 40e3, 120e3)), 3)
  b
}

sf_answers <- function(i, Xrow) {
  r <- lapply(D$sf_items, function(it) num_or_pna(Xrow[[it]])); names(r) <- D$sf_items
  secs <- if (plan$flag[i] == "speeder") runif(1, 1, 1.8) else runif(1, 3, 9)
  r$`_timing` <- timing(round(secs * 20 * 1000), 3); r
}
cand_answers <- function(Xrow) {
  r <- lapply(D$cand_items, function(it) num_or_pna(Xrow[[it]])); names(r) <- D$cand_items
  r$`_order` <- list(method = "interleave", version = 1, positions = as.list(setNames(sample(24), D$cand_items)))
  r$`_timing` <- timing(round(runif(1, 60e3, 200e3)), 3); r
}
validity_answers <- function(i, slug) {
  k <- names(V)[vapply(V, `[[`, "", "slug") == slug]
  if (slug == "sf2-past-week") {
    pw <- D$past_week; r <- list()
    r$crit_pause_days <- if (runif(1) < .006) "pna" else mc(names(pw$crit_pause_days)[match(VI$crit_pause_days[i], pw$crit_pause_days)])
    r$crit_notice_first <- if (runif(1) < .006) "pna" else if (runif(1) < .08) mc("nothing_unexpected") else
      mc(names(pw$crit_notice_first)[match(VI$crit_notice_first[i], pw$crit_notice_first)])
    r$crit_outdoors <- if (runif(1) < .006) "pna" else mc(names(pw$crit_outdoors)[match(VI$crit_outdoors[i], pw$crit_outdoors)])
    r$`_timing` <- timing(round(runif(1, 15e3, 40e3)), 1)
    return(r)
  }
  M <- VI[[k]]; r <- lapply(colnames(M), function(it) num_or_pna(M[i, it])); names(r) <- colnames(M)
  if (!is.null(V[[k]]$check)) {
    correct <- SF2$checks[[V[[k]]$check]]
    fail <- plan$attention[i] == "fail_both" || (plan$attention[i] %in% c("fail_one", "fail_only_reached") && slug == plan$fail_slug[i])
    r[[V[[k]]$check]] <- if (fail) sample(setdiff(seq(V[[k]]$range[1], V[[k]]$range[2]), correct), 1) else correct
  }
  r$`_timing` <- timing(round(runif(1, 20e3, 150e3)), 2)
  r
}
# Which check a one-check failure fails: for completers, either; for partial people, the one they reached.
plan[, fail_slug := NA_character_]
for (i in which(plan$attention == "fail_one")) plan$fail_slug[i] <- sample(checks, 1)
for (i in which(plan$attention == "fail_only_reached")) plan$fail_slug[i] <- intersect(ORDER[i, seq_len(plan$last_step[i] - 2L)], checks)

rows <- list(); add <- function(x) rows[[length(rows) + 1]] <<- x
main_done <- rep(as.POSIXct(NA), N)
for (i in seq_len(N)) {
  start <- t0 + runif(1, 0, 10 * 86400)
  e <- list(participant = plan$participant[i], retest_code = NULL, n_links = NULL, is_test = plan$role[i] == "staff_test",
            prolific_study_id = posting_id[[plan$posting[i]]], repository_consent = runif(1) < .8, enrollment_status = "active",
            consent_date = ts(start), enrolled_at = ts(start - 30), schedule_id = uuid(),
            session_status = if (plan$last_step[i] >= 14) "completed" else "unlocked", session_completed_at = NULL, step_crashes = 0L)
  steps <- unname(c("sf-background", "sf-short-20", "sf2-sfs-candidates", ORDER[i, ], "debrief"))
  t <- start; recs <- list()
  for (s in seq_len(plan$last_step[i] + 1L) - 1L) {
    entered <- t; t <- t + runif(1, 20, 150); exited <- t
    slug <- steps[s + 1]
    sub <- if (slug == "VAS") "vas_life_satisfaction" else if (slug == "debrief") "debrief" else slug
    recs[[length(recs) + 1]] <- list(kind = "step", s = s, sub = sub, entered = entered, exited = exited)
    if (slug == "debrief") next
    if (slug == "VAS") { recs[[length(recs) + 1]] <- list(kind = "vas", s = s, at = exited, value = VI$lifesat[i]); next }
    resp <- switch(slug, `sf-background` = background(i), `sf-short-20` = sf_answers(i, X[i, ]),
                   `sf2-sfs-candidates` = cand_answers(X[i, ]), validity_answers(i, slug))
    recs[[length(recs) + 1]] <- list(kind = "q", s = s, slug = slug, at = exited, resp = resp)
  }
  if (plan$last_step[i] >= 14) e$session_completed_at <- ts(t)
  if (plan$last_step[i] >= 1) main_done[i] <- recs[[which(vapply(recs, function(r) identical(r$slug, "sf-short-20"), TRUE))[1]]]$at
  stopifnot(plan$last_step[i] < 1 || !is.na(main_done[i]))
  add(row("enrollment", "main", e))
  for (r in recs) {
    if (r$kind == "step") add(row("step", "main", e, slug = r$sub, step_index = r$s, entered_at = ts(r$entered), exited_at = ts(r$exited)))
    if (r$kind == "vas") add(row("vas", "main", e, response_id = uuid(), slug = "life_satisfaction", step_index = r$s,
                                 received_at = ts(r$at), value = as.integer(r$value)))
    if (r$kind == "q") {
      id <- uuid()
      add(row("questionnaire", "main", e, response_id = id, slug = r$slug, step_index = r$s, received_at = ts(r$at),
              presentation_format = "stacked", responses = r$resp))
      if (r$slug == "sf-short-20" && plan$resubmit[i])            # byte-identical copy 2 s later: flagged by the platform
        add(row("questionnaire", "main", e, response_id = uuid(), slug = r$slug, step_index = r$s, received_at = ts(r$at + 2),
                resubmission_of = id, presentation_format = "stacked", responses = r$resp))
      if (r$slug == "sf2-sfs-candidates" && plan$reentry[i]) {   # re-entered 10 minutes later with different answers
        Xn <- X[i, ]; Xn[D$cand_items] <- as.character(sample(1:6, 24, TRUE))
        add(row("questionnaire", "main", e, response_id = uuid(), slug = r$slug, step_index = r$s, received_at = ts(r$at + 600),
                presentation_format = "stacked", responses = cand_answers(Xn)))
      }
    }
  }
}

# ---- the retest -----------------------------------------------------------------------------
eligible <- which(plan$role == "completer" & plan$attention != "fail_both")
invited <- sample(eligible, 400)          # Section 2: about 400 invited for about 300 completions
resp <- sample(invited, 300)
days <- c(runif(278, 14, 21.9), runif(16, 22, 25.9), runif(6, 12, 13.9))[sample(300)]
eta2 <- sweep(P$eta, 2, RHO, `*`) + sweep(mvrnorm(N, rep(0, 5), Phi), 2, sqrt(1 - RHO^2), `*`)
ri2 <- .8 * P$ri + .6 * rnorm(N, 0, sqrt(RI_VAR))
X2 <- gen_sfc(P, eta2, ri2)$X; storage.mode(X2) <- "character"
for (it in D$sf_items) X2[runif(N) < PNA_RATE[[it]], it] <- "pna"
PW2 <- gen_pastweek2(P, eta2, ri2)
unlinked <- resp[1]
rt <- data.table(i = resp, days = days, retest_code = code(300, "P-"))
rt <- rbind(rt, data.table(i = NA_integer_, days = 15, retest_code = sprintf("test-%08d", sample(1e7, 1))))
for (j in seq_len(nrow(rt))) {
  i <- rt$i[j]; is_test <- is.na(i)
  at <- if (is_test) t0 + 20 * 86400 else main_done[i] + rt$days[j] * 86400
  e <- list(participant = if (is_test || identical(i, unlinked)) NULL else plan$participant[i], retest_code = rt$retest_code[j],
            n_links = if (is_test || identical(i, unlinked)) 0L else 1L, is_test = is_test, prolific_study_id = "sim-retest-posting",
            repository_consent = NULL, enrollment_status = "active", consent_date = ts(at - 60), enrolled_at = ts(at - 90),
            schedule_id = uuid(), session_status = "completed", session_completed_at = ts(at + 30), step_crashes = 0L)
  Xrow <- if (is_test) setNames(as.character(sample(1:6, 20, TRUE)), D$sf_items) else X2[i, ]
  r <- lapply(D$sf_items, function(it) num_or_pna(Xrow[[it]])); names(r) <- D$sf_items
  r$`_timing` <- timing(round(runif(1, 60e3, 160e3)), 3)
  add(row("enrollment", "retest", e))
  add(row("step", "retest", e, slug = "sf-short-20", step_index = 0L, entered_at = ts(at - 100), exited_at = ts(at)))
  add(row("questionnaire", "retest", e, response_id = uuid(), slug = "sf-short-20", step_index = 0L, received_at = ts(at),
          presentation_format = "stacked", responses = r))
  pw <- PW2[if (is_test) 1 else i, ]; dpw <- D$past_week
  q2 <- lapply(names(pw), function(v) if (runif(1) < .006) "pna" else if (v == "crit_notice_first" && runif(1) < .08) mc("nothing_unexpected")
               else mc(names(dpw[[v]])[match(pw[[v]], dpw[[v]])]))
  names(q2) <- names(pw); q2$`_timing` <- timing(round(runif(1, 15e3, 40e3)), 1)
  add(row("step", "retest", e, slug = "sf2-past-week", step_index = 1L, entered_at = ts(at), exited_at = ts(at + 30)))
  add(row("questionnaire", "retest", e, response_id = uuid(), slug = "sf2-past-week", step_index = 1L, received_at = ts(at + 30),
          presentation_format = "stacked", responses = q2))
  add(row("step", "retest", e, slug = "debrief", step_index = 2L, entered_at = ts(at + 30), exited_at = ts(at + 60)))
}

out <- list(boundary = "simulated", rows = rows,
            warning = "SIMULATED DATA (study2/analysis/00_simulate.R): no participant answered any of this.")
dir.create(sf2_path("raw"), recursive = TRUE, showWarnings = FALSE)
write_json(out, sf2_path("raw", "simulated_extract.json"), auto_unbox = TRUE, null = "null", na = "null", digits = NA)
fwrite(data.table(participant = plan$participant[invited]), sf2_path("raw", "simulated_retest_invited.csv"))

truth$planted <- plan[, .(participant, role, posting, practice_group, attention, fail_slug, flag, resubmit, reentry, last_step, reached_end)]
truth$retest <- data.table(participant = plan$participant[rt$i], retest_code = rt$retest_code, days = rt$days,
                           linked = !is.na(rt$i) & !(rt$i %in% unlinked), is_test = is.na(rt$i))
truth$invited <- plan$participant[invited]
saveRDS(truth, sf2_derived("simulation_truth"))
sf2_log("simulated extraction: ", length(rows), " rows, ", N, " main enrollments, ", nrow(rt), " retest enrollments -> ",
        sf2_path("raw", "simulated_extract.json"))
print(plan[, .N, by = .(role)]); print(plan[role == "completer", .N, by = practice_group])
