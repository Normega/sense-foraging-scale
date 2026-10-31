# 06: Letting be against wording (Section 5, step 5; decision L).
#
# The primary model plus a wording factor W: loadings free on all 19 reverse-keyed
# statements, uncorrelated with the content and agreement factors. If that model
# does not converge or gives an improper solution, W replaces the random
# intercept, and the reason is reported.
# Decision L: Letting be stands if, in this model, at least three of its 7 forward
# statements load >= .40 on the Letting be factor. Otherwise it is reported as a
# wording factor, and Study 2b keeps it only as an exploratory score.
# PREREG: with k != 5 the Letting be factor is the one 05 matched to Letting be
# (most of its defining statements were written for it); if none matches, L
# cannot stand.

source("study2a/analysis/config.R")
source("study2a/analysis/sf2a_models.R")
A <- readRDS(sf2a_derived("sf2a"))
def <- A$definition; X <- A$items_raw; k <- SF2A$k_primary
P <- readRDS(sf2a_derived(sprintf("primary_k%02d", k)))
member <- P$member; fn <- P$fnames

sf2a_log("primary model + wording factor")
m_w <- tryCatch(sf2a_fit_esem(X, member, ri = TRUE, wording = def$reversed), error = function(e) e)
variant <- "primary + wording factor"; reason <- ""
if (sf2a_failed(m_w) || sf2a_improper(m_w)) {
  reason <- if (sf2a_failed(m_w)) sf2a_msg(m_w) else "improper solution"
  sf2a_log("  ", reason, "; the wording factor replaces the random intercept")
  m_w <- sf2a_fit_esem(X, member, ri = FALSE, wording = def$reversed)
  variant <- "wording factor in place of the random intercept"
}
O <- sf2a_orient(sf2a_loadings(m_w, fn, def$items), member, def$reversed)
Kw <- O$keyed

lb <- names(P$labels)[which(P$labels == "Letting be")]
lb <- if (k == length(def$factors)) sf2a_fname("Letting be") else if (length(lb)) lb[1] else NA_character_
fwd <- def$items[def$factor[def$items] == "Letting be" & !(def$items %in% def$reversed)]
own <- if (!is.na(lb)) Kw[fwd, lb] else rep(NA_real_, length(fwd))
n40 <- sum(own >= SF2A$load_min, na.rm = TRUE)
dec <- data.table(model = variant, why_variant = reason, converged = !sf2a_failed(m_w), improper = sf2a_improper(m_w),
                  letting_be_factor = lb, n_forward = length(fwd), n_forward_ge_40 = n40,
                  letting_be_stands = !is.na(lb) && n40 >= 3, as.data.table(sf2a_fit_stats(m_w)))
fwrite(dec, sf2a_out(sprintf("06_decision_L_k%02d.csv", k)))

Wl <- lavInspect(m_w, "std")$lambda[def$items, "W"]
tab <- data.table(item = def$items, source = def$source[def$items], hypothesized = def$factor[def$items],
                  keying = ifelse(def$items %in% def$reversed, "reverse", "forward"), round(Kw, 3),
                  W_std = round(Wl, 3),
                  letting_be_without_W = if (!is.na(lb)) round(P$keyed[def$items, lb], 3) else NA_real_)
fwrite(tab, sf2a_out(sprintf("06_wording_loadings_k%02d.csv", k)))
if (!is.null(sf2a_ri(m_w))) fwrite(sf2a_ri(m_w), sf2a_out(sprintf("06_agreement_factor_with_W_k%02d.csv", k)))
saveRDS(list(fit = m_w, variant = variant, keyed = Kw, letting_be = lb), sf2a_derived(sprintf("wording_k%02d", k)))

sf2a_log("decision L: ", n40, " of ", length(fwd), " forward Letting be statements >= .40 on ", lb,
         " (", variant, "); stands = ", dec$letting_be_stands)
if (!is.na(lb)) print(tab[item %in% fwd, .(item, source, own = get(lb), letting_be_without_W)])
