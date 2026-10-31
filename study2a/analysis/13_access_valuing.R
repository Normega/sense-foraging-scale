# 13: exploratory, not preregistered (2026-10-31). Does Study 1's Access / Valuing split
# (the trial's AV10, latent r .60 in Study 1, one factor CFI .86) hold in Study 2a, and how
# does the AV10 score compare with 12's strongest 8-statement set? Four AV10 statements are
# in the scale in Study 1 wording (sf2a_01, 02, 03, 27); six come at the end.
# Run from the repo root with SF2A_MODE=real, after 12.
source("study2a/analysis/config.R"); source("study2a/analysis/sf2a_models.R")
A <- readRDS(sf2a_derived("sf2a")); def <- A$definition
X <- cbind(A$items_raw[, c("sf2a_01", "sf2a_02", "sf2a_03", "sf2a_27")], A$end)
colnames(X)[1:4] <- c("sf6_9", "sf6_36", "sf6_35", "sf6_44")
acc <- c("sf6_35", "sf6_36", "sf6_9", "sf6_8", "sf6_7"); val <- c("sf6_12", "sf6_18", "sf6_11", "sf6_27", "sf6_44")
it <- c(acc, val); X <- as.data.frame(X[, it])
ri <- c(paste0("RI =~ ", paste0("1*", it, collapse = " + ")), "RI ~~ RI")
m2 <- c(paste0("Acc =~ NA*", paste(acc, collapse = " + ")), paste0("Val =~ NA*", paste(val, collapse = " + ")),
        "Acc ~~ 1*Acc", "Val ~~ 1*Val", "Acc ~~ Val", ri, "RI ~~ 0*Acc + 0*Val")
m1 <- c(paste0("SF =~ NA*", paste(it, collapse = " + ")), "SF ~~ 1*SF", ri, "RI ~~ 0*SF")
f2 <- sf2a_sem(m2, X, ordered = it); f1 <- sf2a_sem(m1, X, ordered = it)
fs <- function(f) round(fitMeasures(f, c("chisq.scaled", "df", "cfi.scaled", "rmsea.scaled", "srmr")), 3)
print(rbind(two = fs(f2), one = fs(f1)))
pe <- parameterEstimates(f2, standardized = TRUE)
print(pe[pe$lhs == "Acc" & pe$rhs == "Val", c("est", "ci.lower", "ci.upper", "std.all")])
print(round(lavInspect(f2, "std")$lambda, 2))
print(lavTestLRT(f1, f2))
# observed subscale correlation, keyed
K <- X; for (v in c("sf6_35", "sf6_36", "sf6_44")) K[[v]] <- 7 - K[[v]]
sa <- rowMeans(K[, acc], na.rm = TRUE); sv <- rowMeans(K[, val], na.rm = TRUE)
cat("observed r Access-Valuing:", round(cor(sa, sv, use = "pair"), 3), "\n")
cat("alpha Access", round(psych::alpha(K[, acc], warnings = FALSE)$total$raw_alpha, 3),
    "Valuing", round(psych::alpha(K[, val], warnings = FALSE)$total$raw_alpha, 3), "\n")
# against the best-8
best <- fread(sf2a_out("12_strongest_scale.csv"))
s8 <- strsplit(best[n == 8 & min_reversed == 2, ids], " ")[[1]]
if (!length(s8)) s8 <- strsplit(best[n == 8][1, ids], " ")[[1]]
print(best[n == 8, .(min_reversed, n_reversed, omega_content, statements)])
sc8 <- rowMeans(A$items_keyed[, s8], na.rm = TRUE)
cat("r best8 with AV10", round(cor(sc8, (sa + sv) / 2, use = "pair"), 3), " with Access", round(cor(sc8, sa, use = "pair"), 3),
    " with Valuing", round(cor(sc8, sv, use = "pair"), 3), "\n")
grp <- A$meta$practice_group
dval <- function(y) { a <- y[grp %in% "current"]; b <- y[grp %in% "never"]
  (mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE)) / sqrt((var(a, na.rm = TRUE) + var(b, na.rm = TRUE)) / 2) }
cat("d Access", round(dval(sa), 2), "Valuing", round(dval(sv), 2), "AV10", round(dval((sa + sv) / 2), 2), "best8", round(dval(sc8), 2), "\n")
cat("r IUS Access", round(cor(sa, A$meta$ius_total, use = "pair"), 2), "Valuing", round(cor(sv, A$meta$ius_total, use = "pair"), 2), "\n")
cat("r self-rating Access", round(cor(sa, A$meta$vig_self, use = "pair"), 2), "Valuing", round(cor(sv, A$meta$vig_self, use = "pair"), 2),
    "best8", round(cor(sc8, A$meta$vig_self, use = "pair"), 2), "\n")
