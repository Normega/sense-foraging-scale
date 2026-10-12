# Review check 5: (a) does 11_followup's "agreement index" track the model's random intercept?
# (b) known-groups differences on the primary model's factor scores, with and without
# spirituality/religiosity as covariates (exploratory); (c) SEs of the primary ESEM loadings, to
# project precision at n = 250 (Study 2's selection half).
suppressPackageStartupMessages({library(lavaan); library(data.table); library(jsonlite)})
A <- readRDS("F:/radlab_data/sf_study1/derived/sf1_real.rds")
P <- readRDS("F:/radlab_data/sf_study1/derived/primary_real_k04.rds")
m <- P$with_ri
FS <- as.data.table(lavPredict(m)); idx <- lavInspect(m, "case.idx")
meta <- as.data.table(A$meta)[idx]
raw <- A$items_raw[idx, ]; pairs <- A$definition$pairs
acq <- rowMeans(raw[, c(pairs[, 1], pairs[, 2])], na.rm = TRUE)
allmean <- rowMeans(raw, na.rm = TRUE)
cat("(a) r(agreement index, RI factor score) =", round(cor(acq, FS$RI), 3), "\n")
cat("    r(mean of all 48 raw items, RI factor score) =", round(cor(allmean, FS$RI), 3), "\n")
cat("    r(agreement index, content factor scores F1..F4) =", round(cor(acq, as.matrix(FS[, .(F1, F2, F3, F4)])), 3), "\n")
# spirituality / religiosity from the raw extraction (same rules as 11)
raw_file <- tail(sort(list.files("F:/radlab_data/sf_study1/raw", "^extract_.*[.]json$", full.names = TRUE)), 1)
rows <- fromJSON(raw_file, simplifyVector = FALSE)$rows
bg <- Filter(function(x) identical(x$slug, "sf-background") && is.null(x$resubmission_of), rows)
bg <- bg[order(vapply(bg, `[[`, "", "received_at"))]; bg <- bg[!duplicated(vapply(bg, `[[`, "", "participant"))]
num1 <- function(v) if (is.numeric(v)) v else NA_real_
B <- rbindlist(lapply(bg, function(x) data.table(participant = x$participant, spiritual = num1(x$responses$bg_spiritual), religious = num1(x$responses$bg_religious),
   types = paste(unlist(lapply(x$responses$bg_practice_types, function(o) if (is.list(o)) o$option_id else o)), collapse = ","))))
d <- cbind(meta[, .(participant, practice_group, posting)], FS)
d <- merge(d, B, by = "participant", all.x = TRUE, sort = FALSE)
orient <- c(1, 1, 1, -1)
for (j in 1:4) set(d, j = paste0("F", j), value = d[[paste0("F", j)]] * orient[j])
dd <- d[practice_group %in% c("current", "never")]; dd[, cur := as.integer(practice_group == "current")]
cat("\n(b) current vs never, standardized difference on factor scores (pooled within-group SD)\n")
for (f in paste0("F", 1:4)) {
  sdp <- sqrt(mean(tapply(dd[[f]], dd$cur, var)))
  b0 <- coef(lm(reformulate(c("cur", "posting"), f), dd))["cur"]
  b1 <- coef(lm(reformulate(c("cur", "posting", "spiritual", "religious"), f), dd))["cur"]
  cat(f, " d (posting) =", round(b0 / sdp, 2), "  d (+ spirituality, religiosity) =", round(b1 / sdp, 2), "\n")
}
cat("\n    practice types among 'current' (share naming each):\n")
tt <- unlist(strsplit(dd[cur == 1, types], ",")); print(round(sort(table(tt), decreasing = TRUE) / sum(dd$cur), 2))
cat("\n(c) SE of primary-model standardized content loadings (N =", length(idx), ")\n")
s <- standardizedSolution(m); s <- s[s$op == "=~" & s$lhs %in% P$fnames, ]
cat("    median SE", round(median(s$se), 3), "; projected at n = 250:", round(median(s$se) * sqrt(length(idx) / 250), 3), "\n")
