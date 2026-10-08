# =====================================================================================
# rq3_run_all.R
# Run this ONCE on a laptop that has the raw ABS files. It writes small AGGREGATE tables and
# figures (no respondent rows) to the folder rq3_results/. The shared Rmd only reads those,
# so teammates can knit without the raw data.
#
#   PART A  primary models   (fibre x activity; BMI, waist, obesity)   -> Rosetta's task
#   PART B  robustness       (outliers, multiple imputation)            -> Amelie's task
#   PART C  four-group obesity rates (fibre high/low x activity)        -> needs the fibre cut-off
#
# Before running:
#   1. RQ3_RAW_DATA_DIR set in ~/.Renviron (folder holding AHSnpa11bp.csv, inp13bsp.csv, inp13bhh.csv)
#   2. Working directory = project folder with the two cleaned workbooks, rq3_spec.R and
#      rq3-rosetta-results.csv  (open the .Rproj, or setwd())
#   3. Source this file in RStudio, or:  Rscript rq3_run_all.R
# Re-run whenever rq3_spec.R changes (confounders, obesity definition, cut-offs).
# =====================================================================================

suppressPackageStartupMessages({library(mice); library(robustbase)})
source("rq3_spec.R")

m_imputations <- 10            # number of imputed datasets
out_dir <- "rq3_results"
dir.create(out_dir, showWarnings = FALSE)
groups <- SPEC$groups
outcomes <- names(outcome_info)

inputs <- setNames(lapply(groups, load_inputs), groups)
frames <- setNames(lapply(groups, function(g) prep_cohort(inputs[[g]], g)), groups)

# ---- 0. Does the data pipeline still reproduce Rosetta's published primary model? ----
her_file <- "rq3-rosetta-results.csv"
if (file.exists(her_file)) {
  her <- subset(read.csv(her_file, stringsAsFactors = FALSE),
                analysis_id == "primary_two_day" & grepl("ratio of odds ratios", estimate))
  chk <- do.call(rbind, lapply(groups, function(g) {
    z <- prep_cohort(inputs[[g]], g, need_waist = FALSE)
    d <- z$data[complete.cases(z$data[, all.vars(legacy_formula)]), ]
    fit <- svyglm(legacy_formula, make_design(d, z$reps), family = quasibinomial())
    h <- her[her$group == g, ]
    data.frame(cohort = g, n_here = nrow(d), n_rosetta = h$n,
               or_here = exp(coef(fit)[["fibre5:active"]]), or_rosetta = h$odds_ratio)
  }))
  print(chk)
  if (!all(chk$n_here == chk$n_rosetta) || max(abs(chk$or_here - chk$or_rosetta)) > 0.001)
    stop("The data pipeline does not reproduce Rosetta's primary model. Fix this before using any output.")
  write.csv(chk, file.path(out_dir, "reproduce-check.csv"), row.names = FALSE)
} else message("rq3-rosetta-results.csv not found: skipping the reproduce check")

# ---- sample flow ----
flow <- do.call(rbind, lapply(groups, function(g) {
  d <- frames[[g]]$data; cc <- complete_rows(d, "obese")
  data.frame(cohort = g, eligible = frames[[g]]$n_eligible, known = nrow(d), complete = nrow(cc),
             cases = sum(cc$obese), mean_bmi = weighted.mean(cc$bmi, cc$mainwt),
             mean_waist = weighted.mean(cc$waist, cc$mainwt),
             obesity_pct = 100 * weighted.mean(cc$obese, cc$mainwt))
}))
write.csv(flow, file.path(out_dir, "sample-flow.csv"), row.names = FALSE)

# =====================================================================================
# PART A: PRIMARY MODELS (survey-weighted; replicate-weight jackknife as in Rosetta's model)
# =====================================================================================
message("PART A: primary models")
weighted_auc <- function(y, p, w) {
  o <- order(p[y == 0]); pn <- p[y == 0][o]; wn <- w[y == 0][o]; cw <- c(0, cumsum(wn))
  pp <- p[y == 1]; wp <- w[y == 1]
  lt <- cw[findInterval(pp, pn, left.open = TRUE) + 1]; le <- cw[findInterval(pp, pn) + 1]
  sum(wp * (lt + 0.5 * (le - lt))) / (sum(wp) * sum(wn))
}
calib_table <- function(y, p, w, g = 10) {
  grp <- cut(p, quantile(p, seq(0, 1, length.out = g + 1)), include.lowest = TRUE, labels = FALSE)
  data.frame(decile = 1:g, predicted = as.numeric(tapply(p * w, grp, sum) / tapply(w, grp, sum)),
             observed = as.numeric(tapply(y * w, grp, sum) / tapply(w, grp, sum)))
}

primary <- list(); key_tab <- list(); coef_tab <- list(); perf_tab <- list(); calib <- list(); vif_tab <- list()
for (g in groups) {
  z <- frames[[g]]
  for (y in outcomes) {
    dcc <- complete_rows(z$data, y); logistic <- outcome_info[[y]]$family == "binomial"
    fit <- fit_svy(y, dcc, z$reps); primary[[paste(g, y)]] <- list(fit = fit, dcc = dcc, reps = z$reps)
    ke <- key_effects(fit, logistic)
    key_tab[[paste(g, y)]] <- cbind(cohort = g, outcome = y, n = nrow(dcc), ke)
    b <- coef(fit); se <- sqrt(diag(vcov(fit))); f <- if (logistic) exp else identity
    coef_tab[[paste(g, y)]] <- data.frame(cohort = g, outcome = y, term = pretty_term(names(b)),
      estimate = f(b), lower = f(b - 1.96 * se), upper = f(b + 1.96 * se), p = 2 * pnorm(-abs(b / se)), row.names = NULL)
    if (logistic) {
      p <- as.numeric(predict(fit, type = "response")); pbar <- weighted.mean(dcc$obese, dcc$mainwt)
      perf_tab[[g]] <- data.frame(cohort = g, n = nrow(dcc), auc = weighted_auc(dcc$obese, p, dcc$mainwt),
        brier = weighted.mean((dcc$obese - p)^2, dcc$mainwt), brier_null = pbar * (1 - pbar),
        condition_number = kappa(model.matrix(make_formula(y), dcc), exact = TRUE))
      calib[[g]] <- cbind(cohort = g, calib_table(dcc$obese, p, dcc$mainwt))
    }
  }
  # multicollinearity (variance inflation factors on the main-effects design, unweighted)
  d0 <- complete_rows(z$data, "obese")
  X <- model.matrix(reformulate(c("fibre5", "active", confounders)), d0)[, -1, drop = FALSE]
  vif_tab[[g]] <- data.frame(cohort = g, term = pretty_term(colnames(X)), vif = diag(solve(cor(X))), row.names = NULL)
}
write.csv(do.call(rbind, key_tab),  file.path(out_dir, "primary-key-effects.csv"), row.names = FALSE)
write.csv(do.call(rbind, coef_tab), file.path(out_dir, "primary-coefficients.csv"), row.names = FALSE)
write.csv(do.call(rbind, perf_tab), file.path(out_dir, "performance.csv"), row.names = FALSE)
write.csv(do.call(rbind, calib),    file.path(out_dir, "calibration.csv"), row.names = FALSE)
write.csv(do.call(rbind, vif_tab),  file.path(out_dir, "vif.csv"), row.names = FALSE)

# residual diagnostics for the two continuous outcomes (one figure per cohort)
for (k in seq_along(groups)) {
  g <- groups[k]
  png(file.path(out_dir, paste0("diagnostics-", k, ".png")), width = 1500, height = 1000, res = 150)
  par(mfrow = c(2, 3), mar = c(4.2, 4.2, 3, 1))
  for (y in c("bmi", "waist")) {
    pr <- primary[[paste(g, y)]]; r <- as.numeric(resid(pr$fit)); fv <- as.numeric(fitted(pr$fit))
    lab <- c(bmi = "BMI", waist = "Waist")[[y]]
    plot(fv, r, pch = 16, cex = 0.4, col = rgb(0, 0, 0, 0.35), xlab = "Fitted value", ylab = "Residual",
         main = paste0(lab, ": residuals vs fitted")); abline(h = 0, lty = 2); lines(lowess(fv, r), col = "red", lwd = 2)
    qqnorm(r / sd(r), pch = 16, cex = 0.4, main = paste0(lab, ": normal Q-Q")); qqline(r / sd(r), col = "red")
    lmf <- lm(make_formula(y), data = pr$dcc)
    plot(cooks.distance(lmf), type = "h", xlab = "Observation", ylab = "Cook's distance",
         main = paste0(lab, ": influence (unweighted)")); abline(h = 4 / nrow(pr$dcc), lty = 2, col = "red")
  }
  mtext(g, outer = TRUE, line = -1.2, cex = 0.9); dev.off()
}

# =====================================================================================
# PART B: ROBUSTNESS (outliers, multiple imputation)
# =====================================================================================
message("PART B: robustness checks (a few minutes)")
outlier_checks <- function(y, dcc, reps) {
  logistic <- outcome_info[[y]]$family == "binomial"; f <- make_formula(y); n <- nrow(dcc)
  fit0 <- fit_svy(y, dcc, reps)
  g0  <- if (logistic) glm(f, data = dcc, family = binomial()) else lm(f, data = dcc)
  gr  <- tryCatch(suppressWarnings(if (logistic) glmrob(f, data = dcc, family = binomial(), weights.on.x = "hat")
                  else lmrob(f, data = dcc)), error = function(e) NULL)
  cd   <- cooks.distance(g0); keep <- cd <= quantile(cd, 0.99)          # drop the 1% most influential
  dtrim <- dcc[keep, ]
  fit_inf <- fit_svy(y, dtrim, reps)
  rbind(int_from_fit(fit0, "Weighted, complete case (primary model)", n, logistic),
        int_from_fit(fit_inf, "Weighted, most influential 1% removed (Cook's distance)", sum(keep), logistic),
        int_from_summary(g0, "Unweighted standard model", n, logistic),
        if (!is.null(gr)) int_from_summary(gr, "Unweighted robust model (robustbase)", n, logistic)
        else data.frame(analysis = "Unweighted robust model (robustbase): did not converge", n = n,
                        estimate = NA, lower = NA, upper = NA, p = NA))
}

mi_run <- function(z, m = m_imputations, seed = 123) {
  d <- z$data
  imp_dat <- d[, c("bmi", "waist", "obese", "fibre5", "active", "age_cat", "female", "energy_mj",
                   "income", "year12", "smoking", "alcohol", "protein10")]
  imp_dat$fibre5_x_active <- imp_dat$fibre5 * imp_dat$active   # keeps the interaction in the imputation model
  imp_dat$logwt <- log(d$mainwt)                                # design information used as a predictor
  imp_dat$year12 <- factor(imp_dat$year12)
  imp <- mice(imp_dat, m = m, maxit = 10, seed = seed, printFlag = FALSE)
  other <- d[, c("mainwt", z$reps), drop = FALSE]
  fits <- lapply(seq_len(m), function(i) {
    di <- complete(imp, i)[, c("bmi", "waist", "obese", "fibre5", "active", "age_cat", "female", "energy_mj",
                               "income", "year12", "smoking", "alcohol", "protein10")]
    di$year12 <- as.integer(as.character(di$year12))
    dd <- cbind(di, other)
    setNames(lapply(outcomes, function(y) fit_svy(y, dd, z$reps)), outcomes)
  })
  list(imp = imp, fits = fits)
}

rob_rows <- list(); mi_store <- list(); missing_tab <- list()
for (g in groups) {
  z <- frames[[g]]
  mi <- mi_run(z); mi_store[[g]] <- mi
  mv <-c("age_cat", "female", "energy_mj", "income", "year12", "smoking", "alcohol", "protein10")
  missing_tab[[g]] <- data.frame(cohort = g, variable = mv, n_missing = colSums(is.na(z$data[, mv])),
                                 pct_missing = round(100 * colMeans(is.na(z$data[, mv])), 1), row.names = NULL)
  for (y in outcomes) {
    message("  ", g, " / ", y)
    logistic <- outcome_info[[y]]$family == "binomial"
    oc <- outlier_checks(y, complete_rows(z$data, y), z$reps)
    p <- pool_contrast(lapply(mi$fits, `[[`, y))
    mrow <- est_row(paste0("Multiple imputation (m = ", m_imputations, ")"), p[["est"]], p[["se"]], nrow(z$data), logistic)
    rob_rows[[paste(g, y)]] <- cbind(cohort = g, outcome = y, rbind(oc, mrow))
  }
}
write.csv(do.call(rbind, rob_rows), file.path(out_dir, "robustness.csv"), row.names = FALSE)
write.csv(do.call(rbind, missing_tab), file.path(out_dir, "missing-covariates.csv"), row.names = FALSE)
png(file.path(out_dir, "imputation-chains.png"), width = 1400, height = 1000, res = 150); print(plot(mi_store[[1]]$imp)); dev.off()
png(file.path(out_dir, "imputation-density.png"), width = 1400, height = 1000, res = 150); print(densityplot(mi_store[[1]]$imp)); dev.off()

# =====================================================================================
# PART C: OBESITY RATES BY FIBRE (high/low) x ACTIVITY (>=150 min/week or not)
# The fibre cut-off is still undecided. The main models above use continuous fibre and do NOT
# depend on it. Pick the rule BEFORE looking at the group results, then keep the others as
# sensitivity analyses. Edit SPEC$fibre_rule_primary in rq3_spec.R to change the main rule.
# =====================================================================================
message("PART C: four-group obesity rates")
wtd_median <- function(x, w) { o <- order(x); cw <- cumsum(w[o]) / sum(w); x[o][which(cw >= 0.5)[1]] }
fibre_rules <- list(
  "NHMRC adequate intake (men 30 g, women 25 g)" = function(d) as.integer(d$fibre >= ifelse(d$female == 1, 25, 30)),
  "Cohort median, by sex (weighted)" = function(d) {
    med <- c(`0` = wtd_median(d$fibre[d$female == 0], d$mainwt[d$female == 0]),
             `1` = wtd_median(d$fibre[d$female == 1], d$mainwt[d$female == 1]))
    as.integer(d$fibre >= med[as.character(d$female)])
  },
  "Fixed 25 g/day for everyone" = function(d) as.integer(d$fibre >= 25)
)
grp_levels <- c("Low fibre, low activity", "High fibre, low activity", "Low fibre, high activity", "High fibre, high activity")
stopifnot(SPEC$fibre_rule_primary %in% names(fibre_rules))

group_analysis <- function(z, rule_name) {
  dcc <- complete_rows(z$data, "obese")
  high <- fibre_rules[[rule_name]](dcc)
  dcc$grp <- factor(grp_levels[1 + high + 2 * dcc$active], levels = grp_levels)
  des <- make_design(dcc, z$reps)
  rhs <- c("grp", setdiff(confounders, character(0)))
  f <- reformulate(rhs, response = "obese")
  fit <- svyglm(f, des, family = quasibinomial())
  b <- coef(fit); se <- sqrt(diag(vcov(fit)))
  # standardised (covariate-adjusted) probability of obesity in each group, with replicate-weight SEs
  theta <- function(w, data) {
    data$.w <- w / mean(w)
    m <- glm(f, data = data, weights = .w, family = quasibinomial())
    sapply(grp_levels, function(l) { nd <- data; nd$grp <- factor(l, levels = grp_levels)
                                     weighted.mean(predict(m, nd, type = "response"), w) })
  }
  st <- withReplicates(des, theta)
  do.call(rbind, lapply(seq_along(grp_levels), function(i) {
    l <- grp_levels[i]; s <- subset(des, grp == l)
    um <- svymean(~obese, s); uest <- as.numeric(coef(um)); use <- as.numeric(SE(um))
    nm <- paste0("grp", l)
    data.frame(rule = rule_name, group = l, n = sum(dcc$grp == l), cases = sum(dcc$obese[dcc$grp == l]),
      unadj_pct = 100 * uest, unadj_low = 100 * max(0, uest - 1.96 * use), unadj_high = 100 * min(1, uest + 1.96 * use),
      adj_pct = 100 * coef(st)[i], adj_low = 100 * (coef(st)[i] - 1.96 * SE(st)[i]), adj_high = 100 * (coef(st)[i] + 1.96 * SE(st)[i]),
      or = if (i == 1) NA else exp(b[nm]), or_low = if (i == 1) NA else exp(b[nm] - 1.96 * se[nm]),
      or_high = if (i == 1) NA else exp(b[nm] + 1.96 * se[nm]),
      or_p = if (i == 1) NA else 2 * pnorm(-abs(b[nm] / se[nm])), row.names = NULL)
  }))
}
grp_res <- do.call(rbind, lapply(groups, function(g)
  do.call(rbind, lapply(names(fibre_rules), function(r) cbind(cohort = g, group_analysis(frames[[g]], r))))))
grp_res$primary_rule <- grp_res$rule == SPEC$fibre_rule_primary
write.csv(grp_res, file.path(out_dir, "four-group-rates.csv"), row.names = FALSE)

message("Done. Commit the files in ", out_dir, "/ (aggregate only), rq3_spec.R and rq3_run_all.R, then push.")
