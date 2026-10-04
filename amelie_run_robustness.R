# =====================================================================================
# amelie_run_robustness.R
# Run this ONCE on your own laptop (it needs the raw ABS files). It writes small AGGREGATE
# tables to the folder amelie_results/ (no respondent rows). The shared Rmd only reads those
# tables, so nobody else needs the raw data to knit.
#
# Before running:
#   1. RQ3_RAW_DATA_DIR set in your .Renviron (folder with AHSnpa11bp.csv and inp13bsp.csv)
#   2. Working directory = the project folder that holds the two cleaned workbooks and
#      rq3-rosetta-results.csv (open the .Rproj, or setwd() to that folder)
#   3. Run in RStudio with Source, or:  Rscript amelie_run_robustness.R
# Re-run it whenever Rosetta's model, the formula or the waist cut-offs change.
# =====================================================================================

suppressPackageStartupMessages({library(readxl); library(survey); library(mice); library(robustbase)})
options(survey.lonely.psu = "adjust")

m_imputations <- 10          # use 5 for a quick test run
out_dir <- "amelie_results"
dir.create(out_dir, showWarnings = FALSE)

# ---- helper functions (copy of the tested analysis code) ----

`%||%` <- function(a, b) if (is.null(a)) b else a
as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
valid  <- function(x, ok) { x <- as_num(x); x[!(ok(x) %in% TRUE)] <- NA; x }   # invalid codes -> NA

full_formula <- obese ~ fibre5 * active + age10 + I(age10^2) + female +
  energy_mj + ses + income + year12 + smoking + alcohol
model_vars <- all.vars(full_formula)

load_raw <- function(group) {
  raw_dir <- Sys.getenv("RQ3_RAW_DATA_DIR", unset = "")
  if (!nzchar(raw_dir)) stop("Set RQ3_RAW_DATA_DIR to the folder with the authorised ABS Basic CURF files.")
  read.csv(file.path(raw_dir, if (group == "General NNPAS") "AHSnpa11bp.csv" else "inp13bsp.csv"),
           stringsAsFactors = FALSE)
}

waist_cut <- c(male = 94, female = 80)

prep_cohort <- function(clean, raw, group, outcome = c("bmi", "waist", "bmi_and_waist")) {
  outcome <- match.arg(outcome)
  general <- group == "General NNPAS"
  clean <- as.data.frame(clean)
  key <- if (general) as.character(raw$ABSPID) else paste(raw$ABSHID, raw$ABSGID, sep = "-")
  idx <- match(as.character(clean$participant_id), key)
  stopifnot(!anyNA(idx), !anyDuplicated(key))
  raw <- raw[idx, , drop = FALSE]
  reps <- if (general) sprintf("W2M01%02d", 1:60) else sprintf("WPM1%03d", 1:250)
  stopifnot(all(reps %in% names(raw)))
  wt_raw <- if (general) as_num(raw$NPAD2WGT) else as_num(raw$IPAFINWT)
  stopifnot(all(abs(wt_raw - as_num(clean$selected_survey_weight)) < 1e-6))

  preg <- as_num(raw$SABDYMS) == 4
  if (general) preg <- preg | as_num(raw$FEMLSBC) == 2
  keep <- is.na(preg) | !preg
  clean <- clean[keep, , drop = FALSE]; raw <- raw[keep, , drop = FALSE]

  fibre   <- valid(clean$fibre_2day_avg_g,   function(x) is.finite(x) & x >= 0)
  energy  <- valid(clean$energy_2day_avg_kj, function(x) is.finite(x) & x > 0)
  alcohol <- valid(clean$alcohol_2day_avg_g, function(x) is.finite(x) & x >= 0)
  age     <- valid(clean$age,                function(x) is.finite(x) & x >= 18 & x <= 65)
  sex     <- valid(clean$sex_raw,            function(x) x %in% 1:2)
  ses     <- valid(clean$seifa_disadvantage, function(x) x %in% 1:5)
  income  <- valid(clean$income_decile,      function(x) x %in% 1:10)
  school  <- valid(clean$highest_school_year_raw, function(x) x %in% 1:5)
  smoking <- valid(clean$smoking_status_raw, function(x) x %in% 1:5)
  pa      <- valid(raw$EXREGUIC,             function(x) x %in% 1:2)
  wt      <- valid(clean$selected_survey_weight, function(x) is.finite(x) & x > 0)

  bmi_flag <- valid(clean$obese_bmi, function(x) x %in% 0:1)
  waist_flag <- function() {
    if (!"waist_cm" %in% names(clean)) stop("No waist_cm column in this workbook")
    w <- valid(clean$waist_cm, function(x) is.finite(x) & x > 0)
    as.integer(w >= ifelse(sex == 2, waist_cut[["female"]], waist_cut[["male"]]))
  }
  obese <- switch(outcome,
    bmi           = bmi_flag,
    waist         = waist_flag(),
    bmi_and_waist = as.integer(bmi_flag == 1 & waist_flag() == 1))   # BOTH criteria (NA if undecidable)

  d <- data.frame(
    obese = obese, fibre5 = (fibre - 20) / 5, active = as.integer(pa == 1),
    age10 = (age - 45) / 10, female = as.integer(sex == 2),
    energy_mj = (energy - 8000) / 1000, ses = factor(ses, levels = 1:5),
    income = income - 5.5, year12 = as.integer(school == 1),
    smoking = factor(ifelse(smoking %in% 1:3, "Current", ifelse(smoking == 4, "Former", "Never")),
                     levels = c("Never", "Former", "Current")),
    alcohol = log1p(alcohol / 10), mainwt = wt)
  # Analysis frame: outcome, exposure, modifier, age, sex and weight must be known.
  # Missing confounders stay NA so they can be imputed.
  frame <- with(d, !is.na(obese) & !is.na(fibre5) & !is.na(active) & !is.na(age10) &
                   !is.na(female) & !is.na(mainwt))
  d <- cbind(d, raw[, reps, drop = FALSE])[frame, , drop = FALSE]
  rownames(d) <- NULL
  stopifnot(!anyNA(d[, reps]), all(as.matrix(d[, reps]) >= 0))
  list(data = d, reps = reps, n_eligible = nrow(clean))
}

make_design <- function(d, reps) {
  G <- length(reps)
  svrepdesign(data = d, weights = ~mainwt, repweights = as.matrix(d[, reps, drop = FALSE]),
              type = "JK1", scale = (G - 1) / G, rscales = rep(1, G),
              combined.weights = TRUE, mse = TRUE)
}

int_row <- function(label, fit, n, term = "fibre5:active") {
  b <- unname(coef(fit)[term]); se <- unname(sqrt(vcov(fit)[term, term]))
  data.frame(Analysis = label, n = n, OR = exp(b), Lower = exp(b - 1.96 * se),
             Upper = exp(b + 1.96 * se), p = 2 * pnorm(-abs(b / se)))
}

outlier_checks <- function(dcc, reps) {
  fit0 <- svyglm(full_formula, make_design(dcc, reps), family = quasibinomial())   # Rosetta's model
  g    <- glm(full_formula, data = dcc, family = binomial())                      # unweighted
  cd   <- cooks.distance(g)
  keep <- cd <= quantile(cd, 0.99)                                                # drop the 1% most influential
  fit_inf <- svyglm(full_formula, make_design(dcc[keep, ], reps), family = quasibinomial())
  gr <- tryCatch(glmrob(full_formula, data = dcc, family = binomial(), weights.on.x = "hat"),
                 error = function(e) NULL)
  rows <- list(int_row("Weighted, complete case (primary model)", fit0, nrow(dcc)),
               int_row("Weighted, most influential 1% removed (Cook's distance)", fit_inf, sum(keep)),
               int_row("Unweighted logistic", g, nrow(dcc)))
  rows <- c(rows, list(if (!is.null(gr)) int_row("Unweighted robust logistic (robustbase)", gr, nrow(dcc))
                       else data.frame(Analysis = "Unweighted robust logistic (robustbase): did not converge",
                                       n = nrow(dcc), OR = NA, Lower = NA, Upper = NA, p = NA)))
  list(fit = fit0, table = do.call(rbind, rows))
}

mi_fit <- function(d, reps, m = 10, seed = 123) {
  imp_dat <- d[, model_vars]
  imp_dat$fibre5_x_active <- imp_dat$fibre5 * imp_dat$active   # keeps the interaction in the imputation model
  imp_dat$year12 <- factor(imp_dat$year12)
  imp <- mice(imp_dat, m = m, maxit = 10, seed = seed, printFlag = FALSE)
  other <- d[, c("mainwt", reps), drop = FALSE]
  fits <- lapply(seq_len(m), function(i) {
    di <- complete(imp, i)
    di$fibre5_x_active <- NULL
    di$year12 <- as.integer(as.character(di$year12))
    svyglm(full_formula, make_design(cbind(di, other), reps), family = quasibinomial())
  })
  list(imp = imp, fits = fits)
}

pool_contrast <- function(fits, terms = "fibre5:active") {   # Rubin's rules for a sum of coefficients
  m <- length(fits)
  est <- sapply(fits, function(f) sum(coef(f)[terms]))
  u   <- sapply(fits, function(f) sum(vcov(f)[terms, terms, drop = FALSE]))
  c(est = mean(est), se = sqrt(mean(u) + (1 + 1 / m) * var(est)))
}
mi_row <- function(fits, n, label = NULL) {
  p <- pool_contrast(fits)
  data.frame(Analysis = label %||% paste0("Multiple imputation (m = ", length(fits), ")"), n = n,
             OR = exp(p[["est"]]), Lower = exp(p[["est"]] - 1.96 * p[["se"]]),
             Upper = exp(p[["est"]] + 1.96 * p[["se"]]), p = 2 * pnorm(-abs(p[["est"]] / p[["se"]])))
}

weighted_auc <- function(y, p, w) {
  o  <- order(p[y == 0]); pn <- p[y == 0][o]; wn <- w[y == 0][o]
  cw <- c(0, cumsum(wn))
  pp <- p[y == 1]; wp <- w[y == 1]
  lt <- cw[findInterval(pp, pn, left.open = TRUE) + 1]
  le <- cw[findInterval(pp, pn) + 1]
  sum(wp * (lt + 0.5 * (le - lt))) / (sum(wp) * sum(wn))
}
calib_table <- function(y, p, w, g = 10) {
  grp <- cut(p, quantile(p, seq(0, 1, length.out = g + 1)), include.lowest = TRUE, labels = FALSE)
  data.frame(decile = 1:g,
             predicted = as.numeric(tapply(p * w, grp, sum) / tapply(w, grp, sum)),
             observed  = as.numeric(tapply(y * w, grp, sum) / tapply(w, grp, sum)))
}
abs_scale <- function(dcc, reps, low = 10, high = 30) {
  f5 <- (c(low, high) - 20) / 5
  theta <- function(w, data) {
    data$.w <- w / mean(w)   # glm treats weights of binomial families as trial counts, so rescale
    fit <- glm(full_formula, data = data, weights = .w, family = quasibinomial())
    pr <- function(f, a) { nd <- data; nd$fibre5 <- f; nd$active <- a
                           weighted.mean(predict(fit, nd, type = "response"), w) }
    p <- c(inactive_low = pr(f5[1], 0), inactive_high = pr(f5[2], 0),
           active_low = pr(f5[1], 1),   active_high = pr(f5[2], 1))
    c(p, rd_inactive = unname(p[2] - p[1]), rd_active = unname(p[4] - p[3]),
      diff_in_diff = unname((p[4] - p[3]) - (p[2] - p[1])))
  }
  withReplicates(make_design(dcc, reps), theta)
}


# ---- inputs ----
groups <- c("General NNPAS", "Indigenous NATSINPAS (non-remote)")
books  <- c("rq3_general_nnpas_cleaned.xlsx", "rq3_indigenous_natsinpas_cleaned.xlsx")
for (f in c(books, "rq3-rosetta-results.csv"))
  if (!file.exists(f)) stop("Can't find ", f, " in ", getwd(), ". Set the working directory to the project folder.")
clean_books <- setNames(lapply(books, read_excel, sheet = "eligible_sample"), groups)
raw_files   <- setNames(lapply(groups, load_raw), groups)
build <- function(outcome) setNames(Map(function(cl, rw, g) prep_cohort(cl, rw, g, outcome),
                                        clean_books, raw_files, groups), groups)

# ---- 1. must reproduce Rosetta's primary model before anything else is trusted ----
s_chk <- build("bmi")
her <- subset(read.csv("rq3-rosetta-results.csv", stringsAsFactors = FALSE),
              analysis_id == "primary_two_day" & grepl("ratio of odds ratios", estimate))
chk <- do.call(rbind, lapply(groups, function(g) {
  z <- s_chk[[g]]; dcc <- z$data[complete.cases(z$data[, model_vars]), ]
  fit <- svyglm(full_formula, make_design(dcc, z$reps), family = quasibinomial())
  h <- her[her$group == g, ]
  data.frame(cohort = g, n_here = nrow(dcc), n_rosetta = h$n,
             or_here = int_row("", fit, nrow(dcc))$OR, or_rosetta = h$odds_ratio)
}))
print(chk)
if (!all(chk$n_here == chk$n_rosetta) || max(abs(chk$or_here - chk$or_rosetta)) > 0.001)
  stop("Does not reproduce Rosetta's primary model. Fix this before using any output.")
write.csv(chk, file.path(out_dir, "reproduce-check.csv"), row.names = FALSE)

# ---- 2. robustness checks for both outcome definitions ----
run_robustness <- function(outcome, m = m_imputations) {
  s <- build(outcome)
  out <- lapply(groups, function(g) {
    message("  ", outcome, " / ", g)
    z <- s[[g]]; dcc <- z$data[complete.cases(z$data[, model_vars]), ]
    oc <- outlier_checks(dcc, z$reps)
    mi <- mi_fit(z$data, z$reps, m = m)
    list(table = rbind(oc$table, mi_row(mi$fits, nrow(z$data))), fit = oc$fit, dcc = dcc,
         reps = z$reps, imp = mi$imp, data = z$data,
         flow = data.frame(eligible = z$n_eligible, known = nrow(z$data),
                           complete = nrow(dcc), cases = sum(dcc$obese)))
  })
  setNames(out, groups)
}
message("Running robustness checks (a few minutes)...")
rob_bmi  <- run_robustness("bmi")             # Set 1: BMI >= 30 (Rosetta's outcome)
rob_both <- run_robustness("bmi_and_waist")   # Set 2: BMI >= 30 AND waist above cut-off (provisional)

tab_of  <- function(rob, set) do.call(rbind, lapply(groups, function(g) cbind(outcome_set = set, cohort = g, rob[[g]]$table)))
flow_of <- function(rob, set) do.call(rbind, lapply(groups, function(g) cbind(outcome_set = set, cohort = g, rob[[g]]$flow)))
write.csv(rbind(tab_of(rob_bmi, "bmi"), tab_of(rob_both, "bmi_and_waist")),
          file.path(out_dir, "robustness.csv"), row.names = FALSE)
write.csv(rbind(flow_of(rob_bmi, "bmi"), flow_of(rob_both, "bmi_and_waist")),
          file.path(out_dir, "sample-flow.csv"), row.names = FALSE)

# ---- 3. performance, calibration, absolute-scale results (BMI-only primary model) ----
perf <- do.call(rbind, lapply(groups, function(g) {
  r <- rob_bmi[[g]]; d <- r$dcc; p <- as.numeric(predict(r$fit, type = "response")); pbar <- weighted.mean(d$obese, d$mainwt)
  data.frame(cohort = g, n = nrow(d), auc = weighted_auc(d$obese, p, d$mainwt),
             brier = weighted.mean((d$obese - p)^2, d$mainwt), brier_null = pbar * (1 - pbar))
}))
write.csv(perf, file.path(out_dir, "performance.csv"), row.names = FALSE)

calib <- do.call(rbind, lapply(groups, function(g) {
  r <- rob_bmi[[g]]; p <- as.numeric(predict(r$fit, type = "response"))
  cbind(cohort = g, calib_table(r$dcc$obese, p, r$dcc$mainwt))
}))
write.csv(calib, file.path(out_dir, "calibration.csv"), row.names = FALSE)

message("Absolute-scale results...")
abs_res <- do.call(rbind, lapply(groups, function(g) {
  r <- rob_bmi[[g]]; ab <- abs_scale(r$dcc, r$reps)
  data.frame(cohort = g, measure = names(coef(ab)), estimate = as.numeric(coef(ab)), se = as.numeric(SE(ab)))
}))
abs_res$lower <- abs_res$estimate - 1.96 * abs_res$se
abs_res$upper <- abs_res$estimate + 1.96 * abs_res$se
write.csv(abs_res, file.path(out_dir, "absolute-scale.csv"), row.names = FALSE)

# ---- 4. imputation diagnostics (General cohort, BMI set): look at these before trusting the imputation ----
png(file.path(out_dir, "imputation-chains.png"), width = 1400, height = 1000, res = 150); print(plot(rob_bmi[[1]]$imp)); dev.off()
png(file.path(out_dir, "imputation-density.png"), width = 1400, height = 1000, res = 150); print(densityplot(rob_bmi[[1]]$imp)); dev.off()
print(colSums(is.na(rob_bmi[[1]]$data[, model_vars])))   # which covariates were missing

message("Done. Commit the files in ", out_dir, "/ (aggregate only) and push.")
