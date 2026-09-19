# Team 5 RQ3 main analysis using the two cleaned eligible_sample sheets.
# Requires readxl and survey. Restricted respondent-level data stay local.
# Set RQ3_RAW_DATA_DIR to the folder containing AHSnpa11bp.csv and inp13bsp.csv.
suppressPackageStartupMessages({
  library(readxl)
  library(survey)
})
options(survey.lonely.psu = "adjust")

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
repo <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
raw_dir <- Sys.getenv("RQ3_RAW_DATA_DIR", unset = "")
if (!nzchar(raw_dir)) stop("Set RQ3_RAW_DATA_DIR to the authorised ABS Basic CURF folder.")
out_dir <- Sys.getenv("RQ3_OUTPUT_DIR", unset = repo)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

prepare <- function(group) {
  general <- group == "General NNPAS"
  book <- file.path(repo, if (general) "rq3_general_nnpas_cleaned.xlsx" else "rq3_indigenous_natsinpas_cleaned.xlsx")
  raw_file <- file.path(raw_dir, if (general) "AHSnpa11bp.csv" else "inp13bsp.csv")
  if (!file.exists(book) || !file.exists(raw_file)) stop("Missing input for ", group)
  clean <- as.data.frame(read_excel(book, sheet = "eligible_sample"))
  raw <- read.csv(raw_file, stringsAsFactors = FALSE)
  key <- if (general) as.character(raw$ABSPID) else paste(raw$ABSHID, raw$ABSGID, sep = "-")
  stopifnot(!anyDuplicated(key), !anyDuplicated(clean$participant_id))
  idx <- match(as.character(clean$participant_id), key)
  stopifnot(all(!is.na(idx)))
  reps <- if (general) sprintf("W2M01%02d", 1:60) else sprintf("WPM1%03d", 1:250)
  stopifnot(all(reps %in% names(raw)))
  raw <- raw[idx, , drop = FALSE]
  expected_weight <- if (general) as_num(raw$NPAD2WGT) else as_num(raw$IPAFINWT)
  stopifnot(all(abs(expected_weight - as_num(clean$selected_survey_weight)) < 1e-6))
  pregnancy <- as_num(raw$SABDYMS) == 4
  if (general) pregnancy <- pregnancy | as_num(raw$FEMLSBC) == 2
  keep <- is.na(pregnancy) | !pregnancy
  clean <- clean[keep, , drop = FALSE]
  raw <- raw[keep, , drop = FALSE]
  n_eligible <- nrow(clean)

  pa <- as.character(clean$pa_guideline)
  stopifnot(all(pa %in% c("Met 150-minute guideline", "Did not meet guideline")))
  fibre <- as_num(clean$fibre_2day_avg_g)
  energy <- as_num(clean$energy_2day_avg_kj)
  alcohol <- as_num(clean$alcohol_2day_avg_g)
  age <- as_num(clean$age)
  sex <- as_num(clean$sex_raw)
  ses <- as_num(clean$seifa_disadvantage)
  income <- as_num(clean$income_decile)
  school <- as_num(clean$highest_school_year_raw)
  smoking <- as_num(clean$smoking_status_raw)
  obese <- as_num(clean$obese_bmi)
  wt <- as_num(clean$selected_survey_weight)
  complete <- is.finite(fibre) & fibre >= 0 & is.finite(energy) & energy > 0 &
    is.finite(alcohol) & alcohol >= 0 & is.finite(age) & age >= 18 & age <= 65 &
    sex %in% 1:2 & ses %in% 1:5 & income %in% 1:10 & school %in% 1:5 &
    smoking %in% 1:5 & obese %in% 0:1 & is.finite(wt) & wt > 0
  n_complete <- sum(complete)
  clean <- clean[complete, , drop = FALSE]
  raw <- raw[complete, , drop = FALSE]
  pa <- pa[complete]; fibre <- fibre[complete]; energy <- energy[complete]
  alcohol <- alcohol[complete]; age <- age[complete]; sex <- sex[complete]
  ses <- ses[complete]; income <- income[complete]; school <- school[complete]
  smoking <- smoking[complete]; obese <- obese[complete]; wt <- wt[complete]
  d <- data.frame(
    obese = obese,
    fibre = fibre,
    fibre5 = (fibre - 20) / 5,
    active = as.integer(pa == "Met 150-minute guideline"),
    age10 = (age - 45) / 10,
    female = as.integer(sex == 2),
    energy_mj = (energy - 8000) / 1000,
    ses = factor(ses, levels = 1:5),
    income = income - 5.5,
    year12 = as.integer(school == 1),
    smoking = factor(ifelse(smoking %in% 1:3, "Current",
                      ifelse(smoking == 4, "Former", "Never")),
                     levels = c("Never", "Former", "Current")),
    alcohol = log1p(alcohol / 10),
    mainwt = wt
  )
  d <- cbind(d, raw[, reps, drop = FALSE])
  stopifnot(all(complete.cases(d)), all(as.matrix(d[, reps, drop = FALSE]) >= 0))
  flow <- data.frame(group = group,
                     stage = c("Cleaned eligible sample", "After pregnancy exclusion", "Complete main-model covariates"),
                     n = c(length(keep), n_eligible, n_complete))
  list(data = d, reps = reps, flow = flow)
}

formula_main <- obese ~ fibre5 * active + age10 + I(age10^2) + female +
  energy_mj + ses + income + year12 + smoking + alcohol

contrast <- function(fit, terms, label) {
  b <- coef(fit); v <- vcov(fit)
  c <- setNames(rep(0, length(b)), names(b))
  c[terms] <- 1
  est <- sum(c * b); se <- sqrt(drop(t(c) %*% v %*% c))
  data.frame(estimate = label, odds_ratio = exp(est),
             ci_low = exp(est - 1.96 * se), ci_high = exp(est + 1.96 * se),
             p_value = 2 * pnorm(-abs(est / se)))
}

analyse <- function(group) {
  z <- prepare(group); d <- z$data; G <- length(z$reps)
  design <- svrepdesign(data = d, weights = ~mainwt,
                       repweights = as.matrix(d[, z$reps, drop = FALSE]),
                       type = "JK1", scale = (G - 1) / G,
                       rscales = rep(1, G), combined.weights = TRUE, mse = TRUE)
  fit <- svyglm(formula_main, design = design, family = quasibinomial())
  stopifnot(all(is.finite(coef(fit))))
  rows <- rbind(
    contrast(fit, "fibre5", "Fibre +5 g/day, did not meet PA guideline"),
    contrast(fit, c("fibre5", "fibre5:active"), "Fibre +5 g/day, met PA guideline"),
    contrast(fit, "fibre5:active", "Fibre x PA interaction: ratio of odds ratios"),
    contrast(fit, "active", "Met vs did not meet PA guideline at 20 g fibre/day")
  )
  rows$group <- group; rows$n <- nrow(d); rows$obesity_cases <- sum(d$obese)
  rows <- rows[, c("group", "n", "obesity_cases", "estimate", "odds_ratio", "ci_low", "ci_high", "p_value")]
  pred <- as.numeric(predict(fit, type = "response"))
  q99 <- as.numeric(quantile(d$fibre, .99))
  fit99 <- svyglm(formula_main, design = subset(design, fibre <= q99), family = quasibinomial())
  cut <- contrast(fit99, "fibre5:active", "Interaction excluding top 1% fibre")
  fit_quad <- svyglm(update(formula_main, . ~ . + I(fibre5^2) + I(fibre5^2):active),
                     design = design, family = quasibinomial())
  qt <- grep("fibre5\\^2", names(coef(fit_quad)), value = TRUE)
  stopifnot(length(qt) == 2L)
  qb <- coef(fit_quad)[qt]; qv <- vcov(fit_quad)[qt, qt, drop = FALSE]
  qstat <- drop(t(qb) %*% solve(qv, qb))
  diagnostics <- data.frame(
    group = group, n = nrow(d), obesity_cases = sum(d$obese), replicate_weights = G,
    weighted_obesity_percent = 100 * weighted.mean(d$obese, d$mainwt),
    weighted_brier = weighted.mean((d$obese - pred)^2, d$mainwt),
    design_matrix_condition_number = kappa(model.matrix(formula_main, d), exact = TRUE),
    fibre_99th_percentile_g = q99, top_1_percent_excluded_n = sum(d$fibre > q99),
    top_1_percent_interaction_or = cut$odds_ratio,
    top_1_percent_interaction_ci_low = cut$ci_low,
    top_1_percent_interaction_ci_high = cut$ci_high,
    quadratic_fibre_joint_p = pchisq(qstat, df = 2, lower.tail = FALSE)
  )
  by_activity <- do.call(rbind, lapply(0:1, function(a) {
    sub <- d[d$active == a, , drop = FALSE]
    data.frame(group = group, activity = if (a) "Met 150-minute guideline" else "Did not meet 150-minute guideline",
               n = nrow(sub), obesity_cases = sum(sub$obese),
               weighted_obesity_percent = 100 * weighted.mean(sub$obese, sub$mainwt))
  }))
  list(flow = z$flow, model = rows, diagnostics = diagnostics, by_activity = by_activity)
}

results <- lapply(c("General NNPAS", "Indigenous NATSINPAS"), analyse)
for (name in c("flow", "model", "diagnostics", "by_activity")) {
  write.csv(do.call(rbind, lapply(results, "[[", name)),
            file.path(out_dir, paste0("rq3-rosetta-", gsub("_", "-", name), ".csv")), row.names = FALSE)
}
cat("Wrote four aggregate tables to", out_dir, "\n")
