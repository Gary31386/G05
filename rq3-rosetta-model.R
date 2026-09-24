# STAT3888 Team 5 — RQ3 fibre × physical activity model
#
# QUESTION
#   Among adults aged 18–65, does physical activity modify the association
#   between reported fibre intake and measured obesity (BMI >= 30 kg/m^2)?
#   NNPAS and non-remote NATSINPAS are fitted separately; this is not a direct
#   Indigenous versus non-Indigenous comparison.
#
# PRIMARY MODEL
#   Exposure: mean Day 1 + Day 2 fibre, per 5 g/day.
#   Modifier: >=150 minutes of activity in the previous week (EXREGUIC).
#   Covariates: age, age^2, sex, two-day mean energy, SEIFA disadvantage,
#   income decile, Year 12, smoking, and two-day mean alcohol.
#   Weights: NNPAS Day 2 weight + 60 W2M replicates; NATSINPAS person weight
#   + 250 WPM replicates because the Basic CURF has no Day 2-specific weight.
#   Therefore the NATSINPAS primary result is a weighted Day 2-completer
#   analysis, not an estimate representative of all non-remote adults.
#
# SENSITIVITY MODELS IN THIS SAME FILE
#   1. Two-day exposure with activity also requiring five sessions (EXREGUID).
#   2. Day 1 exposure in the exact primary-model cohort and weights.
#   3. Day 1 exposure among all eligible Day 1 respondents.
#   4. Day 1 exposure with the core confounder set (income and Year 12 omitted).
#
# INPUTS
#   Repository: rq3_general_nnpas_cleaned.xlsx,
#               rq3_indigenous_natsinpas_cleaned.xlsx
#   RQ3_RAW_DATA_DIR: AHSnpa11bp.csv, inp13bsp.csv, inp13bhh.csv
#
# OUTPUT
#   rq3-rosetta-results.csv — one aggregate table for every model above.
#   No respondent-level rows are written.
#
# RUN
#   RQ3_RAW_DATA_DIR=/authorised/abs/folder Rscript rq3-rosetta-model.R

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
read_abs <- function(file) read.csv(file.path(raw_dir, file), stringsAsFactors = FALSE)

groups <- c("General NNPAS", "Indigenous NATSINPAS (non-remote)")

full_formula <- obese ~ fibre5 * active + age10 + I(age10^2) + female +
  energy_mj + ses + income + year12 + smoking + alcohol
core_formula <- obese ~ fibre5 * active + age10 + I(age10^2) + female +
  energy_mj + ses + smoking + alcohol

prepare_day2_cohort <- function(group, dietary_measure, pa_field) {
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
  raw <- raw[idx, , drop = FALSE]

  reps <- if (general) sprintf("W2M01%02d", 1:60) else sprintf("WPM1%03d", 1:250)
  stopifnot(all(reps %in% names(raw)), pa_field %in% c("EXREGUIC", "EXREGUID"))
  expected_weight <- if (general) as_num(raw$NPAD2WGT) else as_num(raw$IPAFINWT)
  stopifnot(all(abs(expected_weight - as_num(clean$selected_survey_weight)) < 1e-6))

  pregnancy <- as_num(raw$SABDYMS) == 4
  if (general) pregnancy <- pregnancy | as_num(raw$FEMLSBC) == 2
  keep <- is.na(pregnancy) | !pregnancy
  clean <- clean[keep, , drop = FALSE]
  raw <- raw[keep, , drop = FALSE]

  pa_label <- as.character(clean$pa_guideline)
  stopifnot(all(pa_label %in% c("Met 150-minute guideline", "Did not meet guideline")),
            all(as_num(raw$EXREGUIC) %in% 1:2),
            identical(as_num(raw$EXREGUIC) == 1, pa_label == "Met 150-minute guideline"),
            all(as_num(raw$EXREGUID) == as_num(clean$pa_150min_5sessions_raw)))

  if (dietary_measure == "two_day_mean") {
    fibre <- as_num(clean$fibre_2day_avg_g)
    energy <- as_num(clean$energy_2day_avg_kj)
    alcohol <- as_num(clean$alcohol_2day_avg_g)
  } else if (dietary_measure == "day1") {
    fibre <- as_num(clean$fibre_day1_g)
    energy <- as_num(clean$energy_day1_kj)
    alcohol <- as_num(clean$alcohol_day1_g)
  } else {
    stop("Unknown dietary_measure: ", dietary_measure)
  }

  # Both versions use the exact primary-model cohort. This makes the Day 1
  # comparison a change in recall definition, not a change in respondents.
  both_days <- is.finite(as_num(clean$fibre_2day_avg_g)) &
    is.finite(as_num(clean$energy_2day_avg_kj)) & as_num(clean$energy_2day_avg_kj) > 0 &
    is.finite(as_num(clean$alcohol_2day_avg_g)) & as_num(clean$alcohol_2day_avg_g) >= 0
  pa <- as_num(raw[[pa_field]])
  age <- as_num(clean$age)
  sex <- as_num(clean$sex_raw)
  ses <- as_num(clean$seifa_disadvantage)
  income <- as_num(clean$income_decile)
  school <- as_num(clean$highest_school_year_raw)
  smoking <- as_num(clean$smoking_status_raw)
  obese <- as_num(clean$obese_bmi)
  wt <- as_num(clean$selected_survey_weight)
  complete <- both_days & pa %in% 1:2 & is.finite(fibre) & fibre >= 0 &
    is.finite(energy) & energy > 0 & is.finite(alcohol) & alcohol >= 0 &
    is.finite(age) & age >= 18 & age <= 65 & sex %in% 1:2 & ses %in% 1:5 &
    income %in% 1:10 & school %in% 1:5 & smoking %in% 1:5 & obese %in% 0:1 &
    is.finite(wt) & wt > 0

  raw <- raw[complete, , drop = FALSE]
  d <- data.frame(
    obese = obese[complete], fibre = fibre[complete],
    fibre5 = (fibre[complete] - 20) / 5,
    active = as.integer(pa[complete] == 1),
    age10 = (age[complete] - 45) / 10,
    female = as.integer(sex[complete] == 2),
    energy_mj = (energy[complete] - 8000) / 1000,
    ses = factor(ses[complete], levels = 1:5),
    income = income[complete] - 5.5,
    year12 = as.integer(school[complete] == 1),
    smoking = factor(ifelse(smoking[complete] %in% 1:3, "Current",
                      ifelse(smoking[complete] == 4, "Former", "Never")),
                     levels = c("Never", "Former", "Current")),
    alcohol = log1p(alcohol[complete] / 10), mainwt = wt[complete]
  )
  d <- cbind(d, raw[, reps, drop = FALSE])
  stopifnot(all(complete.cases(d)), all(as.matrix(d[, reps, drop = FALSE]) >= 0))
  list(data = d, reps = reps)
}

prepare_all_day1 <- function(group, pa_field, confounder_set) {
  indigenous <- group == "Indigenous NATSINPAS (non-remote)"
  p <- read_abs(if (indigenous) "inp13bsp.csv" else "AHSnpa11bp.csv")
  if (indigenous) {
    h <- read_abs("inp13bhh.csv")
    stopifnot(!anyDuplicated(h$ABSHID))
    p <- merge(p, h[, c("ABSHID", "RAECURF", "SF2SA1DB", "INCDECE")],
               by = "ABSHID", all.x = TRUE, sort = FALSE)
    p <- p[!is.na(p$RAECURF) & as_num(p$RAECURF) == 1, , drop = FALSE]
  }

  age_field <- if (indigenous) "AGEEC" else "AGEC"
  age <- as_num(p[[age_field]])
  p <- p[is.finite(age) & age >= 18 & age <= 65, , drop = FALSE]
  bmi <- as_num(p$BMISC)
  p <- p[is.finite(bmi) & bmi > 0 & bmi < 97, , drop = FALSE]
  pregnancy <- as_num(p$SABDYMS) == 4
  if (!indigenous) pregnancy <- pregnancy | as_num(p$FEMLSBC) == 2
  p <- p[is.na(pregnancy) | !pregnancy, , drop = FALSE]

  reps <- if (indigenous) sprintf("WPM1%03d", 1:250) else sprintf("WPM01%02d", 1:60)
  weight_field <- if (indigenous) "IPAFINWT" else "NPAFINWT"
  stopifnot(all(reps %in% names(p)), pa_field %in% names(p))
  pa <- as_num(p[[pa_field]])
  fibre <- as_num(p$FIBRET1)
  energy <- as_num(p$ENERGYT1)
  alcohol <- as_num(p$ALCT1)
  wt <- as_num(p[[weight_field]])
  sex <- as_num(p$SEX)
  ses <- as_num(p[[if (indigenous) "SF2SA1DB" else "SF2SA1QN"]])
  income <- as_num(p[[if (indigenous) "INCDECE" else "INCDEC"]])
  school <- as_num(p$HYSCHCBC)
  smoking <- as_num(p$SMKSTAT)
  age <- as_num(p[[age_field]])

  complete <- pa %in% 1:2 & is.finite(fibre) & fibre >= 0 &
    is.finite(energy) & energy > 0 & is.finite(alcohol) & alcohol >= 0 &
    is.finite(wt) & wt > 0 & sex %in% 1:2 & ses %in% 1:5 & smoking %in% 1:5
  if (confounder_set == "full") complete <- complete & income %in% 1:10 & school %in% 1:5

  p <- p[complete, , drop = FALSE]
  d <- data.frame(
    obese = as.integer(as_num(p$BMISC) >= 30), fibre = as_num(p$FIBRET1),
    fibre5 = (as_num(p$FIBRET1) - 20) / 5,
    active = as.integer(as_num(p[[pa_field]]) == 1),
    age10 = (as_num(p[[age_field]]) - 45) / 10,
    female = as.integer(as_num(p$SEX) == 2),
    energy_mj = (as_num(p$ENERGYT1) - 8000) / 1000,
    ses = factor(as_num(p[[if (indigenous) "SF2SA1DB" else "SF2SA1QN"]]), levels = 1:5),
    smoking = factor(ifelse(as_num(p$SMKSTAT) %in% 1:3, "Current",
                      ifelse(as_num(p$SMKSTAT) == 4, "Former", "Never")),
                     levels = c("Never", "Former", "Current")),
    alcohol = log1p(as_num(p$ALCT1) / 10), mainwt = as_num(p[[weight_field]])
  )
  if (confounder_set == "full") {
    d$income <- as_num(p[[if (indigenous) "INCDECE" else "INCDEC"]]) - 5.5
    d$year12 <- as.integer(as_num(p$HYSCHCBC) == 1)
  }
  d <- cbind(d, p[, reps, drop = FALSE])
  stopifnot(all(complete.cases(d)), all(d$mainwt > 0),
            all(as.matrix(d[, reps, drop = FALSE]) >= 0))
  list(data = d, reps = reps)
}

contrast <- function(fit, terms, label) {
  b <- coef(fit)
  v <- vcov(fit)
  cvec <- setNames(rep(0, length(b)), names(b))
  cvec[terms] <- 1
  est <- sum(cvec * b)
  se <- sqrt(drop(t(cvec) %*% v %*% cvec))
  data.frame(estimate = label, odds_ratio = exp(est),
             ci_low = exp(est - 1.96 * se), ci_high = exp(est + 1.96 * se),
             p_value = 2 * pnorm(-abs(est / se)))
}

analyse <- function(spec, group) {
  if (spec$cohort == "day2_completers") {
    z <- prepare_day2_cohort(group, spec$dietary_measure, spec$pa_definition)
  } else {
    z <- prepare_all_day1(group, spec$pa_definition, spec$confounder_set)
  }
  d <- z$data
  G <- length(z$reps)
  model_formula <- if (spec$confounder_set == "core") core_formula else full_formula
  design <- svrepdesign(data = d, weights = ~mainwt,
                       repweights = as.matrix(d[, z$reps, drop = FALSE]),
                       type = "JK1", scale = (G - 1) / G, rscales = rep(1, G),
                       combined.weights = TRUE, mse = TRUE)
  fit <- svyglm(model_formula, design = design, family = quasibinomial())
  stopifnot(all(is.finite(coef(fit))))

  threshold <- if (spec$pa_definition == "EXREGUID") {
    "150 min/week and five sessions"
  } else {
    "150 min/week"
  }
  rows <- rbind(
    contrast(fit, "fibre5", paste0("Fibre +5 g/day, below ", threshold)),
    contrast(fit, c("fibre5", "fibre5:active"), paste0("Fibre +5 g/day, met ", threshold)),
    contrast(fit, "fibre5:active", "Fibre x PA interaction: ratio of odds ratios"),
    contrast(fit, "active", paste0("Met vs below ", threshold, " at 20 g fibre/day"))
  )

  pred <- as.numeric(predict(fit, type = "response"))
  q99 <- as.numeric(quantile(d$fibre, 0.99))
  fit99 <- svyglm(model_formula, design = subset(design, fibre <= q99), family = quasibinomial())
  cut <- contrast(fit99, "fibre5:active", "Interaction excluding top 1% fibre")
  fit_quad <- svyglm(update(model_formula, . ~ . + I(fibre5^2) + I(fibre5^2):active),
                     design = design, family = quasibinomial())
  quad_terms <- grep("fibre5\\^2", names(coef(fit_quad)), value = TRUE)
  stopifnot(length(quad_terms) == 2L)
  quad_b <- coef(fit_quad)[quad_terms]
  quad_v <- vcov(fit_quad)[quad_terms, quad_terms, drop = FALSE]
  quad_stat <- drop(t(quad_b) %*% solve(quad_v, quad_b))

  rows$analysis_id <- spec$analysis_id
  rows$analysis_role <- spec$analysis_role
  rows$dietary_measure <- spec$dietary_measure
  rows$cohort <- spec$cohort
  rows$pa_definition <- spec$pa_definition
  rows$confounder_set <- spec$confounder_set
  rows$group <- group
  rows$n <- nrow(d)
  rows$obesity_cases <- sum(d$obese)
  rows$replicate_weights <- G
  rows$weighted_obesity_percent <- 100 * weighted.mean(d$obese, d$mainwt)
  rows$weighted_brier <- weighted.mean((d$obese - pred)^2, d$mainwt)
  rows$design_matrix_condition_number <- kappa(model.matrix(model_formula, d), exact = TRUE)
  rows$fibre_99th_percentile_g <- q99
  rows$top_1_percent_excluded_n <- sum(d$fibre > q99)
  rows$top_1_percent_interaction_or <- cut$odds_ratio
  rows$top_1_percent_interaction_ci_low <- cut$ci_low
  rows$top_1_percent_interaction_ci_high <- cut$ci_high
  rows$quadratic_fibre_joint_p <- pchisq(quad_stat, df = 2, lower.tail = FALSE)
  rows[, c("analysis_id", "analysis_role", "dietary_measure", "cohort",
           "pa_definition", "confounder_set", "group", "n", "obesity_cases",
           "estimate", "odds_ratio", "ci_low", "ci_high", "p_value",
           "replicate_weights", "weighted_obesity_percent", "weighted_brier",
           "design_matrix_condition_number", "fibre_99th_percentile_g",
           "top_1_percent_excluded_n", "top_1_percent_interaction_or",
           "top_1_percent_interaction_ci_low", "top_1_percent_interaction_ci_high",
           "quadratic_fibre_joint_p")]
}

specs <- data.frame(
  analysis_id = c("primary_two_day", "two_day_five_sessions", "day1_same_cohort",
                  "day1_all", "day1_core"),
  analysis_role = c("Primary", rep("Sensitivity", 4)),
  dietary_measure = c("two_day_mean", "two_day_mean", "day1", "day1", "day1"),
  cohort = c("day2_completers", "day2_completers", "day2_completers",
             "all_day1_eligible", "all_day1_eligible"),
  pa_definition = c("EXREGUIC", "EXREGUID", "EXREGUIC", "EXREGUIC", "EXREGUIC"),
  confounder_set = c("full", "full", "full", "full", "core"),
  stringsAsFactors = FALSE
)

results <- do.call(rbind, lapply(seq_len(nrow(specs)), function(i) {
  spec <- as.list(specs[i, , drop = FALSE])
  do.call(rbind, lapply(groups, function(group) analyse(spec, group)))
}))
write.csv(results, file.path(out_dir, "rq3-rosetta-results.csv"), row.names = FALSE)
cat("Wrote one aggregate RQ3 results table to", out_dir, "\n")
