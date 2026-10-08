# =====================================================================================
# rq3_spec.R  --  SHARED definitions for the fibre x physical activity analysis
#
# Everything the group might change (confounders, obesity definition, age bands,
# fibre cut-offs) lives in THIS file. Rosetta's models, Amelie's robustness checks and
# anyone else's code should source() it so everyone uses exactly the same variables.
#
# Needs: the two cleaned workbooks in the working directory, plus the authorised ABS Basic
# CURF files (AHSnpa11bp.csv, inp13bsp.csv, inp13bhh.csv) in the folder named by the
# environment variable RQ3_RAW_DATA_DIR (set it in ~/.Renviron; never commit those files).
# =====================================================================================

suppressPackageStartupMessages({library(readxl); library(survey)})
options(survey.lonely.psu = "adjust")

`%||%` <- function(a, b) if (is.null(a)) b else a
as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
valid  <- function(x, ok) { x <- as_num(x); x[!(ok(x) %in% TRUE)] <- NA; x }   # invalid codes -> NA

# ------------------------------- GROUP DECISIONS -------------------------------------
SPEC <- list(
  # obesity = BMI >= bmi_cut AND waist >= waist_cut (both must be met)
  bmi_cut    = 30,
  waist_cut  = c(male = 94, female = 80),   # group decision: 94 cm men, 80 cm women (health.gov.au increased-risk cut-offs)
                                            # ("increased risk") and 102/88 ("high risk") by sex:
                                            # change the female value here to run that version.
  age_breaks = c(18, 30, 40, 50, 60, 66),   # bands 18-29, 30-39, 40-49, 50-59, 60-65
  fibre_centre = 20, fibre_unit = 5,        # fibre5 = (g/day - 20) / 5  -> effects are per 5 g/day
  protein_centre = 80, protein_unit = 10,   # protein10 = (g/day - 80) / 10
  groups = c("General NNPAS", "Indigenous NATSINPAS (non-remote)"),
  # high / low fibre rule for the 4-group obesity-rate analysis (decided LAST, see run script)
  fibre_rule_primary = "NHMRC adequate intake (men 30 g, women 25 g)"
)

# Confounders, exactly as agreed. Coding of each variable is in prep_cohort() below.
confounders <- c("age_cat", "female", "energy_mj", "income", "year12", "smoking", "alcohol", "protein10")
make_formula <- function(y) reformulate(c("fibre5 * active", confounders), response = y)

# Rosetta's original model (used only to confirm the data pipeline still reproduces her results)
legacy_formula <- obese_bmi ~ fibre5 * active + age10 + I(age10^2) + female +
  energy_mj + ses + income + year12 + smoking + alcohol

# --------------------------------- DATA PREPARATION ----------------------------------
cohort_info <- list(
  "General NNPAS" = list(book = "rq3_general_nnpas_cleaned.xlsx", raw = "AHSnpa11bp.csv",
                         wt = "NPAD2WGT", reps = sprintf("W2M01%02d", 1:60)),
  "Indigenous NATSINPAS (non-remote)" = list(book = "rq3_indigenous_natsinpas_cleaned.xlsx", raw = "inp13bsp.csv",
                         wt = "IPAFINWT", reps = sprintf("WPM1%03d", 1:250))
)

load_inputs <- function(group) {
  raw_dir <- Sys.getenv("RQ3_RAW_DATA_DIR", unset = "")
  if (!nzchar(raw_dir)) stop("Set RQ3_RAW_DATA_DIR to the folder with the authorised ABS Basic CURF files.")
  ci <- cohort_info[[group]]
  if (!file.exists(ci$book)) stop("Can't find ", ci$book, " in ", getwd(), ". Set the working directory to the project folder.")
  list(clean = as.data.frame(read_excel(ci$book, sheet = "eligible_sample")),
       raw = read.csv(file.path(raw_dir, ci$raw), stringsAsFactors = FALSE))
}

# need_waist = TRUE  : analysis frame requires measured BMI AND waist (the agreed obesity definition)
# need_waist = FALSE : only used for the legacy check against Rosetta's BMI-only model
prep_cohort <- function(inp, group, need_waist = TRUE) {
  clean <- inp$clean; raw <- inp$raw; ci <- cohort_info[[group]]
  general <- group == "General NNPAS"
  key <- if (general) as.character(raw$ABSPID) else paste(raw$ABSHID, raw$ABSGID, sep = "-")
  stopifnot(!anyDuplicated(key), !anyDuplicated(clean$participant_id))
  idx <- match(as.character(clean$participant_id), key)
  stopifnot(!anyNA(idx))
  raw <- raw[idx, , drop = FALSE]
  stopifnot(all(ci$reps %in% names(raw)))
  stopifnot(all(abs(as_num(raw[[ci$wt]]) - as_num(clean$selected_survey_weight)) < 1e-6))   # weights line up

  preg <- as_num(raw$SABDYMS) == 4
  if (general) preg <- preg | as_num(raw$FEMLSBC) == 2
  keep <- is.na(preg) | !preg
  clean <- clean[keep, , drop = FALSE]; raw <- raw[keep, , drop = FALSE]

  fibre   <- valid(clean$fibre_2day_avg_g,   function(x) is.finite(x) & x >= 0)
  energy  <- valid(clean$energy_2day_avg_kj, function(x) is.finite(x) & x > 0)
  alcohol <- valid(clean$alcohol_2day_avg_g, function(x) is.finite(x) & x >= 0)
  # protein: use the workbook column if Selina added one, otherwise the two-day mean from the raw file
  protein <- if ("protein_2day_avg_g" %in% names(clean)) {
    valid(clean$protein_2day_avg_g, function(x) is.finite(x) & x >= 0)
  } else {
    (valid(raw$PROTT1, function(x) is.finite(x) & x >= 0) + valid(raw$PROTT2, function(x) is.finite(x) & x >= 0)) / 2
  }
  age     <- valid(clean$age,                function(x) is.finite(x) & x >= 18 & x <= 65)
  sex     <- valid(clean$sex_raw,            function(x) x %in% 1:2)
  ses     <- valid(clean$seifa_disadvantage, function(x) x %in% 1:5)
  income  <- valid(clean$income_decile,      function(x) x %in% 1:10)
  school  <- valid(clean$highest_school_year_raw, function(x) x %in% 1:5)
  smoking <- valid(clean$smoking_status_raw, function(x) x %in% 1:5)
  pa      <- valid(raw$EXREGUIC,             function(x) x %in% 1:2)
  wt      <- valid(clean$selected_survey_weight, function(x) is.finite(x) & x > 0)
  bmi     <- valid(clean$bmi,                function(x) is.finite(x) & x > 0 & x < 97)
  waist   <- valid(clean$waist_cm,           function(x) is.finite(x) & x > 0 & x < 997)

  wcut    <- ifelse(sex == 2, SPEC$waist_cut[["female"]], SPEC$waist_cut[["male"]])
  obese   <- as.integer(bmi >= SPEC$bmi_cut & waist >= wcut)        # both criteria; NA if either missing

  age_lab <- paste0(head(SPEC$age_breaks, -1), "–", SPEC$age_breaks[-1] - 1)
  d <- data.frame(
    bmi = bmi, waist = waist, obese = obese, obese_bmi = as.integer(bmi >= SPEC$bmi_cut),
    fibre = fibre, fibre5 = (fibre - SPEC$fibre_centre) / SPEC$fibre_unit,
    active = as.integer(pa == 1),
    age = age, age_cat = cut(age, SPEC$age_breaks, right = FALSE, labels = age_lab),
    age10 = (age - 45) / 10,
    female = as.integer(sex == 2),
    energy_mj = (energy - 8000) / 1000,
    income = income - 5.5,
    year12 = as.integer(school == 1),
    smoking = factor(ifelse(smoking %in% 1:3, "Current", ifelse(smoking == 4, "Former", "Never")),
                     levels = c("Never", "Former", "Current")),
    alcohol = log1p(alcohol / 10),
    protein10 = (protein - SPEC$protein_centre) / SPEC$protein_unit,
    ses = factor(ses, levels = 1:5),
    mainwt = wt)
  # Analysis frame: outcome, exposure, modifier, age, sex and weight must be known.
  # Missing confounders stay NA here so they can be imputed in the robustness analysis.
  frame <- !is.na(d$bmi) & (!need_waist | !is.na(d$waist)) & !is.na(d$fibre) & !is.na(d$active) &
           !is.na(d$age) & !is.na(d$female) & !is.na(d$mainwt)
  d <- cbind(d, raw[, ci$reps, drop = FALSE])[frame, , drop = FALSE]
  rownames(d) <- NULL
  stopifnot(!anyNA(d[, ci$reps]), all(as.matrix(d[, ci$reps]) >= 0))
  list(data = d, reps = ci$reps, n_eligible = nrow(clean), group = group)
}

make_design <- function(d, reps) {
  G <- length(reps)
  svrepdesign(data = d, weights = ~mainwt, repweights = as.matrix(d[, reps, drop = FALSE]),
              type = "JK1", scale = (G - 1) / G, rscales = rep(1, G), combined.weights = TRUE, mse = TRUE)
}

outcome_info <- list(
  bmi   = list(label = "BMI (kg/m²)", family = "gaussian",
               effect = "Difference in fibre slope (kg/m² per 5 g/day)"),
  waist = list(label = "Waist circumference (cm)", family = "gaussian",
               effect = "Difference in fibre slope (cm per 5 g/day)"),
  obese = list(label = "Obesity (BMI ≥ 30 and waist ≥ 94 cm men, ≥ 80 cm women)", family = "binomial",
               effect = "Ratio of fibre odds ratios")
)
complete_rows <- function(d, y) d[complete.cases(d[, all.vars(make_formula(y))]), , drop = FALSE]

fit_svy <- function(y, dcc, reps) {
  fam <- if (outcome_info[[y]]$family == "binomial") quasibinomial() else gaussian()
  # do.call embeds the objects in the call, which svyglm needs (it re-evaluates its own call)
  do.call(svyglm, list(formula = make_formula(y), design = make_design(dcc, reps), family = fam))
}

# ------------------------------- ESTIMATE EXTRACTION ---------------------------------
# estimate of a sum of coefficients (a contrast); exp() for the logistic model (odds ratio)
contrast_row <- function(fit, terms, label, logistic) {
  b <- coef(fit); v <- vcov(fit)
  est <- sum(b[terms]); se <- sqrt(sum(v[terms, terms, drop = FALSE]))
  f <- if (logistic) exp else identity
  data.frame(effect = label, estimate = f(est), lower = f(est - 1.96 * se), upper = f(est + 1.96 * se),
             p = 2 * pnorm(-abs(est / se)), row.names = NULL)
}
key_effects <- function(fit, logistic) {
  rbind(contrast_row(fit, "fibre5", "Fibre, per 5 g/day: below 150 min/week", logistic),
        contrast_row(fit, c("fibre5", "fibre5:active"), "Fibre, per 5 g/day: at least 150 min/week", logistic),
        contrast_row(fit, "fibre5:active", "Fibre x activity interaction", logistic),
        contrast_row(fit, "active", "Activity (at least vs below 150 min/week) at 20 g/day fibre", logistic))
}

pretty_term <- function(x) {
  map <- c("(Intercept)" = "Intercept", fibre5 = "Fibre (per 5 g/day)",
           active = "Activity ≥150 min/week (vs below)", "fibre5:active" = "Fibre x activity interaction",
           female = "Female (vs male)", energy_mj = "Energy (per MJ/day)",
           income = "Household income (per decile)", year12 = "Completed Year 12 (vs not)",
           smokingFormer = "Former smoker (vs never)", smokingCurrent = "Current smoker (vs never)",
           alcohol = "Alcohol (log scale)", protein10 = "Protein (per 10 g/day)")
  out <- ifelse(x %in% names(map), map[x], x)
  ifelse(grepl("^age_cat", x), paste0("Age ", sub("^age_cat", "", x), " years (vs 18–29)"), out)
}

# ---------------------------- ROBUSTNESS / MI HELPERS --------------------------------
pool_contrast <- function(fits, terms = "fibre5:active") {   # Rubin's rules
  m <- length(fits)
  est <- sapply(fits, function(f) sum(coef(f)[terms]))
  u   <- sapply(fits, function(f) sum(vcov(f)[terms, terms, drop = FALSE]))
  c(est = mean(est), se = sqrt(mean(u) + (1 + 1 / m) * var(est)))
}
est_row <- function(label, est, se, n, logistic) {
  f <- if (logistic) exp else identity
  data.frame(analysis = label, n = n, estimate = f(est), lower = f(est - 1.96 * se),
             upper = f(est + 1.96 * se), p = 2 * pnorm(-abs(est / se)))
}
int_from_fit <- function(fit, label, n, logistic, term = "fibre5:active") {
  b <- coef(fit); se <- sqrt(diag(vcov(fit)))
  est_row(label, unname(b[term]), unname(se[term]), n, logistic)
}
# glm / lm / glmrob / lmrob summaries differ slightly; this works for all of them
int_from_summary <- function(fit, label, n, logistic, term = "fibre5:active") {
  cf <- coef(summary(fit)); est_row(label, cf[term, 1], cf[term, 2], n, logistic)
}
