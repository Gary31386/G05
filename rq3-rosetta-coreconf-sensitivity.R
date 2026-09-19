# Rosetta's core-confounder sensitivity uses Day 1 intake and person weights.
# Exports aggregate tables only; two-day intake is a separately reported sensitivity.
# Requires the authorised 2011-12 NNPAS and 2012-13 NATSINPAS Basic CURFs.
suppressPackageStartupMessages(library(survey))
options(survey.lonely.psu = "adjust")

raw_dir <- Sys.getenv("RQ3_RAW_DATA_DIR", unset = "")
if (!nzchar(raw_dir)) stop("Set RQ3_RAW_DATA_DIR to the authorised ABS Basic CURF folder.")
args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
repo <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
out_dir <- Sys.getenv("RQ3_OUTPUT_DIR", unset = repo)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
read_abs <- function(file) read.csv(file.path(raw_dir, file), stringsAsFactors = FALSE)

prepare <- function(group, pa_field) {
  indigenous <- group == "Indigenous NATSINPAS (non-remote)"
  p <- read_abs(if (indigenous) "inp13bsp.csv" else "AHSnpa11bp.csv")
  if (indigenous) {
    h <- read_abs("inp13bhh.csv")
    stopifnot(!anyDuplicated(h$ABSHID))
    p <- merge(p, h[, c("ABSHID", "RAECURF", "SF2SA1DB", "INCDECE")],
               by = "ABSHID", all.x = TRUE, sort = FALSE)
    p <- p[!is.na(p$RAECURF) & num(p$RAECURF) == 1, , drop = FALSE]
  }
  age_field <- if (indigenous) "AGEEC" else "AGEC"
  age <- num(p[[age_field]])
  p <- p[is.finite(age) & age >= 18 & age <= 65, , drop = FALSE]
  n_age <- nrow(p)
  bmi <- num(p$BMISC)
  p <- p[is.finite(bmi) & bmi > 0 & bmi < 97, , drop = FALSE]
  n_bmi <- nrow(p)
  pregnant <- num(p$SABDYMS) == 4
  if (!indigenous) pregnant <- pregnant | num(p$FEMLSBC) == 2
  p <- p[is.na(pregnant) | !pregnant, , drop = FALSE]
  n_pregnancy <- nrow(p)

  reps <- if (indigenous) sprintf("WPM1%03d", 1:250) else sprintf("WPM01%02d", 1:60)
  weight_field <- if (indigenous) "IPAFINWT" else "NPAFINWT"
  stopifnot(all(reps %in% names(p)), pa_field %in% names(p))
  pa <- num(p[[pa_field]])
  fibre <- num(p$FIBRET1)
  energy <- num(p$ENERGYT1)
  alcohol <- num(p$ALCT1)
  wt <- num(p[[weight_field]])
  sex <- num(p$SEX)
  ses <- num(p[[if (indigenous) "SF2SA1DB" else "SF2SA1QN"]])
  income <- num(p[[if (indigenous) "INCDECE" else "INCDEC"]])
  school <- num(p$HYSCHCBC)
  smoking <- num(p$SMKSTAT)
  age <- num(p[[age_field]])
  complete <- pa %in% 1:2 & is.finite(fibre) & fibre >= 0 &
    is.finite(energy) & energy > 0 & is.finite(alcohol) & alcohol >= 0 &
    is.finite(wt) & wt > 0 & sex %in% 1:2 & ses %in% 1:5 &
    smoking %in% 1:5
  p <- p[complete, , drop = FALSE]
  d <- data.frame(
    obese = as.integer(num(p$BMISC) >= 30),
    fibre = num(p$FIBRET1),
    fibre5 = (num(p$FIBRET1) - 20) / 5,
    active = as.integer(num(p[[pa_field]]) == 1),
    age10 = (num(p[[age_field]]) - 45) / 10,
    female = as.integer(num(p$SEX) == 2),
    energy_mj = (num(p$ENERGYT1) - 8000) / 1000,
    ses = factor(num(p[[if (indigenous) "SF2SA1DB" else "SF2SA1QN"]]), levels = 1:5),
    smoking = factor(ifelse(num(p$SMKSTAT) %in% 1:3, "Current",
                            ifelse(num(p$SMKSTAT) == 4, "Former", "Never")),
                     levels = c("Never", "Former", "Current")),
    alcohol = log1p(num(p$ALCT1) / 10),
    mainwt = num(p[[weight_field]])
  )
  d <- cbind(d, p[, reps, drop = FALSE])
  stopifnot(all(complete.cases(d)), all(d$mainwt > 0),
            all(as.matrix(d[, reps, drop = FALSE]) >= 0))
  flow <- data.frame(group = group, pa_definition = pa_field,
                     stage = c("Aged 18-65", "Valid measured BMI", "After pregnancy exclusion",
                               "Complete model covariates"),
                     n = c(n_age, n_bmi, n_pregnancy, nrow(d)))
  list(data = d, reps = reps, flow = flow)
}

formula_main <- obese ~ fibre5 * active + age10 + I(age10^2) + female +
  energy_mj + ses + smoking + alcohol

contrast <- function(fit, terms, label) {
  b <- coef(fit)
  v <- vcov(fit)
  c <- setNames(rep(0, length(b)), names(b))
  c[terms] <- 1
  est <- sum(c * b)
  se <- sqrt(drop(t(c) %*% v %*% c))
  data.frame(estimate = label, odds_ratio = exp(est),
             ci_low = exp(est - 1.96 * se), ci_high = exp(est + 1.96 * se),
             p_value = 2 * pnorm(-abs(est / se)))
}

analyse <- function(group, pa_field) {
  z <- prepare(group, pa_field)
  d <- z$data
  G <- length(z$reps)
  design <- svrepdesign(data = d, weights = ~mainwt,
                       repweights = as.matrix(d[, z$reps, drop = FALSE]),
                       type = "JK1", scale = (G - 1) / G,
                       rscales = rep(1, G), combined.weights = TRUE, mse = TRUE)
  fit <- svyglm(formula_main, design = design, family = quasibinomial())
  stopifnot(all(is.finite(coef(fit))))
  rows <- rbind(
    contrast(fit, "fibre5", "Fibre +5 g/day, below activity threshold"),
    contrast(fit, c("fibre5", "fibre5:active"), "Fibre +5 g/day, met activity threshold"),
    contrast(fit, "fibre5:active", "Fibre x activity: ratio of odds ratios"),
    contrast(fit, "active", "Met vs below threshold at 20 g fibre/day")
  )
  rows$group <- group
  rows$pa_definition <- pa_field
  rows$n <- nrow(d)
  rows$obesity_cases <- sum(d$obese)
  rows <- rows[, c("group", "pa_definition", "n", "obesity_cases", "estimate",
                   "odds_ratio", "ci_low", "ci_high", "p_value")]
  pred <- as.numeric(predict(fit, type = "response"))
  q99 <- as.numeric(quantile(d$fibre, 0.99))
  fit99 <- svyglm(formula_main, design = subset(design, fibre <= q99),
                  family = quasibinomial())
  cut <- contrast(fit99, "fibre5:active", "Interaction after excluding top 1% fibre")
  fit_quad <- svyglm(update(formula_main, . ~ . + I(fibre5^2) + I(fibre5^2):active),
                     design = design, family = quasibinomial())
  qt <- grep("fibre5\\^2", names(coef(fit_quad)), value = TRUE)
  stopifnot(length(qt) == 2L)
  qb <- coef(fit_quad)[qt]
  qv <- vcov(fit_quad)[qt, qt, drop = FALSE]
  qstat <- drop(t(qb) %*% solve(qv, qb))
  diagnostics <- data.frame(
    group = group, pa_definition = pa_field, n = nrow(d),
    obesity_cases = sum(d$obese), replicate_weights = G,
    weighted_obesity_percent = 100 * weighted.mean(d$obese, d$mainwt),
    weighted_brier = weighted.mean((d$obese - pred)^2, d$mainwt),
    design_matrix_condition_number = kappa(model.matrix(formula_main, d), exact = TRUE),
    fibre_99th_percentile_g = q99,
    top_1_percent_excluded_n = sum(d$fibre > q99),
    top_1_percent_interaction_or = cut$odds_ratio,
    top_1_percent_interaction_ci_low = cut$ci_low,
    top_1_percent_interaction_ci_high = cut$ci_high,
    quadratic_fibre_joint_p = pchisq(qstat, df = 2, lower.tail = FALSE)
  )
  by_activity <- do.call(rbind, lapply(0:1, function(a) {
    sub <- d[d$active == a, , drop = FALSE]
    data.frame(group = group, pa_definition = pa_field,
               activity = if (a) "Met activity threshold" else "Below activity threshold",
               n = nrow(sub), obesity_cases = sum(sub$obese),
               weighted_obesity_percent = 100 * weighted.mean(sub$obese, sub$mainwt))
  }))
  list(flow = z$flow, model = rows, diagnostics = diagnostics,
       by_activity = by_activity)
}

groups <- c("General NNPAS", "Indigenous NATSINPAS (non-remote)")
pa_fields <- c("EXREGUIC", "EXREGUID") # 150 min; 150 min and five sessions
results <- lapply(groups, function(g) lapply(pa_fields, function(pa) analyse(g, pa)))
results <- unlist(results, recursive = FALSE)
for (name in c("flow", "model", "diagnostics", "by_activity")) {
  write.csv(do.call(rbind, lapply(results, `[[`, name)),
            file.path(out_dir, paste0("rq3-rosetta-coreconf-", gsub("_", "-", name), ".csv")),
            row.names = FALSE)
}
cat("Wrote four aggregate Day 1 tables to", out_dir, "\n")
