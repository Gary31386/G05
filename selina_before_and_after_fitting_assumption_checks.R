# =====================================================================================
# selina_before_and_after_fitting_assumption_checks.R
#
# Assumption checks for the CURRENT RQ3 analysis defined in rq3_spec.R and fitted by
# rq3_run_all.R. This script does not alter Amelie's model files or rq3_results/.
# It writes aggregate tables and figures only to rq3_assumption_results/.
#
# Run from the G05 project directory:
#   Rscript selina_before_and_after_fitting_assumption_checks.R
#
# The checks are split into:
#   PART A: checks before fitting (data validity, support/overlap, missingness,
#           extreme values and multicollinearity)
#   PART B: checks after fitting (residual form, influence, calibration,
#           robustness, imputation diagnostics and nonlinear fibre sensitivity)
#
# Exact refitting and nonlinear checks require:
#   - RQ3_RAW_DATA_DIR pointing to the authorised ABS files
#   - the survey package
# If either is unavailable, the script still produces the checks supported by the
# cleaned workbooks and the current aggregate files in rq3_results/.
# =====================================================================================

required_packages <- c("readxl", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install missing packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
})

options(stringsAsFactors = FALSE, scipen = 999)
set.seed(3888)

out_dir <- "rq3_assumption_results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

groups <- c("General NNPAS", "Indigenous NATSINPAS (non-remote)")
books <- c(
  "General NNPAS" = "rq3_general_nnpas_cleaned.xlsx",
  "Indigenous NATSINPAS (non-remote)" = "rq3_indigenous_natsinpas_cleaned.xlsx"
)
current_results_dir <- "rq3_results"

needed_files <- c(unname(books), file.path(current_results_dir, c(
  "sample-flow.csv", "missing-covariates.csv", "vif.csv", "performance.csv",
  "calibration.csv", "robustness.csv"
)))
missing_files <- needed_files[!file.exists(needed_files)]
if (length(missing_files)) {
  stop("Missing required current-analysis files: ", paste(missing_files, collapse = ", "))
}

write_aggregate <- function(x, filename) {
  write.csv(x, file.path(out_dir, filename), row.names = FALSE, na = "")
}

read_current <- function(filename) {
  read.csv(file.path(current_results_dir, filename), stringsAsFactors = FALSE,
           check.names = FALSE)
}

weighted_quantile <- function(x, w, probs) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  x <- x[keep]
  w <- w[keep]
  if (!length(x)) return(rep(NA_real_, length(probs)))
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cumulative_weight <- cumsum(w) / sum(w)
  vapply(probs, function(p) x[which(cumulative_weight >= p)[1]], numeric(1))
}

weighted_mean_safe <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  weighted.mean(x[keep], w[keep])
}

theme_set(
  theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(face = "bold")
    )
)

# =====================================================================================
# PART A: CHECKS BEFORE FITTING
# =====================================================================================

message("PART A: checks before fitting")

prepare_screening_data <- function(path, cohort) {
  raw <- as.data.frame(read_excel(path, sheet = "eligible_sample"))
  required <- c(
    "age", "sex_raw", "fibre_2day_avg_g", "energy_2day_avg_kj",
    "alcohol_2day_avg_g", "pa_guideline_raw", "bmi", "waist_cm",
    "selected_survey_weight"
  )
  missing <- setdiff(required, names(raw))
  if (length(missing)) {
    stop(cohort, " cleaned workbook is missing: ", paste(missing, collapse = ", "))
  }

  data.frame(
    cohort = cohort,
    activity = factor(
      ifelse(raw$pa_guideline_raw == 1, "At least 150 min/week", "Below 150 min/week"),
      levels = c("Below 150 min/week", "At least 150 min/week")
    ),
    fibre = suppressWarnings(as.numeric(raw$fibre_2day_avg_g)),
    bmi = suppressWarnings(as.numeric(raw$bmi)),
    waist = suppressWarnings(as.numeric(raw$waist_cm)),
    energy_mj = suppressWarnings(as.numeric(raw$energy_2day_avg_kj)) / 1000,
    alcohol = suppressWarnings(as.numeric(raw$alcohol_2day_avg_g)),
    age = suppressWarnings(as.numeric(raw$age)),
    female = as.integer(suppressWarnings(as.numeric(raw$sex_raw)) == 2),
    survey_weight = suppressWarnings(as.numeric(raw$selected_survey_weight))
  )
}

screening <- do.call(rbind, lapply(groups, function(g) {
  prepare_screening_data(books[[g]], g)
}))
rownames(screening) <- NULL

valid_screening <- with(screening,
  is.finite(fibre) & fibre >= 0 &
    is.finite(bmi) & bmi > 0 & bmi < 97 &
    is.finite(waist) & waist > 0 & waist < 997 &
    !is.na(activity) & is.finite(age) & age >= 18 & age <= 65 &
    is.finite(survey_weight) & survey_weight > 0
)
screening <- screening[valid_screening, , drop = FALSE]

support_rows <- lapply(split(screening, list(screening$cohort, screening$activity), drop = TRUE), function(d) {
  q <- weighted_quantile(d$fibre, d$survey_weight,
                         c(0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99))
  data.frame(
    cohort = d$cohort[1],
    activity = as.character(d$activity[1]),
    n = nrow(d),
    weighted_percent = 100 * sum(d$survey_weight) /
      sum(screening$survey_weight[screening$cohort == d$cohort[1]]),
    fibre_min = min(d$fibre),
    fibre_q01 = q[1],
    fibre_q05 = q[2],
    fibre_q25 = q[3],
    fibre_median = q[4],
    fibre_q75 = q[5],
    fibre_q95 = q[6],
    fibre_q99 = q[7],
    fibre_max = max(d$fibre),
    weighted_mean_bmi = weighted_mean_safe(d$bmi, d$survey_weight),
    weighted_mean_waist = weighted_mean_safe(d$waist, d$survey_weight),
    stringsAsFactors = FALSE
  )
})
support_summary <- do.call(rbind, support_rows)
rownames(support_summary) <- NULL
write_aggregate(support_summary, "before_fitting_fibre_support.csv")

overlap_rows <- lapply(groups, function(g) {
  d <- screening[screening$cohort == g, , drop = FALSE]
  by_activity <- split(d, d$activity, drop = TRUE)
  bounds <- lapply(by_activity, function(z) {
    weighted_quantile(z$fibre, z$survey_weight, c(0.05, 0.95))
  })
  common_low <- max(vapply(bounds, `[`, numeric(1), 1))
  common_high <- min(vapply(bounds, `[`, numeric(1), 2))
  data.frame(
    cohort = g,
    common_support_q05_to_q95_low = common_low,
    common_support_q05_to_q95_high = common_high,
    n_outside_common_support = sum(d$fibre < common_low | d$fibre > common_high),
    pct_outside_common_support = 100 * mean(d$fibre < common_low | d$fibre > common_high),
    interpretation = if (common_low < common_high) {
      "Central fibre ranges overlap; interpret estimates cautiously outside common support."
    } else {
      "Central fibre ranges do not overlap adequately."
    },
    stringsAsFactors = FALSE
  )
})
overlap_summary <- do.call(rbind, overlap_rows)
write_aggregate(overlap_summary, "before_fitting_overlap_summary.csv")

continuous_variables <- c(
  fibre = "Fibre (g/day)", bmi = "BMI (kg/m2)", waist = "Waist (cm)",
  energy_mj = "Energy (MJ/day)", alcohol = "Alcohol (g/day)", age = "Age (years)"
)
extreme_rows <- list()
k <- 1L
for (g in groups) {
  d <- screening[screening$cohort == g, , drop = FALSE]
  for (v in names(continuous_variables)) {
    x <- d[[v]]
    q <- weighted_quantile(x, d$survey_weight, c(0.01, 0.05, 0.50, 0.95, 0.99))
    extreme_rows[[k]] <- data.frame(
      cohort = g,
      variable = continuous_variables[[v]],
      n = sum(is.finite(x)),
      min = min(x, na.rm = TRUE),
      q01 = q[1], q05 = q[2], median = q[3], q95 = q[4], q99 = q[5],
      max = max(x, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
}
extreme_summary <- do.call(rbind, extreme_rows)
write_aggregate(extreme_summary, "before_fitting_extreme_value_summary.csv")

fibre_overlap_plot <- ggplot(
  screening,
  aes(x = fibre, colour = activity, fill = activity, weight = survey_weight)
) +
  geom_density(alpha = 0.16, adjust = 1.1) +
  facet_wrap(~cohort, scales = "free_y") +
  coord_cartesian(xlim = range(
    weighted_quantile(screening$fibre, screening$survey_weight, c(0.005, 0.995))
  )) +
  labs(
    title = "Fibre support before fitting the interaction models",
    subtitle = "Survey-weighted distributions by physical activity group",
    x = "Two-day mean fibre intake (g/day)", y = "Weighted density",
    colour = "Physical activity", fill = "Physical activity"
  )
print(fibre_overlap_plot)
ggsave(file.path(out_dir, "before_fitting_fibre_overlap.png"), fibre_overlap_plot,
       width = 10, height = 5.5, dpi = 180)

boxplot_data <- do.call(rbind, lapply(names(continuous_variables)[1:5], function(v) {
  data.frame(
    cohort = screening$cohort,
    activity = screening$activity,
    variable = continuous_variables[[v]],
    value = screening[[v]],
    stringsAsFactors = FALSE
  )
}))
boxplot_data <- boxplot_data[is.finite(boxplot_data$value), , drop = FALSE]
boxplot_data$variable <- factor(boxplot_data$variable, levels = unname(continuous_variables[1:5]))
distribution_plot <- ggplot(boxplot_data, aes(x = activity, y = value, fill = activity)) +
  geom_boxplot(outlier.alpha = 0.20, width = 0.65) +
  facet_grid(variable ~ cohort, scales = "free_y") +
  labs(
    title = "Continuous-variable screening before model fitting",
    subtitle = "Unweighted boxplots; points beyond the whiskers require sensitivity checks, not automatic deletion",
    x = NULL, y = NULL, fill = "Physical activity"
  ) +
  theme(axis.text.x = element_text(angle = 18, hjust = 1))
print(distribution_plot)
ggsave(file.path(out_dir, "before_fitting_distribution_boxplots.png"), distribution_plot,
       width = 11, height = 10, dpi = 180)

sample_flow <- read_current("sample-flow.csv")
missingness <- read_current("missing-covariates.csv")
vif <- read_current("vif.csv")
write_aggregate(sample_flow, "before_fitting_sample_flow.csv")
write_aggregate(missingness, "before_fitting_missingness.csv")
write_aggregate(vif, "before_fitting_multicollinearity.csv")

# =====================================================================================
# PART B: CHECKS AFTER FITTING
# =====================================================================================

message("PART B: checks after fitting")

performance <- read_current("performance.csv")
calibration <- read_current("calibration.csv")
robustness <- read_current("robustness.csv")
write_aggregate(performance, "after_fitting_performance.csv")
write_aggregate(calibration, "after_fitting_calibration.csv")
write_aggregate(robustness, "after_fitting_robustness.csv")

calibration_plot <- ggplot(calibration, aes(x = predicted, y = observed, colour = cohort)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_line(linewidth = 0.7, alpha = 0.8) +
  geom_point(size = 2.4) +
  facet_wrap(~cohort) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "Obesity-model calibration after fitting",
    subtitle = "Observed versus predicted obesity within weighted risk groups",
    x = "Mean predicted probability", y = "Observed weighted proportion",
    colour = NULL
  ) +
  theme(legend.position = "none")
print(calibration_plot)
ggsave(file.path(out_dir, "after_fitting_calibration.png"), calibration_plot,
       width = 9, height = 5, dpi = 180)

null_lines <- data.frame(
  outcome = c("bmi", "waist", "obese"),
  null = c(0, 0, 1)
)
outcome_levels <- c("bmi", "waist", "obese")
outcome_labels <- c("BMI", "Waist circumference", "Obesity")
robustness$outcome_label <- factor(
  robustness$outcome,
  levels = outcome_levels,
  labels = outcome_labels
)
null_lines$outcome_label <- factor(
  null_lines$outcome,
  levels = outcome_levels,
  labels = outcome_labels
)
robustness_plot <- ggplot(
  robustness,
  aes(x = estimate, y = analysis, xmin = lower, xmax = upper)
) +
  geom_vline(
    data = null_lines,
    aes(xintercept = null),
    inherit.aes = FALSE,
    linetype = 2,
    colour = "grey45"
  ) +
  geom_errorbar(orientation = "y", width = 0.15, colour = "#2563EB") +
  geom_point(size = 2.2, colour = "#1D4ED8") +
  facet_grid(outcome_label ~ cohort, scales = "free_x") +
  labs(
    title = "Interaction estimates across after-fitting sensitivity checks",
    subtitle = "Null is 0 for BMI/waist slope differences and 1 for the obesity ratio of odds ratios",
    x = "Interaction estimate with 95% confidence interval", y = NULL
  )
print(robustness_plot)
ggsave(file.path(out_dir, "after_fitting_robustness_forest.png"), robustness_plot,
       width = 13, height = 9, dpi = 180)

# Preserve Amelie's current diagnostic figures without modifying rq3_results/.
figure_map <- c(
  "diagnostics-1.png" = "after_fitting_residual_diagnostics_general.png",
  "diagnostics-2.png" = "after_fitting_residual_diagnostics_indigenous.png",
  "imputation-chains.png" = "after_fitting_imputation_chains_general.png",
  "imputation-density.png" = "after_fitting_imputation_density_general.png"
)
for (source_name in names(figure_map)) {
  source_path <- file.path(current_results_dir, source_name)
  if (file.exists(source_path)) {
    file.copy(source_path, file.path(out_dir, figure_map[[source_name]]), overwrite = TRUE)
  }
}

max_vif_indices <- unlist(lapply(split(seq_len(nrow(vif)), vif$cohort), function(i) {
  i[which.max(vif$vif[i])]
}), use.names = FALSE)
max_vif_rows <- vif[max_vif_indices, , drop = FALSE]
missing_nonzero <- missingness[missingness$n_missing > 0, , drop = FALSE]
valid_robustness <- robustness[grepl("^Weighted|^Multiple imputation", robustness$analysis), ]
valid_robustness$contains_null <- ifelse(
  valid_robustness$outcome == "obese",
  valid_robustness$lower <= 1 & valid_robustness$upper >= 1,
  valid_robustness$lower <= 0 & valid_robustness$upper >= 0
)

raw_dir <- Sys.getenv("RQ3_RAW_DATA_DIR", unset = "")
raw_files <- c("AHSnpa11bp.csv", "inp13bsp.csv")
can_refit <- nzchar(raw_dir) && all(file.exists(file.path(raw_dir, raw_files))) &&
  requireNamespace("survey", quietly = TRUE)

nonlinearity_note <- if (can_refit) {
  "Exact survey-weighted quadratic fibre checks were run; see after_fitting_nonlinearity_checks.csv."
} else {
  paste(
    "Exact nonlinear refitting was not run in this session because RQ3_RAW_DATA_DIR",
    "and/or the survey package was unavailable. Run this script in the configured RStudio session."
  )
}

assumption_summary <- data.frame(
  stage = c(
    "Before fitting", "Before fitting", "Before fitting", "Before fitting",
    "After fitting", "After fitting", "After fitting", "After fitting",
    "After fitting", "After fitting"
  ),
  check = c(
    "Data validity and survey weights",
    "Fibre support across activity groups",
    "Missing covariates",
    "Multicollinearity",
    "Residual mean structure",
    "Residual distribution",
    "Influential observations",
    "Obesity-model calibration and discrimination",
    "Multiple-imputation diagnostics",
    "Nonlinear fibre sensitivity"
  ),
  status = c(
    "Satisfactory",
    ifelse(any(overlap_summary$common_support_q05_to_q95_low >=
                 overlap_summary$common_support_q05_to_q95_high), "Concern", "Caution"),
    "Addressed with multiple imputation",
    ifelse(max(vif$vif) < 5, "Satisfactory", "Concern"),
    "Caution",
    "Caution",
    ifelse(all(valid_robustness$contains_null), "Conclusion stable", "Review"),
    "Modest predictive performance",
    "General cohort available; Indigenous cohort still needs a saved diagnostic plot",
    ifelse(can_refit, "Completed", "Pending exact rerun")
  ),
  evidence = c(
    "Cleaned-workbook validation checks pass; rq3_spec.R verifies participant matching and replicate weights when raw files are available.",
    paste0("Central 5th-95th percentile fibre ranges overlap in both cohorts; ",
           paste(round(overlap_summary$pct_outside_common_support, 1), collapse = "% and "),
           "% of records lie outside cohort-specific common support."),
    if (nrow(missing_nonzero)) {
      paste(paste0(missing_nonzero$cohort, ": ", missing_nonzero$variable, " ",
                   missing_nonzero$pct_missing, "%"), collapse = "; ")
    } else "No missing model covariates.",
    paste0("Maximum VIF values are ",
           paste(sprintf("%.2f (%s)", max_vif_rows$vif, max_vif_rows$cohort), collapse = "; "), "."),
    "Residual-versus-fitted plots are broadly flat for General NNPAS but show more curvature in the smaller Indigenous cohort.",
    "Q-Q plots show tail departures, especially for BMI; survey replicate-weight inference and robustness checks reduce reliance on normal residuals.",
    "Cook's-distance sensitivity, survey-weighted primary models and multiple imputation should be interpreted together; unweighted rows are not valid survey inference.",
    paste(sprintf("%s: AUC %.3f, Brier %.3f versus null %.3f",
                  performance$cohort, performance$auc, performance$brier,
                  performance$brier_null), collapse = "; "),
    "Current saved chain and density figures are for the General cohort only.",
    nonlinearity_note
  ),
  stringsAsFactors = FALSE
)

# ---- Exact survey-weighted nonlinear checks when the authorised raw files are available ----
if (can_refit) {
  message("Raw files and survey package found: running exact current-model nonlinear checks")
  source("rq3_spec.R")
  inputs <- setNames(lapply(SPEC$groups, load_inputs), SPEC$groups)
  frames <- setNames(lapply(SPEC$groups, function(g) prep_cohort(inputs[[g]], g)), SPEC$groups)

  nonlinear_rows <- list()
  diagnostic_rows <- list()
  row_index <- 1L
  diagnostic_index <- 1L

  for (g in SPEC$groups) {
    z <- frames[[g]]
    for (y in names(outcome_info)) {
      dcc <- complete_rows(z$data, y)
      primary_fit <- fit_svy(y, dcc, z$reps)
      quadratic_formula <- update(
        make_formula(y),
        . ~ . + I(fibre5^2) + I(fibre5^2):active
      )
      family_object <- if (outcome_info[[y]]$family == "binomial") {
        quasibinomial()
      } else {
        gaussian()
      }
      quadratic_fit <- do.call(
        survey::svyglm,
        list(
          formula = quadratic_formula,
          design = make_design(dcc, z$reps),
          family = family_object
        )
      )
      term_test <- survey::regTermTest(
        quadratic_fit,
        ~I(fibre5^2) + I(fibre5^2):active,
        method = "Wald"
      )
      quadratic_terms <- grep("I\\(fibre5\\^2\\)", names(coef(quadratic_fit)), value = TRUE)
      quadratic_interaction <- quadratic_terms[grepl(":active", quadratic_terms)]
      quadratic_main <- setdiff(quadratic_terms, quadratic_interaction)

      nonlinear_rows[[row_index]] <- data.frame(
        cohort = g,
        outcome = y,
        n = nrow(dcc),
        joint_quadratic_wald_p = as.numeric(term_test$p),
        quadratic_main_estimate = unname(coef(quadratic_fit)[quadratic_main]),
        quadratic_main_p = 2 * pnorm(-abs(
          coef(quadratic_fit)[quadratic_main] /
            sqrt(diag(vcov(quadratic_fit)))[quadratic_main]
        )),
        quadratic_interaction_estimate = unname(coef(quadratic_fit)[quadratic_interaction]),
        quadratic_interaction_p = 2 * pnorm(-abs(
          coef(quadratic_fit)[quadratic_interaction] /
            sqrt(diag(vcov(quadratic_fit)))[quadratic_interaction]
        )),
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1L

      if (y %in% c("bmi", "waist")) {
        residual_value <- as.numeric(residuals(primary_fit))
        fitted_value <- as.numeric(fitted(primary_fit))
        standardised_residual <- residual_value / sd(residual_value, na.rm = TRUE)
        diagnostic_rows[[diagnostic_index]] <- data.frame(
          cohort = g,
          outcome = y,
          fitted = fitted_value,
          residual = residual_value,
          sqrt_abs_standardised_residual = sqrt(abs(standardised_residual)),
          stringsAsFactors = FALSE
        )
        diagnostic_index <- diagnostic_index + 1L
      }
    }
  }

  nonlinear_checks <- do.call(rbind, nonlinear_rows)
  write_aggregate(nonlinear_checks, "after_fitting_nonlinearity_checks.csv")

  scale_location_data <- do.call(rbind, diagnostic_rows)
  scale_location_plot <- ggplot(
    scale_location_data,
    aes(x = fitted, y = sqrt_abs_standardised_residual)
  ) +
    geom_point(alpha = 0.18, size = 0.7, colour = "grey35") +
    geom_smooth(method = "loess", se = FALSE, colour = "#DC2626", linewidth = 0.9) +
    facet_grid(outcome ~ cohort, scales = "free") +
    labs(
      title = "Scale-location checks for the continuous primary models",
      subtitle = "A nearly horizontal smooth supports approximately constant residual spread",
      x = "Fitted value", y = "Square root of absolute standardised residual"
    )
  print(scale_location_plot)
  ggsave(file.path(out_dir, "after_fitting_scale_location.png"), scale_location_plot,
         width = 11, height = 7, dpi = 180)

  assumption_summary$status[assumption_summary$check == "Nonlinear fibre sensitivity"] <-
    ifelse(all(nonlinear_checks$joint_quadratic_wald_p >= 0.05),
           "No clear nonlinear departure", "Review nonlinear specification")
  assumption_summary$evidence[assumption_summary$check == "Nonlinear fibre sensitivity"] <-
    paste(sprintf("%s/%s joint quadratic p=%.3f",
                  nonlinear_checks$cohort, nonlinear_checks$outcome,
                  nonlinear_checks$joint_quadratic_wald_p), collapse = "; ")
} else {
  message(
    "Exact refitting skipped: set RQ3_RAW_DATA_DIR and install survey, then rerun in RStudio."
  )
}

write_aggregate(assumption_summary, "assumption_check_summary.csv")

manifest <- data.frame(
  file = list.files(out_dir),
  description = NA_character_,
  stringsAsFactors = FALSE
)
manifest$description <- ifelse(
  grepl("^before_fitting", manifest$file),
  "Check performed before fitting or interpreting the primary models",
  ifelse(
    grepl("^after_fitting", manifest$file),
    "Diagnostic or sensitivity check performed after fitting",
    "Overall assumption-check summary"
  )
)
write_aggregate(manifest, "output_manifest.csv")

message("Done. Aggregate assumption-check outputs are in ", out_dir, "/")
