# RQ3 model assumption checks for the Group Rmd
# ============================================================
# Research question:
# Among Australian Indigenous and general-population adults, how does
# physical activity modify the association between fibre intake and obesity?
#
# This file is deliberately organised so Sections 5-8 can be copied into the
# Group Rmd one section at a time. It creates no CSV, PNG or PDF files.
# Every plot is printed immediately after its corresponding model.
#
# Important scope note: General NNPAS is not verified non-Indigenous.

# ============================================================
# 1. Packages
# ============================================================

required_packages <- c(
  "readxl", "dplyr", "ggplot2", "caret", "rpart",
  "rpart.plot", "randomForest", "mgcv", "viridis"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages)) {
  stop("Install missing packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(caret)
  library(rpart)
  library(rpart.plot)
  library(randomForest)
  library(mgcv)
  library(viridis)
})

set.seed(3888)

# ============================================================
# 2. Use the data objects already created by the Group Rmd
# ============================================================

# The existing load-data chunk already runs:
# library(readxl)
# general <- read_excel("rq3_general_nnpas_cleaned.xlsx", sheet = "eligible_sample")
# indigenous <- read_excel("rq3_indigenous_natsinpas_cleaned.xlsx", sheet = "eligible_sample")
# Do not read the workbooks again here.

if (!exists("general") || !is.data.frame(general)) {
  stop("Run the existing load-data chunk first: object 'general' is unavailable.")
}

if (!exists("indigenous") || !is.data.frame(indigenous)) {
  stop("Run the existing load-data chunk first: object 'indigenous' is unavailable.")
}

# ============================================================
# 3. Create common analysis variables
# ============================================================

prepare_model_data <- function(raw_data, population_name) {
  required_variables <- c(
    "age", "sex", "fibre_2day_avg_g", "energy_2day_avg_kj",
    "alcohol_2day_avg_g", "pa_guideline_raw", "obese_bmi",
    "smoking_status_raw", "highest_school_year_raw",
    "seifa_disadvantage", "income_decile", "selected_survey_weight",
    "include_primary_analysis"
  )

  missing_variables <- setdiff(required_variables, names(raw_data))
  if (length(missing_variables)) {
    stop(
      population_name, " is missing: ",
      paste(missing_variables, collapse = ", ")
    )
  }

  raw_data <- raw_data[
    raw_data$include_primary_analysis %in% TRUE,
    ,
    drop = FALSE
  ]

  model_data <- raw_data %>%
    transmute(
      Population = population_name,

      # Binary outcome: 1 = BMI-defined obesity
      obese = as.integer(obese_bmi),
      obesity = factor(
        obese_bmi,
        levels = c(0, 1),
        labels = c("Not obese", "Obese")
      ),

      # Two-day mean fibre, centred at 20 g/day and scaled per 5 g/day
      fibre5 = (fibre_2day_avg_g - 20) / 5,

      # PA raw code 1 = met 150-minute guideline
      active = as.integer(pa_guideline_raw == 1),
      activity = factor(
        active,
        levels = c(0, 1),
        labels = c("Below guideline", "Met guideline")
      ),

      # Person 1-style transformations
      age10 = (age - 45) / 10,
      age10_sq = age10^2,
      female = as.integer(sex == "Women"),
      sex_factor = factor(
        female,
        levels = c(0, 1),
        labels = c("Men", "Women")
      ),
      energy_mj = (energy_2day_avg_kj - 8000) / 1000,
      ses = factor(seifa_disadvantage, levels = 1:5),
      income = as.numeric(income_decile) - 5.5,
      year12 = as.integer(highest_school_year_raw == 1),
      year12_factor = factor(
        year12,
        levels = c(0, 1),
        labels = c("No", "Yes")
      ),
      smoking = factor(
        case_when(
          smoking_status_raw %in% 1:3 ~ "Current",
          smoking_status_raw == 4 ~ "Former",
          smoking_status_raw == 5 ~ "Never",
          TRUE ~ NA_character_
        ),
        levels = c("Never", "Former", "Current")
      ),
      alcohol = log1p(alcohol_2day_avg_g / 10),
      survey_weight = as.numeric(selected_survey_weight)
    )

  model_data <- model_data[
    complete.cases(model_data) &
      is.finite(model_data$survey_weight) &
      model_data$survey_weight > 0,
    ,
    drop = FALSE
  ]

  # rpart/model.frame is most reliable with a base data.frame rather than a
  # tibble returned by read_excel/dplyr.
  model_data <- as.data.frame(model_data)

  if (length(unique(model_data$obese)) != 2L) {
    stop(population_name, " does not contain both obesity outcomes.")
  }

  model_data
}

general_data <- prepare_model_data(general, "General NNPAS")
indigenous_data <- prepare_model_data(
  indigenous,
  "Indigenous NATSINPAS"
)

cat(
  "General NNPAS complete cases:", nrow(general_data), "\n",
  "Indigenous NATSINPAS complete cases:", nrow(indigenous_data), "\n"
)

normalise_weights <- function(weights) {
  weights / mean(weights)
}

# ============================================================
# 4. Person 1-style logistic regression baseline
# ============================================================

baseline_formula <- obese ~
  fibre5 * active +
  age10 +
  I(age10^2) +
  female +
  energy_mj +
  ses +
  income +
  year12 +
  smoking +
  alcohol

general_glm <- glm(
  baseline_formula,
  data = general_data,
  weights = normalise_weights(general_data$survey_weight),
  family = quasibinomial()
)

indigenous_glm <- glm(
  baseline_formula,
  data = indigenous_data,
  weights = normalise_weights(indigenous_data$survey_weight),
  family = quasibinomial()
)

extract_interaction <- function(model, population_name) {
  coefficient_table <- coef(summary(model))
  interaction_term <- "fibre5:active"

  data.frame(
    Population = population_name,
    Estimate = coefficient_table[interaction_term, "Estimate"],
    Standard_error = coefficient_table[interaction_term, "Std. Error"],
    P_value = coefficient_table[interaction_term, "Pr(>|t|)"],
    Ratio_of_odds_ratios = exp(
      coefficient_table[interaction_term, "Estimate"]
    )
  )
}

glm_interaction_results <- bind_rows(
  extract_interaction(general_glm, "General NNPAS"),
  extract_interaction(indigenous_glm, "Indigenous NATSINPAS")
)

print(glm_interaction_results)

# The cleaned workbooks contain main survey weights but not the replicate
# weights used in Person 1's primary svyglm analysis. Therefore, this GLM is
# the common predictive baseline for Sections 5-8, not a replacement for the
# primary replicate-weighted inference.

# ============================================================
# 5. Decision tree and plots
# ============================================================

tree_formula <- obesity ~
  fibre5 +
  activity +
  age10 +
  age10_sq +
  sex_factor +
  energy_mj +
  ses +
  income +
  year12_factor +
  smoking +
  alcohol

fit_pruned_tree <- function(df) {
  if (!is.data.frame(df)) {
    stop("fit_pruned_tree() requires a data.frame. Run Sections 1-3 first.")
  }

  required_tree_variables <- all.vars(tree_formula)
  missing_tree_variables <- setdiff(required_tree_variables, names(df))
  if (length(missing_tree_variables)) {
    stop(
      "Tree data are missing: ",
      paste(missing_tree_variables, collapse = ", "),
      ". Run Sections 1-3 before Section 5."
    )
  }

  # Keep the tree unweighted. rpart evaluates weights through model.frame,
  # which can fail after this section is copied into an R Markdown chunk.
  # Survey weights remain in the GLM and held-out CV metrics.
  large_tree <- rpart(
    tree_formula,
    data = df,
    method = "class",
    parms = list(split = "information"),
    control = rpart.control(
      cp = 0.0001,
      minsplit = 20,
      xval = 10
    )
  )

  cp_table <- large_tree$cptable
  best_cp <- cp_table[which.min(cp_table[, "xerror"]), "CP"]

  list(
    large_tree = large_tree,
    pruned_tree = prune(large_tree, cp = best_cp),
    best_cp = best_cp
  )
}

# 5.1 General NNPAS decision tree
general_tree <- fit_pruned_tree(general_data)

plotcp(
  general_tree$large_tree,
  main = "General NNPAS: cross-validation error by tree size"
)

rpart.plot(
  general_tree$pruned_tree,
  type = 4,
  extra = 104,
  fallen.leaves = TRUE,
  main = "General NNPAS: pruned classification tree"
)

# Simple lecture-style scatter plot. The black lines are descriptive reference
# thresholds (20 and 30 g/day fibre; ages 45 and 65), not estimated tree splits.
general_tree_scatter <- ggplot(
  general_data,
  aes(
    x = 20 + 5 * fibre5,
    y = 45 + 10 * age10,
    colour = obese
  )
) +
  geom_point(alpha = 0.45, size = 1.6) +
  scale_colour_viridis(
    limits = c(0, 1),
    breaks = c(0, 1),
    labels = c("Not obese", "Obese"),
    name = "Observed outcome"
  ) +
  geom_vline(xintercept = 20) +
  annotate(
    "segment",
    x = 20,
    xend = max(20 + 5 * general_data$fibre5),
    y = 45,
    yend = 45
  ) +
  annotate(
    "segment",
    x = 20,
    xend = max(20 + 5 * general_data$fibre5),
    y = 65,
    yend = 65
  ) +
  annotate(
    "segment",
    x = 30,
    xend = 30,
    y = 45,
    yend = 65
  ) +
  facet_wrap(~activity) +
  labs(
    title = "General NNPAS: fibre and age by obesity status",
    x = "Two-day mean fibre intake (g/day)",
    y = "Age (years)"
  ) +
  theme_minimal()

print(general_tree_scatter)

# 5.2 Indigenous NATSINPAS decision tree
indigenous_tree <- fit_pruned_tree(indigenous_data)

plotcp(
  indigenous_tree$large_tree,
  main = "Indigenous NATSINPAS: cross-validation error by tree size"
)

rpart.plot(
  indigenous_tree$pruned_tree,
  type = 4,
  extra = 104,
  fallen.leaves = TRUE,
  main = "Indigenous NATSINPAS: pruned classification tree"
)

indigenous_tree_scatter <- ggplot(
  indigenous_data,
  aes(
    x = 20 + 5 * fibre5,
    y = 45 + 10 * age10,
    colour = obese
  )
) +
  geom_point(alpha = 0.55, size = 2) +
  scale_colour_viridis(
    limits = c(0, 1),
    breaks = c(0, 1),
    labels = c("Not obese", "Obese"),
    name = "Observed outcome"
  ) +
  geom_vline(xintercept = 20) +
  annotate(
    "segment",
    x = 20,
    xend = max(20 + 5 * indigenous_data$fibre5),
    y = 45,
    yend = 45
  ) +
  annotate(
    "segment",
    x = 20,
    xend = max(20 + 5 * indigenous_data$fibre5),
    y = 65,
    yend = 65
  ) +
  annotate(
    "segment",
    x = 30,
    xend = 30,
    y = 45,
    yend = 65
  ) +
  facet_wrap(~activity) +
  labs(
    title = "Indigenous NATSINPAS: fibre and age by obesity status",
    x = "Two-day mean fibre intake (g/day)",
    y = "Age (years)"
  ) +
  theme_minimal()

print(indigenous_tree_scatter)

# ============================================================
# 6. Random forest and plots
# ============================================================

fit_random_forest <- function(df, seed) {
  if (!is.data.frame(df)) {
    stop("fit_random_forest() requires a data.frame. Run Sections 1-3 first.")
  }

  set.seed(seed)

  predictor_matrix <- model.matrix(
    delete.response(terms(tree_formula)),
    data = df
  )

  number_predictors <- max(1L, ncol(predictor_matrix) - 1L)

  randomForest(
    tree_formula,
    data = df,
    ntree = 500,
    mtry = max(1L, floor(sqrt(number_predictors))),
    importance = TRUE,
    proximity = FALSE,
    na.action = na.omit
  )
}

# 6.1 General NNPAS random forest
general_rf <- fit_random_forest(general_data, seed = 6101)

plot(
  general_rf,
  main = "General NNPAS: random forest OOB error"
)

varImpPlot(
  general_rf,
  type = 1,
  main = "General NNPAS: permutation variable importance"
)

# 6.2 Indigenous NATSINPAS random forest
indigenous_rf <- fit_random_forest(indigenous_data, seed = 6102)

plot(
  indigenous_rf,
  main = "Indigenous NATSINPAS: random forest OOB error"
)

varImpPlot(
  indigenous_rf,
  type = 1,
  main = "Indigenous NATSINPAS: permutation variable importance"
)

# randomForest does not use the ABS survey weights in the same way as svyglm.
# Treat these fits as predictive assumption checks, not population estimates.

# ============================================================
# 7. Generalized additive model and plots
# ============================================================

gam_formula <- obese ~
  activity +
  s(fibre5, by = activity, k = 6, bs = "cr") +
  s(age10, k = 6, bs = "cr") +
  s(energy_mj, k = 6, bs = "cr") +
  sex_factor +
  ses +
  income +
  year12_factor +
  smoking +
  alcohol

general_gam <- gam(
  gam_formula,
  data = general_data,
  weights = normalise_weights(general_data$survey_weight),
  family = quasibinomial(),
  method = "REML"
)

indigenous_gam <- gam(
  gam_formula,
  data = indigenous_data,
  weights = normalise_weights(indigenous_data$survey_weight),
  family = quasibinomial(),
  method = "REML"
)

# 7.1 Lecture-style GAM smooths and diagnostics: General NNPAS
plot(
  general_gam,
  pages = 1,
  shade = TRUE,
  residuals = FALSE,
  scale = 0
)

gam.check(general_gam)

# 7.2 Lecture-style GAM smooths and diagnostics: Indigenous NATSINPAS
plot(
  indigenous_gam,
  pages = 1,
  shade = TRUE,
  residuals = FALSE,
  scale = 0
)

gam.check(indigenous_gam)

# 7.3 Adjusted GAM probability curves
modal_value <- function(x) {
  names(which.max(table(x)))
}

make_gam_prediction_data <- function(model, data, population_name) {
  fibre_sequence <- seq(
    quantile(data$fibre5, 0.02),
    quantile(data$fibre5, 0.98),
    length.out = 120
  )

  prediction_data <- expand.grid(
    fibre5 = fibre_sequence,
    activity = factor(
      c("Below guideline", "Met guideline"),
      levels = levels(data$activity)
    )
  )

  prediction_data$age10 <- weighted.mean(
    data$age10,
    data$survey_weight
  )

  prediction_data$energy_mj <- weighted.mean(
    data$energy_mj,
    data$survey_weight
  )

  prediction_data$sex_factor <- factor(
    modal_value(data$sex_factor),
    levels = levels(data$sex_factor)
  )

  prediction_data$ses <- factor(
    modal_value(data$ses),
    levels = levels(data$ses)
  )

  prediction_data$income <- weighted.mean(
    data$income,
    data$survey_weight
  )

  prediction_data$year12_factor <- factor(
    modal_value(data$year12_factor),
    levels = levels(data$year12_factor)
  )

  prediction_data$smoking <- factor(
    modal_value(data$smoking),
    levels = levels(data$smoking)
  )

  prediction_data$alcohol <- weighted.mean(
    data$alcohol,
    data$survey_weight
  )

  prediction <- predict(
    model,
    newdata = prediction_data,
    type = "link",
    se.fit = TRUE
  )

  prediction_data$Predicted_probability <- plogis(prediction$fit)
  prediction_data$Lower <- plogis(
    prediction$fit - 1.96 * prediction$se.fit
  )
  prediction_data$Upper <- plogis(
    prediction$fit + 1.96 * prediction$se.fit
  )
  prediction_data$Fibre_g_day <- 20 + 5 * prediction_data$fibre5
  prediction_data$Population <- population_name

  prediction_data
}

gam_prediction_data <- bind_rows(
  make_gam_prediction_data(
    general_gam,
    general_data,
    "General NNPAS"
  ),
  make_gam_prediction_data(
    indigenous_gam,
    indigenous_data,
    "Indigenous NATSINPAS"
  )
)

gam_probability_plot <- ggplot(
  gam_prediction_data,
  aes(
    x = Fibre_g_day,
    y = Predicted_probability,
    colour = activity,
    fill = activity
  )
) +
  geom_ribbon(
    aes(ymin = Lower, ymax = Upper),
    alpha = 0.15,
    colour = NA
  ) +
  geom_line(linewidth = 1) +
  facet_wrap(~Population, scales = "free_x") +
  labs(
    title = "GAM: fibre intake, physical activity and obesity",
    subtitle = "Adjusted nonlinear predicted probability curves",
    x = "Two-day mean fibre intake (g/day)",
    y = "Predicted probability of obesity",
    colour = "Physical activity",
    fill = "Physical activity"
  ) +
  theme_minimal()

print(gam_probability_plot)

# ============================================================
# 8. Stratified ten-fold cross-validation and plot
# ============================================================

clip_probability <- function(probability, epsilon = 1e-6) {
  pmin(pmax(as.numeric(probability), epsilon), 1 - epsilon)
}

weighted_auc <- function(outcome, probability, weights) {
  outcome <- as.integer(outcome)
  probability <- as.numeric(probability)
  weights <- as.numeric(weights)

  if (sum(weights[outcome == 1]) == 0 ||
      sum(weights[outcome == 0]) == 0) {
    return(NA_real_)
  }

  probability_table <- data.frame(
    probability = probability,
    positive_weight = weights * outcome,
    negative_weight = weights * (1 - outcome)
  ) %>%
    group_by(probability) %>%
    summarise(
      positive_weight = sum(positive_weight),
      negative_weight = sum(negative_weight),
      .groups = "drop"
    ) %>%
    arrange(probability) %>%
    mutate(
      negative_below = cumsum(negative_weight) - negative_weight
    )

  numerator <- sum(
    probability_table$positive_weight *
      (
        probability_table$negative_below +
          0.5 * probability_table$negative_weight
      )
  )

  denominator <-
    sum(probability_table$positive_weight) *
    sum(probability_table$negative_weight)

  numerator / denominator
}

calculate_metrics <- function(outcome, probability, weights) {
  probability <- clip_probability(probability)

  data.frame(
    AUC = weighted_auc(outcome, probability, weights),
    Brier = weighted.mean(
      (outcome - probability)^2,
      weights
    ),
    Log_loss = weighted.mean(
      -(
        outcome * log(probability) +
          (1 - outcome) * log(1 - probability)
      ),
      weights
    )
  )
}

run_ten_fold_cv <- function(data, population_name, seed) {
  set.seed(seed)

  test_folds <- createFolds(
    data$obesity,
    k = 10,
    list = TRUE,
    returnTrain = FALSE
  )

  fold_results <- vector("list", length(test_folds))

  for (fold_number in seq_along(test_folds)) {
    test_index <- test_folds[[fold_number]]
    train_index <- setdiff(seq_len(nrow(data)), test_index)

    training_data <- data[train_index, , drop = FALSE]
    testing_data <- data[test_index, , drop = FALSE]

    # Store the normalised weights inside the training data. glm() and gam()
    # evaluate weights through model.frame; a separate local object such as
    # training_weights may not be visible in that evaluation environment.
    training_data$cv_weight <- normalise_weights(
      training_data$survey_weight
    )

    # GLM baseline
    cv_glm <- glm(
      baseline_formula,
      data = training_data,
      weights = cv_weight,
      family = quasibinomial()
    )

    glm_probability <- predict(
      cv_glm,
      newdata = testing_data,
      type = "response"
    )

    # Decision tree
    cv_tree <- fit_pruned_tree(training_data)$pruned_tree

    tree_probability <- predict(
      cv_tree,
      newdata = testing_data,
      type = "prob"
    )[, "Obese"]

    # Random forest
    cv_rf <- fit_random_forest(
      training_data,
      seed = seed + fold_number
    )

    rf_probability <- predict(
      cv_rf,
      newdata = testing_data,
      type = "prob"
    )[, "Obese"]

    # GAM
    cv_gam <- gam(
      gam_formula,
      data = training_data,
      weights = cv_weight,
      family = quasibinomial(),
      method = "REML"
    )

    gam_probability <- predict(
      cv_gam,
      newdata = testing_data,
      type = "response"
    )

    model_probabilities <- list(
      GLM = glm_probability,
      `Decision tree` = tree_probability,
      `Random forest` = rf_probability,
      GAM = gam_probability
    )

    fold_results[[fold_number]] <- bind_rows(
      lapply(names(model_probabilities), function(model_name) {
        cbind(
          data.frame(
            Population = population_name,
            Fold = fold_number,
            Model = model_name
          ),
          calculate_metrics(
            outcome = testing_data$obese,
            probability = model_probabilities[[model_name]],
            weights = testing_data$survey_weight
          )
        )
      })
    )
  }

  bind_rows(fold_results)
}

general_cv_results <- run_ten_fold_cv(
  general_data,
  population_name = "General NNPAS",
  seed = 8101
)

indigenous_cv_results <- run_ten_fold_cv(
  indigenous_data,
  population_name = "Indigenous NATSINPAS",
  seed = 8102
)

cv_results <- bind_rows(
  general_cv_results,
  indigenous_cv_results
)

cv_summary <- cv_results %>%
  group_by(Population, Model) %>%
  summarise(
    Mean_AUC = mean(AUC, na.rm = TRUE),
    SD_AUC = sd(AUC, na.rm = TRUE),
    Mean_Brier = mean(Brier, na.rm = TRUE),
    SD_Brier = sd(Brier, na.rm = TRUE),
    Mean_log_loss = mean(Log_loss, na.rm = TRUE),
    SD_log_loss = sd(Log_loss, na.rm = TRUE),
    .groups = "drop"
  )

print(cv_summary)

cv_plot_data <- bind_rows(
  cv_results %>%
    transmute(
      Population,
      Fold,
      Model,
      Metric = "AUC (higher is better)",
      Score = AUC
    ),
  cv_results %>%
    transmute(
      Population,
      Fold,
      Model,
      Metric = "Brier score (lower is better)",
      Score = Brier
    ),
  cv_results %>%
    transmute(
      Population,
      Fold,
      Model,
      Metric = "Log loss (lower is better)",
      Score = Log_loss
    )
)

cv_comparison_plot <- ggplot(
  cv_plot_data,
  aes(x = Model, y = Score, fill = Model)
) +
  geom_boxplot(outlier.alpha = 0.45) +
  facet_grid(
    Population ~ Metric,
    scales = "free_y"
  ) +
  labs(
    title = "Ten-fold cross-validation model comparison",
    x = NULL,
    y = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

print(cv_comparison_plot)

# Interpretation guide:
# 1. Higher AUC indicates better discrimination.
# 2. Lower Brier score and log loss indicate better probability prediction.
# 3. If Tree, RF or GAM do not consistently improve on GLM, the simpler
#    Person 1-style model is reasonably supported as a predictive form.
# 4. GAM EDF values clearly above 1 suggest nonlinearity.
# 5. RF importance is predictive importance, not a causal effect estimate.
