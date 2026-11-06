# STAT3888 Group 5

## Does physical activity modify the association between dietary fibre intake and adiposity in Australian adults?

This repository contains the R code, analysis documents and results for our University of Sydney STAT3888 group project.

We investigated whether the association between dietary fibre intake and adiposity differs according to physical activity level.

## Data

We used data from the Australian Health Survey, analysing two cohorts separately.

The general cohort came from the 2011 to 2012 National Nutrition and Physical Activity Survey (NNPAS). The Indigenous cohort came from the 2012 to 2013 National Aboriginal and Torres Strait Islander Nutrition and Physical Activity Survey (NATSINPAS), restricted to non remote areas.

The final complete case analysis included 4,286 adults aged 18 to 65, consisting of 3,897 general cohort participants and 389 Indigenous cohort participants.

The general cohort was not restricted to non Indigenous adults. The two cohorts therefore do not represent a direct Indigenous versus non Indigenous comparison.

## Variables

The outcomes were BMI, waist circumference and obesity.

Obesity was defined as meeting both the BMI threshold of 30 kg/m² and the waist circumference threshold of 94 cm for men or 80 cm for women.

Fibre intake was averaged across two dietary recall days and analysed per 5 g/day. Physical activity was classified as at least 150 minutes per week or below 150 minutes per week.

Models adjusted for age group, sex, energy intake, household income, Year 12 completion, smoking status, alcohol intake and protein intake.

## Methods

We fitted survey weighted linear regression models for BMI and waist circumference, and logistic regression models for obesity. Each model included an interaction between fibre intake and physical activity. Jackknife replicate weights were used to estimate uncertainty.

Assumption checks examined fibre intake overlap between activity groups, missing data, multicollinearity, residual patterns and influential observations.

Robustness checks included removing influential observations, multiple imputation for missing covariates, and comparisons with unweighted and robust regression models.

We also compared obesity across four groups defined by fibre intake and physical activity. The primary fibre thresholds were 30 g/day for men and 25 g/day for women.

## Main findings

In the general cohort, higher fibre intake was associated with lower BMI, smaller waist circumference and lower odds of obesity.

The survey weighted models found no clear evidence that physical activity modified these associations in either cohort. Estimates in the Indigenous cohort were less precise because of the smaller sample.

Adults meeting both fibre and physical activity recommendations had the lowest estimated obesity prevalence. Some unweighted analyses gave different interaction results, so a small modifying effect cannot be ruled out.

These findings describe associations and do not establish causation.

## Files used in the final report

Our final report uses Amelie’s primary models and robustness analyses, together with Selina’s assumption checks.

| File or folder | Description |
| --- | --- |
| `amelie_run_robustness.R` | Amelie’s code for modelling and robustness analyses |
| `amelie_results/` | Saved results from Amelie’s analysis |
| `selina_before_and_after_fitting_assumption_checks.R` | Selina’s code for assumption checks and model diagnostics |
| `selina_assumption_results/` | Saved assumption check tables and figures |

Other files in the repository contain earlier work and were not used in the final report.

## Using the code

Download the repository and open `STAT3888 Groupwork.Rproj` in RStudio. Keep the repository as the working directory so that relative file paths work correctly.

To view the saved analysis, open `STAT3888-Project-Group-5.html` in a browser.

To render the shared analysis document, open `STAT3888 Project Group 5.Rmd` and click Knit. Install the packages loaded by the document before running it.

Refitting the survey models requires authorised access to the original ABS data files. Set their folder location in R before running scripts that require them:

```r
Sys.setenv(RQ3_RAW_DATA_DIR = "/path/to/authorised/ABS/files")
```

The original ABS files should be kept outside the repository.

## Group contributions

Anika developed the research aims, limitations and conclusion. Varty contributed the nutrition background, interpretation, abstract and final editing.

Amelie conducted the primary statistical analysis. Rosetta conducted the robustness and sensitivity analyses. Keran performed the group comparison. Selina conducted assumption checks and prepared the report formatting and editing.

## Limitations

The study used observational data collected at one point in time, so causal conclusions cannot be drawn.

Dietary recall errors, missing data and residual confounding may affect the estimates. The smaller Indigenous sample limited precision, and differences between the surveys restrict direct comparisons between cohorts.
