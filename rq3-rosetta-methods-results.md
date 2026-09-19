# Rosetta: RQ3 main model and results

Prepared 19 September 2026. The R script and four CSV tables in this repository are Rosetta's contribution to the main statistical analysis. The results were regenerated from the two `eligible_sample` sheets already in this repository. Only aggregate model output was added; the original ABS person files remain outside the repository.

## Population and variables

The cleaned datasets include adults aged **18–65** with two dietary recalls, measured BMI, known physical-activity guideline status and a positive survey weight. The model uses two-day mean total fibre intake as a continuous exposure, in **5 g/day** increments. Obesity means measured BMI **≥30 kg/m²**. Physical activity is the cleaned `pa_guideline` indicator: met versus did not meet the survey's 150-minute guideline. This is the effect modifier. The model is fitted separately in each survey.

The complete-case adjustment set is age and age squared, sex, two-day mean total energy intake, SEIFA disadvantage quintile, household income decile, Year 12 completion, smoking status and two-day mean alcohol intake. Physical activity also has a main effect in the interaction model. Broader diet, medication and disease history are not comparably available as defensible pre-exposure confounders in both cleaned files. An observational association from these data should not be presented as a causal effect.

The model uses the survey weights and jackknife replicate weights that match each sample: **NNPAS `NPAD2WGT` with 60 Day 2 replicates `W2M0101`–`W2M0160`**; **NATSINPAS `IPAFINWT` with 250 person replicates `WPM1001`–`WPM1250`**. The script checks every cleaned record against its original person record and verifies the full-sample weight. The ABS [NNPAS Basic CURF data item list](https://www.abs.gov.au/statistics/microdata-tablebuilder/available-microdata-tablebuilder/australian-health-survey-nutrition-and-physical-activity) identifies the Day 2 weight; the local data item lists identify the matching replicate fields. Jackknife variance uses the factor `(G−1)/G`; 95% confidence intervals use estimate ±1.96 standard errors on the log-odds scale.

## Sample flow

| Stage | General NNPAS | Indigenous NATSINPAS |
| --- | ---: | ---: |
| Cleaned eligible sample | 4,239 | 448 |
| After pregnancy exclusion | 4,239 | 448 |
| Complete main-model covariates | **3,960** | **408** |
| Obesity cases in model | **1,055** | **190** |

Most additional exclusions for the full model come from missing household income: 275 in NNPAS and 38 in NATSINPAS. Some other covariate missingness overlaps. These complete-case results need Amelie's planned missing-data sensitivity analysis.

## Adjusted results

An odds ratio (OR) below 1 means lower odds of obesity per 5 g/day higher reported fibre intake within the listed activity group. The interaction OR compares the fibre ORs between activity groups. An interaction OR of 1 means no multiplicative difference.

| Estimate | General NNPAS, n=3,960 | Indigenous NATSINPAS, n=408 |
| --- | ---: | ---: |
| Fibre +5 g/day, did not meet PA guideline | **0.84** (95% CI 0.77–0.93) | **1.06** (0.75–1.50) |
| Fibre +5 g/day, met PA guideline | **0.81** (0.73–0.91) | **0.84** (0.59–1.20) |
| **Fibre × PA interaction: ratio of ORs** | **0.97** (0.86–1.09), p=0.559 | **0.79** (0.54–1.16), p=0.232 |

Neither sample provides clear evidence that the fibre–obesity association differs by physical-activity status. The Indigenous estimate is imprecise, so its interval remains compatible with a meaningful difference in either direction. The interaction is the answer to RQ3; the fibre associations within each activity group are supporting estimates.

Weighted obesity prevalence in these complete-case samples is **25.6%** in General NNPAS and **42.3%** in NATSINPAS. These are unadjusted within-sample descriptions, not an age-standardised comparison or a causal effect of Indigenous status. The by-activity file gives the corresponding unadjusted subgroup figures.

## Diagnostics

| Check | General NNPAS | Indigenous NATSINPAS |
| --- | ---: | ---: |
| Weighted Brier score (squared prediction error) | 0.178 | 0.211 |
| Model-matrix condition number | 22.7 | 26.2 |
| Joint quadratic-fibre screen | p=0.766 | p=0.282 |
| Excluding the top 1% of fibre: interaction OR | 0.96 (0.85–1.09) | 0.85 (0.55–1.30) |

The outlier check does not materially change the interaction conclusion. The quadratic screen provides no strong evidence against a linear fibre term within these samples, but it cannot prove linearity. The small Indigenous complete-case sample limits precision. Fibre is the average of two recall days, not a direct measure of habitual intake. Physical activity is self-reported.

## Scope corrections for the manuscript

1. The NNPAS Basic CURF available here has **no Indigenous-status field**. Its sample is **General NNPAS**, not a verified non-Indigenous group. Do not label this a direct Indigenous versus non-Indigenous comparison.
2. The cleaned files restrict age to **18–65**, while the current research question says “adults.” State the age restriction in the question/methods, or rebuild the cleaning and models for all adults.
3. The Indigenous eligible sample is **non-remote**, because the activity guideline item used here was collected for non-remote adults. The surveys also have different sampling frames and years. Discuss them side by side, with a brief comparison only.
4. These are complete-case results. Multiple imputation and other robustness analyses are separate tasks and are not claimed here.

## Reproduce

Install or use existing R packages `readxl` and `survey`. Keep the authorised original ABS Basic CURF person files on your own computer. Set `RQ3_RAW_DATA_DIR` to the directory containing `AHSnpa11bp.csv` and `inp13bsp.csv`, then run `Rscript rq3-rosetta-main.R` from this repository. Set `RQ3_OUTPUT_DIR` only if you want output elsewhere. The script writes four aggregate CSVs: flow, model estimates, diagnostics and descriptive results by activity. Do not add original ABS person files or unredacted respondent records to GitHub.
