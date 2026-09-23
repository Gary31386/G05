# Rosetta — RQ3 methods and results

Updated 23 September 2026. This file is the handoff for the RQ3 statistical contribution. The primary analysis uses the two-day mean dietary intake specified in the group proposal and confirmed with the tutor. Day 1 analyses are retained as sensitivity checks.

## Research question and scope

Among adults aged 18–65 in the 2011–12 general-population NNPAS sample and the 2012–13 non-remote Aboriginal and Torres Strait Islander NATSINPAS sample, does activity in the previous week modify the association between reported fibre intake and measured obesity? The two survey samples are fitted **separately**. NNPAS has no Indigenous-status variable in the available Basic CURF, so its participants must not be labelled a verified non-Indigenous group. The surveys also differ in year and geographic coverage. [ABS NNPAS scope](https://www.abs.gov.au/methodologies/food-and-nutrients/2011-12); [ABS NATSINPAS microdata](https://www.abs.gov.au/statistics/microdata-tablebuilder/available-microdata-tablebuilder/australian-aboriginal-and-torres-strait-islander-health-survey-nutrition-and-physical-activity).

## Primary model

The exposure is the mean of Day 1 and Day 2 total dietary fibre (`FIBRET1` and `FIBRET2`), modelled continuously per **5 g/day**. The outcome is measured BMI **≥30 kg/m²** (`BMISC`), with invalid BMI codes excluded during nutrition cleaning. The modifier is `EXREGUIC`: **at least 150 minutes of physical activity in the preceding week**, coded 1 versus 2. It does **not** require five sessions; `EXREGUID`, which does, is tested separately. The script matches every cleaned participant ID back to the authorised ABS person file and verifies the PA codes and survey weight. [ABS NNPAS activity methodology](https://www.abs.gov.au/methodologies/australian-health-survey-physical-activity-methodology/2011-12).

Each survey has its own weighted logistic regression with `fibre5 × activity`, age and age squared, sex, two-day mean total energy, SEIFA disadvantage quintile, household income decile, Year 12 completion, smoking status and two-day mean alcohol. This extended set follows the assigned full-confounder model. Each fibre odds ratio therefore compares people reporting 5 g/day more fibre **at the same reported energy and included covariates**. These are observational associations, not causal effects. Medication, chronic disease and broader dietary pattern cannot be consistently adjusted across the available files; reporting bias and residual confounding remain possible.

For the two-day NNPAS sample, the model uses `NPAD2WGT` with 60 `W2M0101`–`W2M0160` Day 2 replicates. NATSINPAS has no Day 2-specific weight in the Basic CURF, so the model uses `IPAFINWT` with 250 `WPM1001`–`WPM1250` person replicates. The replicate-weight variance uses the group-jackknife factor `(G−1)/G`; intervals and p-values use a normal approximation. Complete cases have no missing model variables; pregnancy is excluded where recorded.

The proposal and tutor direction specify a two-day mean exposure, so it remains the primary analysis. The weight limitation changes the interpretation, not the planned exposure definition: the NATSINPAS estimate is a weighted analysis of non-remote Day 2 completers and is not fully representative of all non-remote adults. ABS reports 771 Day 2 completers across all ages; the cleaned file contains 492 adults aged 18–65 before later BMI, PA and model-covariate exclusions. [ABS NATSINPAS Explanatory Notes](https://www.abs.gov.au/AUSSTATS/abs%40.nsf/Lookup/4727.0.55.008Explanatory%20Notes12012-13).

## Sample flow and primary result

| Stage | General NNPAS | Indigenous NATSINPAS, non-remote |
| --- | ---: | ---: |
| Aged 18–65 with two dietary recalls | 4,824 | 492 |
| Valid measured BMI | 4,271 | 451 |
| Known 150-minute activity status | 4,239 | 448 |
| Complete-case main model | **3,960** | **408** |
| Obesity cases in model | **1,055** | **190** |

The pregnancy rule removed no further records from the cleaned eligible samples. The remaining complete-case exclusions arise mainly from the extended adjustment set and missing two-day energy or alcohol.

| Estimate, per 5 g/day higher reported fibre | General NNPAS, n=3,960 | Indigenous NATSINPAS, n=408 |
| --- | ---: | ---: |
| Below 150 min/week | **0.84** (95% CI 0.77–0.93) | **1.06** (0.75–1.50) |
| At least 150 min/week | **0.81** (0.73–0.91) | **0.84** (0.59–1.20) |
| **Interaction: ratio of fibre odds ratios** | **0.97** (0.86–1.09), p=0.559 | **0.79** (0.54–1.16), p=0.232 |

An interaction ratio of 1 means equal fibre odds ratios across the two activity groups on the multiplicative odds scale. Absolute probability differences were not estimated, so this analysis does not answer whether modification exists on an absolute-risk scale. These intervals provide no clear evidence of modification in either survey. They do not establish that the underlying interaction is exactly zero. The two survey estimates are not a formal Indigenous versus non-Indigenous comparison.

Weighted obesity percentages within the complete-case primary samples are 25.6% and 42.3%. They are unadjusted, drawn from different survey frames, and are not age-standardised group contrasts.

## Sensitivity analyses

All three rows use the same extended adjustment set and the 150-minute activity variable. Rows 1 and 2 use the **exact same Day 2 completers and replicate weights**, so the comparison changes fibre, energy and alcohol recall definitions together. Row 3 expands to all eligible Day 1 respondents and, for NNPAS, changes from Day 2 weights to person weights. The comparisons show where estimates move; they do not isolate a causal source.

| Dietary sample and recall definition | NNPAS: n; interaction OR (95% CI); p | NATSINPAS: n; interaction OR (95% CI); p |
| --- | ---: | ---: |
| **Day 2 completers, mean of both days — primary** | 3,960; **0.965** (0.857–1.087); 0.559 | 408; **0.790** (0.536–1.163); 0.232 |
| Same Day 2 completers, Day 1 intake | 3,960; **0.981** (0.883–1.091); 0.724 | 408; **0.850** (0.581–1.242); 0.400 |
| All eligible Day 1 respondents — sensitivity | 5,863; **0.985** (0.912–1.063); 0.699 | 862; **0.967** (0.787–1.189); 0.751 |

The NATSINPAS estimate moves more across recall definitions than the NNPAS estimate, while every interval still includes 1. Its person weight does not correct for selection into the second interview, so the primary NATSINPAS line is a completer analysis. The Day 1 sensitivity has more day-to-day measurement noise than a two-day mean; neither recall definition measures long-term habitual fibre intake directly.

| Additional check | General NNPAS | Indigenous NATSINPAS |
| --- | ---: | ---: |
| Two-day model: activity requires 150 minutes **and five sessions** (`EXREGUID`) | n=3,953; 0.944 (0.846–1.052), p=0.297 | n=408; 0.886 (0.629–1.249), p=0.489 |
| All-respondent Day 1 model, core confounders only | n=6,364; 0.994 (0.926–1.068), p=0.877 | n=974; 0.951 (0.789–1.146), p=0.598 |

The Day 1 core model keeps age, age squared, sex, energy, SEIFA, smoking and alcohol and retains people whose income or schooling is missing. Neither sensitivity changes the conclusion about lack of clear statistical evidence for modification. Amelie's planned imputation analysis remains a separate task.

## Primary-model checks and limits

| Check | General NNPAS | Indigenous NATSINPAS |
| --- | ---: | ---: |
| Weighted Brier score | 0.178 | 0.211 |
| Model-matrix condition number | 22.7 | 26.2 |
| Joint quadratic fibre screen | p=0.766 | p=0.282 |
| Exclude top 1% of fibre: interaction OR | 0.961 (0.851–1.085) | 0.848 (0.553–1.301) |

The upper-tail exclusion and quadratic screen do not alter the main interaction conclusion, but neither proves the linear form correct. Self-reported activity, measurement and reporting error, complete-case selection and unmeasured confounding limit interpretation. The NATSINPAS Day 2-completer selection and lack of a matching Day 2 weight are the main additional limitations.

## Reproduce and data handling

Use an R environment with `survey` and `readxl`. Keep the original authorised ABS Basic CURF files outside the repository. Set `RQ3_RAW_DATA_DIR` to the folder containing `AHSnpa11bp.csv`, `inp13bsp.csv` and `inp13bhh.csv`. Run `Rscript rq3-rosetta-main.R`, then `rq3-rosetta-day1-in-day2-sensitivity.R`, `rq3-rosetta-day1-sensitivity.R` and `rq3-rosetta-coreconf-sensitivity.R`. `RQ3_OUTPUT_DIR` can redirect outputs. All scripts write aggregate tables only. The two pre-existing cleaned workbooks in this repository contain respondent-level records; keep repository access restricted to authorised course collaborators and do not export those rows in the knitted report.
