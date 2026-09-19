# Rosetta — RQ3 methods and results

Updated 19 September 2026. This file is the handoff for Rosetta's statistical contribution. The primary analysis now uses Day 1 dietary intake in both surveys. The earlier two-day analysis remains fully reported as a sensitivity analysis.

## Research question and scope

Among adults aged 18–65 in the 2011–12 general-population NNPAS sample and the 2012–13 non-remote Aboriginal and Torres Strait Islander NATSINPAS sample, does activity in the previous week modify the association between reported fibre intake and measured obesity? The two survey samples are fitted **separately**. NNPAS has no Indigenous-status variable in the available Basic CURF, so its participants must not be labelled a verified non-Indigenous group. The surveys also differ in year and geographic coverage. [ABS NNPAS scope](https://www.abs.gov.au/methodologies/food-and-nutrients/2011-12); [ABS NATSINPAS microdata](https://www.abs.gov.au/statistics/microdata-tablebuilder/available-microdata-tablebuilder/australian-aboriginal-and-torres-strait-islander-health-survey-nutrition-and-physical-activity).

## Primary model

The exposure is Day 1 total dietary fibre (`FIBRET1`), modelled continuously per **5 g/day**. The outcome is measured BMI **≥30 kg/m²** (`BMISC`), with codes 0, 98 and 99 excluded. The modifier is `EXREGUIC`: **at least 150 minutes of physical activity in the preceding week**, coded 1 versus 2 in both ABS data item lists. It does **not** require five sessions; `EXREGUID`, which does, is tested separately. The primary model does not use a generic “met physical activity guidelines” label. The original codebooks are `AHSnpaBasicCURFdataItemList.xls` and `AHSinpBasicCURFdataItemList.xls`; the two-day script additionally verifies the cleaned labels against raw codes. [ABS NNPAS activity methodology](https://www.abs.gov.au/methodologies/australian-health-survey-physical-activity-methodology/2011-12).

Each survey has its own weighted logistic regression with `fibre5 × activity`, age and age squared, sex, Day 1 total energy, SEIFA disadvantage quintile, household income decile, Year 12 completion, smoking status and Day 1 alcohol. This extended set follows Rosetta's assigned full-confounder model; the core set below checks what changes when income and schooling are omitted and more people remain in the analysis. Thus each fibre odds ratio compares people reporting 5 g/day more fibre **at the same reported energy and included covariates**. These are observational associations, not causal effects. Medication, chronic disease and broader dietary pattern cannot be consistently adjusted across the available files; reporting bias and residual confounding remain possible.

For all Day 1 respondents, NNPAS uses `NPAFINWT` with 60 `WPM0101`–`WPM0160` person replicates; NATSINPAS uses `IPAFINWT` with 250 `WPM1001`–`WPM1250` person replicates. The replicate-weight variance uses the group-jackknife factor `(G−1)/G`; intervals and p-values use a normal approximation (1.96 standard errors). The NNPAS fitted design has 43 residual degrees of freedom: using a t interval instead would increase the critical value from 1.960 to 2.017, about 3%, without changing the interaction conclusion. This matches the person-weight definitions and replicate method in the ABS data item lists and [NATSINPAS technical note](https://www.abs.gov.au/ausstats/abs%40.nsf/Lookup/4727.0.55.002Technical%2BNote12012-13). Complete cases have no missing model variables; pregnancy is excluded where recorded.

The switch from the earlier two-day main analysis was prompted by the absence of a Day 2-specific weight in the NATSINPAS Basic CURF. Preliminary Day 1 estimates were already available when the revised plan was finalised; the choice of primary analysis follows weight and population coverage, and both earlier results remain visible below. NNPAS has an identified Day 2 weight (`NPAD2WGT`) and 60 `W2M0101`–`W2M0160` replicates. NATSINPAS has only its person weight for the smaller Day 2 completer subset; its two-day estimate therefore has weaker population interpretation. ABS published NATSINPAS nutrition estimates from Day 1 and reported that only 771 respondents completed Day 2. [ABS NATSINPAS nutrition notes](https://www.abs.gov.au/ausstats/abs@.nsf/Latestproducts/4727.0.55.005Explanatory%20Notes12012-13).

## Sample flow and primary result

| Stage | General NNPAS | Indigenous NATSINPAS, non-remote |
| --- | ---: | ---: |
| Aged 18–65 | 7,528 | 1,119 |
| Valid measured BMI | 6,415 | 980 |
| After pregnancy exclusion | 6,415 | 980 |
| Complete-case main model | **5,863** | **862** |
| Obesity cases in model | **1,585** | **374** |

The pregnancy rule removed no further records after requiring measured BMI. Most remaining complete-case exclusions arise from the extended adjustment set. The core-covariate sensitivity below assesses its effect.

| Estimate, per 5 g/day higher reported fibre | General NNPAS, n=5,863 | Indigenous NATSINPAS, n=862 |
| --- | ---: | ---: |
| Below 150 min/week | **0.91** (95% CI 0.85–0.96) | **1.01** (0.85–1.18) |
| At least 150 min/week | **0.89** (0.84–0.95) | **0.97** (0.82–1.15) |
| **Interaction: ratio of fibre odds ratios** | **0.99** (0.91–1.06), p=0.699 | **0.97** (0.79–1.19), p=0.751 |

An interaction ratio of 1 means equal fibre odds ratios across the two activity groups on the multiplicative odds scale. Absolute probability differences were not estimated, so this analysis does not answer whether modification exists on an absolute-risk scale. These intervals provide no clear evidence of modification in either survey. They do not establish that the underlying interaction is exactly zero. The two survey estimates are not a formal Indigenous versus non-Indigenous comparison.

Weighted obesity percentages within the complete-case primary samples are 25.8% and 41.9%. They are unadjusted, drawn from different survey frames, and are not age-standardised group contrasts.

## Sensitivity analyses

The table keeps the original two-day result in full. All three rows use the same extended adjustment set and the 150-minute activity variable. Row 2 restricts to the **exact** Row 3 respondents, so the two rows have identical sample size and replicate weights; their dietary fibre, energy and alcohol measurements change together. Row 1 versus Row 2 also changes sample and, in NNPAS, the applicable person versus Day 2 weighting. The comparisons identify where estimates move; they do not isolate a causal source of movement.

| Dietary sample and recall definition | NNPAS: n; interaction OR (95% CI); p | NATSINPAS: n; interaction OR (95% CI); p |
| --- | ---: | ---: |
| **All eligible Day 1 respondents — primary** | 5,863; **0.985** (0.912–1.063); 0.699 | 862; **0.967** (0.787–1.189); 0.751 |
| Day 2 completers, Day 1 intake | 3,960; **0.981** (0.883–1.091); 0.724 | 408; **0.850** (0.581–1.242); 0.400 |
| Day 2 completers, mean of both days | 3,960; **0.965** (0.857–1.087); 0.559 | 408; **0.790** (0.536–1.163); 0.232 |

The NATSINPAS two-day subset's point estimate moves more than the NNPAS estimate, while every interval still includes 1. The NATSINPAS Day 2 person weight does not correct for selection into the second interview, so its two-day line should be interpreted as a completer analysis. The Day 1 primary measure has more day-to-day measurement noise than a two-day mean; this can attenuate fibre associations toward 1, and differing noise by activity group could also affect the interaction. Neither recall definition measures long-term habitual fibre intake directly.

| Additional check on Day 1 model | General NNPAS | Indigenous NATSINPAS |
| --- | ---: | ---: |
| Activity requires 150 minutes **and five sessions** (`EXREGUID`): interaction OR | n=5,848; 0.982 (0.905–1.065), p=0.664 | n=860; 1.031 (0.835–1.273), p=0.778 |
| Core confounders only (exclude income and Year 12): interaction OR | n=6,364; 0.994 (0.926–1.068), p=0.877 | n=974; 0.951 (0.789–1.146), p=0.598 |

The core model keeps age, age squared, sex, energy, SEIFA, smoking and alcohol. It retains people whose income or schooling is missing. Neither sensitivity changes the conclusion about lack of clear statistical evidence for modification. Amelie's planned imputation analysis remains a separate task.

## Primary-model checks and limits

| Check | General NNPAS | Indigenous NATSINPAS |
| --- | ---: | ---: |
| Weighted Brier score | 0.181 | 0.220 |
| Model-matrix condition number | 26.6 | 26.9 |
| Joint quadratic fibre screen | p=0.847 | p=0.906 |
| Exclude top 1% of fibre: interaction OR | 0.980 (0.909–1.057) | 1.016 (0.808–1.276) |

The upper-tail exclusion and quadratic screen do not alter the main interaction conclusion, but neither proves the linear form correct. Single-day diet recall, self-reported activity, measurement and reporting error, complete-case selection and unmeasured confounding limit interpretation. The old two-day diagnostic tables are retained separately and must not be copied into the new primary results.

## Reproduce and data handling

Use an R environment with `survey` and `readxl`. Keep the original authorised ABS Basic CURF files outside the repository. Set `RQ3_RAW_DATA_DIR` to the folder containing `AHSnpa11bp.csv`, `inp13bsp.csv` and `inp13bhh.csv`. Run `Rscript rq3-rosetta-main.R`, then `rq3-rosetta-day1-in-day2-sensitivity.R`, `rq3-rosetta-two-day-sensitivity.R` and `rq3-rosetta-coreconf-sensitivity.R`. `RQ3_OUTPUT_DIR` can redirect outputs. All scripts write aggregate tables only. The two pre-existing cleaned workbooks in this repository contain respondent-level records; keep repository access restricted to authorised course collaborators and do not export those rows in the knitted report.
