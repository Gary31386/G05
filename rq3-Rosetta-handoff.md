# Rosetta — STAT3888 RQ3 handoff (current)

Updated 19 September 2026. Rosetta's assigned model, diagnostics and formatted results are complete in G05 `main`. The analysis revision is commit `caa24ea`; the full specification and tables are in [rq3-rosetta-methods-results.md](rq3-rosetta-methods-results.md). The shared R Markdown report reads the aggregate result tables and renders successfully.

## Use this analysis

- **Primary:** [rq3-rosetta-main.R](rq3-rosetta-main.R). Separately fit survey-weighted logistic regressions in adults aged 18–65. Outcome: measured BMI ≥30. Exposure: Day 1 total fibre per 5 g/day. Modifier: at least 150 minutes of activity in the preceding week (`EXREGUIC`). Adjust for age, age squared, sex, Day 1 energy, SEIFA, income decile, Year 12, smoking and Day 1 alcohol.
- **Weights:** NNPAS `NPAFINWT` + 60 `WPM01` replicates; NATSINPAS `IPAFINWT` + 250 `WPM1` replicates. The NATSINPAS activity sample is non-remote. NNPAS is **general population**, not a verified non-Indigenous sample.
- **Primary interaction:** NNPAS n=5,863, OR **0.985** (95% CI 0.912–1.063), p=0.699; NATSINPAS n=862, OR **0.967** (0.787–1.189), p=0.751. Both intervals include 1, so there is no clear evidence of modification on the odds-ratio scale.
- **Why Day 1 is primary:** the NATSINPAS Basic CURF has no Day 2-specific weight. A Day 1 recall is noisier than a two-day average, so both approaches and their limitations are reported. Do not describe the two-day NATSINPAS result as representative of all non-remote adults.

## Sensitivities retained

| Dietary definition | NNPAS interaction OR (95% CI), n | NATSINPAS interaction OR (95% CI), n |
| --- | --- | --- |
| All Day 1 respondents — primary | 0.985 (0.912–1.063), 5,863 | 0.967 (0.787–1.189), 862 |
| Day 2 completers, using their Day 1 intake | 0.981 (0.883–1.091), 3,960 | 0.850 (0.581–1.242), 408 |
| Same Day 2 completers, two-day mean intake | 0.965 (0.857–1.087), 3,960 | 0.790 (0.536–1.163), 408 |

The first comparison changes sample and, for NNPAS, the applicable weights. The last two rows use the same people and weights, changing fibre, energy and alcohol recall measures together. The five-session activity definition (`EXREGUID`), reduced-confounder model, upper-tail exclusion and quadratic screen are in the full methods file and matching CSVs.

## Reproduce and hand over

Use authorised local Basic CURF files `AHSnpa11bp.csv`, `inp13bsp.csv`, `inp13bhh.csv`; set `RQ3_RAW_DATA_DIR` to their folder. Run `Rscript rq3-rosetta-main.R`, then the `rq3-rosetta-day1-in-day2-sensitivity.R`, `rq3-rosetta-two-day-sensitivity.R` and `rq3-rosetta-coreconf-sensitivity.R` scripts. Set `RQ3_OUTPUT_DIR` to redirect output. Only aggregate CSVs are generated. The two existing cleaned workbooks in G05 contain respondent-level records; keep access to authorised collaborators and do not print their rows in reports.

Rosetta's older local `rq3-main-two-surveys.R` and the previous version of this handoff are **historical**: they used ages 18+, fibre per 10 g/day and different counts. Do not mix their results into the current manuscript.

**Remaining group decisions:** correct the research-question wording to name NNPAS general population and NATSINPAS non-remote Aboriginal and Torres Strait Islander adults aged 18–65; incorporate Amelie's missing-data checks and Gary/Selina's assumption checks when available. Those checks may change the final interpretation but are outside Rosetta's assigned model delivery.
