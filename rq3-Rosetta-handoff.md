# Rosetta — STAT3888 RQ3 handoff (current)

Updated 23 September 2026. The RQ3 model is aligned with the proposal and tutor direction: two-day mean fibre is primary and Day 1 is sensitivity. The full specification and tables are in [rq3-rosetta-methods-results.md](rq3-rosetta-methods-results.md). The shared R Markdown report reads aggregate result tables only.

## Use this analysis

- **Primary:** [rq3-rosetta-main.R](rq3-rosetta-main.R). Separately fit survey-weighted logistic regressions in adults aged 18–65 who completed both dietary recalls. Outcome: measured BMI ≥30. Exposure: two-day mean total fibre per 5 g/day. Modifier: at least 150 minutes of activity in the preceding week (`EXREGUIC`). Adjust for age, age squared, sex, two-day mean energy, SEIFA, income decile, Year 12, smoking and two-day mean alcohol.
- **Weights:** NNPAS `NPAD2WGT` + 60 `W2M01` Day 2 replicates; NATSINPAS `IPAFINWT` + 250 `WPM1` person replicates because the Basic CURF has no Day 2-specific weight.
- **Primary interaction:** NNPAS n=3,960, OR **0.965** (95% CI 0.857–1.087), p=0.559; NATSINPAS n=408, OR **0.790** (0.536–1.163), p=0.232.
- **Interpretation:** both intervals include 1, so there is no clear evidence of modification on the multiplicative odds scale. The NATSINPAS result is a weighted Day 2-completer analysis and must not be described as representative of all non-remote adults. NNPAS is **general population**, not a verified non-Indigenous sample.

## Sensitivities retained

| Dietary definition | NNPAS interaction OR (95% CI), n | NATSINPAS interaction OR (95% CI), n |
| --- | --- | --- |
| Day 2 completers, two-day mean intake — primary | 0.965 (0.857–1.087), 3,960 | 0.790 (0.536–1.163), 408 |
| Same Day 2 completers, using Day 1 intake | 0.981 (0.883–1.091), 3,960 | 0.850 (0.581–1.242), 408 |
| All Day 1 respondents — sensitivity | 0.985 (0.912–1.063), 5,863 | 0.967 (0.787–1.189), 862 |

The first two rows use the same people and weights, changing fibre, energy and alcohol recall measures together. The third row changes sample and, for NNPAS, the applicable weights. The five-session activity definition (`EXREGUID`), Day 1 reduced-confounder model, upper-tail exclusion and quadratic screen are in the full methods file and matching CSVs.

## Reproduce and hand over

Use authorised local Basic CURF files `AHSnpa11bp.csv`, `inp13bsp.csv`, `inp13bhh.csv`; set `RQ3_RAW_DATA_DIR` to their folder. Run `Rscript rq3-rosetta-main.R`, then `rq3-rosetta-day1-in-day2-sensitivity.R`, `rq3-rosetta-day1-sensitivity.R` and `rq3-rosetta-coreconf-sensitivity.R`. Set `RQ3_OUTPUT_DIR` to redirect output. Only aggregate CSVs are generated. The two existing cleaned workbooks in G05 contain respondent-level records; keep access to authorised collaborators and do not print their rows in reports.

Rosetta's older local `rq3-main-two-surveys.R` and the previous version of this handoff are **historical**: they used ages 18+, fibre per 10 g/day and different counts. Do not mix their results into the current manuscript.

**Remaining group decisions:** correct the research-question wording to name NNPAS general population and NATSINPAS non-remote Aboriginal and Torres Strait Islander adults aged 18–65; incorporate Amelie's missing-data checks and Gary/Selina's assumption checks when available. Those checks may change the final interpretation but are outside Rosetta's assigned model delivery.
