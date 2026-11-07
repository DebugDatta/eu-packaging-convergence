# EU Packaging Waste Convergence Analysis

This project analyzes packaging waste per capita across the 27 EU member states
from 2014 to 2023. It includes annual panel convergence models, country-level
cross-sections, dispersion and direction measures, Phillips-Sul tests, timing
tests, and diagnostics.

## Project files

- [`ecosem_analysis.py`](./ecosem_analysis.py): analysis script.
- [`ecosem_data.xlsx`](./ecosem_data.xlsx): input workbook and required data
  sheets.
- [`results/`](./results/): generated CSV output files.

## Data sources

| Sheet | Variable | Source | Notes |
|---|---|---|---|
| `packaging_waste` | Packaging waste generated, kg per capita, 27 EU countries, 2014 to 2023 | Eurostat, table `env_waspac`, Packaging waste by waste management operations (CC BY 4.0) | 8 of 270 country-years have no value. `obs_flag` keeps Eurostat's own flags (16 country-years). |
| `eu27_aggregate` | EU27 packaging waste, kg per capita, 2014 to 2023 | Eurostat, `env_waspac`, row "European Union - 27 countries (from 2020)" | Population-weighted Eurostat aggregate. Used only as a cross-check on direction. |
| `gdp_pps` | GDP per capita in purchasing power standards, EU27 = 100, 2014 to 2023 | Eustat (Basque Statistics Office), based on Eurostat Economic Accounts. [Source page](https://en.eustat.eus/elementos/ele0012300/ti_pib-per-capita-ppc-de-la-ca-de-euskadi-y-de-los-paises-de-la-union-europea-ue-27100-2010-2023/tbl0012365_i.html) | Complete for all 27 countries and years. The series is a republication, not Eurostat's own table. |
| `wtp`, column `wtp_2014_qa10` | Share who totally or tend to agree they would buy environmentally friendly products even if they cost a little more | Special Eurobarometer 416, item QA10, fieldwork April to May 2014. "Total Agree" from the Volume A results file `OP_VolumeAEB813ENVenvironment_v2.xls`, [EU Open Data Portal](https://data.europa.eu/88u/dataset/s2008_81_3_416). [Report](https://doi.org/10.2779/25662) | Whole-percent values stored as shares. The UK is excluded. |
| `wtp`, column `wtp_2024_qb8` | Share answering yes to paying more for products that are easier to repair, recyclable or sustainably produced (furniture, textiles, electronic devices) | Special Eurobarometer 550, item QB8, fieldwork March to April 2024. [Report](https://doi.org/10.2779/07854) | Different scale and goods from QA10, so the two columns are never compared over time. |

Flag codes in `obs_flag`: `b` break in time series, `d` definition differs,
`e` estimated, and combinations such as `be` and `bd`.

Data vintage: the Eurostat waste and EU27 values were checked against the live
`env_waspac` table on 7 November 2026 and match it, including the flags. Eurostat
revises this table, so later downloads can differ. The Eustat GDP table was last
updated by Eustat on 10 January 2025.

Country names follow Eurostat (for example Czechia). `iso2` uses the
Eurobarometer convention (Greece is EL).

## Sample

- 27 EU member states, 2014 to 2023.
- Possible yearly changes: 27 x 9 = 243. A change needs both consecutive
  years to be present.
- Usable changes: 232. Missing values: Bulgaria 2020 and 2023, Greece 2020,
  Czechia 2023, Hungary 2023, Poland 2020, 2021 and 2023.
- Cross-section: 27 countries. The end year is the last year with data (2023
  for 23 countries, 2022 for Bulgaria, Czechia, Hungary and Poland).
- Balanced panel (all ten years present): 22 countries. It is used for the
  Phillips-Sul test and as a check on sigma convergence.

## Variables

- W: packaging waste per capita (kg).
- Δln W(it) = ln W(it) − ln W(i,t−1), the annual growth in country i in year t.
- ln W(i,t−1): lagged log level.
- Y: GDP per capita in PPS (EU27 = 100), entered as ln Y(i,t−1) in the panel
  and ln Y(i,2014) in the cross-section.
- WTP: 2014 QA10 agreement share (main measure). QB8 2024 is a cross-sectional
  robustness check only.
- τ(t): year effects. Post(t) equals 1 for t at or after the break year, where
  t is the year of the change.

## Methodology

All panel regressions use ordinary least squares with standard errors clustered
by country (27 clusters). Cross-section regressions use classical standard
errors, with HC3 as a check.

P-values for clustered standard errors use the normal distribution (statsmodels
default).

1. Absolute beta convergence:
   `Δln W(it) = α + β ln W(i,t−1) + ε(it)`
2. Conditional convergence:
   `Δln W(it) = α + β ln W(i,t−1) + γ ln Y(i,t−1) + τ(t) + ε(it)`
3. Speed and half-life: `λ = −ln(1 + β)`,
   `half-life = ln 2 / λ`
4. Cross-section:
   `[ln W(i,T) − ln W(i,2014)] / (T − 2014) = α + β1 ln W(i,2014) + β2 ln Y(i,2014) + β3 WTP(i) + ε(i)`,
   where T is the last year with data for country i.
5. Nested F-tests:
   `F = [(RSS_r − RSS_u) / q] / [RSS_u / (n − k)]`, for income, WTP given
   income, income and WTP together, and QB8 given income.
6. Sigma convergence: coefficient of variation
   `CV(t) = SD(t) / mean(t)` and `σ(t)` = standard deviation of ln W, using
   the population formula (divisor N). Each is regressed on a linear time
   trend. The mean is the simple cross-country mean. Reported for all
   available countries each year and for the balanced panel.
7. Direction: dispersion (σ and CV) and the mean of W are compared between
   2014 and 2022 (all 27 countries have both years). The four cells are
   desirable convergence (dispersion falls, mean falls), convergence through
   deterioration (dispersion falls, mean rises), improvement with uneven
   progress (dispersion rises, mean falls) and divergence with deterioration
   (both rise). Countries are classified the same way: the gap
   `|ln W(i) − ln mean|` narrows or widens, and the country's own waste falls
   or rises. The Spearman correlation between the 2014 level and the
   percentage change is also reported.
8. Phillips-Sul log-t test on the balanced panel:
   `h(it) = W(it) / mean(W(t))`,
   `H(t) = mean over i of (h(it) − 1)²`, and
   `log(H(1)/H(t)) − 2 log L(t) = a + b log t + u(t)` with
   `L(t) = log(t + 1)`. HAC standard errors with one lag. Convergence is
   rejected if the t-statistic of b is below −1.65. The test runs on all ten
   years and again after dropping the first two years.
9. Convergence clubs: a simplified Phillips-Sul club search. Countries are
   ordered by their mean waste in the last three years. The core group is the
   first k countries with the highest t-statistic above −1.65. Each remaining
   country joins if the t-statistic stays above −1.65, and the search repeats
   on the rest. Club membership is indicative only because T is 10.
10. Timing: for each candidate break year from 2016 to 2022, the Chow test
    `F = [(RSS_p − RSS_1 − RSS_2) / k] / [(RSS_1 + RSS_2) / (n − 2k)]`
    with k = 2 on the model in equation 1, and the interaction
    `Δln W(it) = α + β ln W(i,t−1) + θ [ln W(i,t−1) x Post(t)] + τ(t) + ε(it)`
    with year effects. The Post main effect is left out because the year
    effects absorb it. The Bonferroni threshold for 7 candidate years is
    0.05 / 7. The pre-specified dates are 2018 and 2020.
11. Diagnostics: Breusch-Pagan, Jarque-Bera, variance inflation factors, the
    condition number of standardised regressors, Cook's distance, and a
    leave-one-country-out re-estimate of equation 1. The cross-section model
    gets the same checks.
12. Robustness: two-way fixed effects (country and year) and country fixed
    effects, which are biased in a short panel with a lagged outcome and are
    read as an upper bound on speed. A re-run of equations 1 and 2 without any
    transition that touches a flagged value. Drops of Estonia and Luxembourg
    from the cross-section.

## Run the analysis

Install the required Python packages:

```shell
pip install pandas numpy scipy statsmodels openpyxl
```

From the project directory, run:

```shell
python ecosem_analysis.py --out results
```

The script uses `ecosem_data.xlsx` by default, prints each result table, and
writes a CSV with the same table name into the output directory. `--data` and
`--out` can be set to other paths.

The script returns exit code `0` on success, `1` for input or file errors, and
`2` if one or more analysis stages fail. A failed stage may mean that some
other output files were still produced; check the error message printed for
the failed stage.

Tested with Python 3.13, pandas 3.0.5, numpy 2.5.3, scipy 1.18.1,
statsmodels 0.15.0, and openpyxl 3.1.5.

## Output files

All files below are CSVs in the output directory. Blank values can occur where
a statistic does not apply (for example, speed cannot be calculated for an
unstable coefficient, or a singleton convergence club has no club test).

| File | Meaning and columns |
|---|---|
| `beta_convergence.csv` | Panel convergence model results. One row per specification: absolute convergence, income with and without year effects, and two-way or country fixed effects. `beta` is the coefficient on lagged log waste (`level`); `se`, `p`, `ci_low`, and `ci_high` describe its uncertainty. `lambda` and `half_life_years` are implied convergence speed and half-life. `r2` and `n` are fit and sample size. `gamma`, `gamma_se`, and `gamma_p` are the income coefficient and its uncertainty where income is included. Panel standard errors are clustered by country. |
| `cross_section_models.csv` | Country-level regressions of average annual log waste growth from 2014 to each country's last available year. `model` identifies the specification (`A_level`, `B_income`, `C_income_wtp2014`, or `D_income_qb8`); `term` identifies each coefficient (including the intercept). `coef`, `se`, `p`, `r2`, `adj_r2`, and `n` report coefficient, classical standard error, p-value, model fit, and sample size. |
| `cross_section_f_tests.csv` | Nested F-tests for adding income, 2014 willingness-to-pay (WTP), both income and WTP, or 2024 QB8 WTP to the relevant country-level model. `F`, `p`, `df1`, and `df2` are the test statistic, p-value, and numerator and denominator degrees of freedom. |
| `cross_section_extras.csv` | Additional cross-section checks in `measure`/`value` rows: the 2014 WTP coefficient's confidence interval and HC3 p-value; Pearson and Spearman correlations between WTP measures; correlations of 2014 WTP with income and initial waste; and selected coefficient/p-value checks after dropping Estonia or Luxembourg. |
| `dispersion_all_available.csv` | Yearly waste distribution using all countries with data in that year. `countries` is the number observed; `mean_kg` and `sd_kg` are the mean and population standard deviation of waste per capita; `cv` is the coefficient of variation; `sd_ln` is the population standard deviation of log waste. |
| `dispersion_balanced.csv` | The same yearly distribution measures (`mean_kg`, `cv`, and `sd_ln`) for the balanced sample with data in every year. `countries` gives its fixed sample size. |
| `dispersion_trends.csv` | Linear time-trend regressions for `cv`, `sd_ln`, and mean waste (`mean_kg`) for both all-available and balanced samples. `slope_per_year`, `se`, `p`, and `r2` describe each trend; `countries_last_year` reports the number of countries in the final year of that sample. |
| `direction_summary.csv` | Summary comparing 2014 with 2022 for countries observed in both years: number of countries, mean waste at each endpoint and its percent change, and the dispersion (`sd_ln` and `cv`) at each endpoint. |
| `direction_countries.csv` | Country-by-country 2014–2022 comparison. `start_kg`, `end_kg`, and `change_pct` give waste levels and percentage change; `gap_narrowed` indicates whether the country's absolute log-waste gap from the cross-country mean narrowed; `own_waste_fell` indicates whether its own waste decreased. |
| `direction_cells.csv` | Counts countries across the two direction classifications. The row index `gap_narrowed` and the `True`/`False` columns indicate whether own waste fell (`True`) or did not fall (`False`). |
| `direction_spearman.csv` | Spearman rank correlation (`spearman_initial_level_vs_change`) between 2014 waste level and percentage change through 2022, plus its p-value. |
| `eu27_aggregate.csv` | Eurostat EU27 aggregate packaging waste per capita by `year`, included as a directional cross-check. It is not the country-level simple mean. |
| `phillips_sul_tests.csv` | Phillips-Sul log-t convergence tests for the balanced panel over all ten years and with the first two years dropped. `b`, `t`, and `one_sided_p` give the slope, HAC t-statistic, and one-sided p-value; `reject_convergence` applies the project's −1.65 threshold. `h_first` and `h_last` are relative-transition dispersion at the endpoints. |
| `phillips_sul_clubs.csv` | Indicative convergence-club assignments for each Phillips-Sul window. `club` and `size` identify each group; `members` lists countries; `club_b` and `club_t` give the group's log-t slope and t-statistic when the club has enough members for the test. |
| `timing_tests.csv` | Candidate structural breaks for change years 2016–2022. `chow_F` and `chow_p` are the Chow test; `theta` and `theta_p` describe the post-break interaction in the year-effects model; `beta_before` and `beta_after` are separate unadjusted slope estimates. The remaining columns report observations and mean growth on each side and the Bonferroni-adjusted significance threshold. |
| `diagnostics.csv` | Regression diagnostics in `diagnostic`/`value` rows for panel and country-level models: Breusch-Pagan and Jarque-Bera tests, skewness/kurtosis, VIFs, standardized condition numbers, Cook's distance summaries, and leave-one-country-out panel slope bounds. Names prefixed with `endpoint_` refer to the country-level model. |
| `flagged_check.csv` | Sensitivity check excluding panel transitions where either endpoint year has a Eurostat observation flag. `dropped_transitions` and `n` give the excluded and retained sample sizes; `beta_eq1`/`p_eq1` and `beta_eq2`/`p_eq2` give the two convergence slopes and p-values; `gamma_p` is the conditional model's income p-value. |

## Limitations

- The design has no untreated comparison group, so the results describe
  convergence during the rollout and do not identify a causal effect of the
  Directive.
- There are only ten annual points for sigma and Phillips-Sul, and 27 clusters
  for the panel standard errors. Treat those tests as indicative.
- The yearly sigma series on all available countries changes composition in
  2020, 2021 and 2023. The balanced series is the cleaner comparison.
- QA10 and QB8 differ in scale and in the goods covered, and neither concerns
  packaging.
- GDP is the Eustat republication of Eurostat accounts.
- The recycling-rate robustness check described in the paper is not run,
  because the recycling rate is not in the input file.
