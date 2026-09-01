# Optimization Methods for Bioprocess Engineering — R Course

Synthetic experimental data for a hands-on **R** course: a bioprocessing student
learns to code by finding the operating conditions that maximize bacterial
cellulose (BC) yield.

> R equivalent of `python-optimization-methods/` (the original Python course).

## Files

```
data/
  bc_full_factorial_dataset.csv   # 9375 rows — the student-facing dataset
  hidden_truth.json               # TRUE model + true optimum (teacher only!)
notebooks/
  Lesson01_data_exploration.Rmd   # Lesson 1 (student version, exercises blank)
  Lesson01_solutions.Rmd          # Lesson 1 solutions (teacher only)
  ...
  Lesson11_capstone.Rmd           # Lesson 11 capstone (student)
  Lesson11_solutions.Rmd          # Lesson 11 solutions (teacher only)
scripts/
  generate_dataset.R             # reproducible data generator (seed 42)
  build_lesson1.R                # rebuilds the Lesson 1 notebooks
  ...
  build_lesson11.R               # rebuilds the Lesson 11 notebooks
```

## Package mapping (Python → R)

| Python package | R equivalent |
|---|---|
| `pandas` | `dplyr` + `tidyr` + `readr` |
| `numpy` | Base R vectors/matrices + `matrixStats` |
| `matplotlib` / `seaborn` | `ggplot2` + `scales` + `viridis` |
| `scipy.stats` | `stats` (base) |
| `statsmodels` | `stats::lm()` + `car` |
| `scikit-learn` | `caret` + `e1071` + `randomForest` + `glmnet` |
| `scipy.optimize` | `stats::optim()` + `BB::optim()` |
| `simulated annealing` | `stats::optim(method="SANN")` |

## R packages required

```r
# Core
library(tidyverse)   # dplyr, ggplot2, readr, tidyr, purrr, stringr
library(scales)
library(viridis)

# Statistics
library(car)         # Type III ANOVA, vif()

# Machine Learning
library(caret)       # Unified ML interface
library(e1071)       # SVR
library(randomForest)
library(glmnet)      # Regularized regression
library(rpart)       # Decision trees

# Notebooks
library(rmarkdown)
```

## Dataset

Full Factorial Design: **5 factors × 5 levels each = 5⁵ = 3125 conditions**,
each run in **triplicate** (9375 rows total).

| Factor | Levels |
|---|---|
| pH | 4.0 / 4.5 / 5.0 / 5.5 / 6.0 |
| Temperature (°C) | 25 / 27.5 / 30 / 32.5 / 35 |
| RPM | 50 / 87.5 / 125 / 162.5 / 200 |
| Inoculum volume (%) | 5 / 7.5 / 10 / 12.5 / 15 |
| Sugar concentration (g/L) | 20 / 30 / 40 / 50 / 60 |
| **Response: BC yield (g/L)** | observed, with measurement noise |

Columns: `run_number`, `condition_id`, `replicate`, `pH`, `Temp_C`, `RPM`,
`Inoculum_pct`, `Sugar_gL`, `BC_yield_gL`.

## Ground truth (teacher only)

```
True optimum:  pH 5.53,  27.4 °C, 113 RPM,  12.1% inoculum,  49.3 g/L sugar
               -> 8.08 g/L BC  (best observed row: ~8.2 g/L)
```

## Course roadmap

1. **[DONE] R basics + tidyverse + ggplot2** — variables, vectors, functions,
   data frames; load CSV, filter, `group_by`, `summarise`; scatter/histogram/bar/
   facets/themes. Includes 6 coding exercises.
2. **[DONE] Noise-aware statistics** — variance/std/SE, simulated proof of why
   triplicates beat single runs, t-based 95% CIs, honest comparison of conditions
   (Welch t-test, p-value), std vs SE error bars, reproducible workflow (seed +
   `write_csv`).
3. **[DONE] Brute-force search** — formalise the maximisation problem, the
   cost explosion of grids (`L^F`, curse of dimensionality), a lookup-table
   `evaluate()` with an evaluation counter, `which.max` over all 3125 conditions,
   and a heatmap view of the response surface.
4. **[DONE] Random search & grid refinement** — sample N random points, refine
   around best; introduce the "noise-aware" idea (compare means, not single runs).
5. **[DONE] Local search: hill climbing & coordinate ascent** — from a random
   start, step factor-by-factor; discuss local optima and why smooth RSM helps.
6. **[DONE] Global search: simulated annealing** — accepts worse moves with
   decreasing probability; use `stats::optim(method="SANN")`; visualize
   acceptance over iterations.
7. **[DONE] RSM: multiple linear regression** — code a quadratic model with
   `lm()`, check R², ANOVA significance of terms, canonical analysis, predict
   optimum.
8. **[DONE] RSM vs truth** — compare each method's prediction to
   `hidden_truth.json`; report bias and distance to the true optimum.
9. **[DONE] Support vector regression (SVR)** — a second non-parametric
   surrogate. Why scaling matters, the epsilon-tube and the penalty `C`,
   linear vs RBF kernels, CV against RSM, and maximising the surrogate over
   the continuous box.
10. **[DONE] Machine learning surrogates** — Random Forest / Gradient Boosting
    on (factors → yield), feature importance, use the surrogate to maximize
    yield; compare ML vs RSM fit quality.
11. **[DONE] Capstone** — a single notebook comparing all approaches
    (brute force, random, hill climb, annealing, RSM, SVR, ML) on accuracy
    vs cost, a truth reveal, and a written conclusion in "lab report" style.

## Regenerating the data

```bash
/mnt/BackUp/Applications/.miniconda/envs/r_env/bin/Rscript scripts/generate_dataset.R
```

Same seed = same dataset.
