# Replication package

This repository contains the replication codes for Cosma, Kostyrka, Tripathi (2026) Missing endogenous variables in conditional-moment-restriction models. Journal of Business and Economic Statistics.

# Contents

The files should be processed in the following order.

1a. Download the data files (`pums80-twoa.rds`, `pums90-onea.rds`);

1b. Alternatively, run the pre-processor script to get the aforementioned two as RData from SAS (`rep-AE98-preprocess-rev2.R`);

2a. Run the main script (`rep-AE98-smooth-rev2.R`).

2b. The tables should be printed in the console; the images should be generated in the working directory.

The simulation is resource-intensive and is best run in a super-computer cluster.

# What is under the hood of `rep-AE98-preprocess-rev2.R`

This file requires `m_d_806.sas7bdat` from the [Angrist Data Archive](https://economics.mit.edu/people/faculty/josh-angrist/angrist-data-archive) (specifically, the [AngEv98.zip](https://economics.mit.edu/sites/default/files/publications/AngEv98.zip) file with all the data).
Put `m_d_806.sas7bdat` into the same folder as `rep-AE98-preprocess-rev2.R` and run the latter in R (no user input needed) to get the outputs `pums80-twoa.rds` and `pums90-onea.rds`.
The already-processed files `pums80-twoa.rds` and `pums90-onea.rds` are provided in this repository for research purposes.

# What is under the hood of `rep-AE98-smooth-rev2.R`

In this numerical exercise, we evaluate the impact of missingness on the accuracy of 3 estimators: Generalised Method of Moments (GMM) and two Smoothed Empirical Likelihood (SEL) variants: validation-sample-only (SELg) and full-sample efficient (SELrho); the latter is our main contribution to the literature on missing data.

The script replicates and extends the Angrist--Evans (AER, 1998) same-sex IV design to estimate the causal effect of multi-child parenting (`MOREKIDS`) on mother’s labour-market outcomes (currently restricted to mother’s income, `INCOMEM`).

We evaluate the three estimators in a **Monte-Carlo study** with simulated missingness/selection.

## Data construction and sample selection

Goal: prepare AE98-style dataset.

Load 1980 (or 1990) Census PUMS extracts. Restrict to married women, whites only, 1980 baseline specification with full controls. Create income in $1,000 units.

Excluded instruments: `BOYS2`, `GIRLS2`

Included instruments: age, age at first birth, first child sex.

## Baseline IV / GMM estimation

Goal: provide benchmark IV estimates and over-identification diagnostics (columns 2-4 of Table 1).

Using `momentfit`, estimates the baseline AE98 IV model. Estimate richer instrument expansions (interactions, squares) as a robustness check. Run a top-coded age-robustness version.

Write summaries of coefficients, SEs, and GMM diagnostics to the console.


## Smoothed Empirical Likelihood (SEL)

Goal: provide columns 5-7 of Table 1.

The `optMixSELfull()` function:
* Treats discrete and continuous variable differently for conditioning;
* Uses kernel smoothing for weighted empirical-likelihood computation;
* Uses blocking, de-duplication, and sparse kernel weights for speed;
* Performs BFGS optimisation and Hessian-based variance estimation.

This likelihood-style semiparametric estimator is designed to be most efficient.

We save the SEL coefficient vector, the variance matrix, and timing / optimization diagnostics.

## Monte Carlo simulation with missingness

For each pair across:
* Many random seeds (`doOneSeed <- function(s)`),
* Missingness strength grid (`k in seq_along(svec)`),

-- do the following:

1. Generate missingness
2. Re-estimates with multiple estimators (`GMM`, `SELg`, `SELrho`, each with 4 different bandwidths)

Save per-seed results. It is computationally faster to gradually increase missingness and use the previous optimum as the starting value. The smoothing bandwidth does not have any substantial effect on the efficient estimator (`SELrho`) owing to the double-robustness property.

We save the output into many per-seed `.RData` result files.

## Aggregation and figure generation

After loading simulation results, we produce:
* Missingness-rate curves;
* Estimator success-rate diagnostics;
* Distribution summaries of estimates and CI behaviour;
* Publication-ready plots (PDF + TikZ/LaTeX).

This section also generates Table 2 and Figures 1 and 2.
