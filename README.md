# intertest

An R package for **interaction testing of causal effects**. Given a binary
treatment X, outcome Y, pre-treatment confounders Z, and mediators W, `intertest`
decomposes the total variation in Y into a direct effect (DE), indirect effect
(IE), and spurious effect (SE), and tests whether these pathway-specific effects
interact with one another.

The five interaction measures are:

| Measure | Question |
|---|---|
| TE × SE | Does the total effect interact with the spurious effect? |
| DE × IE | Does the direct effect interact with the indirect effect? |
| DE × SE | Does the direct effect interact with the spurious effect? |
| IE × SE | Does the indirect effect interact with the spurious effect? |
| DE × IE × SE | Higher-order three-way interaction |

Estimation uses a semiparametric one-step debiased estimator with 10-fold
cross-fitting via gradient boosting (XGBoost). Three scales are supported:
difference, log-risk, and log-odds.

## Installation

```r
# install.packages("devtools")
devtools::install_github("dplecko/intertest")
```

## Quick start

```r
library(intertest)

set.seed(1)
n <- 500
Z <- rnorm(n)
X <- rbinom(n, 1, plogis(Z))
W <- X + Z + rnorm(n)
Y <- X * W + rnorm(n)          # X*W interaction -> DE x IE signal
dat <- data.frame(Z = Z, X = X, W = W, Y = Y)

fit <- inter_test(dat, X = "X", Z = "Z", W = "W", Y = "Y")
summary(fit)
ggplot2::autoplot(fit)
```

`inter_test()` returns an `"intertest"` object with a `results` data frame
(columns `measure`, `value`, `sd`, `p_value`, `scale`) plus the five
pathway-specific estimands (TV, Ctf-DE, Ctf-IE, Ctf-SE, ETT) accessible via
`summary()`.

`one_step_debias()` is the lower-level workhorse and can be called directly when
you need finer control (e.g. to supply a pre-computed cross-fit object or to
request pseudo-outcomes via `save_pso = TRUE`).

## Reproducing paper results

All the code used for the paper "Interaction Testing in Variation Analysis" 
can be found in `scripts/`, and should be run from the project root.

### Data

Download the necessary datasets from the
[data link](https://www.dropbox.com/scl/fo/5d4a2tqjy5rn73xesalir/AKHTU_oieip9-XBaFRP9954?rlkey=s0tt1oymkhil2pql36wvaf5td&st=vwlg0gon&dl=0)
and place them in the `data/` folder.

### Dependencies

Install the package and its dependencies:

```r
devtools::install_github("dplecko/ia-testing")
install.packages(c("xgboost", "ggplot2", "data.table", "rprojroot",
                   "faircause", "readr", "scales"))
```

### Scripts

Below is the overview of the files that perform different analyses.

| Script | Purpose |
|---|---|
| `ia-synthetic.R` | Synthetic SCM experiments (p-value distributions, type I/II error rates) |
| `ia-mimic.R` | MIMIC-IV analysis: interaction table + conditional direct effect by admission type |
| `ia-in-the-wild.R` | Real-data interaction testing across ten fairness datasets |

Results are written to `results/`.
