
<!-- README.md is generated from README.Rmd. Please edit that file -->

# plac

<!-- badges: start -->

[![R-CMD-check](https://github.com/942kid/plac/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/942kid/plac/actions/workflows/R-CMD-check.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/plac)](https://CRAN.R-project.org/package=plac)
<!-- badges: end -->

This R package implements a semi-parametric estimation method for the
Cox model introduced in the paper *[A Pairwise Likelihood Augmented Cox
Estimator for Left-truncated data by Wu et
al. (2018)](https://doi.org/10.1111/biom.12746)*. It gives more
efficient estimate for left-truncated survival data using the marginal
survival information up to the start of follow-up (when the subject
enters the risk set). The independence between the underlying truncation
time distribution and the covariates is the only additional assumption,
which holds true for most applications of length-biased sampling problem
and beyond.

## Installation

The package can be installed from CRAN:

``` r
install.packages("plac")
```

You can also install the development version of it from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("942kid/plac")
```

## Example

The main wrapper function `PLAC()` calls the appropriate working
function according to the covariate types in the dataset. For example,

``` r
library(plac)
#> Loading required package: survival

# When only time-invariant covariates are involved
dat1 <- sim.ltrc(n = 50)$dat
PLAC(
  ltrc.formula = Surv(As, Ys, Ds) ~ Z1 + Z2,
  ltrc.data = dat1,
  td.type = "none"
)
#> Calling PLAC_TI()...
#> 12 Iterations
#> Coefficient Estimates:
#>    est.Cox se.Cox p.Cox est.PLAC se.PLAC p.PLAC
#> Z1   2.055  0.431 0.000    1.804   0.357  0.000
#> Z2   0.919  0.347 0.008    0.804   0.259  0.002

# When there is a time-dependent covariate that is independent of the truncation time
dat2 <- sim.ltrc(n = 50, time.dep = TRUE, distr.A = "binomial", p.A = 0.8, Cmax = 5)$dat
PLAC(
  ltrc.formula = Surv(As, Ys, Ds) ~ Z,
  ltrc.data = dat2, td.type = "independent",
  td.var = "Zv", t.jump = "zeta"
)
#> Calling PLAC_TD()...
#> 100 Iterations
#> 
#> Coefficient Estimates:
#>    est.Cox se.Cox p.Cox est.PLAC se.PLAC p.PLAC
#> Z    0.866  0.330 0.009    0.795   0.224      0
#> Zv   0.877  0.355 0.014    0.864   0.214      0

# When there is a time-dependent covariate that depends on the truncation time
dat3 <- sim.ltrc(n = 50, time.dep = TRUE, Zv.depA = TRUE, Cmax = 5)$dat
PLAC(
  ltrc.formula = Surv(As, Ys, Ds) ~ Z,
  ltrc.data = dat3, td.type = "post-trunc",
  td.var = "Zv", t.jump = "zeta"
)
#> Calling PLAC_TDR()...
#> 8 Iterations
#> 
#> Coefficient Estimates:
#>    est.Cox se.Cox p.Cox est.PLAC se.PLAC p.PLAC
#> Z    0.668  0.301 0.027    0.487   0.246  0.047
#> Zv   0.915  0.327 0.005    0.938   0.301  0.002
```

For computation details, please refer to the document of the main
wrapper function:

``` r
help(PLAC)
```

## References

Wu, F., Kim, S., Qin, J., Saran, R., & Li, Y. (2018). A pairwise
likelihood augmented Cox estimator for left‐truncated data.
*Biometrics*, 74(1), 100-108.
