# R package 'plac' version 0.1.1

[![Travis-CI Build Status](https://travis-ci.org/942kid/plac.svg?branch=master)](https://travis-ci.org/942kid/plac)
[![Coverage Status](https://img.shields.io/codecov/c/github/942kid/plac/master.svg)](https://codecov.io/github/942kid/plac?branch=master)
[![CRAN version](http://www.r-pkg.org/badges/version/plac)](https://cran.r-project.org/package=plac)
![](http://cranlogs.r-pkg.org/badges/grand-total/plac)


## A Pairwise Likelihood Augmented Estimator for the Cox Model under Left-Truncation 

This R package implements a semi-parametric estimation method for the Cox model
introduced in the paper *[A Pairwise Likelihood Augmented Cox Estimator for
Left-truncated data by Wu et al.  (2018)](https://doi.org/10.1111/biom.12746)*.
It gives more efficient estimate for left-truncated survival data using the
marginal survival information upto the start of follow-up (when the subject
enters the risk set).  The independence between the underlying truncation time
distribution and the covariates is the only additional assumption, which holds
true for most applications of length-biased sampling problem and beyond.

##  Installation 

This package can be installed from CRAN:

```R
install.packages("plac")
```

or from github:

```R
# install.packages("devtools")
devtools::install_github("942kid/plac")
```

## Examples

The main wrapper function `PLAC()` calls the appropriate working function
according to the covariate types in the dataset. For example,

```R
library(plac)

# When only time-invariant covariates are involved
dat1 = sim.ltrc(n = 50)$dat
PLAC(ltrc.formula = Surv(As, Ys, Ds) ~ Z1 + Z2,
     ltrc.data = dat1, td.type = "none")

# When there is a time-dependent covariate that is independent of the truncation time
dat2 = sim.ltrc(n = 50, time.dep = TRUE,
               distr.A = "binomial", p.A = 0.8, Cmax = 5)$dat
PLAC(ltrc.formula = Surv(As, Ys, Ds) ~ Z,
     ltrc.data = dat2, td.type = "independent",
     td.var = "Zv", t.jump = "zeta")

# When there is a time-dependent covariate that depends on the truncation time
dat3 = sim.ltrc(n = 50, time.dep = TRUE, Zv.depA = TRUE, Cmax = 5)$dat
PLAC(ltrc.formula = Surv(As, Ys, Ds) ~ Z,
     ltrc.data = dat3, td.type = "post-trunc",
     td.var = "Zv", t.jump = "zeta")
```

For details, please refer to the document

```R
help(PLAC)
```
