## R package 'plac' version 0.1.1

[![Travis-CI Build Status](https://travis-ci.org/942kid/plac.svg?branch=master)](https://travis-ci.org/942kid/plac)
[![Coverage Status](https://img.shields.io/codecov/c/github/942kid/plac/master.svg)](https://codecov.io/github/942kid/plac?branch=master)
[![CRAN version](http://www.r-pkg.org/badges/version/plac)](https://cran.r-project.org/package=plac)
![](http://cranlogs.r-pkg.org/badges/grand-total/plac)


### A Pairwise Likelihood Augmented Estimator for the Cox Model under Left-Truncation 

This is an R package to implement the semi-parametric estimation method for the Cox model introduced in the paper [<it>A Pairwise Likelihood Augmented Estimator for the Cox Model under Left-Truncation</it> by Wu et al. (2015)](http://biostats.bepress.com/umichbiostat/paper118/).
It gives more efficient estimates for left-truncated survival data under the Cox model by making use of the marginal survival information upto the entry of the subjects.
The method will be most helpful when the sample size or number of observed events are not large and estimation efficiency is more of a concern.

####  Installation 

This package can be installed from github with the help of `devtools`:

```R
# install.packages("devtools")
devtools::install_github("942kid/plac")
```

#### Examples

The main wrapper function `PLAC()` calls the appropriate working function according to the covariate types in the dataset. For example,

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

For more details on the arguments of the function, please run

```R
help(PLAC)
```
