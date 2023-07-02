# plac 0.1.3

* Reduce sample sizes in data generation examples

# plac version 0.1.2

* Update package documentation.
* Update NAMESPACE 
  - Add `useDynLib(plac, .registration = TRUE)` to avoid notes from `R CMD check`
* Update DESCRIPTION
  - Remove Date field
  - Remove LazyData field
  - Add UTF-8 encoding
  - Increase roxygen2 version number
  - Change license field with `usethis::use_gpl_license()`
* Update README
  - Use `README.Rmd`
  - Update badges
  - Correct bibliography

# plac version 0.1.1

* Add
```R 
importFrom("stats", "as.formula", "coef", "model.frame",
          "model.matrix", "model.response", "plnorm", "pnorm",
          "qlnorm", "quantile", "rbinom", "rexp", "runif", "rweibull",
          "stepfun")
```
to NAMESPACE.

# plac version 0.1.0

* Add `sim.ltrc()` used to simulate left-truncated and right censored data
* Add `plr()` to test stationarity (uniform) assumption of the truncation times
* Add all the model fitting functions `PLAC_TXX()` which can handle multiple time-invariant and time-dependent covariates
* Add the main wrapper function `PLAC()` as the interface to the source

