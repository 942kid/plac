## 'plac' version 0.1.0

* Adding `sim.ltrc()` used to simulate left-truncated and right censored data
* Adding `plr()` to test stationarity (uniform) assumption of the trunation times
* Adding all the model fitting functions `PLAC_TXX()` which can handle multiple time-invariant and time-dependent covariates
* Adding the main wrapper function `PLAC()` as the interface to the source

## 'plac' version 0.1.1

* Adding 
```R 
importFrom("stats", "as.formula", "coef", "model.frame",
          "model.matrix", "model.response", "plnorm", "pnorm",
          "qlnorm", "quantile", "rbinom", "rexp", "runif", "rweibull",
          "stepfun")
```
to the NAMESPACE file.
