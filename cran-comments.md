## Resubmission
This is a resubmission. In this version I have:

* Added 
```R 
importFrom("stats", "as.formula", "coef", "model.frame",
          "model.matrix", "model.response", "plnorm", "pnorm",
          "qlnorm", "quantile", "rbinom", "rexp", "runif", "rweibull",
          "stepfun")
```
to the NAMESPACE file.

## Test environments
* local ubuntu 14.04, R 3.2.5
* ubuntu 12.04 (on travis-ci), R 3.2.4
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE:

* checking installed package size ... NOTE
  installed size is 30.0Mb
  sub-directories of 1Mb or more:
    libs  29.9Mb 

  This because of the compiled source code.
  
## Downstream dependencies

None.
