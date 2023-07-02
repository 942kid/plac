## Resubmission

This is a resubmission for package **plac**.

Changes this version (v0.1.2):

* Update package documentation.
* Update NAMESPACE with `useDynLib(plac, .registration = TRUE)` to avoid notes from `R CMD check`.
* Update DESCRIPTION
  - Remove Date field
  - Add UTF-8 encoding
  - Increase roxygen2 version number
  - Change license field with `usethis::use_gpl_license()`
* Update README
  - Use `README.Rmd`
  - Update badges
  - Correct bibliography

## R CMD check Results

We used GitHub Actions to run `R CMD check` on different operating systems and R version.

There were no ERROR or WARNING.

There was 1 NOTE on Linux:

* checking installed package size ... NOTE
  installed size is 31.4Mb
  sub-directories of 1Mb or more:
  libs  31.2Mb

This because of the compiled source code.

## Reverse Dependencies

None.
