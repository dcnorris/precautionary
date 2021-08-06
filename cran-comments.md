## Test environments
* local OS X install, R 4.1.0
* R-CMD-check GitHub matrix action
  - macOS-latest (release)
  - ubuntu-18.04 (release & devel)
  - See https://github.com/dcnorris/precautionary/actions/runs/1105715875

## R CMD check --as-cran results

There was 1 NOTE about sub-directory sizes:

> checking installed package size ... NOTE
>   installed size is  6.1Mb
>   sub-directories of 1Mb or more:
>     doc    2.9Mb
>     libs   2.2Mb

This reflects that the package includes compiled code, and that
its vignettes are substantial.
