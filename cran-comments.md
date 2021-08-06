## Test environments
* local OS X install, R 4.0.5
* R-CMD-check GitHub matrix action
  - macOS-latest (release)
  - windows-latest (release)
  - ubuntu-18.04 (release & devel)
  - See https://github.com/dcnorris/precautionary/actions/runs/1001668187

## Rust library

Thank you for your close analysis of v0.2-2 (submitted 2021-04-15),
which I understand you found in conformance with CRAN policies.

My v0.2.3 submission (per Uwe Ligges's emailed request of 24 Apr)
never got acted on, I think. Uwe's plan at that time had been to
"disable the runtime tests so that we do not need to install [Rust
package manager] cargo on the Windows machine" while avoiding the
need for an 'OS_type: unix' constraint in my DESCRIPTION.

## R CMD check --as-cran results
There was 1 NOTE about sub-directory sizes:

> checking installed package size
>   installed size is  5.9Mb
>   sub-directories of 1Mb or more:
>     doc    2.9Mb
>     libs   2.2Mb

This reflects that the package includes compiled code, and that
its vignettes are substantial.
