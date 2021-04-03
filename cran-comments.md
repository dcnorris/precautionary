## Test environments
* local OS X install, R 4.0.4
* rhub::check(platform="linux-x86_64-rocker-gcc-san")

## R CMD check --as-cran results
There was one NOTE regarding sub-directory sizes:

> N  checking installed package size
>      installed size is  5.9Mb
>      sub-directories of 1Mb or more:
>        doc    2.9Mb
>        libs   2.2Mb

This reflects that the package now uses some
compiled code, and its vignettes are substantial.
