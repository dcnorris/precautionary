## Test environments
* local OS X install, R 4.0.5
* rhub::check(platform="debian-gcc-devel")

## Rust library

Thank you for your close analysis of v0.2-2
(submitted 2021-04-15), which I understand
you find in conformance with CRAN policies.

Per Uwe Ligges's emailed request of 24 Apr,
I now submit v0.2-3 without any OS_type
restriction.

## R CMD check --as-cran results
There was 1 NOTE re sub-dir sizes:

> checking installed package size
>   installed size is  5.9Mb
>   sub-directories of 1Mb or more:
>     doc    2.9Mb
>     libs   2.2Mb

This reflects that the package includes
compiled code, and that its vignettes are
substantial.
