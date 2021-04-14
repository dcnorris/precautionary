## Test environments
* local OS X install, R 4.0.4
* rhub::check(platform="debian-gcc-devel")

## Re: Windows and Rust

In order to progress step-wise with
CRAN release of new 'precautionary'
that includes its new Rust library,
I have set OS_type: unix in the
DESCRIPTION file.

This will allow me to deal with the
Windows build issues in a separate,
non-blocking thread of activity.

## R CMD check --as-cran results
There was 1 NOTE re sub-dir sizes:

> checking installed package size
>   installed size is  5.9Mb
>   sub-directories of 1Mb or more:
>     doc    2.9Mb
>     libs   2.2Mb

This reflects that the package now
uses some compiled code, and that
its vignettes are substantial.
