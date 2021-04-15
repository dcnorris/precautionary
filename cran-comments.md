## Test environments
* local OS X install, R 4.0.5
* rhub::check(platform="debian-gcc-devel")

## Avoiding write to ~/.cargo

Digging deep in CRAN package 'salso', I see
that it employs a custom R build script in
tools/. I have adapted this for use in package
'precautionary', and also changed the Makevars
LIBDIR variable from ./rust/target/release
to rust/target/release.

## Re: Windows and Rust

In order to progress step-wise with
CRAN release of new 'precautionary'
that includes its new Rust library,
I have set OS_type: unix in the
DESCRIPTION file.

I have a new vignette under development,
however, which will create the opportunity
for a v0.2-3 release in the near future,
allowing an attempt to undo the OS_type
constraint in keeping with CRAN spirit
of maximum cross-platform compatibility.

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
