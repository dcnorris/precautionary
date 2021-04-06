## Test environments
* local OS X install, R 4.0.4
* rhub::check(platform="linux-x86_64-rocker-gcc-san")

## Re: Windows and Rust

Pursuant to our recent exchange after failed
incoming check on Win-builder, I have adapted
the INSTALL file from package 'baseflow'.
This explains in detail how to install Rust's
compilation manager 'cargo' on Windows.

By "LinkingTo: cargo (>= 0.1.28)" in the
DESCRIPTION file, I believe I should be
pointing Win-builder toward a CRAN-friendly
verson of the Rust 'run' command, which
avoids writing to ~/.cargo in accord with
CRAN prohibitions against writing to user
filesystem.


## R CMD check --as-cran results
There was one NOTE regarding sub-directory sizes:

> N  checking installed package size
>      installed size is  5.9Mb
>      sub-directories of 1Mb or more:
>        doc    2.9Mb
>        libs   2.2Mb

This reflects that the package now uses some
compiled code, and its vignettes are substantial.
