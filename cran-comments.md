## Test environments
* local OS X install, R 4.1.1
* R-CMD-check GitHub matrix action
  - macOS-latest (release)
  - ubuntu-18.04 (release & devel)
  - window-latest (release)
  - See https://github.com/dcnorris/precautionary/actions/runs/1170449271

## R CMD check --as-cran results

There was 1 NOTE about sub-directory sizes:

> checking installed package size ... NOTE
>   installed size is  6.7Mb
>   sub-directories of 1Mb or more:
>     doc    3.4Mb
>     libs   2.2Mb

This reflects that the package includes compiled code, and that
its vignettes are substantial.

## Actions pursuant to comments on v0.2.6, addressed below:

> "Package authors should make all reasonable efforts to provide
> cross-platform portable code."
>
> has not been complied with.
>
> You have added a requirement for cargo but not declared it nor tested
> for it.  And your Makevars is using GNU make extensions without checking
> for them.

* Removed `OS_type: unix` from DESCRIPTION
* Added `SystemRequirements: Cargo` to DESCRIPTION
* Added a `configure` script that tests for cargo and rustc
* Rectified non-portable usage `export CARGO_HOME=...` in Makevars,
  expressing the export on a separate line from the var definition.

> Using rust is not portable, as 'Writing R Extensions' makes clear (as it
> does about make extensions).  It is not something that should be added
> to an existing package, and most definitely not in a patch-level update.

* I have used Rust here to speed up computations of certain integrands,
  an enhancement that is incidental to the conception of this package
  (and so most definitely does not indicate starting a fresh package),
  but is essential for feasible realization of Complete Path Enumeration
  (CPE) for the important Continual Reassessment Method (CRM) class of
  dose-escalation designs.

* Rust's widely celebrated support for safe concurrency in fact creates
  new technical opportunities for me to extend multicore CPE to Windows
  users in the future, bypassing disadvantages R suffers on that platform.
  Thus, my introduction of Rust will serve the interests of a pragmatic,
  relevant cross-platform support for valuable package features, on major
  platforms represented in the academic & pharmaceutical contexts where
  `precautionary` will be of greatest interest.

* The versioning scheme I have adopted is a conservative one, in step with
  typical academic culture in this field, such that several related CRAN
  packages have 0.1.x versions. This serves to emphasize the experimental
  and provisional status of the underlying ideas, as well as the lengthy
  anticipated course of their future development. I am using a `-<patch>`
  suffix to indicate patches, such as the present 'v0.2.6-1' addressing
  the above comments.
