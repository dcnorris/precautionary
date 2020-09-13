## Test environments
* local OS X install, R 4.0.2
* win-builder (devel)

## R CMD check --as-cran results
There were no ERRORs or WARNINGs

There was 1 **spurious** NOTE relating to a
known "unable to verify current time" issue;
see: https://www.mail-archive.com/r-package-devel@r-project.org/msg05869.html.

## Additional comments
I have reduced by 25% the size of the simulation in the
example for 'format.safetytab' (num_sims = 60 vs prev. 80),
since it consumed 10.95s elapsed time on Win arch 'x64'.
