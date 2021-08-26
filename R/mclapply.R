## This is a 'passthru' adapter that stands in for a progress-reporting
## enhancement on the main branch, which however would not pass CRAN checks
## because it uses unexported `parallel` functionality.
## The effect of including this should be that the extended signature may
## be invoked, with the progress reporting omitted by 'silent failure'.

mclapply <- function(X, FUN, ..., mc.preschedule = TRUE, mc.set.seed = TRUE,
                     mc.silent = FALSE, mc.cores = getOption("mc.cores", 2L),
                     mc.cleanup = TRUE, mc.allow.recursive = TRUE,
                     progreport = NULL,
                     progmetric = length,
                     proginit = 0L,
                     affinity.list = NULL)
{
    if (.Platform$OS.type == "windows") # On Windows, which lacks fork(),
        return(lapply(X, FUN, ...))     # mclapply is just lapply.

    parallel::mclapply(X, FUN, ..., mc.preschedule = mc.preschedule,
                       mc.set.seed = mc.set.seed, mc.silent = mc.silent,
                       mc.cores = mc.cores, mc.cleanup = mc.cleanup,
                       mc.allow.recursive = mc.allow.recursive,
                       affinity.list = affinity.list)
}
