## This is a 'passthru' adapter which stands in for an enhanced `precautionary:::mclapply`
## that reports progress on parallelized computations. This is used on Windows, which lacks
## the fork() system call needed for `mclapply`-based parallelism. It is also used on CRAN
## releases, because the enhanced `precautionary::mclapply` uses unexported `parallel`
## functionality, and so would not pass CRAN checks.
## The effect of including this code should be that the extended signature may be invoked,
## with the 'nicety' of progress reporting simply omitted as a 'silent failure'.

mclapply <- function(X, FUN, ..., mc.preschedule = TRUE, mc.set.seed = TRUE,
                     mc.silent = FALSE, mc.cores = getOption("mc.cores", 2L),
                     mc.cleanup = TRUE, mc.allow.recursive = TRUE,
                     progreport = NULL,
                     progmetric = length,
                     proginit = 0L,
                     affinity.list = NULL)
    parallel::mclapply(X, FUN, ..., mc.preschedule = mc.preschedule,
                       mc.set.seed = mc.set.seed, mc.silent = mc.silent,
                       mc.cores = mc.cores, mc.cleanup = mc.cleanup,
                       mc.allow.recursive = mc.allow.recursive,
                       affinity.list = affinity.list)
