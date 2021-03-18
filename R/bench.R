## This seems to be converging on a reasonable generic pattern.
## The idea is to write an expression 'expr' that includes '...'
## which can be filled in using (most appropriately!) this very
## function's '...' argument!)
##' @importFrom microbenchmark microbenchmark
##' @importFrom stats median
benchtab <- function(expr, ...) {
  ## In general, this function should evaluate <expr> with various
  ## substitutions as provided in the ... args.
  ## Perhaps the ... will consist of a small number (1--3) of
  ## arguments, each of which is given as a short vector of
  ## alternative values. The resulting performance comparison
  ## ought to tabulate the performance against the full outer
  ## product of the possible values.
  fun <- eval(substitute(function(...) expr))
  perf <- CJ(..., sorted = FALSE)
  mb <- list()
  for (i in seq(nrow(perf))) {
    mb[[i]] <- microbenchmark(do.call(fun, perf[i, ]))
  }
  perf$milliseconds <- sapply(mb, function(.) median(.$time)/1e6)
  perf$evals.per.sec <- sapply(mb, function(.) mean(1e9/.$time))
  perf
}
