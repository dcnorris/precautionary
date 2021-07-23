#' Plot duty cycle from performance report on a parallel computation
#'
#' This is intended for application to the tables which get attached
#' to `Cpe` instances after invocation the `$trace_paths()` method.
#' But any table with columns named `pid`, `t1` and `t2` suffices.
#' @param perftab A performance table with columns:
#' * `pid` An integer, character or factor with process ids
#' * `t1`, `t2` Task start and end times, in milliseconds
#' @param ... Optional parameters passed along to `lattice::xyplot`
#' @return An `xyplot` with a duty-cycle panel for each worker,
#' plus an overall average
#' @examples
#' if (interactive()) {
#'   ## Example from Braun2020
#'   d1_maxn <- 5
#'   cum_maxn <- 10
#'   mod <- Crm$new(skeleton = c(0.03, 0.11, 0.25, 0.42, 0.58, 0.71),
#'                  scale = 0.85, # aka 'sigma'
#'                  target = 0.25)$
#'     no_skip_esc(TRUE)$    # compare Braun's 'restrict = T'
#'     no_skip_deesc(FALSE)$
#'       stop_func(function(x) {
#'         enrolled <- tabulate(x$level, nbins = length(x$prior))
#'         x$stop <- enrolled[1] >= d1_maxn || max(enrolled) >= cum_maxn
#'         x
#'       })
#'   mod$trace_paths(1, rep(2, 13), unroll = 4
#'                 , mc.cores = parallelly::availableCores(omit=2))
#'   print(mod$performance)
#'   plot_dutycycle(mod$performance)
#' }
#' @export
plot_dutycycle <- function(perftab, ...) {
  jpt <- perftab
  jpt$job <- seq(nrow(jpt)) # TODO: Better if perftab supplies names

  jpt$pid <- as.factor(jpt$pid) # TODO: Get 'pid' colname from a formula

  T <- 0:ceiling(max(jpt$t2))
  m <- matrix(FALSE, nrow=length(T), ncol=nlevels(jpt$pid))
  dimnames(m) <- list(t = T, pid = levels(jpt$pid))

  for (j in 1:nrow(jpt)) {
    pid <- jpt$pid[j]
    m[,pid] <- m[,pid] | (jpt$t1[j] < T & T < jpt$t2[j])
  }

  ## Add an 'avg' column
  m <- addmargins(m, margin = 2, FUN = c(avg = mean))

  ## Reshape this matrix to 'tall' form yielding series to plot
  dt <- data.table(t = T
                 , active = 1.0 * as.vector(m)
                 , pid = factor(rep(colnames(m), each=nrow(m)))
                   )

  lattice::xyplot(active ~ t | pid, data = dt, type='l'
       , layout=c(1, NA), as.table=TRUE
       , ylab = ""
       , xlab = "Sys.time() since parallelization [ms]"
       , scales = list(y = list(at = NULL))
       , ...
         )
}
