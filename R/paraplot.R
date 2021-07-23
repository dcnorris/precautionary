## A pair of plots for parallel performance

plot.dutycycle <- function(perftab, ...) {
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
