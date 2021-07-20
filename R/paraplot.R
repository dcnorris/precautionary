## A pair of plots for parallel performance

plot.dutycycle <- function(perftab) {
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

  ## Reshape this matrix to 'tall' form yielding series to plot
  active <- as.vector(m)
  pid <- rep(levels(jpt$pid), each=nrow(jpt))

  dt <- data.table(t = T
                 , active = 1.0 * as.vector(m)
                 , pid = factor(rep(levels(jpt$pid), each=nrow(m)))
                   )

  ## Here's one option for the plot
  xyplot(active ~ t, groups = pid, data = dt, type='l'
       , auto.key=list(points=FALSE, lines=TRUE) #, columns=6)
         )
}
