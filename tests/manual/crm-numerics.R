## Compare the objective functions visually
comp_objfn <- function(s = 500) {

  skeleton <- seq(0.1, 0.6, 0.1)

  ## cohort-wise set-up
  level <- c(1L, 2L, 3L, 2L, 3L, 3L)
  n_tox <- c(0L, 0L, 1L, 0L, 0L, 1L)
  obs <- encode_cohorts(enr = tabulate(level, nbins=8)
                       ,tox = xtabs(n_tox ~ factor(level, levels=1:8),
                                    data=data.frame(n_tox = n_tox,
                                                    level = level))
                        )

  ## The 'dfcrm' & 'rusti' implementations take patient-wise vectors x and y,
  ## and patient-wise weights w that 'break exchangeability'.
  x <- skeleton[rep(level, each=3)]
  coh <- function(tox, n=3) c(rep(0L,n-tox), rep(1L,tox))
  y <- as.integer(sapply(y, coh))
  w <- rep(1,length(y))
  w[y == 1] <- 0.0 # encode y in w for rusti's benefit (dfcrm is unaffected)

  a <- seq(-10, 10, 0.1)
  objectives <- data.table(
    u = a
  , dfcrm = dfcrm::crmh(a, x, y, w, s)
  , rusti = crmh(a, log(x), w, s)
  , ruste = crmh_ev(a, obs, log(skeleton), s, 0)
  )

  ## Calculate deltas vs dfcrm::crmh
  deltas <- objectives
  deltas[, `:=`(
    rusti = rusti - dfcrm
   ,ruste = ruste - dfcrm
 )]

  ## Arrange this data frame *tall*
  talldata <- melt(deltas[,c('u','rusti','ruste')], measure.vars = c("rusti","ruste"),
                   variable.name = "impl", value.name = "delta.f")

  Hmisc::xYplot(delta.f ~ u, groups = impl, data = talldata,
                type='l', auto.key=TRUE)
}
