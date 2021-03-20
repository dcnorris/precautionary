## Compare the objective functions visually
comp_objfn <- function(s = 500) {
  x <- c(1,2,3,2,3,3)*0.1
  y <- c(0L,0L,1L,0L,0L,1L)
  w <- rep(1,length(y))

  a <- seq(-10, 10, 0.1)
  objectives <- data.table(
    u = a
  , dfcrm = dfcrm::crmh(a, x, y, w, s)
  , Ri = precautionary:::crmh(a, x, y, w, s)
  , rusti = rcrmh(a, x, y, w, s)
  )

  ## Calculate deltas vs dfcrm::crmh
  deltas <- objectives
  deltas[, `:=`(
    Ri = Ri - dfcrm,
    rusti = rusti - dfcrm
  )]

  ## Arrange this data frame *tall*
  talldata <- melt(deltas[,c('u','Ri','rusti')], measure.vars = c("Ri","rusti"),
                   variable.name = "impl", value.name = "delta.f")

  Hmisc::xYplot(delta.f ~ u, groups = impl, data = talldata,
                type='l', auto.key=TRUE)
}
