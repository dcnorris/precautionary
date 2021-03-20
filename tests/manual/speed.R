## To avoid a NOTE on package check, a dummy '...' formal argument is written.
## But this should not be taken as an invitation to provide actual arguments!
bt_crms <- function(...) {
  if (length(list(...))) stop("Don't invoke bt_crms with arguments; it's a trap!")
  benchtab(crm(prior=c(0.05, 0.10, 0.20, 0.35, 0.50, 0.70),
               target=0.2,
               tox=c(0L, 0L, 1L, 0L, 0L, 1L, 1L, 0L, 0L, 0L),
               level=c(3, 4, 4, 3, 3, 4, 3, 2, 2, 2),
               ...),
           impl=c("dfcrm","rusti","rustq"))
}

## This set-up is verbatim from the example in dfcrm::crm
comp_crms <- function(prior = c(0.05, 0.10, 0.20, 0.35, 0.50, 0.70),
                      target = 0.2,
                      level = c(3, 4, 4, 3, 3, 4, 3, 2, 2, 2),
                      y = c(0, 0, 1, 0, 0, 1, 1, 0, 0, 0),
                      impl = c("dfcrm","rusti","rustq")) {
  ## TODO: Get this working:
  ##comp_impl(quote(crm(prior, target, y, level, impl)), impl=impl)
  ## Until the above works, do it 'by hand':
  perf <- CJ(impl = impl, sorted = FALSE)
  for (i in seq(nrow(perf))) {
    crm_call <- substitute(crm(prior, target, y, level, impl=impl), perf[i,])
    perf$median.ms[i] <- median(microbenchmark::microbenchmark({
      eval(crm_call, perf[i,])
    })$time/1000000) # convert ns -> ms
    perf[i,'call'] <- deparse(crm_call)
  }
  perf
}

comp_icrm <- function(x=c(1,2,3,2,3,3)*0.1,
                      y=c(0L,0L,1L,0L,0L,1L),
                      w=rep(1,length(y)),
                      s=500) {
  mb_old <- microbenchmark::microbenchmark(
                              old0 <- integrate(dfcrm::crmh, -Inf, Inf, x, y, w, s)[[1]])
  old1 <- integrate(dfcrm::crmht, -Inf, Inf, x, y, w, s)[[1]]
  old2 <- integrate(dfcrm::crmht2, -Inf, Inf, x, y, w, s)[[1]]
  mb_rust <- microbenchmark::microbenchmark(rust0 <- icrm(x, y, w, s, 0))
  rust1 <- icrm(x, y, w, s, 1)
  rust2 <- icrm(x, y, w, s, 2)

  speedup <- mean(mb_old$time) / mean(mb_rust$time)
  cat(paste0("speedup: ", format(speedup, digits=2), "x\n"))

  ## TODO: Undertake a more comprehensive examination of convergence
  data.frame(old0 = old0, rust0 = rust0,
             old1 = old1, rust1 = rust1,
             old2 = old2, rust2 = rust2)
}

## TODO: Test in case where w not identically 1.
comp_crmh <- compiler::cmpfun(function(a=seq(-0.5, 0.5, 0.05),
                      x=c(1,2,3,2,3,3)*0.1,
                      y=c(0L,0L,1L,0L,0L,1L),
                      w=rep(1,length(y)), s=500) {
  old <- dfcrm::crmh(a,x,y,w,s)
  rust <- crmh(a,x,y,w,s)
  deltar <- max(abs(rust - old))
  if (deltar > .Machine$double.eps)
    cat("crmh |rust - old| =", deltar, "\n")
  ##
  old <- dfcrm::crmht(a,x,y,w,s)
  rust <- crmht(a,x,y,w,s)
  tdeltar <- max(abs(rust - old))
  if (tdeltar > .Machine$double.eps)
    cat("crmht |rust - old| =", tdeltar, "\n")
  ##
  old <- dfcrm::crmht(a,x,y,w,s)
  rust <- crmht(a,x,y,w,s)
  t2deltar <- max(abs(rust - old))
  if (t2deltar > .Machine$double.eps)
    cat("crmht2 |rust - old| =", t2deltar, "\n")

  t_old <- microbenchmark::microbenchmark(integrate(dfcrm::crmh, -Inf, Inf, x, y, w, s))
  t_rust <- microbenchmark::microbenchmark(integrate(crmh, -Inf, Inf, x, y, w, s))
  print(t_old)
  print(t_rust)
  speedup_message(t_old, t_rust)
  invisible(list(old=t_old, rust=t_rust))
})
