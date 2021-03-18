## This file is copied directly from the 'dfcrm.R' source file
## of Ken Cheung's package 'dfcrm', and key objective functions
## have been optimized, with microbenchmarked speedups noted.

## Parameters of objective functions are as follows:
## a: The single parameter of 1-parameter CRM. This function is
##    vectorized over this parameter, as required by 'integrate'.
## x: Vector of prior probabilities of tox at administered doses
## y: Vector of 0/1 indicators of toxicity, patient-wise
## w: A vector of weights?
## s: A (scalar) scale factor


## These 3 functions reimplement the 'empiric' (aka 'power') model
## integrands, achieving a 2x speedup for crm(impl="Ri") vs "dfcrm".
crmh = compiler::cmpfun(function(a,x,y,w,s) {  ## posterior
  exp_a_ <- exp(a)
  v = exp( -0.5*(a/s)^2  +  sum(log(x[y==1])) * exp_a_ )
  for (i in which(y==0)) {
    v = v * (1 - w[i] * x[i]^exp_a_)
  }
  v[!is.finite(v)] <- 0.0
  return(v)
})

crmht = compiler::cmpfun(function(a,x,y,w,s) { ## posterior times x
  exp_a_ <- exp(a)
  v = a * exp( -0.5*(a/s)^2  +  sum(log(x[y==1])) * exp_a_ )
  for (i in which(y==0)) {
    v = v * (1 - w[i] * x[i]^exp_a_)
  }
  v[!is.finite(v)] <- 0.0
  return(v)
})

crmht2 = compiler::cmpfun(function(a,x,y,w,s) { ## posterior times x^2
  exp_a_ <- exp(a)
  v = a^2 * exp( -0.5*(a/s)^2  +  sum(log(x[y==1])) * exp_a_ )
  for (i in which(y==0)) {
    v = v * (1 - w[i] * x[i]^exp_a_)
  }
  v[!is.finite(v)] <- 0.0
  return(v)
})

## To avoid a NOTE on package check, a dummy '...' formal argument is written.
## But this should not be taken as an invitation to provide actual arguments!
bt_crms <- function(...) {
  if (length(list(...))) stop("Don't invoke bt_crms with arguments; it's a trap!")
  benchtab(crm(prior=c(0.05, 0.10, 0.20, 0.35, 0.50, 0.70),
               target=0.2,
               tox=c(0L, 0L, 1L, 0L, 0L, 1L, 1L, 0L, 0L, 0L),
               level=c(3, 4, 4, 3, 3, 4, 3, 2, 2, 2),
               ...),
           impl=c("dfcrm","Ri","rusti","rustq"))
}

## This set-up is verbatim from the example in dfcrm::crm
comp_crms <- function(prior = c(0.05, 0.10, 0.20, 0.35, 0.50, 0.70),
                      target = 0.2,
                      level = c(3, 4, 4, 3, 3, 4, 3, 2, 2, 2),
                      y = c(0, 0, 1, 0, 0, 1, 1, 0, 0, 0),
                      impl = c("dfcrm","Ri","rusti","rustq")) {
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

## Figure out why I'm having trouble invoking crm(impl="rustq") ...
qcrm <- function(prior = c(0.05, 0.10, 0.20, 0.35, 0.50, 0.70),
                 target = 0.2,
                 level = c(3, 4, 4, 3, 3, 4, 3, 2, 2, 2),
                 y = c(0, 0, 1, 0, 0, 1, 1, 0, 0, 0),
                 impl = "dfcrm") {
  crm(prior, target, y, level, impl=impl)
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
  new <- crmh(a,x,y,w,s)
  rust <- rcrmh(a,x,y,w,s)
  delta <- max(abs(new - old))
  if (delta > .Machine$double.eps)
    cat("crmh |new - old| =", delta, "\n")
  deltar <- max(abs(rust - old))
  if (deltar > .Machine$double.eps)
    cat("crmh |rust - old| =", deltar, "\n")
  ##
  old <- dfcrm::crmht(a,x,y,w,s)
  new <- crmht(a,x,y,w,s)
  rust <- rcrmht(a,x,y,w,s)
  tdelta <- max(abs(new - old))
  if (tdelta > .Machine$double.eps)
    cat("crmht |new - old| =", tdelta, "\n")
  tdeltar <- max(abs(rust - old))
  if (tdeltar > .Machine$double.eps)
    cat("crmht |rust - old| =", tdeltar, "\n")
  ##
  old <- dfcrm::crmht(a,x,y,w,s)
  new <- crmht(a,x,y,w,s)
  rust <- rcrmht(a,x,y,w,s)
  t2delta <- max(abs(new - old))
  if (t2delta > .Machine$double.eps)
    cat("crmht2 |new - old| =", t2delta, "\n")
  t2deltar <- max(abs(rust - old))
  if (t2deltar > .Machine$double.eps)
    cat("crmht2 |rust - old| =", t2deltar, "\n")

  t_old <- microbenchmark::microbenchmark(integrate(dfcrm::crmh, -Inf, Inf, x, y, w, s))
  t_new <- microbenchmark::microbenchmark(integrate(crmh, -Inf, Inf, x, y, w, s))
  t_rust <- microbenchmark::microbenchmark(integrate(rcrmh, -Inf, Inf, x, y, w, s))
  print(t_old)
  print(t_new)
  print(t_rust)
  speedup <- median(t_old$time)/median(t_new$time)
  message("speedup (new): ", round(speedup,2), "x")
  speedupr <- median(t_old$time)/median(t_rust$time)
  message("speedup (rust): ", round(speedupr,2), "x")
  invisible(list(old=t_old, new=t_new, rust=t_rust))
})

## TODO: Re-implement and test these objective functions:

## crmhlgt <- function(a,x,y,w,s,alp0)  { ## posterior logit model
##   v = exp(-a^2/2/s^2)
##   for (i in 1:length(x)) {
##  	PSI <- (1 + exp(-alp0-exp(a)*x[i]))^{-1}
##     	v <- v * (PSI^y[i]) * (1-w[i]*PSI)^(1-y[i])
##   }
##   return(v)
## }
## crmhtlgt <- function(a,x,y,w,s,alp0)  { ## posterior times x
##   v = a * exp(-a^2/2/s^2)
##   for (i in 1:length(x)) {
##  	PSI <- (1 + exp(-alp0-exp(a)*x[i]))^{-1}
##     	v <- v * (PSI^y[i]) * (1-w[i]*PSI)^(1-y[i])
##   }
##   return(v)
## }
## crmht2lgt <- function(a,x,y,w,s,alp0)  { ## posterior times x^2
##   v = a^2 * exp(-a^2/2/s^2)
##   for (i in 1:length(x)) {
##  	PSI <- (1 + exp(-alp0-exp(a)*x[i]))^{-1}
##     	v <- v * (PSI^y[i]) * (1-w[i]*PSI)^(1-y[i])
##   }
##   return(v)
## }
## lcrm <- function(a,x,y,w) { #loglikelihood of empiric function
##   v <- 0
##   for (i in 1:length(x))
##     v <- v + y[i]*log(x[i])*exp(a) + (1-y[i])*log(1 - w[i]*x[i]^exp(a))
##   return(v)
## }
## lcrmlgt <- function(a,x,y,w,alp0) { #loglikelihood of logit function
##   v <- 0
##   for (i in 1:length(x)) {
##     PSI <- (1 + exp(-alp0-exp(a)*x[i]))^{-1}
##     v <- v + y[i]*log(PSI) + (1-y[i])*log(1-w[i]*PSI)
##   }
##   return(v)
## }

## This local (as-yet, unexported) copy of dfcrm:crm serves as a test harness
## for various performance tuning experiments. The 'impl' arg allows selection
## of various alternative implementations:
## - Ri uses the improved (roughly 2x speedup) R integrands [crmh,crmht,crmh2] above
## - rusti substitutes integrands [rcrmh,rcrmht,rcrmht2] written in Rust
## - rustq uses Rust for the quadrature itself
##' @importFrom stats integrate optimize qnorm
##' @importFrom dfcrm crmhlgt crmhtlgt crmht2lgt
##' @importFrom dfcrm lcrm lcrmlgt
crm <- compiler::cmpfun(function(prior, target, tox, level, n=length(level),
                dosename=NULL, include=1:n, pid=1:n, conf.level=0.90,
                method="bayes", model="empiric", intcpt=3,
                scale=sqrt(1.34), model.detail=TRUE, patient.detail=TRUE, var.est=TRUE,
                impl=c("Ri","rusti","rustq","dfcrm")) { # implementation switch
  if (impl[1]=="dfcrm")
    return(dfcrm::crm(prior, target, tox, level, n, dosename, include, pid, conf.level,
                      method, model, intcpt, scale, model.detail, patient.detail, var.est))
  if (method != "bayes")
    warning("CRM method '", method, "' has no '", impl, "' implementation as yet.")
  if (model != "empiric")
    warning("The '", model, "' model has no '", impl, "' implementation as yet.")
  ## NB: If we reach this point without warnings, then we're estimating
  ##     the empiric model using Bayes method. This is the test case for
  ##     my alternative implementations.
  y1p <- tox[include]
  y1p <- as.integer(y1p) # Rust methods require this
  w1p <- rep(1,length(include))
  if (model=="empiric") {
    dosescaled <- prior

    x1p <- prior[level[include]]
    if (method=="mle") {
      if (sum(y1p)==0 | sum(y1p)==length(y1p)) stop(" mle does not exist!")
      est <- optimize(lcrm,c(-10,10),x1p,y1p,w1p,tol=0.0001,maximum=TRUE)$max
      if (var.est) { e2 <- integrate(crmht2,-Inf,Inf,x1p,y1p,w1p,500,abs.tol=0)[[1]] / integrate(crmh,-Inf,Inf,x1p,y1p,w1p,500,abs.tol=0)[[1]]; }
    }
    else if (method=="bayes") {
      switch(impl[1],
             Ri = {
               den <- integrate(crmh,-Inf,Inf,x1p,y1p,w1p,scale,abs.tol=0)[[1]]
               est <- integrate(crmht,-Inf,Inf,x1p,y1p,w1p,scale,abs.tol=0)[[1]] / den
               if (var.est)
                 e2 <- integrate(crmht2,-Inf,Inf,x1p,y1p,w1p,scale,abs.tol=0)[[1]] / den
             },
             rusti = {
               den <- integrate(rcrmh,-Inf,Inf,x1p,y1p,w1p,scale,abs.tol=0)[[1]]
               est <- integrate(rcrmht,-Inf,Inf,x1p,y1p,w1p,scale,abs.tol=0)[[1]] / den
               if (var.est)
                 e2 <- integrate(rcrmht2,-Inf,Inf,x1p,y1p,w1p,scale,abs.tol=0)[[1]] / den
             },
             rustq = {
               den <- icrm(x1p, y1p, w1p, scale, 0)
               est <- icrm(x1p, y1p, w1p, scale, 1) / den
               if (var.est)
                 e2 <- icrm(x1p, y1p, w1p, scale, 2) / den
             })
    }
    else { stop(" unknown estimation method"); }
    ptox <- prior^exp(est)
    if (var.est) {
      post.var <- e2-est^2
      crit <- qnorm(0.5+conf.level/2)
      lb <- est - crit*sqrt(post.var)
      ub <- est + crit*sqrt(post.var)
      ptoxL <- prior^exp(ub)
      ptoxU <- prior^exp(lb)
    }
  }
  else if (model=="logistic") {
    dosescaled <- log(prior/(1-prior)) - intcpt
    if (!all(dosescaled<0)) {
      stop( "Intercept parameter in logit model is too small: scaled doses > 0!")
    }

    x1p <- dosescaled[level[include]]

    if (method=="mle") {
      if (sum(y1p)==0 | sum(y1p)==length(y1p)) stop(" mle does not exist!")
      est <- optimize(lcrmlgt,c(-10,10),x1p,y1p,w1p,intcpt,tol=0.0001,maximum=TRUE)$max
      if (var.est) { e2 <- integrate(crmht2lgt,-Inf,Inf,x1p,y1p,w1p,500,intcpt,abs.tol=0)[[1]] / integrate(crmhlgt,-Inf,Inf,x1p,y1p,w1p,500,intcpt,abs.tol=0)[[1]]; }
    }
    else if (method=="bayes") {
      den <- integrate(crmhlgt,-Inf,Inf,x1p,y1p,w1p,scale,intcpt,abs.tol=0)[[1]]
      est <- integrate(crmhtlgt,-Inf,Inf,x1p,y1p,w1p,scale,intcpt,abs.tol=0)[[1]] / den
      if (var.est) { e2 <- integrate(crmht2lgt,-Inf,Inf,x1p,y1p,w1p,scale,intcpt,abs.tol=0)[[1]] / den; }
    }
    else { stop(" unknown estimation method"); }
    ptox <- (1 + exp(-intcpt-exp(est)*dosescaled))^{-1}
    if (var.est) {
      post.var <- e2-est^2
      crit <- qnorm(0.5+conf.level/2)
      lb <- est - crit*sqrt(post.var)
      ub <- est + crit*sqrt(post.var)
      ptoxL <- (1 + exp(-intcpt-exp(ub)*dosescaled))^{-1}
      ptoxU <- (1 + exp(-intcpt-exp(lb)*dosescaled))^{-1}
    }
  }
  else { stop(" model specified not available."); }

  if (all(ptox<=target)) { rec <- length(prior); }
  else if (all(ptox>=target)) { rec <- 1; }
  else { rec <- order(abs(ptox-target))[1]; }
  if (!var.est) { post.var <- ptoxU <- ptoxL <- NA; }
  foo <- list(prior=prior, target=target, tox=tox, level=level,
              dosename=dosename, subset=pid[include], estimate=est,
              model=model, prior.var=scale^2, post.var=post.var,method=method,
              mtd=rec, include=include, pid=pid, model.detail=model.detail,intcpt=intcpt,
              ptox=ptox, ptoxL=ptoxL, ptoxU=ptoxU, conf.level=conf.level,
              patient.detail=patient.detail,tite=FALSE,dosescaled=dosescaled,var.est=var.est)
  class(foo) <- "mtd"
  foo
})

