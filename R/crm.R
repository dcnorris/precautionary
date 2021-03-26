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


##' A package-local (as-yet, unexported) test harness adapted from dfcrm::crm().
##'
##' for various performance tuning experiments. The 'impl' arg allows selection
##' of various alternative implementations:
##' - rusti substitutes integrands crmh, crmht, crmht2 written in Rust
##' - dfcrm is the original as implemented in package \code{dfcrm}.
##' @param prior The CRM skeleton: dose-wise prior probabilities of toxicity
##' @param target Target toxicity rate
##' @param tox A patient-wise vector of toxicity counts
##' @param level A patient-wise vector of dose level assignments
##' @param n The number of patients enrolled
##' @param dosename Optional designators for the doses
##' @param include Index of patients to include
##' @param pid Vector of patient ID labels
##' @param conf.level TODO
##' @param method Estimation method:
##' @param model Presently, only the \sQuote{empiric} (or \sQuote{power}) model
##' has a Rust likelihood implementation.
##' @param intcpt Intercept for \sQuote{logistic} model
##' @param scale TODO: Investgate the history and intent of this scale factor
##' @param model.detail TODO
##' @param patient.detail TODO
##' @param var.est TODO: Appreciate the history and usage of this
##' @param impl Switch between \code{'rusti'} and \code{'dfcrm'} implementations.
##' Currently the \code{'rusti'} option is implemented only for the Bayes method
##' of the empirical (\sQuote{power}) model. An experimental \code{'ruste'}
##' implementaton is in the works.
##' @importFrom stats integrate optimize qnorm
##' @importFrom dfcrm crmhlgt crmhtlgt crmht2lgt
##' @importFrom dfcrm lcrm lcrmlgt
crm <- function(prior, target, tox, level, n=length(level),
                dosename=NULL, include=1:n, pid=1:n, conf.level=0.90,
                method="bayes", model="empiric", intcpt=3,
                scale=sqrt(1.34), model.detail=TRUE, patient.detail=TRUE, var.est=TRUE,
                impl=c("rusti","ruste","dfcrm")) { # implementation switch
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
  w1p[y1p == 1] <- 0; # encode y in w1p
  if (model=="empiric") {
    dosescaled <- prior

    x1p <- prior[level[include]]
    ln_x1p <- log(x1p) # Rust integrands {crmh,crmht,crmht2} take a log-skeleton to spare txops
    if (method=="mle") {
      if (sum(y1p)==0 | sum(y1p)==length(y1p)) stop(" mle does not exist!")
      est <- optimize(lcrm,c(-10,10),x1p,y1p,w1p,tol=0.0001,maximum=TRUE)$max
      if (var.est) { e2 <- integrate(crmht2,-Inf,Inf,ln_x1p,w1p,500,abs.tol=0)[[1]] / integrate(crmh,-Inf,Inf,x1p,y1p,w1p,500,abs.tol=0)[[1]]; }
    }
    else if (method=="bayes") {
      switch(impl[1],
             rusti = {
               den <- integrate(crmh,-Inf,Inf,ln_x1p,w1p,scale,abs.tol=0)[[1]]
               est <- integrate(crmht,-Inf,Inf,ln_x1p,w1p,scale,abs.tol=0)[[1]] / den
               if (var.est)
                 e2 <- integrate(crmht2,-Inf,Inf,ln_x1p,w1p,scale,abs.tol=0)[[1]] / den
             },
             ruste = {
               obs <- encode_cohorts(enr = tabulate(level, nbins=8)/3
                                    ,tox = xtabs(tox ~ factor(level, levels=1:8))
                                     )
               ln_prior <- log(prior)
               den <- integrate(crmh_ev,-Inf,Inf,obs,ln_prior,scale,0,abs.tol=0)[[1]]
               est <- integrate(crmh_ev,-Inf,Inf,obs,ln_prior,scale,1,abs.tol=0)[[1]] / den
               if (var.est)
                 e2 <- integrate(crmh_ev,-Inf,Inf,obs,ln_prior,scale,2,abs.tol=0)[[1]] / den
             },
             stop(paste("impl =", impl, "not recognized.")))
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
}

