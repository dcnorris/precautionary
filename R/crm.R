## Introducing an R6 class for CRM models, capable of supporting
## efficient caching maneuvers.

## TODO: Organize a class hierarchy of model types (empiric, logistic),
##       exploiting inheritance to clarify what is special vs shared.

##' @name Crm-class
##' @title An R6 class encapsulating CRM models
##'
##' This class wraps the functionality of package \CRANpkg{dfcrm}, adding
##' efficient Rust implementations of some numerical routines.
##'
##' @details
##' Syntactically, the method chaining supported by R6 classes makes the
##' invocation of CRM models more transparent. The mutability conferred
##' by reference semantics enables memoization (caching) of results, which
##' can speed up DTP calculations significantly.
##'
##' Presently, this class supports only the 'empiric' (aka 'power') model.
##' But it is hoped that inheritance will assist in rendering other models
##' implemented in package \CRANpkg{dfcrm} clearly, with code reuse.
##' @export
Crm <- R6::R6Class("Crm",
               public = list(
                 ##' @details
                 ##' Create a new `Crm` object.
                 ##'
                 ##' @param skeleton CRM skeleton
                 ##' @param scale Sigma parameter of prior on beta parameter
                 ##' @param target Target toxicity rate
                 ##' @param cohort.size Size of enrolled cohorts
                 ##' @return A Crm object.
                 ##'
                 ##' @examples
                 ##' ## TODO.
                 initialize = function(skeleton, scale = sqrt(1.34), target, cohort.size = 3) {
                   private$ln_skel = log(skeleton)
                   private$scale = scale
                   private$cohort.size = 3
                   private$target = target
                 }
                 ##' @details
                 ##' Set the stopping function
                 ##'
                 ##' @param sfunc A function taking \code{mtd} objects to \code{mtd} objects,
                 ##' attaching suitable stopping information
                 ##' @return Self, invisibly
                ,stop_func = function(sfunc){
                  private$.stop_func <- sfunc
                  ## TODO: Does this invalidate CRM model cache? (I don't believe so.)
                  invisible(self)
                }
                 ##' @details
                 ##' Set the \code{no_skip_esc} behavior
                 ##'
                 ##' @param tf An atomic logical value, TRUE or FALSE
                 ##' @return Self, invisibly
                ,no_skip_esc = function(tf){
                  stopifnot(is.logical(tf))
                  stopifnot(length(tf) == 1)
                  private$.no_skip_esc = tf
                  invisible(self)
                }
                 ##' @details
                 ##' Set the \code{no_skip_deesc} behavior
                 ##'
                 ##' @param tf An atomic logical value, TRUE or FALSE
                 ##' @return Self, invisibly
                ,no_skip_deesc = function(tf){
                  stopifnot(is.logical(tf))
                  stopifnot(length(tf) == 1)
                  private$.no_skip_deesc = tf
                  invisible(self)
                }
                 ##' @details
                 ##' Set the \code{global_coherent_esc} behavior
                 ##'
                 ##' @param tf An atomic logical value, TRUE or FALSE
                 ##' @return Self, invisibly
                ,global_coherent_esc = function(tf){
                  stopifnot(is.logical(tf))
                  stopifnot(length(tf) == 1)
                  private$.global_coherent_esc = tf
                  invisible(self)
                }
                 ##' @details
                 ##' Set the required confidence level for escalation decisions
                 ##'
                 ##' @param conf A numeric confidence less than 1.0
                 ##' @return Self, invisibly
                ,conf_level = function(conf){
                  private$conf.level <- conf
                  invisible(self)
                }
                 ##' @details
                 ##' Set patient-wise toxicity observations
                 ##'
                 ##' @param level A patient-wise vector of dose assignments
                 ##' @param tox A patient-wise vector of 0/1 toxicity assessments
                 ##' @return Self, invisibly
                ,observe = function(level, tox){ # TODO: Preserve order for dfcrm's sake?
                  stopifnot(length(level) == length(tox))
                  stopifnot(all(tox %in% c(0,1)))
                  D <- length(private$ln_skel)
                  stopifnot(all(level %in% 1:D))
                  level <- factor(level, levels=1:D)
                  self$tally(x = xtabs(tox ~ level),
                             o = xtabs(!tox ~ level))
                }
                 ##' @details
                 ##' Set dose-wise toxicity observations
                 ##'
                 ##' @param x A dose-wise vector of toxicity counts
                 ##' @param o A dose-wise vector of non-toxicity counts
                 ##' @return Self, invisibly
                ,tally = function(x, o){
                  D <- length(private$ln_skel)
                  ## As a special case, we allow 'o' to be missing, and 'x' to be
                  ## a full x/o tally encoded in a numeric(1) -- i.e. a single f64.
                  if (missing(o) && is.numeric(x) && length(x) == 1) {
                    if (private$cohort.size == 3) {
                      enr <- sapply(1:8, function(d) (x %% 6^d) %/% 6^(d-1))
                      tox <- sapply(-(0:7), function(d) ((x %% 1) %% 16^d) %/% 16^(d-1))
                      nos <- 3 * enr - tox
                      ## Ideally, the tox and no-tox individuals are contiguous.
                      ## This is especially helpful (numerically) for the no-tox (w>0) cases.
                      private$w <- c(rep(1, sum(nos)),
                                     rep(0, sum(tox)))
                      ## Now, we need to lay out a vector of dose-level assignments
                      ## consistent with the given data.
                      private$level <- c(rep(1:D, nos[1:D]),
                                         rep(1:D, tox[1:D]))
                    } else {
                      stop("unimplemented")
                      ## TODO: implement n=1 case:
                      ## - enroll at most 15 per dose ==> 8*4 = 32 bits
                      ## - suppose up to 8 tox may occur at each dose
                      ## - then 8*log2(8+1) = 27 bits needed
                      ## (32+27=59 bits just fits in a float if I use 9 bits of exponent)
                    }
                    invisible(self)
                  }

                  ## Otherwise ...
                  ## we expect 'x' to be a dosewise vector of exchangeable toxicity counts
                  if (length(x) != D) {
                    cat("x =", deparse(x), "\n")
                  }
                  stopifnot(length(x) == D)
                  stopifnot(is.numeric(x))
                  stopifnot(all(round(x) == x))
                  private$x <- as.integer(x)

                  ## And we expect a dosewise 'o'
                  stopifnot(length(o) == D)
                  ## that is either a numeric vector of exchangeable nontox counts ...
                  if (is.numeric(o)) {
                    private$o <- as.integer(o)
                    private$level <- c(rep(1:D, o),
                                       rep(1:D, x))
                    private$w <- c(rep(1, sum(o)),
                                   rep(0, sum(x)))
                  } else { # .. or a _list_ of weight vectors breaking exchangeability.
                    ## TODO: Consider a more functional construction of $w and $level:
                    private$w <- numeric(length(o) +
                                         sum(x)) # initialized to 0, thus encoding y
                    ix <- 0
                    for (level in 1:length(o)) {
                      len <- length(o[[level]])
                      private$w[ix + 1:len] <- o[[level]]
                      private$level[ix + 1:len] <- level
                      ix <- ix + len
                    }
                    private$o <- lapply(o, length)
                  }
                  invisible(self)
                }
                ##' @details
                ##' Set patient-wise toxicity observations
                ##'
                ##' @param impl A string choosing the low-level implementation to use.
                ##' Possible values include \code{"dfcrm"}, \code{"rusti"} and \code{"ruste"}.
                ##' @param abbrev Logical; if TRUE (the default), an abbreviated \code{mtd}
                ##' object is returned to save execution time. If FALSE, a complete object is
                ##' returned, suitable for regression testing against package \CRANpkg{dfcrm}.
                ##' @return An object of class \code{mtd} as per package \CRANpkg{dfcrm}
               ,est = function(impl, abbrev=TRUE){
                 scale <- private$scale
                 model <- "empiric"
                 method <- "bayes"
                 include <- seq_along(private$w)
                 pid <- include
                 switch(impl,
                        dfcrm = {
                          return(
                            dfcrm::crm(prior = exp(private$ln_skel), target = private$target,
                                       tox=(private$w==0), level=private$level,
                                       n=length(private$w), dosename=NULL,
                                       include=include, pid=pid, conf.level=private$conf.level,
                                       method=method, model=model, intcpt=3,
                                       scale=scale, model.detail=TRUE, patient.detail=TRUE,
                                       var.est=TRUE)
                          )
                        }
                       ,rusti = {
                         ln_x1p <- private$ln_skel[private$level]
                         w <- private$w
                         m0 <- integrate(crmh  ,-Inf,Inf, ln_x1p, w, scale, abs.tol=0)[[1]]
                         m1 <- integrate(crmht ,-Inf,Inf, ln_x1p, w, scale, abs.tol=0)[[1]]
                         m2 <- integrate(crmht2,-Inf,Inf, ln_x1p, w, scale, abs.tol=0)[[1]]
                       }
                      ,ruste = {
                        ln_p <- private$ln_skel
                        tox <- private$x
                        nos <- private$o
                        m0 <- integrate(crmh_xo,-Inf,Inf,ln_p,tox,nos,scale,0,abs.tol=0)[[1]]
                        m1 <- integrate(crmh_xo,-Inf,Inf,ln_p,tox,nos,scale,1,abs.tol=0)[[1]]
                        m2 <- integrate(crmh_xo,-Inf,Inf,ln_p,tox,nos,scale,2,abs.tol=0)[[1]]
                      }
                      ,stop("must specify impl in Crm$est()")
                      )
                 estimate <- m1/m0
                 post.var <- m2/m0 - estimate^2
                 prior <- exp(private$ln_skel)
                 ptox <- prior^exp(estimate)
                 ans <- list(prior = prior,       # used by stop_for_excess_toxicity_empiric
                             estimate = estimate, # used by stop_for_excess_toxicity_empiric
                             post.var = post.var, # used by stop_for_excess_toxicity_empiric
                             level = private$level, # used by stop_for_consensus_reached
                             mtd = which.min(abs(ptox - private$target)))
                 if (abbrev)
                   return(ans)
                 ## Otherwise we build out the whole return value of class "mtd"
                 crit <- qnorm(0.5 + private$conf.level/2)
                 lb <- estimate - crit*sqrt(post.var)
                 ub <- estimate + crit*sqrt(post.var)
                 ## TODO: Might 'within(.)' itself be costly? Try flattening
                 ##       this code to a series of straight component assignments.
                 ans <- within(ans, {
                   target <- private$target
                   conf.level <- private$conf.level
                   ptox <- ptox
                   ptoxL <- prior^exp(ub)
                   ptoxU <- prior^exp(lb)
                   tox <- 1.0*(private$w==0) # NB: dfcrm's "mtd" holds tox as double
                   dosename <- NULL
                   subset <- pid[include]
                   model <- model
                   prior.var <- scale^2
                   method <- method
                   include <- include
                   pid <- pid
                   model.detail <- TRUE
                   intcpt <- 3 # TODO: Un-hard-code this
                   ptoxL <- ptoxL
                   ptoxU <- ptoxU
                   patient.detail <- TRUE
                   tite <- FALSE
                   dosescaled <- prior # TODO: this is correct for EMPIRIC MODEL ONLY
                   var.est <- TRUE
                 })

                 ## Order components to enable testthat::expect_equal()
                 ans <- ans[c("prior","target","tox","level","dosename","subset",
                              "estimate","model","prior.var","post.var","method",
                              "mtd","include","pid","model.detail","intcpt",
                              "ptox","ptoxL","ptoxU","conf.level","patient.detail",
                              "tite","dosescaled","var.est")]
                 class(ans) <- "mtd"
                 return(ans)
               } #</est()>
               ##' @details
               ##' Return dose recommendation for given tox/no-tox tallies.
               ##'
               ##' @param x A dose-wise vector of toxicity counts
               ##' @param o A dose-wise vector of non-toxicity counts
               ##' @param last_dose The most recently given dose, as required to implement
               ##' the \code{global_coherent_esc=TRUE} behavior
               ##' @param ... Optional parameters passed to \code{Crm$esc()}, enabling
               ##' the passthru of the \code{impl} parameter.
               ##' @return An object of class \code{mtd} as per package \CRANpkg{dfcrm}
              ,applied = function(x, o, last_dose = NA, ...){
                level <- which((x+o) > 0) # vector of levels that have been tried
                est <- self$tally(x, o)$est(...)
                if (private$.no_skip_esc) {
                  est$mtd <- min(est$mtd, max(level) + 1)
                }
                if (private$.no_skip_deesc) {
                  est$mtd <- max(est$mtd, min(level) - 1)
                }
                if (private$.global_coherent_esc) {
                  if(is.na(last_dose)) stop("last_dose required if global_coherent_esc = TRUE")
                  if (x[last_dose] > private$target * (x+o)[last_dose]) {
                    est$mtd <- min(est$mtd, last_dose)
                  }
                }
                if (!is.null(private$.stop_func)) {
                  est = private$.stop_func(est)
                }
                return(est)
              }
              ) # </public>
              ,private = list(
                 ln_skel = NA
               , scale = NA
               , target = NA
               , cohort.size = 3 # TODO: Factor this out? It seems like a 'wart'.
               , .stop_func = NULL
               , .no_skip_esc = TRUE     # These would seem to be the
               , .no_skip_deesc = FALSE  # 'safest' defaults.
               , .global_coherent_esc = TRUE  # This sounds like a good thing.
               , x = NA
               , o = NA
               , w = NA # in general, this will encode y as well
               , level = NA # patient-wise dose level, such that ln_skel[l] aligns with w
               , obs = NaN
               , cached = list()
               , conf.level = 0.90
               )
               )

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
##' @author Adapted by David C. Norris, from Ken Cheung's \CRANpkg{dfcrm}
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
               ln_skel <- log(prior)
               den <- integrate(crmh_ev,-Inf,Inf,obs,ln_skel,scale,0,abs.tol=0)[[1]]
               est <- integrate(crmh_ev,-Inf,Inf,obs,ln_skel,scale,1,abs.tol=0)[[1]] / den
               if (var.est)
                 e2 <- integrate(crmh_ev,-Inf,Inf,obs,ln_skel,scale,2,abs.tol=0)[[1]] / den
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

