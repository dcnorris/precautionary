#' Simulate trials defined via package \code{escalation}
#' 
#' @docType methods
#' @param selector_factory An object of S3 class \code{selector_factory}
#' @param num_sims A vector of simulation dimensions
#' @param true_prob_tox A generator of toxicity distributions
#' @param ... Passed to subroutines
#' 
#' @importFrom escalation simulate_trials selector_factory
#' @export
setGeneric("simulate_trials")

setOldClass("ordtox")
setOldClass(c("three_plus_three_selector_factory",
              "tox_selector_factory",
              "selector_factory"))
setOldClass(c('boin_selector_factory',
              'tox_selector_factory',
              'selector_factory'))
setOldClass(c('dfcrm_selector_factory',
              'tox_selector_factory',
              'selector_factory'))
# TODO: Register also various hierarchies for derived_dose_selector_factory?

# Interestingly, package 'escalation' needs to implement methods for
# 'num_doses' and 'dose_indices' only for the (post-fit) selectors,
# and not for the selector_factory. But it is characteristic of my
# 'ordtox' class that it involves itself with the number of doses
# under investigation.

#' @export
num_doses.three_plus_three_selector_factory <- function(x, ...) {
  return(x$num_doses)
}

#' @export
num_doses.dfcrm_selector_factory <- function(x, ...) {
  return(length(x$skeleton))
}

#' @export
num_doses.boin_selector_factory <- function(x, ...) {
  return(x$num_doses)
}

#' @export
dose_indices.ordtox <- function(x, ...) {
  n <- num_doses(x)
  if(n > 0) {
    return(1:n)
  } else {
    return(integer(length = 0))
  }
}

#' @export
dose_indices.default <- function(x, ...) {
  n <- num_doses(x)
  if(n > 0) {
    return(1:n)
  } else {
    return(integer(length = 0))
  }
}

setMethod("simulate_trials", c(selector_factory="ordtox"),
  function(selector_factory, num_sims, true_prob_tox, ...){
    protocol <- selector_factory # separate naming from implementation details
    lambda_CV <- 3
    median_mtd <- num_doses(protocol) - 1
    median_sd <- median_mtd/3
    r0 <- seq(0.5, 2.5, 0.5)
    K <- num_sims
    M <- num_sims
    # TODO: Consider making true_prob_tox argument optionally
    #       a data table (or data frame) of the following form,
    #       _OR_ (more usefully) a *generative model* for same.
    true_prob_tox <- data.table(CV = rexp(n=K, rate=lambda_CV))
    true_prob_tox[, `:=`(
      sigma = sqrt(log(CV^2 + 1)) # NB: CV>sigma, w/ near-identity below 0.4
      , mu = rnorm(n=K, mean=median_mtd, sd=median_sd)
    )]
    P_ <- paste0("P", dose_indices(protocol))
    true_prob_tox[, c(P_) := lapply(dose_indices(protocol),
                                    function(q) pnorm(q, mean=mu, sd=sigma))]
    ensembles <<- list()
    toxicities <<- list()
    for(k in 1:K){
      cat("k =", k, "\n")
      # Now invoke the default method from package 'escalation':
      sims <- callNextMethod(protocol, num_sims,
                             true_prob_tox = as.numeric(true_prob_tox[k, ..P_]))
      ensemble <- rbindlist(lapply(sims[[1]], function(.) .[[1]]$fit$outcomes)
                            , idcol = "rep")
      ensemble <- merge(data.frame(r0=r0), ensemble) # cartesian product
      ensemble <- as.data.table(ensemble) # restore data.table after cartesian product
      ensemble[, MTDi3 := qnorm(p = u_i
                                , mean = true_prob_tox[k]$mu
                                , sd = true_prob_tox[k]$sigma
      )]
      ensemble[, `:=`(
        MTDi1 = MTDi3 - 2*r0
        , MTDi2 = MTDi3 - r0
        , MTDi4 = MTDi3 + r0
        , MTDi5 = MTDi3 + 2*r0
      )]
      ensemble[, toxgrade := (dose>MTDi1) + (dose>MTDi2) + (dose>MTDi3) +
                 (dose>MTDi4) + (dose>MTDi5)]
      stopifnot(with(ensemble, all(tox == as.integer(dose >= MTDi3))))
      # N.B.: Any trial realization ('rep') is just as likely as any other.
      #       Therefore, I have to aggregate with reps as the denominator.
      #       Specifically, pooling all the patients would fail to weight
      #       them inversely to trial size.
      #       Another way to see this is, we are asking about the probability
      #       of a fatality in *this* trial. Alternatively, we can ask for
      #       the expected number of each grade of toxicity in *this* trial.
      #       Note that this should partition the expected size of this trial!
      ensemble[, Tox := factor(toxgrade, levels=0:5, labels=paste0("Gr",0:5))]
      ensembles[[k]] <<- ensemble
      # See https://stackoverflow.com/a/16519612/3338147 explaining below syntax
      toxicities[[k]] <<- ensemble[, .(Tox=ordered(levels(Tox)), N=c(table(Tox)))
                                   , by=r0]
      # Importantly, because these protocols run considering only DLT = Gr>=3,
      # they are r0-agnostic. But I ought to allow for the general case where
      # ordinal toxicities (MTDig for g != 3) affect trial conduct.
      #
      # Initially, let me go simply for *counts* and worry secondarily about
      # obtaining (and dividing by) the denominators.
      # Want a tabulation with several r0 values defining rows, and columns for
      # the counts of toxicity grades. A final column may show total enrollments,
      # which as noted above will all be equal unless ordinal toxicities affect
      # escalation or termination decisions.
    }
    # Indeed, I like the idea of abstracting the calculation of high-level
    # summary statistics into a separate function, possibly even a 'summary'
    # method for a suitably defined S3 class. But for now, let me implement
    # these summaries here, as additional components of the returned list.
    toxdt <- rbindlist(toxicities, idcol = "k")
    counts <- toxdt[, .(n=sum(N)), by=.(r0,Tox)]
    toxtab <- ftable(xtabs(n ~ Tox + r0, counts))
    # It now seems to me there are 2 perspectives on the probabilities here.
    # One POV is the per-trial perspective, which might ask for expected
    # numbers of each grade of toxicity.
    # The second POV is that of the enrolling patient, who might ask what
    # is the probability of each grade of toxicity. But this POV ought to
    # be conditioned on the prior results seen in the trial, and so demands
    # much more sophisticated modeling -- indeed, modeling of the kind I
    # employed in the AFM11 paper.
    # Note in fact that, as soon as such a perspective BEGINS to be acknowledged,
    # I have already 'won' the argument about dose individualization!
    list(tpt = true_prob_tox
        ,toxdt = toxdt
        ,toxtab = toxtab
        ,expect = toxtab/(K*M)
        ,K = K
        ,M = M
        )
  }
)


#' Syntactic sugar for injecting safety analysis into \code{simulate_trials}
#'
#' @param x A selector_factory, required as first argument to \code{simulate_trials}
#'
#' @return x with class \sQuote{ordtox} prepended
#' @export
#'
#' @examples
#' 
#' # Simulate a 3+3 trial
#' protocol <- get_three_plus_three(num_doses = 6)
#' mu <- 5              # Values estimated in arXiv:2004.12755 [stat, q-bio];
#' sigma <- 1/sqrt(1.3) # see subsection 'Model estimates' and Table 2, pp. 2-3.
#' dose_levels <- c(2, 6, 20, 60, 180, 400)
#' MTDi <- Lognormal$new(meanlog = mu, sdlog = sigma)
#' DLT_probs <- MTDi$cdf(dose_levels)
#' sims <- simulate_trials(protocol, num_sims = 20, true_prob_tox = DLT_probs)
#' 
#' # Do a richer simulation, examining safety qua avoidance of Grade 4 & 5
#' # toxicities, under a default generative model for ordinal toxicities:
#' safety_sims <- simulate_trials(check_safety(protocol))
check_safety <- function(x){
  stopifnot("x must be a selector_factory" = is(x, "selector_factory"))
  class(x) <- c("ordtox", class(x))
  x
}
