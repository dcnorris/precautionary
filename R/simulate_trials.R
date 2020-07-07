#' @include toxicity_generators.R
NULL

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

# This is a simple generalization of escalation::simulate_trials,
# to the case where true_prob_tox is specified implicitly through
# a generative model ('prior') rather than explicitly as a vector.
# My aim here is to EXAMINE and ELABORATE the MEANING of these
# generative models as introduced in this source file.
# A also wish to EXPLORE and understand more fully the behavior
# and intent of the existing escalation::simulate_trials function,
# so that I can offer up a more focused extension and/or correction
# of its functionality. REMEMBER: to achieve 'depth', this package
# ought to 'correct escalation while explaining it'!

#' @examples
#' options(dose_levels = c(0.5, 1, 2, 4, 6, 8))
#' mtdi_gen <- hyper_mtdi_lognormal(lambda_CV = 3
#'                                  , median_mtd = 6, median_sd = 2
#'                                  , units="mg/kg")
#' hsims <- get_three_plus_three(num_doses = 6) %>%
#'   simulate_trials(num_sims = c(30, 10), true_prob_tox = mtdi_gen)
#' summary(hsims)
#' # Now attach a proper ordinalizer to 'mtdi_gen':
#' mtdi_gen@ordinalizer <- function(dose, r0) {
#'   c(Gr1=dose*r0^2, Gr2=dose*r0, Gr3=dose, Gr4=dose/r0, Gr5=dose/r0^2)
#' }
#' mtdi_gen@ordinalizer
#' hsimsOT <- get_three_plus_three(num_doses = 6) %>%
#'   simulate_trials(num_sims = c(30, 10), true_prob_tox = mtdi_gen, r0 = 1.5)
#' cat("class(hsimsOT) = ", class(hsimsOT), "\n")
#' summary(hsimsOT)
#' @rdname simulate_trials
setMethod(
  "simulate_trials"
  , c(selector_factory="selector_factory",
      num_sims="numeric",
      true_prob_tox="hyper_mtdi_distribution"),
  function(selector_factory, num_sims, true_prob_tox, ...){
    protocol <- u_i(selector_factory) # separate naming from implementation details
    stopifnot("num_sims must be of length 1 or 2" = length(num_sims) %in% 1:2)
    M <- num_sims[1]
    K <- num_sims[2]; if(is.na(K)) K <- num_sims[1]
    # The most parsimonious generalization of the default function
    # will substitute a *matrix* for the default result's vector
    # attribute 'true_prob_tox'.
    dose_levels <- getOption("dose_levels", default = stop(
      "simulate_trials methods require option(dose_levels)."))
    tpt_matrix <- tox_probs_at(true_prob_tox, dose_levels, K, ...)
    P_ <- seq_along(dose_levels)
    fits <- list()
    pb <- txtProgressBar(max=K, style=3)
    for(k in 1:K){
      fits <- c(fits,
                simulate_trials(protocol
                                ,num_sims = as.integer(M)
                                ,true_prob_tox = tpt_matrix[k, P_]
                )[[1]]
      )
      setTxtProgressBar(pb, k)
    }
    sims <- list(
      fits = fits
      , mean_prob_tox = colMeans(tpt_matrix[, P_])
      , true_prob_tox_matrix = tpt_matrix
      , hyper_mtdi_distribution = true_prob_tox
    )
    sims$dose_levels <- dose_levels
    sims$dose_units <- true_prob_tox@units
    class(sims) <- c("precautionary","simulations")
    return(sims)
  }
)

# Now another incremental generalization of escalation::simulate_trials,
# only this time in the direction of a greater PHARMACOLOGIC REALISM.
# The true_prob_tox argument of class "mtdi_distribution" specifies
# real dosing (including dose units), and deals directly with latent,
# individual-level characteristics of toxicity thresholds/sensitivity.
#
# What I may learn from this is that the (real) dose levels must be
# specified in package-specific options whenever mtdi_distributions
# are employed.

#' @examples
#' # TODO: In keeping with spirit of REALISM, attach UNITS to 'dose_levels'.
#' options(dose_levels = c(2, 6, 20, 60, 180, 400))
#' mtdi_dist <- mtdi_lognormal(CV = 0.5
#'                            ,median = 140
#'                            ,units = "ng/kg/week")
#' sims <- get_three_plus_three(num_doses = 6) %>%
#'   simulate_trials(num_sims = 50, true_prob_tox = mtdi_dist)
#' # Now attach a proper ordinalizer to 'mtdi_dist':
#' mtdi_dist@ordinalizer <- function(dose, r0) {
#'   c(Gr1=dose*r0^2, Gr2=dose*r0, Gr3=dose, Gr4=dose/r0, Gr5=dose/r0^2)
#' }
#' simsOT <- get_three_plus_three(num_doses = 6) %>%
#'   simulate_trials(num_sims = 50, true_prob_tox = mtdi_dist, r0 = 1.5)
#' summary(simsOT)
#' @rdname simulate_trials
setMethod(
  "simulate_trials"
  , c(selector_factory="selector_factory",
      num_sims="numeric",
      true_prob_tox="mtdi_distribution"),
  function(selector_factory, num_sims, true_prob_tox, ...){
    # TODO: Try moving this next line into the argument defaults
    dose_levels <- getOption("dose_levels", default = stop(
      "simulate_trials methods require option(dose_levels)."))
    tpt_vector <- tox_probs_at(true_prob_tox, dose_levels, ...)
    sims <- simulate_trials(selector_factory = u_i(selector_factory)
                            , num_sims = num_sims
                            , true_prob_tox = tpt_vector
    )
    sims$dose_levels <- dose_levels
    sims$dose_units <- true_prob_tox@units
    prependClass("precautionary", sims)
  }
)

# setMethod("simulate_trials", c(selector_factory="ordtox"),
#   function(selector_factory, num_sims, true_prob_tox, ...){
#     protocol <- selector_factory # separate naming from implementation details
#     lambda_CV <- 3
#     median_mtd <- num_doses(protocol) - 1
#     median_sd <- median_mtd/3
#     r0 <- seq(0.5, 2.5, 0.5)
#     K <- num_sims
#     M <- num_sims
#     # TODO: Consider making true_prob_tox argument optionally
#     #       a data table (or data frame) of the following form,
#     #       _OR_ (more usefully) a *generative model* for same.
#     true_prob_tox <- data.table(CV = rexp(n=K, rate=lambda_CV))
#     true_prob_tox[, `:=`(
#       sigma = sqrt(log(CV^2 + 1)) # NB: CV>sigma, w/ near-identity below 0.4
#       , mu = rnorm(n=K, mean=median_mtd, sd=median_sd)
#     )]
#     P_ <- paste0("P", dose_indices(protocol))
#     true_prob_tox[, c(P_) := lapply(dose_indices(protocol),
#                                     function(q) pnorm(q, mean=mu, sd=sigma))]
#     ensembles <<- list()
#     toxicities <<- list()
#     for(k in 1:K){
#       cat("k =", k, "\n")
#       # Now invoke the default method from package 'escalation':
#       sims <- callNextMethod(protocol, num_sims,
#                              true_prob_tox = as.numeric(true_prob_tox[k, ..P_]))
#       ensemble <- rbindlist(lapply(sims[[1]], function(.) .[[1]]$fit$outcomes)
#                             , idcol = "rep")
#       ensemble <- merge(data.frame(r0=r0), ensemble) # cartesian product
#       ensemble <- as.data.table(ensemble) # restore data.table after cartesian product
#       ensemble[, MTDi3 := qnorm(p = u_i
#                                 , mean = true_prob_tox[k]$mu
#                                 , sd = true_prob_tox[k]$sigma
#       )]
#       ensemble[, `:=`(
#         MTDi1 = MTDi3 - 2*r0
#         , MTDi2 = MTDi3 - r0
#         , MTDi4 = MTDi3 + r0
#         , MTDi5 = MTDi3 + 2*r0
#       )]
#       ensemble[, toxgrade := (dose>MTDi1) + (dose>MTDi2) + (dose>MTDi3) +
#                  (dose>MTDi4) + (dose>MTDi5)]
#       stopifnot(with(ensemble, all(tox == as.integer(dose >= MTDi3))))
#       # N.B.: Any trial realization ('rep') is just as likely as any other.
#       #       Therefore, I have to aggregate with reps as the denominator.
#       #       Specifically, pooling all the patients would fail to weight
#       #       them inversely to trial size.
#       #       Another way to see this is, we are asking about the probability
#       #       of a fatality in *this* trial. Alternatively, we can ask for
#       #       the expected number of each grade of toxicity in *this* trial.
#       #       Note that this should partition the expected size of this trial!
#       ensemble[, Tox := factor(toxgrade, levels=0:5, labels=paste0("Gr",0:5))]
#       ensembles[[k]] <<- ensemble
#       # See https://stackoverflow.com/a/16519612/3338147 explaining below syntax
#       toxicities[[k]] <<- ensemble[, .(Tox=ordered(levels(Tox)), N=c(table(Tox)))
#                                    , by=r0]
#       # Importantly, because these protocols run considering only DLT = Gr>=3,
#       # they are r0-agnostic. But I ought to allow for the general case where
#       # ordinal toxicities (MTDig for g != 3) affect trial conduct.
#       #
#       # Initially, let me go simply for *counts* and worry secondarily about
#       # obtaining (and dividing by) the denominators.
#       # Want a tabulation with several r0 values defining rows, and columns for
#       # the counts of toxicity grades. A final column may show total enrollments,
#       # which as noted above will all be equal unless ordinal toxicities affect
#       # escalation or termination decisions.
#     }
#     # Indeed, I like the idea of abstracting the calculation of high-level
#     # summary statistics into a separate function, possibly even a 'summary'
#     # method for a suitably defined S3 class. But for now, let me implement
#     # these summaries here, as additional components of the returned list.
#     toxdt <- rbindlist(toxicities, idcol = "k")
#     counts <- toxdt[, .(n=sum(N)), by=.(r0,Tox)]
#     toxtab <- ftable(xtabs(n ~ Tox + r0, counts))
#     # It now seems to me there are 2 perspectives on the probabilities here.
#     # One POV is the per-trial perspective, which might ask for expected
#     # numbers of each grade of toxicity.
#     # The second POV is that of the enrolling patient, who might ask what
#     # is the probability of each grade of toxicity. But this POV ought to
#     # be conditioned on the prior results seen in the trial, and so demands
#     # much more sophisticated modeling -- indeed, modeling of the kind I
#     # employed in the AFM11 paper.
#     # Note in fact that, as soon as such a perspective BEGINS to be acknowledged,
#     # I have already 'won' the argument about dose individualization!
#     list(tpt = true_prob_tox
#         ,toxdt = toxdt
#         ,toxtab = toxtab
#         ,expect = toxtab/(K*M)
#         ,K = K
#         ,M = M
#         )
#   }
# )

