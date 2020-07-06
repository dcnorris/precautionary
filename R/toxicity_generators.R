# Generative models of binary and ordinal toxicities

# These generative models will abstract from 'true_tox_probs' in several dimensions:
# 1. There is a dimension of *uncertainty*, itself comprising 2 types:
# (a) There is a 'hyperparameter' dimension of (subjective, 'epistemic') uncertainty
#     about the population-level distribution of toxicity responses.
# (b) There is *aleatoric* ('purely statistical') uncertainty arising from the
#     random sampling that study enrollment amounts to.
# 2. There are questions of *pharmacologic realism*, including:
# (a) Deriving binary toxicities from latent u_i (equiv., MTDi)
# (b) Explicit dimensioning of doses
# (c) The scaling of the pharmacological effects of these doses
# (d) Elaboration of pharmacologic effects (e.g., *graded* toxicities)
#
# I have dealt with (1b) in isolation already, by generalizing 'simulate_trials'
# to the case where 'true_tox_probs' is (effectively) a matrix.
#
# A next step may be to handle (2a,b) *together*, in a class which allows for
# specification of an MTDi distribution over real doses.
# In constructing this class, I may wish to anticipate the need to abstract
# the parameter specification to facilitate (1a).
#
#' @importFrom distr6 SDistribution Lognormal
setOldClass(c("Lognormal","SDistribution"))

# The @ordinalizer slot is an optional 'mixin' defaulting to
# the absent case of a NULL body(). This function in general
# is a dose-space transformation that spreads any given dose
# into an ordered vector of doses at which graded toxicities
# occur.
# One appealing manner of specifying an @ordinalizer would be
# to take an MTDi3 as input dose, and output an ascending
# vector of doses c(MTDi1=, MTDi2=, ..., MTDi5=). But using
# this scheme would impose on my code the obligation to invert
# this relation.
# So instead, I will implement @ordinalizers directly through
# the inverse function, one that transforms a given dose into
# a *descending* vector of doses that yield equivalent rates
# of toxicity from the mtdi_generator.
# If I am using 5 toxicity grades (Gr1,..Gr5) and if we take
# MTDi to be the Gr3 threshold, then the @ordinalizer's task
# is to transform an given dose into Gr1..5 *equivalents*
# vis-à-vis the mtdi_generator. This is really so abstract
# that it may need to be hidden from users; consider.
# I will gain little from imposing premature concreteness,
# and indeed may obscure some essential aspects with excessive
# hand-holding. The user should perhaps be forced to specify
# the toxicity grades explicitly!
# Let me keep things simple by assuming this function is of
# a single dose, i.e., is NOT vectorized.
setClass("mtdi_generator",
  slots = list(
    units = "character"
  , ordinalizer = "function" # optional 'mixin', present if !is.null(body(.))
  ),
  contains = "VIRTUAL")

setClass("mtdi_distribution",
  slots = list(
    CV = "numeric"      # Establish CV and median as universal parameters
  , median = "numeric"  # for specifying MTDi distributions.
  , dist = "SDistribution"
  ),
  contains = c("mtdi_generator","VIRTUAL")
)

# This function may offer excellent information-hiding access
# to the optonal @ordinalizer of its first argument.
# The return value may simply have an optional 'ordtox' attribute
# that mirrors the optionality of the @ordinalizer. After invoking
# 'tox_probs_at', client code can check whether this attribute
# exists, and act accordingly.
# The ... argument provides an obvious means to pass further
# information like r0 (or r[1,2,4,5]) to the @ordinalizer.

#' Obtain toxicity probabilities from an mtdi_generator
#'
#' @param mtdi An \code{mtdi_generator}
#' @param doses A vector of (dimensioned) doses at which to calculate
#' the toxicity probabilities.
#' @param ... Arguments passed to \code{mtdi}'s ordinalizer
#'
#' @return A suitable data structure (e.g., vector or matrix) of
#' binary-toxicity probabilities, with an optional \sQuote{ordtox}
#' attribute (in case \code{mtdi} has a non-NULL ordinalizer slot)
#' that spreads this into the further dimension of multiple ordinal
#' toxicities.
#' @export
setGeneric("tox_probs_at", function(mtdi, doses, ...) {
  standardGeneric("tox_probs_at")
})

#' @examples 
#' mtdi_gen <- mtdi_lognormal(CV = 0.5
#'                           ,median = 140
#'                           ,units = "ng/kg/week")
#' doses <- seq(100, 150, 10)
#' mtdi_gen %>% tox_probs_at(doses)
#' # Now attach a proper ordinalizer to 'mtdi_gen':
#' mtdi_gen@ordinalizer <- function(dose, r0) {
#'   c(Gr1=dose*r0^2, Gr2=dose*r0, Gr3=dose, Gr4=dose/r0, Gr5=dose/r0^2)
#' }
#' tpa <- mtdi_gen %>% tox_probs_at(doses, r0=1.5)
#' tpa
#' @rdname tox_probs_at
setMethod(
  "tox_probs_at"
  , c("mtdi_distribution", "numeric"),
  function(mtdi, doses, ...){
    probs <- mtdi@dist$cdf(doses)
    names(probs) <- paste0("Ptox(", doses, ")")
    if( !is.null(body(mtdi@ordinalizer)) ) {
      doses_matrix <- sapply(doses, mtdi@ordinalizer, ...)
      print(doses_matrix)
      probs_matrix <- mtdi@dist$cdf(doses_matrix)
      # restore matrix dims flattened by cdf()
      probs_matrix <- matrix(probs_matrix, ncol=length(doses))
      rownames(probs_matrix) <- rownames(doses_matrix)
      colnames(probs_matrix) <- paste0("Ptox(", doses, ")")
      attr(probs,'ordtox') <- probs_matrix
    }
    probs
  })

#' @export mtdi_lognormal
mtdi_lognormal <- setClass("mtdi_lognormal",
  slots = list(
    dist = "Lognormal"
  ),
  contains = "mtdi_distribution"
)

setMethod("initialize", "mtdi_lognormal",
  function(.Object, CV, median, ...) {
    .Object@dist <- Lognormal$new(meanlog = log(median),
                                  sdlog = sqrt(log(CV^2 + 1)))
    # Initialize 'bottom-up' to avoid a premature validity check:
    .Object <- callNextMethod(.Object, CV=CV, median=median, ...)
  })

setClass("hyper_mtdi_distribution",
  contains = c("mtdi_generator","VIRTUAL")
  )

#' @export hyper_mtdi_lognormal
hyper_mtdi_lognormal <- setClass("hyper_mtdi_lognormal",
  slots = list( # hyperparameters
    lambda_CV = "numeric"
  , median_mtd = "numeric"
  , median_sd = "numeric" # TODO: Cosider whether an argument can be made
                          #       for dispensing with this 2nd uncertainty
                          #       measure. Or is it CV that I might argue
                          #       is superfluous, on grounds that we may
                          #       (intuitively) impute our uncertainty in
                          #       median_mtd to the population variance?
  ),
  contains = "hyper_mtdi_distribution"
  )

setMethod("initialize", "hyper_mtdi_lognormal",
  function(.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
  })

#' @export
setMethod(
  "tox_probs_at"
  , c("hyper_mtdi_lognormal","numeric"),
  function(mtdi, doses, K) {
    if(missing(K)) K <- 1
    # Generate a K-row matrix of sample tox probabilities
    CV <- rexp(n=K, rate=mtdi@lambda_CV)
    sdlog <- sqrt(log(CV^2 + 1))
    meanlog <- rnorm(n=K, mean=mtdi@median_mtd, sd=mtdi@median_sd)
    tox_probs <- matrix(nrow=K, ncol=length(doses)+3)
    colnames(tox_probs) <- c(paste0("P",1:length(doses)), "CV", "meanlog", "sdlog")
    tox_probs[, paste0("P", seq_along(doses))] <- sapply(doses, function(q) pnorm(q, mean=meanlog, sd=sdlog))
    tox_probs[,"CV"] <- CV
    tox_probs[,"meanlog"] <- meanlog
    tox_probs[,"sdlog"] <- sdlog
    tox_probs
  })

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
#' sims <- get_three_plus_three(num_doses = 6) %>%
#'   simulate_trials(num_sims = c(30, 10), true_prob_tox = mtdi_gen)
#' summary(sims)
#' @rdname simulate_trials
setMethod(
  "simulate_trials"
  , c(selector_factory="selector_factory",
      num_sims="numeric",
      true_prob_tox="hyper_mtdi_distribution"),
  function(selector_factory, num_sims, true_prob_tox, ...){
    cat('simulate_trials(true_prob_tox="hyper_mtdi_distribution") method ...\n')
    protocol <- selector_factory # separate naming from implementation details
    stopifnot("num_sims must be of length 1 or 2" = length(num_sims) %in% 1:2)
    M <- num_sims[1]
    K <- num_sims[2]; if(is.na(K)) K <- num_sims[1]
    # The most parsimonious generalization of the default function
    # will substitute a *matrix* for the default result's vector
    # attribute 'true_prob_tox'.
    dose_levels <- getOption("dose_levels"
                             , default = stop("hyper_mtdi_distribution methods require option(dose_levels)."))
    tpt_matrix <- tox_probs_at(true_prob_tox, dose_levels, K)
    P_ <- paste0("P", dose_indices(protocol))
    fits <- list()
    pb <- txtProgressBar(max=K, style=3)
    for(k in 1:K){
      fits <- c(fits,
                simulate_trials(protocol
                                ,num_sims = as.integer(M)
                                ,true_prob_tox = tpt_matrix[k, P_]
                                , ...)[[1]]
                )
      setTxtProgressBar(pb, k)
    }
    sims <- list(
      fits = fits
    , true_prob_tox = colMeans(tpt_matrix[, P_])
    , true_prob_tox_matrix = tpt_matrix
    , hyper_mtdi_distribution = true_prob_tox
    , WARNING = paste("The 'true_prob_tox' component gives the mean of",
                         K, "samples drawn from the hyper_mtdi_distribution.")
    )
    sims$dose_levels <- dose_levels
    sims$dose_units <- true_prob_tox@units
    class(sims) <- c("realdose_simulations","simulations")
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
#' sims <- get_three_plus_three(num_doses = 6) %>%
#'   simulate_trials(num_sims = 50
#'                  ,true_prob_tox = mtdi_lognormal(CV = 0.5
#'                                                 ,median = 140
#'                                                 ,units = "ng/kg/week"))
#' summary(sims)
#' @rdname simulate_trials
setMethod(
  "simulate_trials"
  , c(selector_factory="selector_factory",
      num_sims="numeric",
      true_prob_tox="mtdi_distribution"),
  function(selector_factory, num_sims, true_prob_tox, ...){
    cat('simulate_trials(true_prob_tox="mtdi_distribution") method ...\n')
    # TODO: Try moving this next line into the argument defaults
    dose_levels <- getOption("dose_levels"
                             , default = stop("mtdi_distribution methods require option(dose_levels)."))
    sims <- simulate_trials(selector_factory = selector_factory
                           , num_sims = num_sims
                           , true_prob_tox = tox_probs_at(true_prob_tox, dose_levels)
                           , ...)
    sims$dose_levels <- dose_levels
    sims$dose_units <- true_prob_tox@units
    class(sims) <- c("realdose_simulations", class(sims))
    return(sims)
  }
)

#' @importFrom dplyr rename rename_with mutate select
#' @export
summary.realdose_simulations <- function(x, ...) {
  summary <- NextMethod()
  dose_units <- paste0("dose (", x$dose_units, ")")
  summary %>%
    mutate("real_dose" = c(0, x$dose_levels)[as.integer(summary$dose)]) %>%
    select(dose, real_dose, everything()) %>%
    rename_with(.fn = function(.) dose_units, .cols = real_dose)
}

