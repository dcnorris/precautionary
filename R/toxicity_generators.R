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
setOldClass("SDistribution")
setOldClass("Lognormal")

setClass("mtdi_distribution",
  slots = list(
    CV = "numeric"      # Establish CV and median as universal parameters
  , median = "numeric"  # for specifying MTDi distributions.
  , units = "character"
  , dist = "ANY"
  #, dist = "SDistribution"
  # I would prefer to specify @dist class as "SDistribution",
  # but reconcilePropertiesAndPrototype() called upon loading
  # this package sees this as a "conflict", as if it does not
  # recognize that "Lognormal" inherits from "SDistribution".
  ),
  contains = "VIRTUAL"
)

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

# Initially, let me focus on generating MTDi distributions.
# Because I aim to induce REALISTIC DOSE-RESPONSE THINKING,
# this class will actually demand an explicit dose scale.
# I will also PRESUME from the outset that pharmacologic
# scaling of doses is always LOGARITHMIC. Certainly, any
# relaxation of that assumption can await later versions!
# Furthermore, NO CONCESSIONS can be made here to any fixed
# number or set of 'dose levels'. That is an abstraction
# belonging to dose-escalation trials, NOT to realistic
# thinking about pharmacology.
# One essential requirement for promoting realistic thinking
# is SUPPORTING it with good visualization! Thus, among the
# first functions I should implement are some good graphics
# showing the implications of any given mtdi_generator.
# By assuming logarithmic scaling of the dimensioned doses,
# I will be able to define these graphics without wrestling
# with multiple alternative (or indefinite) scalings.
setClass("mtdi_generator",
  slots = list(doses = "numeric" # TODO: Delegate knowledge of doses to client code?
              ),
  contains = "VIRTUAL"
  )


setGeneric("sample_tox_probs", function(this, K) {
  standardGeneric("sample_tox_probs")
})

# TODO: Maybe this class is too generically named. The key feature,
#       once I've worked out issues of number of parameters, indeed
#       may be that very number. Thus, maybe I have here a class
#       better named (provisionally) "mtdi_lognormal_3hyperparam".
#       Furthermore, the precise details of the argument may rather
#       correspond to a certain paper.
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
  contains = "mtdi_generator"
  )

setMethod("initialize", "hyper_mtdi_lognormal",
  function(.Object, ...) {
    .Object <- callNextMethod(.Object
                              , doses = 1:6 # TODO: Don't hard-code this!
                              , ...)
  })

#' @export
setMethod("sample_tox_probs", signature("hyper_mtdi_lognormal","numeric"),
  function(this, K) {
    if(missing(K)) K <- 1
    # Generate a K-row matrix of sample tox probabilities
    CV <- rexp(n=K, rate=this@lambda_CV)
    sdlog <- sqrt(log(CV^2 + 1))
    meanlog <- rnorm(n=K, mean=this@median_mtd, sd=this@median_sd)
    tox_probs <- matrix(nrow=K, ncol=length(this@doses)+3)
    colnames(tox_probs) <- c(paste0("P",1:length(this@doses)), "CV", "meanlog", "sdlog")
    tox_probs[, paste0("P", seq_along(this@doses))] <- sapply(this@doses, function(q) pnorm(q, mean=meanlog, sd=sdlog))
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
#' mtdi_gen <- hyper_mtdi_lognormal(lambda_CV = 3, median_mtd = 4, median_sd = 1)
#' sims <- get_three_plus_three(num_doses = 6) %>%
#'   simulate_trials(num_sims = c(30, 10), true_prob_tox = mtdi_gen)
#' summary(sims)
#' @rdname simulate_trials
setMethod(
  "simulate_trials"
  , c(selector_factory="three_plus_three_selector_factory",
      # Note that attempts to generalize the class of 'selector_factory' argument
      # result in dispatch of the non-generic default in package escalation.
      # The section 'Specialized Local Generics' in documentation for setGeneric
      # suggests to me that my attempts to 'rewrite' the escalation package from
      # the outside may strain the design intent of S4.
      # In theory, I could repeat this setMethod call for each concrete class of
      # selector_factory, but the main intention here has been to illuminate the
      # connection with escalation::simulate_trials; that aim is already achieved
      # in principle by this method, notwithstanding its overly-specific signature.
      num_sims="numeric",
      true_prob_tox="mtdi_generator"),
  function(selector_factory, num_sims, true_prob_tox, ...){
    cat('simulate_trials(true_prob_tox="mtdi_generator") method ...\n')
    protocol <- selector_factory # separate naming from implementation details
    stopifnot("num_sims must be of length 1 or 2" = length(num_sims) %in% 1:2)
    M <- num_sims[1]
    K <- num_sims[2]; if(is.na(K)) K <- num_sims[1]
    # The most parsimonious generalization of the default function
    # will substitute a *matrix* for the default result's vector
    # attribute 'true_prob_tox'.
    tpt_matrix <- sample_tox_probs(true_prob_tox, K)
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
    , mtdi_generator = true_prob_tox
    , WARNING = paste("The 'true_prob_tox' component gives the mean of",
                         K, "samples drawn from the mtdi_generator.")
    )
    class(sims) <- "simulations"
    return(sims)
  }
)


