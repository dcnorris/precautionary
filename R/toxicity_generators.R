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
# Note that the 'optional' nature of the @ordinalizer carries
# through into the treatment of ordinal toxicity information
# using an attr(sims,'ordtox') which may or may not be NULL.
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
#' @export
setMethod(
  "tox_probs_at"
  , c("mtdi_distribution", "numeric"),
  function(mtdi, doses, ...){
    probs <- mtdi@dist$cdf(doses)
    names(probs) <- paste0(doses, mtdi@units)
    if( !is.null(body(mtdi@ordinalizer)) ) {
      doses_matrix <- sapply(doses, mtdi@ordinalizer, ...)
      # Trick below comes via https://stackoverflow.com/a/8579498/3338147:
      probs_matrix <- doses_matrix # template structure
      probs_matrix[] <- vapply(doses_matrix, mtdi@dist$cdf, numeric(1))
      colnames(probs_matrix) <- paste0(doses, mtdi@units)
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

setGeneric("draw_samples", function(hyper, K, ...) {
  standardGeneric("draw_samples")
})

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

setMethod(
  "draw_samples"
  , c("hyper_mtdi_lognormal","numeric"),
  function(hyper, K=NULL, ...) {
    if(missing(K)) K <- 1
    mapply(mtdi_lognormal
          , CV = rexp(n=K, rate=hyper@lambda_CV)
          , median = rnorm(n=K, mean=hyper@median_mtd, sd=hyper@median_sd)
          , MoreArgs = list(units = hyper@units
                           ,ordinalizer = hyper@ordinalizer)
          , SIMPLIFY = FALSE)
  })

#' @examples 
#' mtdi_gen <- hyper_mtdi_lognormal(lambda_CV = 3
#'                                  , median_mtd = 6, median_sd = 2
#'                                  , units="mg/kg")
#' doses <- c(0.5, 1, 2, 4, 6, 8)
#' mtdi_gen %>% tox_probs_at(doses, K=3)
#' # Now attach a proper ordinalizer to 'mtdi_gen':
#' mtdi_gen@ordinalizer <- function(dose, r0) {
#'   c(Gr1=dose*r0^2, Gr2=dose*r0, Gr3=dose, Gr4=dose/r0, Gr5=dose/r0^2)
#' }
#' tpa <- mtdi_gen %>% tox_probs_at(doses, K=3, r0=1.5)
#' tpa
#' ordtox <- attr(tpa,'ordtox')
#' class(ordtox)
#' dimnames(ordtox)
#' @rdname tox_probs_at
#' @export
setMethod(
  "tox_probs_at"
  , c("hyper_mtdi_lognormal","numeric"),
  function(mtdi, doses, K, ...) {
    # Generate a list of K mtdi_lognormal objects:
    mtdi_list <- draw_samples(hyper = mtdi, K = K)
    ## TEST whether the elements of mtdi_list have 'ordtox' attribute
    cat("attr(mtdi_list[[1]],'ordtox') = ")
    print(attr(mtdi_list[[1]],'ordtox'))
    # Transform this list to a list of tox_probs_at results:
    tpa_list <- lapply(mtdi_list, tox_probs_at, doses, ...)
    # Generate a higher-dimensional version of the result
    # yielded by tox_probs_at("mtdi_distribution") method.
    # This will be a K-row matrix of sample tox probabilities,
    # with an optional 'ordtox' attribute that is an *array*.
    result <- do.call(rbind, tpa_list)
    ordtox <- sapply(tpa_list, attr, which='ordtox', simplify="array")
    if(all(sapply(ordtox, is.null))) ordtox <- NULL
    # TODO: Consider naming the 3rd dimension or ordtox
    #       with labels like 'sample1', 'sample2', ...
    #       or perhaps "k=1", "k=2", ... in case I wish
    #       to settle on 'k' as a sample index.
    attr(result,'ordtox') <- ordtox
    result
  })

