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

setGeneric("plot")

#' @export
setMethod("plot", "mtdi_distribution", function(x, y=NULL, ...) {
  xlab <- paste0("Dose (", x@units, ")")
  ylab <- "CDF"
  title <- paste(class(x@dist)[1], "MTDi Distribution")
  params <- paste0("CV = ", x@CV, "; median = ", x@median, x@units)
  # I will presume most pharmacology should be done in log-dose space...
  CDFs <- seq(0.01, 0.99, 0.01)
  quantiles <- x@dist$quantile(CDFs)
  plot.default(CDFs ~ quantiles, type="l", log="x"
               , xlab = xlab, ylab = ylab, main = title
               , sub = params, font.sub = 3
               , las = 1
               , lab = c(x=20, y=6, len=3)
               )
  # Obtain a reasonable default calculation for the minor tick locations.
  qrange <- range(quantiles)
  minor_step <- 10^floor(mean(log10(qrange)))
  steprange <- qrange/minor_step
  minor_ticks <- minor_step * ceiling(steprange[1]):floor(steprange[2])
  axis(1, at=minor_ticks, tcl=-0.3, labels=NA) # minor ticks
  dose_levels <- getOption('dose_levels')
  if( !is.null(dose_levels) ){
    tox_probs <- x@dist$cdf(dose_levels)
    for( i in seq_along(dose_levels) ) {
      lines(x = c(0.1, rep(dose_levels[i], 2)),
            y = c(rep(tox_probs[i], 2), -1),
            lty = 3)
    }
  }
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
  , median_sd = "numeric" # TODO: Consider whether an argument can be made
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

