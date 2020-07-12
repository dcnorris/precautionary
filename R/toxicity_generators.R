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
#' @importFrom distr6 SDistribution Lognormal Rayleigh
setOldClass(c("Lognormal","SDistribution"))

setClass("mtdi_generator",
  slots = list(
    units = "character"
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
  # Locate old-fashioned, decade-wise logarithmic axis minor ticks
  erange <- floor(log10(range(quantiles)))
  exponents <- erange[1]:erange[2]
  minor_ticks <- as.vector(outer(2:9, exponents, function(x,y) x*10^y))
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
    CV = "numeric"
  , median_mtd = "numeric"
  , median_sdlog = "numeric"
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
          , CV = Rayleigh$new(mode = hyper@CV)$rand(n=K, simplify=TRUE)
          , median = rlnorm(n=K
                            , meanlog=log(hyper@median_mtd)
                            , sdlog=hyper@median_sdlog)
          , MoreArgs = list(units = hyper@units
                           )
          , SIMPLIFY = FALSE)
  })

