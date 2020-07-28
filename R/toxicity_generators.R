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
  , dist = "list" #"SDistribution"
  ),
  contains = c("mtdi_generator","VIRTUAL")
)

setGeneric("plot")

#' Visualize n samples from an \code{mtdi_generator} object
#'
#' @param x An \code{mtdi_generator} object
#' @param y Included for compatibility with generic signature
#' @param n Number of samples to draw from hyperprior for visualization
#' @param col Color of lines used to depict samples
#' @param \dots Additional arguments passed onward to \code{plot}
#'
#' @importFrom graphics abline axis lines plot.default
#' @examples
#' if (interactive()) {
#' mtdi_gen <- hyper_mtdi_lognormal(CV = 1
#'                                 ,median_mtd = 5
#'                                 ,median_sdlog = 0.5
#'                                 ,units="mg/kg")
#' plot(mtdi_gen, n=100, col=adjustcolor("red", alpha=0.5))
#' } 
#' @export
setMethod("plot", "mtdi_generator", function(x, y=NULL, n=20, col="gray", ...) {
  xlab <- paste0("Dose (", x@units, ")")
  ylab <- "CDF"
  params <- paste0("CV ~ Raleigh(mode=", x@CV, ");  median = "
                   , x@median_mtd, x@units, " \u00b1 " # <-- plus/minus character
                   , 100*x@median_sdlog, "%")
  mtdi_samples <- draw_samples(x, n = n)
  title <- paste(n, mtdi_samples[[1]]@dist$name
                 , "MTDi distributions sampled from hyperprior")
  CDFs <- seq(0.01, 0.99, 0.01)
  quantiles <- lapply(mtdi_samples
                      , function(mtdi) mtdi@dist$quantile(CDFs))
  xlim <- range(do.call(c, quantiles))
  plot.default(CDFs ~ quantiles[[1]], type="l", log="x"
               , xlab = xlab, ylab = ylab, main = title
               , xlim = xlim
               , sub = params, font.sub = 3
               , las = 1
               , lab = c(x=20, y=6, len=3)
               , col = col
               , ...
  )
  # Locate old-fashioned, decade-wise logarithmic axis minor ticks
  erange <- floor(log10(range(quantiles)))
  exponents <- erange[1]:erange[2]
  minor_ticks <- as.vector(outer(2:9, exponents, function(x,y) x*10^y))
  axis(1, at=minor_ticks, tcl=-0.3, labels=NA) # minor ticks
  dose_levels <- getOption('dose_levels')
  for(i in 2:n){
    lines(CDFs ~ quantiles[[i]], col = col, ...)
  }
  if( !is.null(dose_levels) ){
    abline(v = dose_levels, lty = 3)
  }
})

#' Visualize an \code{mtdi_distribution} object
#'
#' @param x An \code{mtdi_distribution} object 
#' @param y Included for compatibility with generic signature
#' @param \dots Additional arguments passed onward to \code{plot}
#'
#' @examples 
#' if (interactive()) {
#' mtdi_dist <- mtdi_lognormal(CV = 2
#'                            ,median = 5
#'                            ,units = "mg/kg")
#' # Setting pre-specified dose levels via options() causes
#' # toxicity probabilities to be annotated on the plot.
#' old <- options(dose_levels = c(0.5, 1, 2, 4, 6))
#' plot(mtdi_dist, col = "red")
#' options(old)
#' }
#' @export
setMethod("plot", "mtdi_distribution", function(x, y=NULL, ...) {
  xlab <- paste0("Dose (", x@units, ")")
  ylab <- "CDF"
  title <- paste(x@dist$name, "MTDi Distribution")
  params <- paste0("CV = ", x@CV, "; median = ", x@median, x@units)
  # I will presume most pharmacology should be done in log-dose space...
  CDFs <- seq(0.01, 0.99, 0.01)
  quantiles <- x@dist$quantile(CDFs)
  plot.default(CDFs ~ quantiles, type="l", log="x"
               , xlab = xlab, ylab = ylab, main = title
               , sub = params, font.sub = 3
               , las = 1
               , lab = c(x=20, y=6, len=3)
               , ...
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


#' A lognormal MTDi distribution
#' 
#' @slot dist An object of class \code{distr6::Lognormal}
#'
#' @export mtdi_lognormal
mtdi_lognormal <- setClass("mtdi_lognormal",
  slots = list(
    dist = "list" #"Lognormal"
  ),
  contains = "mtdi_distribution"
)

#' @importFrom stats qlnorm plnorm
setMethod("initialize", "mtdi_lognormal",
  function(.Object, CV, median, ...) {
    # .Object@dist <- Lognormal$new(meanlog = log(median),
    #                               sdlog = sqrt(log(CV^2 + 1)))
    meanlog <- log(median)
    sdlog <- sqrt(log(CV^2 + 1))
    .Object@dist <- list(
      quantile = function(p) qlnorm(p, meanlog = meanlog, sdlog = sdlog)
    , cdf = function(q) plnorm(q, meanlog = meanlog, sdlog = sdlog)
    , name = "Lognormal"
    )
    # Initialize 'bottom-up' to avoid a premature validity check:
    .Object <- callNextMethod(.Object, CV=CV, median=median, ...)
  })

setClass("hyper_mtdi_distribution",
  contains = c("mtdi_generator","VIRTUAL")
  )

setGeneric("draw_samples", function(hyper, n, ...) {
  standardGeneric("draw_samples")
})

#' Hyperprior for lognormal MTDi distributions
#' 
#' This hyperprior generates lognormal MTDi distributions with their
#' coefficient of variation being drawn from a Rayleigh distribution
#' with mode parameter \code{CV}. Because the standard deviation of
#' this distribution is
#' 
#' this implicitly link our uncertainty about \code{CV} to its value.
#' 
#' The medians of the lognormal distributions generated are themselves
#' drawn from a lognormal distribution with \code{meanlog = log(median_mtd)}
#' and \code{sdlog = median_sdlog}. Thus, parameter \code{median_sdlog}
#' represents a proportional uncertainty in \code{median_mtd}.
#'
#' @slot CV Coefficient of variation
#' @slot median_mtd Median MTDi
#' @slot median_sdlog Proportional uncertainty in median MTDi
#'
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

#' @importFrom stats rlnorm
setMethod(
  "draw_samples"
  , c("hyper_mtdi_lognormal","numeric"),
  function(hyper, n=NULL, ...) {
    if(missing(n)) n <- 1
    mapply(mtdi_lognormal
          , CV = Rayleigh$new(mode = hyper@CV)$rand(n=n, simplify=TRUE)
          , median = rlnorm(n=n
                            , meanlog=log(hyper@median_mtd)
                            , sdlog=hyper@median_sdlog)
          , MoreArgs = list(units = hyper@units
                           )
          , SIMPLIFY = FALSE)
  })

#' Test performance of \code{draw_samples} function
#' 
#' @param n Number of samples to draw
#'
#' @export
test_draw_samples <- function(n = 500){
  mtdi_gen <- hyper_mtdi_lognormal(CV = 1
                                   , median_mtd = 6, median_sdlog = 0.5
                                   , units="mg/kg")
  invisible(draw_samples(mtdi_gen, n = n))
}
