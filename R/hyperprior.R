## An R6 class for the MTDi *hyperprior* of the NEW & IMPROVED  EscRisk app
## (and other simulations) build around Complete Path Enumeration (CPE).

## TODO: Factor this class into a class hierarchy that clearly exhibits
##       where the choice of hyperprior distribution does and does not
##       impact the implementation.

#' @name HyperMTDi-class
#' @title An R6 base class for hyperpriors over MTDi distributions
#'
#' @details
#' With CPE liberating `precautionary` from the need for *nested* simulations,
#' the opportunity arises for a more encapsulated treatment of MTDi scenario
#' generators ('hyperpriors') and their sampling. Specifically, individual
#' sampled scenarios need only yield two functions:
#' * A CDF taking a dose vector X = (X1,...,Xd) to probabilities p = (p1,...,pd)
#' * A function F(X,kappa) yielding a fractionation matrix $F = \[0 G; H 0\]$.
#' Crucially, this class actually represents an APPROXIMATION to the hyperprior,
#' via a finite set of samples from it. The use of reference classes enables us
#' to improve this approximation efficiently by in-place updating.
#' @importFrom R6 R6Class
#' @importFrom stats qlnorm plnorm
#' @export
HyperMTDi_lognormal <- R6Class(
  "HyperMTDi_lognormal",
  public = list(
    #' @details
    #' Create a new `HyperMTDi` object.
    #' @note
    #' This class implements a finite approximation to the infinite
    #' set of MTDi scenarios which it describes---an approximation
    #' which may be improved dynamically by expanding the samples.
    #'
    #' @param CV Coefficient of variation of median MTDi
    #' @param median_mtd Median MTDi in the population
    #' @param median_sdlog Uncertainty in median MTDi, on log scale
    #' @param units A short string specifying dose units
    #' @param n Number of samples to draw
    #' @return A `HyperMTDi` object.
    initialize = function(CV, median_mtd, median_sdlog, units, n=100) {
      private$CV <- CV
      private$median_mtd <- median_mtd
      private$median_sdlog <- median_sdlog
      private$units <- units
      private$samples <- data.table(CV = numeric(0), median = numeric(0))
      self$extend(n)
    },
    #' @details
    #' Extend the samples, typically improving the approximation
    #' TODO: Investigate how much variance reduction QRNG yields.
    #' @param n Number of additional MTDi scenarios to sample
    #' @return Self, invisibly
    extend = function(n=1) {
      more <- data.table(
        CV = sqrt(rchisq(n=n, df=2))*private$CV # Rayleigh(1) is chi dist w/ 2 d.f.
      , median = rlnorm(n = n
                      , meanlog = log(private$median_mtd)
                      , sdlog = private$median_sdlog
                        )
      )
      private$samples <- rbind(private$samples, more)
      invisible(self)
    },
    #' @details
    #' Apply a distribution-type function over the sampled realizations
    #' TODO: Consider taking this method private.
    #'
    #' @param f A closure that realizes a distribution-type function (such as
    #' a quantile function or CDF) when evaluated in the environment defined
    #' by any row of the sampled parameters.
    #' @param ... Arguments upon which to evaluate the enclosed function
    #' @return A list of values of f
    apply = function(f, ...) {
      lapply(seq(nrow(private$samples))
           , function(j)
             with(private$samples[j,], f(CV, median)(...)))
    },
    #' @details
    #' Visualize the samples of a `HyperMTDi` object
    #'
    #' @param col Color of lines used to depict samples
    #' @param ... Additional arguments passed onward to `plot`
    #'
    #' @importFrom graphics abline axis lines mtext par plot.default
    #' @examples
    #' if (interactive()) {
    #' mtdi_gen <- HyperMTDi_lognormal$new(CV = 1
    #'                                    ,median_mtd = 5
    #'                                    ,median_sdlog = 0.5
    #'                                    ,units="mg/kg")
    #' mtdi_gen$plot()
    #' }
    #' @export
    plot = function(col="gray", ...) {
      n <- nrow(private$samples)
      xlab <- paste0("Dose (", private$units, ")")
      ylab <- "CDF"
      params <- paste0("CV ~ Raleigh(mode=", private$CV
                     , ");  median = ", private$median_mtd
                     , private$units, " \u00b1 " # <-- plus/minus character
                     , 100*private$median_sdlog, "%")
      title <- paste(n, private$dist$name
                   , "MTDi distributions sampled from hyperprior")
      CDFs <- seq(0.01, 0.99, 0.01)
      quantiles <- self$apply(private$dist$quantile, CDFs)
      xlim <- range(do.call(c, quantiles))
      oldpar <- par(mar = c(5,5,4,1) + 0.1)
      plot.default(CDFs ~ quantiles[[1]], type="l", log="x"
                 , xlab = xlab, ylab = "", main = title
                 , xlim = xlim
                 , sub = params, font.sub = 3
                 , las = 1
                 , lab = c(x=20, y=6, len=3)
                 , col = col
                 , ...
                   )
      mtext(ylab, side = 2, line = 2.5, las = 1) # horiz. axis title as per Tufte
      ## Locate old-fashioned, decade-wise logarithmic axis minor ticks
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
      par(oldpar)
    }
  ), # </public>
  private = list(
    dist = list(name = "Lognormal"
              , quantile = function(CV, median) function(p)
                qlnorm(p, meanlog = log(median), sdlog = sqrt(log(CV^2+1)))
              , cdf = function(CV, median) function(q)
                plnorm(q, meanlog = log(median), sdlog=sqrt(log(CV^2+1)))
                )
  , CV = NA
  , median_mtd = NA
  , median_sdlog = NA
  , units = NA
  , samples = NULL
  )
) # </HyperMTDi_lognormal>
