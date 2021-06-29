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
      self$resample(n)
    },
    #' @details
    #' Draw fresh samples
    #' @param n Number of samples to draw
    #' @return Self, invisibly
    resample = function(n) {
      if (missing(n))
        n <- nrow(private$samples)
      private$samples <- data.table(CV = numeric(0), median = numeric(0))
      self$extend(n)
      invisible(self)
    },
    #' @details
    #' Get number of samples
    #' TODO: Consider a higher-level interface to progress-bar info
    #' @return Number of samples drawn so far
    nsamples = function() {
      nrow(private$samples)
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
    #' Set or query the vector of pre-specified dose levels
    #'
    #' @param x A vector of dose levels
    #' @return Self (invisibly), unless `x` is missing,
    #' in which case the dose vector is returned.
    doses = function(x) {
      if (missing(x))
        return(private$.doses)
      private$.doses <- x
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
    #' Get average toxicity probabilities over the sample
    #' @return Toxicity probabilites at fixed doses, averaged over sample
    avg_tox_probs = function() {
      self$apply(function(CV, median)
        function() {
          private$dist$cdf(CV, median)(private$.doses)
        }) -> probs
      colMeans(do.call('rbind', probs))
    },
    #' @details
    #' Return expected counts of ordinal toxicities
    #'
    #' @param b The pathwise (length-J) vector defined by Eq (4) of Norris2020c
    #' @param U The J*2D matrix defined by Eq (4) of Norris2020c
    #' @param kappa A log-therapeutic index as in Eq (5) of Norris2020c
    #' @return A 6-column matrix, each row being the expected counts for toxicity
    #' grades 0 through 5, at one sampled scenario.
    #' @seealso Documentation for `Cpe-class`
    fractionate = function(b, U, kappa) {
      self$apply(function(CV, median)
        function(b, U, kappa) {
          ## This function should compute pi %*% U %*% F at each sampled pair (CV,median)!
          p <- private$dist$cdf(CV, median)(private$.doses)
          q <- 1 - p
          log_pq <- pmax(c(log(p), log(q)), .Machine$double.min.exp) # avoid -Inf from log(0.0)
          LOG_PQ <<- log_pq
          stopifnot(all(is.finite(log_pq)))
          log_pi <- b + U %*% log_pq
          ## Now I need to construct the fractionation matrix F ...
          ## 1. Obtain a 5*D matrix of shifted doses, marking off MTDi thresholds
          ##    for toxicities of Gr1+, Gr2+, Gr3+, Gr4+ and Gr5, which we will
          ##    denote through strict inequalities as STRICTLY WORSE THAN Gr0..Gr4.
          ##    (Note that this would even be in keeping with the left-continuity
          ##    required to preserve the 'MTD' intuition; we are looking for maximum
          ##    doses that still produce a given grade of toxicity, with dose+epsilon
          ##    resulting in the next-higher grade.)
          gradescale <- c(Gr0=2, Gr1=1, Gr2=0, Gr3=-1, Gr4=-2)
          X <- outer(exp(gradescale*kappa), private$.doses)
          ## 2. Evaluate the COMPLEMENTARY CDF on this matrix, to obtain INCREASING
          ##    vectors with probability of Gr0, Gr1 or less, Gr2 or less, etc.
          F <- 1 - private$dist$cdf(CV, median)(X)
          ## 3. 'Cap' these probability vectors with [0,1] endpoints,
          ##    then difference the columns to obtain a 6*D matrix
          F <- diff(rbind(0, F, Gr5=1)) # NB: prob of Gr5 or less is exactly 1
          ## 4. Transpose F to D*6 matrix
          F <- t(F)
          ## 5. 'Explode' F to a 2D*6 block-antidiagonal matrix
          H <- F[,1:3] # Gr0, Gr1, Gr2
          G <- F[,4:6] # Gr3, Gr4, Gr5
          O <- 0*H
          F <- rbind(cbind(O, G),
                     cbind(H, O))
          ## 6. Normalize the rows of F
          F <- F / rowSums(F)
          ## TODO: Analyze the cases where F can't be normalized, and find a
          ##       principled way of dealing with these. For the time being,
          ##       I'll presume we are dealing with zero-probability events
          ##       in the t(pi) %*% U matrix that render the product with F
          ##       insensitive to those rows of the F matrix.
          ##       I suspect that these cases arise when we draw an extremely
          ##       narrow MTDi distribution (with implausibly small CV) from
          ##       our hyperprior. (One way to deal with these troublesome
          ##       cases would be to avert them entirely by truncating the
          ##       Rayleigh distribution from which CV is drawn, or indeed
          ##       finding another standard distribution altogether. Maybe
          ##       there would be some real benefit from encouraging USERS
          ##       to posit the lower and upper bounds of a uniform density!)
          if (any(is.nan(F))) {
            ## Check that the NaN's are negligible, then zero them
            stopifnot(abs(t(exp(log_pi)) %*% U %*% (is.nan(F))) < 10*.Machine$double.eps)
            F[is.nan(F)] <- 0
          }
          ## Finally, return the matrix product
          t(exp(log_pi)) %*% U %*% F
        }, b, U, kappa) -> expectations
      toxTab <- do.call('rbind', expectations) %>%
        addmargins(margin = 2, FUN = list(Total=sum))
      expectation <- rbind("Expected participants" = colMeans(toxTab)
                          ,"MCSE" = apply(toxTab, MARGIN = 2, FUN = sd) / sqrt(nrow(toxTab))
                           )
      prependClass("safetytab", expectation)
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
  , .doses = numeric(0)
  )
) # </HyperMTDi_lognormal>
