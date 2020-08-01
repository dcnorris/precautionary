#' @include toxicity_generators.R
NULL

#' Simulate trials defined via package \code{escalation}
#' 
#' @docType methods
#' @param selector_factory An object of S3 class \code{selector_factory}
#' @param num_sims Number of simulations to run
#' @param true_prob_tox A generator of toxicity distributions
#' @param ... Passed to subroutines
#' 
#' @details If invoked interactively with \code{num_sims} > 10, then a
#' \code{txtProgressBar} is displayed in the console. The condition on
#' \code{num_sims} has the useful side effect of allowing this function
#' to be invoked iteratively by [extend] (with \code{num_sims} = 10)
#' without the nuisance of nested progress bars.
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
setOldClass(c('stop_at_n_selector_factory',
              'derived_dose_selector_factory',
              'selector_factory'))

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
#' old <- options(dose_levels = c(0.5, 1, 2, 4, 6, 8))
#' mtdi_gen <- hyper_mtdi_lognormal(CV = 1
#'                                  , median_mtd = 6, median_sdlog = 0.5
#'                                  , units="mg/kg")
#' num_sims <- ifelse(interactive()
#' , 300
#' , 15 # avoid taxing CRAN servers
#' )
#' hsims <- get_three_plus_three(num_doses = 6) %>%
#'   simulate_trials(
#'     num_sims = num_sims
#'   , true_prob_tox = mtdi_gen)
#' summary(hsims, ordinalizer=NULL) # vanilla summary with binary toxicity
#' summary(hsims, ordinalizer = function(dose, r0 = sqrt(2))
#'   c(Gr1=dose/r0^2, Gr2=dose/r0, Gr3=dose, Gr4=dose*r0, Gr5=dose*r0^2)
#' )
#' hsims <- hsims %>% extend(num_sims = num_sims)
#' summary(hsims, ordinalizer = function(dose, r0 = sqrt(2))
#'   c(Gr1=dose/r0^2, Gr2=dose/r0, Gr3=dose, Gr4=dose*r0, Gr5=dose*r0^2)
#' )$safety
#' # Set a CRM skeleton from the average probs in above simulation
#' num_sims <- ifelse(interactive()
#' , 16
#' ,  4  # avoid taxing CRAN servers
#' )
#' get_dfcrm(skeleton = hsims$avg_prob_tox
#'          ,target = 0.25
#'          ) %>% stop_at_n(n = 24) %>%
#'   simulate_trials(
#'     num_sims = num_sims
#'   , true_prob_tox = mtdi_gen
#'   ) -> crm_hsims
#' summary(crm_hsims
#' , ordinalizer = function(MTDi, r0 = sqrt(2))
#'     MTDi * r0^c(Gr1=-2, Gr2=-1, Gr3=0, Gr4=1, Gr5=2)
#' )
#' options(old)
#' @rdname simulate_trials
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
setMethod(
  "simulate_trials"
  , c(selector_factory="selector_factory",
      num_sims="numeric",
      true_prob_tox="hyper_mtdi_distribution"),
  function(selector_factory, num_sims, true_prob_tox, ...){
    protocol <- u_i(selector_factory) # separate naming from implementation details
    stopifnot("num_sims must be of length 1" = length(num_sims) == 1)
    # The most parsimonious generalization of the default function
    # will substitute a *matrix* for the default result's vector
    # attribute 'true_prob_tox'.
    dose_levels <- getOption("dose_levels", default = stop(
      "simulate_trials methods require option(dose_levels)."))
    mtdi_samples <- draw_samples(hyper = true_prob_tox, n = num_sims)
    P_ <- seq_along(dose_levels)
    fits <- list()
    tpts <- list()
    if (interactive() && num_sims > 10) pb <- txtProgressBar(max = num_sims, style = 3)
    for(k in 1:num_sims){
      sims_k <- callGeneric(selector_factory = protocol
                           ,num_sims = 1
                           ,true_prob_tox = mtdi_samples[[k]]
                           ,...)
      fits <- c(fits, sims_k[[1]])
      tpts[[k]] <- sims_k$true_prob_tox
      if (exists("pb")) setTxtProgressBar(pb, k)
    }
    if (exists("pb")) close(pb)
    tpt_matrix <- do.call(rbind, tpts)
    colnames(tpt_matrix) <- paste0(dose_levels, true_prob_tox@units)
    sims <- list(
      fits = fits
    # NB: We take the trouble to select the leftmost length(dose_levels) columns
    #     in order to allow for the possibility that the tpt_matrix in general
    #     may also include hyperparameters in columns to the right. This was a
    #     feature (mainly for debugging) of earlier versions of this code, but
    #     seems worth allowing for, pace 'speculative generality'.
    , avg_prob_tox = colMeans(tpt_matrix[, seq_along(dose_levels)])
    , hyper = list(true_prob_tox = tpt_matrix
                  ,mtdi = true_prob_tox
                  ,mtdi_samples = t(sapply(mtdi_samples,function(.)
                    c(CV = .@CV, median = .@median)))
                  )
    , protocol = protocol
    , extra_params = list(...)
    )
    sims$dose_levels <- dose_levels
    sims$dose_units <- true_prob_tox@units
    class(sims) <- "simulations" # impose class (we did not use constructor)
    prependClass(c("hyper","precautionary"), sims)
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
#' old <- options(dose_levels = c(2, 6, 20, 60, 180, 400))
#' mtdi_dist <- mtdi_lognormal(CV = 0.5
#'                            ,median = 140
#'                            ,units = "ng/kg/week")
#' num_sims <- ifelse(interactive()
#' , 100
#' ,  10  # avoid taxing CRAN servers
#' )
#' sims <- get_three_plus_three(num_doses = 6) %>%
#'   simulate_trials(
#'     num_sims = num_sims
#'   , true_prob_tox = mtdi_dist)
#' # Now set a proper ordinalizer via options():
#' options(ordinalizer = function(dose, r0) {
#'   c(Gr1=dose/r0^2, Gr2=dose/r0, Gr3=dose, Gr4=dose*r0, Gr5=dose*r0^2)
#' })
#' summary(sims, r0=2)
#' # Set a CRM skeleton from the average probs in above simulation
#' get_dfcrm(skeleton = sims$true_prob_tox
#'          ,target = 0.25
#'          ) %>% stop_at_n(n = 24) %>%
#'   simulate_trials(
#'     num_sims = 20
#'   , true_prob_tox = mtdi_dist
#'   ) -> crm_sims
#' summary(crm_sims
#' , ordinalizer = function(MTDi, r0 = sqrt(2))
#'     MTDi * r0^c(Gr1=-2, Gr2=-1, Gr3=0, Gr4=1, Gr5=2)
#' )
#' if (interactive()) { # don't overtax CRAN servers
#' crm_sims <- crm_sims %>% extend(target_mcse = 0.1)
#' summary(crm_sims
#' , ordinalizer = function(MTDi, r0 = sqrt(2))
#'     MTDi * r0^c(Gr1=-2, Gr2=-1, Gr3=0, Gr4=1, Gr5=2)
#' )$safety
#' }
#' options(old)
#' @rdname simulate_trials
#' @export
setMethod(
  "simulate_trials"
  , c(selector_factory="selector_factory",
      num_sims="numeric",
      true_prob_tox="mtdi_distribution"),
  function(selector_factory, num_sims, true_prob_tox, ...){
    protocol <- u_i(selector_factory) # separate naming from implementation details
    # TODO: Try moving this next line into the argument defaults
    dose_levels <- getOption("dose_levels", default = stop(
      "simulate_trials methods require option(dose_levels)."))
    tpt_vector <- true_prob_tox@dist$cdf(dose_levels)
    sims <- simulate_trials(selector_factory = protocol
                            , num_sims = num_sims
                            , true_prob_tox = tpt_vector
                            , ... # including, e.g., i_like_big_trials param?
    )
    # In order NOT to have to keep lots of S4 mtdi_distribution objects lying around,
    # especially in the hyperprior-based simulations, we calculate MTDi corresponding
    # to the given u_i, using the current @dist:
    # TODO: Develop accessor methods to abstract away these details of internal structure:
    if( is(sims$fits[[1]][[1]]$fit, "derived_dose_selector") ){
      for(i in seq_along(sims$fits)){
        u_i <- sims$fits[[i]][[1]]$fit$parent$outcomes$u_i
        sims$fits[[i]][[1]]$fit$parent$outcomes$MTDi <- true_prob_tox@dist$quantile(u_i)
      }
    } else {
      for(i in seq_along(sims$fits)){
        u_i <- sims$fits[[i]][[1]]$fit$outcomes$u_i
        sims$fits[[i]][[1]]$fit$outcomes$MTDi <- true_prob_tox@dist$quantile(u_i)
      }
    }
    sims$protocol <- protocol
    sims$extra_params = list(...)
    sims$dose_levels <- dose_levels
    sims$dose_units <- true_prob_tox@units
    sims$mtdi_dist <- true_prob_tox
    prependClass("precautionary", sims)
  }
)

#' Extend an existing simulation, using one of several stopping criteria
#'
#' A trial simulation carried through a predetermined number of replications
#' may not achieve desired precision as judged by Monte Carlo standard errors
#' (MCSE) of estimated toxicity counts. This method enables simulations of
#' class \code{c('precautionary','simulations')} to be extended until a given
#' level of precision has been achieved on expected counts of enrollment and
#' DLTs, or (optionally) for a fixed additional number of simulations.
#'
#' @param sims An existing object of class \code{c('precautionary','simulations')}
#' @param num_sims Optionally, a fixed number of additional replications to accumulate
#' @param target_mcse Optionally, an MCSE constraint to be imposed on expected counts
#'  of DLTs, non-DLTs, and Total enrollment.
#'
#' @return An extended simulation of same class as \code{sims}.
#' 
#' @note The MCSE constraint is imposed during trial \emph{simulation}, at which
#'  point only \emph{binary} toxicities are available. Thus, as a practical matter,
#'  \code{extend} can target MCSEs only for DLTs, non-DLTs and Total enrollment.
#'  The subsequent subdivision of these categories during trial \emph{summary}
#'  (at which point the ordinalizer comes into play along with its parameters) thus
#'  may generate expected counts with MCSEs exceeding \code{target_mcse}.
#'  In practice, however, this tends to affect the estimated counts only for the
#'  \emph{lowest} toxicity grades---those of least concern from a trial-safety
#'  perspective.
#' 
#' @export
#'
#' @seealso
#' See examples under [simulate_trials] and [format.safetytab].
extend <- function(sims, num_sims = NULL, target_mcse = 0.05) {
  UseMethod('extend')
}

toxcount_mcse <- function(sims) {
  se <- function(x) sqrt(var(x)/length(x))
  tox <- notox <- total <- NULL # avoid 'no visible binding for global var' NOTEs 
  mcse <- as.data.table(sims, ordinalizer = NULL)[
    , list(tox=sum(tox), notox=sum(!tox), total=.N), by = rep][
      , list(tox=se(tox), notox=se(notox), total=se(total))
      ]
  # Given that trial duration will be inversely correlated with
  # observed toxicity rate, we may generally expect the MCSE of
  # non-DLT counts to be greater than that of DLT counts.
  if (mcse$tox > mcse$notox)
    warning("MCSE of DLT counts exceeds that of non-DLTs.",
            " This is contrary to expectations, given that",
            " the sum of these 2 counts is trial enrollment,",
            " which will generally be inversely correlated",
            " to observed DLT rate.")
  max(mcse)
}

extension_to_target_mcse <- function(sims, target_mcse) {
  current_mcse <- toxcount_mcse(sims)
  if ( current_mcse < target_mcse )
    return(0)
  extension_factor <- ( current_mcse / target_mcse )^2 - 1
  num_sims <- length(sims$fits) * extension_factor
  num_sims <- ceiling(num_sims) # awkward to do fractional sims ;^)
  num_sims
}

#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
extend.precautionary <- function(sims, num_sims = NULL, target_mcse = 0.05) {
  options(dose_levels = sims$dose_levels)
  if (!missing(num_sims)) { # base case for recursion
    more <- do.call(simulate_trials, c(list(selector_factory = sims$protocol
                                            ,num_sims = num_sims
                                            ,true_prob_tox = sims$mtdi_dist)
                                       ,sims$extra_params
                                       )
                    )
    sims$fits <- c(sims$fits, more$fits)
    return(sims)
  }
  # The more sensible use case: extending sim to target MCSEs
  sims_todo_est <- extension_to_target_mcse(sims, target_mcse = target_mcse)
  if(interactive()) pb <- txtProgressBar(max = 1, style = 3)
  while (sims_todo_est > 0) {
    sims <- extend(sims, num_sims = 10) # NB: num_sims < 11 avoids nested progress bar
    sims_done <- length(sims$fits)
    sims_todo_est <-extension_to_target_mcse(sims, target_mcse = target_mcse)
    fraction_complete <- sims_done / (sims_done + sims_todo_est)
    if (interactive()) setTxtProgressBar(pb, fraction_complete)
  }
  if (interactive()) close(pb)
  sims
}

#' @export
extend.hyper <- function(sims, num_sims = NULL, target_mcse = 0.05) {
  if (!missing(num_sims)) { # base case for recursion
    more <- do.call(simulate_trials, c(list(selector_factory = sims$protocol
                                            ,num_sims = num_sims
                                            ,true_prob_tox = sims$hyper$mtdi)
                                       ,sims$extra_params
                                       )
                    )
    N1 <- length(sims$fits)
    N2 <- length(more$fits)
    sims$avg_prob_tox <- (sims$avg_prob_tox*N1 + more$avg_prob_tox*N2)/(N1+N2)
    sims$fits <- c(sims$fits, more$fits)
    sims$hyper$true_prob_tox <- rbind(sims$hyper$true_prob_tox
                                      ,more$hyper$true_prob_tox)
    sims$hyper$mtdi_samples <- c(sims$hyper$mtdi_sample
                                 ,more$hyper$mtdi_sample)
    return(sims)
  }
  # The more sensible use case: extending sim to target MCSEs
  NextMethod() # recursive case handled fine by extend.precautionary
}

