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
#' options(dose_levels = c(0.5, 1, 2, 4, 6, 8))
#' mtdi_gen <- hyper_mtdi_lognormal(CV = 1
#'                                  , median_mtd = 6, median_sdlog = 0.5
#'                                  , units="mg/kg")
#' num_sims <- ifelse(interactive()
#' , c(15, 20)
#' , c( 3,  5) # avoid taxing CRAN servers
#' )
#' hsims <- get_three_plus_three(num_doses = 6) %>%
#'   simulate_trials(
#'     num_sims = num_sims
#'   , true_prob_tox = mtdi_gen)
#' summary(hsims, ordinalizer=NULL) # vanilla summary with binary toxicity
#' summary(hsims, ordinalizer = function(dose, r0 = sqrt(2))
#'   c(Gr1=dose/r0^2, Gr2=dose/r0, Gr3=dose, Gr4=dose*r0, Gr5=dose*r0^2)
#' )
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
    stopifnot("num_sims must be of length 1 or 2" = length(num_sims) %in% 1:2)
    M <- num_sims[1]
    K <- num_sims[2]; if(is.na(K)) K <- num_sims[1]
    # The most parsimonious generalization of the default function
    # will substitute a *matrix* for the default result's vector
    # attribute 'true_prob_tox'.
    dose_levels <- getOption("dose_levels", default = stop(
      "simulate_trials methods require option(dose_levels)."))
    mtdi_samples <- draw_samples(hyper = true_prob_tox, K = K)
    P_ <- seq_along(dose_levels)
    fits <- list()
    ordtox <- list()
    tpts <- list()
    if(interactive()) pb <- txtProgressBar(max=K, style=3)
    for(k in 1:K){
      sims_k <- callGeneric(selector_factory = protocol
                           ,num_sims = M
                           ,true_prob_tox = mtdi_samples[[k]]
                           ,...)
      fits <- c(fits, sims_k[[1]])
      ordtox[[k]] <- attr(sims_k,'ordtox')
      tpts[[k]] <- sims_k$true_prob_tox
      if(interactive()) setTxtProgressBar(pb, k)
    }
    tpt_matrix <- do.call(rbind, tpts)
    colnames(tpt_matrix) <- paste0(dose_levels, true_prob_tox@units)
    sims <- list(
      fits = fits
    , ordtox_check = ordtox
    # NB: We take the trouble to select the leftmost length(dose_levels) columns
    #     in order to allow for the possibility that the tpt_matrix in general
    #     may also include hyperparameters in columns to the right. This was a
    #     feature (mainly for debugging) of earlier versions of this code, but
    #     seems worth allowing for, pace 'speculative generality'.
    , avg_prob_tox = colMeans(tpt_matrix[, seq_along(dose_levels)])
    , hyper = list(true_prob_tox = tpt_matrix
                  ,mtdi = true_prob_tox
                  ,mtdi_samples = mtdi_samples
                  )
    )
    sims$dose_levels <- dose_levels
    sims$dose_units <- true_prob_tox@units
    if(!is.null(unlist(ordtox))){
      attr(sims,'ordtox') <- rbindlist(ordtox, idcol = "k")
    }
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
#' # TODO: In keeping with spirit of REALISM, attach UNITS to 'dose_levels'.
#' options(dose_levels = c(2, 6, 20, 60, 180, 400))
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
#' # Now attach a proper ordinalizer to 'mtdi_dist':
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

#' @rdname simulate_trials
#' @export
setMethod(
  "simulate_trials"
  , c(selector_factory="selector_factory",
      num_sims="numeric",
      true_prob_tox="mtdi_distribution"),
  function(selector_factory, num_sims, true_prob_tox, ...){
    # TODO: Try moving this next line into the argument defaults
    dose_levels <- getOption("dose_levels", default = stop(
      "simulate_trials methods require option(dose_levels)."))
    tpt_vector <- true_prob_tox@dist$cdf(dose_levels)
    sims <- simulate_trials(selector_factory = u_i(selector_factory)
                            , num_sims = num_sims
                            , true_prob_tox = tpt_vector
                            , ... # including, e.g., i_like_big_trials param?
    )
    sims$dose_levels <- dose_levels
    sims$dose_units <- true_prob_tox@units
    sims$mtdi_dist <- true_prob_tox
    prependClass("precautionary", sims)
  }
)

