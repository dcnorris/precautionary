## Enhancements to functionality from package 'escalation'

prependClass <- function(preclass, object) {
  class(object) <- unique(c(preclass, class(object)))
  object
}

u_i <- function(tox_selector_factory) {
  stopifnot("Class 'u_i' may be applied only to a tox_selector"
            = is(tox_selector_factory,"tox_selector_factory"))
  prependClass("u_i", tox_selector_factory)
}

#' @export
simulation_function.u_i <- function(selector_factory) {
  return(phase1_sim) # returns an override defined in this package
}

#' Override \code{escalation::cohorts_of_n} to include latent toxicity tolerances
#' 
#' The original function in package \code{escalation} recognizes that individual
#' trial participants arrive at distinct times. Building upon this acknowledgement
#' of individuality, this override adds an extra line of code to draw as well a
#' latent toxicity tolerance \code{u_i} for each individual participant.
#' 
#' @seealso [phase1_sim()], which this package also overrides with similarly
#' minute changes to incorporate \code{u_i}.
#'
#' @param n integer, sample arrival times for this many patients.
#' @param mean_time_delta the average gap between patient arrival times. I.e.
#' the reciprocal of the rate parameter in an Exponential distribution.
#' @return \code{data.frame} with columns \code{u_i} and \code{time_delta}
#' containing respectively the uniformly-distributed latent toxicity tolerance
#' and intervals of time between patient arrivals.
#' 
#' @importFrom stats rexp runif
#' @examples
#' cohorts_of_n()
#' cohorts_of_n(n = 10, mean_time_delta = 5)
#' @export
cohorts_of_n <- function(n = 3, mean_time_delta = 1) {
  u_i <- runif(n = n)
  time_delta <- rexp(n = n, rate = 1 / mean_time_delta) %>% round(1)
  data.frame(u_i = u_i, time_delta = time_delta)
}

#' Override \code{escalation::phase1_sim} to incorporate latent toxicity tolerances
#'
#' @param selector_factory 
#' @param true_prob_tox 
#' @param sample_patient_arrivals 
#' @param previous_outcomes 
#' @param next_dose 
#' @param i_like_big_trials 
#' @param return_all_fits 
#'
#' @importFrom magrittr %>%
#' @importFrom utils tail
phase1_sim <- function(
  selector_factory,
  true_prob_tox,
  sample_patient_arrivals = function(df) cohorts_of_n(n=3, mean_time_delta=1),
  previous_outcomes = '',
  next_dose = NULL,
  i_like_big_trials = FALSE, # Safety net if stop_trial_func is mis-specified...
  return_all_fits = FALSE
) {
  if(is.character(previous_outcomes)) {
    base_df <- parse_phase1_outcomes(previous_outcomes, as_list = FALSE)
  } else if(is.data.frame(previous_outcomes)) {
    base_df <- spruce_outcomes_df(previous_outcomes)
  } else{
    base_df <- parse_phase1_outcomes('', as_list = FALSE)
  }
  dose <- base_df$dose
  u_i <- base_df$u_i
  tox <- base_df$tox
  cohort <- base_df$cohort
  next_cohort <- ifelse(length(cohort) > 0, max(cohort) + 1, 1)
  if('time' %in% colnames(base_df)) {
    time <- previous_outcomes$time
  } else {
    time <- rep(0, length(dose))
  }
  
  i <- 1 # loop counter
  max_i <- 30
  time_now <- 0
  fit <- selector_factory %>% fit(base_df)
  if(is.null(next_dose)) next_dose <- fit %>% recommended_dose()
  fits <- list()
  fits[[1]] <- list(.depth = i, time = time_now, fit = fit)
  while(fit %>% continue() & !is.na(next_dose) &
        (i_like_big_trials | i < max_i)) {
    
    current_data = data.frame(
      cohort = cohort,
      patient = seq_along(dose),
      dose = dose,
      tox = tox,
      time = time
    )
    new_pts <- sample_patient_arrivals(current_data)
    arrival_time_deltas <- cumsum(new_pts$time_delta)
    n_new_pts <- nrow(new_pts)
    new_dose <- rep(next_dose, n_new_pts)
    new_tox <- ( new_pts$u_i < true_prob_tox[next_dose] )
    new_cohort <- rep(next_cohort, n_new_pts)
    
    dose <- c(dose, new_dose)
    u_i <- c(u_i, new_pts$u_i)
    tox <- c(tox, new_tox)
    cohort <- c(cohort, new_cohort)
    time <- c(time, time_now + arrival_time_deltas)
    new_data = data.frame(
      cohort = cohort,
      patient = 1:length(dose),
      dose = dose,
      u_i = u_i,
      tox = tox,
      time = time
    )
    
    time_now <- time_now + max(arrival_time_deltas)
    i <- i + 1
    fit <- selector_factory %>% fit(new_data)
    next_cohort <- next_cohort + 1
    fits[[i]] <- list(.depth = i, time = time_now, fit = fit)
    next_dose <- fit %>% recommended_dose()
  }
  
  # Warn about i_like_big_trials if sim stopped because of too big i.
  if(!i_like_big_trials & i >= max_i) {
    warning(paste(
      "Simulation stopped because max depth reached.",
      "Set 'i_like_big_trials = TRUE' to avoid this constraint. "))
  }
  
  if(return_all_fits) {
    return(fits)
  } else {
    return(tail(fits, 1))
  }
}


#' @importFrom dplyr rename rename_with mutate select
#' @export
summary.precautionary <- function(x, ...) {
  summary <- NextMethod()
  if(!is(x,"simulations")) return(summary) # override only summary.simulations
  dose_units <- paste0("dose (", x$dose_units, ")")
  summary <- summary %>%
    mutate("real_dose" = c(0, x$dose_levels)[as.integer(summary$dose)]) %>%
    select(dose, real_dose, everything()) %>%
    rename_with(.fn = function(.) dose_units, .cols = real_dose)
  if(!is.null(attr(x$true_prob_tox,'ordtox'))){
    # TODO: Handle case where attr(x,'ordtox') is present
    attr(summary,'ordtox') <- "TODO: Summarize ORDTOX attribute, too."
  }
  summary
}

# TODO: Check whether I still need these num_doses.* and dose_indices.*
#       methods, after abandoning the 'ordtox' class and its associated
#       'check_safety' syntactical sugar.

#' @export
num_doses.three_plus_three_selector_factory <- function(x, ...) {
  return(x$num_doses)
}

#' @export
num_doses.dfcrm_selector_factory <- function(x, ...) {
  return(length(x$skeleton))
}

#' @export
num_doses.boin_selector_factory <- function(x, ...) {
  return(x$num_doses)
}

#' @export
dose_indices.default <- function(x, ...) {
  n <- num_doses(x)
  if(n > 0) {
    return(1:n)
  } else {
    return(integer(length = 0))
  }
}

