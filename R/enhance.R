## Enhancements to functionality from package 'escalation'

prependClass <- function(preclass, object) {
  class(object) <- unique(c(preclass, class(object)))
  object
}

setOldClass(c("u_i","tox_selector_factory","selector_factory"))

u_i <- function(tox_selector_factory) {
  stopifnot("Class 'u_i' may be applied only to a tox_selector_factory"
            = is(tox_selector_factory,"tox_selector_factory"))
  prependClass("u_i", tox_selector_factory)
}

#' Get a function that simulates dose-escalation trials using latent \code{u_i}
#' 
#' Overrides \code{escalation::simulation_function.tox_selector_factory}
#' to return a [phase1_sim()] that employs latent \code{u_i} in place of
#' the version native to package \code{escalation}, which merely invokes
#' \code{rbinom}.
#' 
#' @param selector_factory Presently, this must be a \code{tox_selector_factory};
#' no equivalent for \code{simulation_function.derived_dose_selector_factory}
#' is yet implemented.
#'
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
#' and arrival-time increment for each trial participant.
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
  suppressWarnings(u_i <- base_df$u_i) # RStudio "Unknown or uninitialised column"
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

#' @export
prob_recommend.precautionary <- function(x, ...) {
  prob_recs <- NextMethod()
  names(prob_recs)[-1] <- paste0(x$dose_levels, x$dose_units)
  prob_recs
}

#' @export
prob_administer.precautionary <- function(x, ...) {
  prob_admin <- NextMethod()
  names(prob_admin) <- paste0(x$dose_levels, x$dose_units)
  prob_admin
}

#' @export
print.precautionary <- function(x, ...) {
  
  cat('Number of iterations:', length(x$fits), '\n')
  cat('\n')
  
  cat('Number of doses:', num_doses(x), '\n')
  cat('\n')
  
  cat('True probability of toxicity:\n')
  print(x$true_prob_tox, digits = 3)
  cat('\n')
  
  cat('Probability of recommendation:\n')
  print(prob_recommend(x), digits = 3)
  cat('\n')
  
  cat('Probability of administration:\n')
  print(prob_administer(x), digits = 3)
  cat('\n')
  
  cat('Sample size:\n')
  print(summary(num_patients(x)))
  cat('\n')
  
  cat('Total toxicities:\n')
  print(summary(num_tox(x)))
  cat('\n')
  
  cat('Trial duration:\n')
  print(summary(trial_duration(x)))
  cat('\n')
}


#' @export
print.hyper <- function(x, ...) {
  
  cat('Number of iterations:', length(x$fits), '\n')
  cat('\n')
  
  cat('Number of doses:', num_doses(x), '\n')
  cat('\n')
  
  cat('Average probability of toxicity:\n')
  print(x$avg_prob_tox, digits = 3)
  cat('\n')
  
  cat('Probability of recommendation:\n')
  print(prob_recommend(x), digits = 3)
  cat('\n')
  
  cat('Probability of administration:\n')
  print(prob_administer(x), digits = 3)
  cat('\n')
  
  cat('Sample size:\n')
  print(summary(num_patients(x)))
  cat('\n')
  
  cat('Total toxicities:\n')
  print(summary(num_tox(x)))
  cat('\n')
  
  cat('Trial duration:\n')
  print(summary(trial_duration(x)))
  cat('\n')
}


#' @export
as.data.table.precautionary <- function(x, ordinalizer = getOption('ordinalizer')
                                        , ...) {
  ensemble <- rbindlist(lapply(x[[1]], function(.) .[[1]]$fit$outcomes)
                        , idcol = "rep")
  # Go 'straight to dose-space' by generating MTDi,g columns
  if( is(x,"hyper") ){ # TODO: Consider handling this via 'as.data.table.hyper'
    K <- length(x$hyper$mtdi_samples)
    M <- max(ensemble$rep)/K  # so now K*M is nrow(ensemble)
    ensemble[, k := 1 + (rep-1) %/% M]
    ensemble[, rep := 1 + (rep-1) %% M]
    setcolorder(ensemble,c("k","rep")) # so key .(k,rep) is on far left
    look <<- ensemble
    ensemble[, MTDi := x$hyper$mtdi_samples[[k]]@dist$quantile(u_i), by = k]
  } else {
    ensemble[, MTDi := x$mtdi_dist@dist$quantile(ensemble$u_i)]
  }
  if( is.null(ordinalizer) )
    return(ensemble)
  # TODO: Do add these columns to 'ensemble', so that the whole table
  #       may later be inspected by user to improve understanding.
  MTDig <- t(sapply(ensemble$MTDi, ordinalizer, ...))
  if( is.null(colnames(MTDig)) ) {
    warning("Ordinalizer returns unnamed vector; using default names for toxicity grades.")
    colnames(MTDig) <- paste("Grade", 1:ncol(MTDig))
  }
  tox_grades <- colnames(MTDig)
  # Compare with actual dose to obtain toxicity grade indicator matrix
  tox_ind <- ( x$dose_levels[ensemble$dose] > MTDig )
  # Tally the thresholds crossed to obtain integer toxgrade
  ensemble$toxgrade <- rowSums(tox_ind)
  # Convert toxgrade to an ordered factor Tox
  ensemble$Tox <- ordered(ensemble$toxgrade+1
                          , levels=seq(1+length(tox_grades))
                          , labels=c('None', tox_grades))
  ensemble
}


#' @importFrom dplyr rename rename_with mutate select
#' @export
summary.precautionary <- function(x, ordinalizer = getOption('ordinalizer'), ...) {
  summary <- NextMethod()
  if(!is(x,"simulations")) return(summary) # override only summary.simulations
  dose_units <- paste0("dose (", x$dose_units, ")")
  summary <- summary %>%
    mutate("real_dose" = c(0, x$dose_levels)[as.integer(summary$dose)]) %>%
    select(dose, real_dose, everything()) %>%
    rename_with(.fn = function(.) dose_units, .cols = real_dose)
  ensemble <- as.data.table(x, ordinalizer = ordinalizer)
  attr(summary,'ensemble') <- ensemble # for debugging
  if( !is.null(ordinalizer) ){
    K <- c(nrow(x$hyper$true_prob_tox), 1)[1] # NB: c(NULL,1) = c(1)
    expectation <- colMeans(xtabs(~ rep + Tox, data=ensemble))/K
    expectation <- c(expectation, All = sum(expectation))
    expectation <- t(as.matrix(expectation))
    rownames(expectation) <- "Expected participants"
    attr(summary,'safety') <- expectation
  }
  summary
}

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

