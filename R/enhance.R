## Enhancements to functionality from package 'escalation'

prependClass <- function(preclass, object) {
  class(object) <- unique(c(preclass, class(object)))
  object
}

setOldClass(c("u_i","tox_selector_factory","selector_factory"))

u_i <- function(selector_factory) {
  stopifnot("Class 'u_i' applies only to a (tox|derived_dose)_selector_factory"
            = is(selector_factory,"tox_selector_factory") || 
              is(selector_factory,"derived_dose_selector_factory")
            )
  prependClass("u_i", selector_factory)
}

#' Get a function that simulates dose-escalation trials using latent \code{u_i}
#' 
#' Overrides \code{escalation::simulation_function.tox_selector_factory}
#' to return a [phase1_sim()] that employs latent \code{u_i} in place of
#' the version native to package \code{escalation}, which merely invokes
#' \code{rbinom}.
#' 
#' This function is exported for the purpose of effecting this override,
#' and is not meant to be invoked directly by the user.
#' 
#' @param selector_factory Presently, this must be a \code{tox_selector_factory};
#' no equivalent for \code{simulation_function.derived_dose_selector_factory}
#' is yet implemented.
#'
#' @importFrom escalation simulation_function
#' @export
simulation_function.u_i <- function(selector_factory) {
  return(phase1_sim) # returns an override defined in this package
}

#' Override \code{escalation::cohorts_of_n} to include latent toxicity tolerances
#' 
#' The original function in package \code{escalation} recognizes that individual
#' trial participants arrive at distinct times. Building upon this acknowledgment
#' of individuality, this override adds an extra line of code to draw as well a
#' latent toxicity tolerance \code{u_i} for each individual participant.
#' 
#' @seealso [phase1_sim()], which this package also overrides with similarly
#' minute changes in order to incorporate \code{u_i}.
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
#' @param selector_factory An \code{escalation::selector_factory} object
#' @param true_prob_tox A vector of toxicity probabilities for the doses
#'  defined in \code{selector_factory} 
#' @param sample_patient_arrivals A function implementing an arrivals process
#'  for trial enrollment
#' @param previous_outcomes This may or may not apply in applications of
#'  package \code{precautionary}
#' @param next_dose Undocumented
#' @param i_like_big_trials I didn't choose this parameter name
#' @param return_all_fits Don't do this
#'
#' @importFrom magrittr %>%
#' @importFrom utils tail
#' @importFrom escalation recommended_dose continue parse_phase1_outcomes
phase1_sim <- function(
  selector_factory,
  true_prob_tox,
  sample_patient_arrivals = function(df) cohorts_of_n(n=3, mean_time_delta=1),
  previous_outcomes = '',
  next_dose = NULL,
  i_like_big_trials = FALSE, # Safety net if stop_trial_func is mis-specified...
  return_all_fits = FALSE
) {
  spruce_outcomes_df <- function(df) { # TODO: Request 'escalation' export this fun
    df$dose <- as.integer(df$dose)
    df$tox <- as.integer(df$tox)
    if('cohort' %in% colnames(df)) df$cohort <- as.integer(df$cohort)
    if('patient' %in% colnames(df)) df$patient <- as.integer(df$patient)
    df
  }
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

#' @importFrom escalation prob_recommend
prob_recommend.precautionary <- function(x, ...) {
  prob_recs <- NextMethod()
  names(prob_recs)[-1] <- paste0(x$dose_levels, x$dose_units)
  prob_recs
}

#' @importFrom escalation prob_administer
prob_administer.precautionary <- function(x, ...) {
  prob_admin <- NextMethod()
  names(prob_admin) <- paste0(x$dose_levels, x$dose_units)
  prob_admin
}

#' Specialize print method for objects of class \code{escalation::simulations}
#'
#' @param x An object of class  c("precautionary","simulations")
#'
#' @param ... Additional arguments; ignored
#'
#' @importFrom escalation num_patients num_tox trial_duration
print.precautionary <- function(x, ...) {
  
  cat('Number of iterations:', length(x$fits), '\n')
  cat('\n')
  
  cat('Number of doses:', length(x$dose_levels), '\n')
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


#' Specialize print method defined for class \code{escalation::simulations}
#'
#' @param x An object of class  c("hyper","precautionary","simulations")
#'
#' @param ... Additional arguments; ignored
#'
#' @importFrom escalation num_patients num_tox trial_duration
print.hyper <- function(x, ...) {
  
  cat('Number of iterations:', length(x$fits), '\n')
  cat('\n')
  
  cat('Number of doses:', length(x$dose_levels), '\n')
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


#' Convert an object of class c('precautionary','simulations') to a data.table
#'
#' @param x An object of class c('precautionary','simulations')
#'
#' @param keep.rownames Unused; retained for S3 generic/method consistency
#' @param ordinalizer If not NULL, this is a function mapping the threshold
#'  dose ('MTDi') at which an individual experiences a binary toxicity (as
#'  recognized by the dose-escalation design) to a named vector giving dose
#'  thresholds for multiple grades of toxicity. The names of this vector will
#'  be taken as designations of the toxicity grades.
#' @param ... Additional parameters passed to the \code{ordinalizer}
#'
#' @export
as.data.table.precautionary <- function(x, keep.rownames = FALSE
                                        , ordinalizer = getOption('ordinalizer')
                                        , ...) {
  extractor <- ifelse(is(x$fits[[1]][[1]]$fit, "derived_dose_selector")
                     ,function(.) .[[1]]$fit$parent$outcomes
                     ,function(.) .[[1]]$fit$outcomes
                     )
  ensemble <- rbindlist(lapply(x$fits, extractor), idcol = "rep")
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

#' Specialize a method defined in package 'escalation' for class 'simulations'
#' 
#' Simulations produced by package \code{precautionary} incorporate a \sQuote{u_i}
#' latent toxicity tolerance that characterizes the toxic dose-response of each
#' simulated individual trial participant. In conjunction with an
#' \sQuote{ordinalizer} function, these extra data enable questions to be asked
#' about trial safety, in terms of the probabilities of high-grade toxicities.
#' This function specializes the \code{escalation::summary.simulations} method
#' accordingly.
#' 
#' @param object An object of class c('precautionary','simulations') 
#'
#' @param ordinalizer An ordinalizer function
#' @param ... Additional parameters passed to the ordinalizer
#'
#' @importFrom dplyr mutate rename_with select everything
#' @importFrom stats xtabs addmargins sd var
#' @importFrom rlang .data
#' @export
summary.precautionary <- function(object, ordinalizer = getOption('ordinalizer'), ...) {
  summary <- NextMethod()
  dose_units <- paste0("dose (", object$dose_units, ")")
  # For an explanation of the .data 'pronoun' used below,
  # see https://dplyr.tidyverse.org/articles/programming.html.
  summary <- summary %>%
    mutate("real_dose" = c(0, object$dose_levels)[as.integer(.data$dose)]) %>%
    select(.data$dose, .data$real_dose, everything()) %>%
    rename_with(.fn = function(.) dose_units, .cols = .data$real_dose)
  ensemble <- as.data.table(object, ordinalizer = ordinalizer, ...)
  if( !is.null(ordinalizer) ){
    summary <- list(escalation = summary, safety = NULL)
    toxTab <- xtabs(~ rep + Tox, data=ensemble) %>%
      addmargins(margin = 2, FUN = list(Total=sum))
    expectation <- rbind("Expected participants" = colMeans(toxTab)
                        ,"MCSE" = apply(toxTab, MARGIN = 2, FUN = sd) / sqrt(nrow(toxTab))
                        )
    summary$safety <- prependClass("safetytab", expectation)
    summary$toxTab <- toxTab # (for DEBUGGING purposes)
  }
  summary
}

#' Format a phase 1 trial safety tabulation to show significant digits only
#'
#' The essential insight of package [precautionary] is distilled into the
#' \emph{safety tabulation} which it generates from trial simulations, reporting
#' the expected number of patients who will experience each grade of toxicity.
#' To render this table for easy interpetation, these expectations are simply
#' displayed with a number of significant digits appropriate to their Monte Carlo
#' standard errors (MCSEs). 
#'
#' @param x A safety tabulation as found in the \code{safety} component of the
#'  list returned by [summary.precautionary].
#'
#' @param ... Unused; included for compatibility with generic signature
#'
#' @note The MCSEs of safety tabulations remain available for inspection
#'  (see example), but are omitted from standard displays because they may lend
#'  themselves to misinterpretation as \emph{confidence bounds} on the number
#'  of patients who will experience each toxicity grade \emph{in any given trial}.
#'
#' @examples 
#' mtdi_gen <- hyper_mtdi_lognormal(CV = 1
#'                                 ,median_mtd = 5
#'                                 ,median_sdlog = 0.5
#'                                 ,units="mg/kg")
#' ordinalizer <- function(MTDi, r0 = 1.5)
#'   MTDi * r0 ^ c(Gr1=-2, Gr2=-1, Gr3=0, Gr4=1, Gr5=2)
#' old <- options(dose_levels = c(0.5, 1, 2, 4, 6)
#'               ,ordinalizer = ordinalizer)
#' get_boin(num_doses = 5, target = 0.25) %>%
#'   stop_at_n(n = 24) %>%
#'   simulate_trials(
#'     num_sims = 80
#'   , true_prob_tox = mtdi_gen) -> boin_hsims
#' safety <- summary(boin_hsims)$safety
#' safety # The print method invokes 'format.safetytab' ..
#' # .. but we can also inspect the underlying matrix by indexing:
#' safety[,]  # indexing strips 'safetytab' class, returning plain matrix
#' # Note that, by extend()ing the simulation we can increase precision:
#' if (interactive()) { # may run a bit too long for CRAN servers' taste
#'   boin_hsims %>% extend(target_mcse = 0.1) -> boin_hsimsX 
#'   summary(boin_hsimsX)$safety
#' }
#' options(old)
#' @export
format.safetytab <- function(x, ...){
  sigtenths <- x['MCSE',] < 0.1
  mapply(function(x, d) format(round(x, digits=d), nsmall=d), x[1,], 0+sigtenths)
}

#' @importFrom stringr str_pad
#' @export
print.safetytab <- function(x, ...) {
  fx <- format(x, ...)
  width <- 2 + nchar(names(fx))
  writeLines(paste(str_pad(names(fx), width=width, "left"), collapse=""))
  writeLines(paste(str_pad(unname(fx), width=width, "left"), collapse=""))
  invisible(x)
}

#' Output a 'kable' for a simulation summary of class 'safetytab'
#'
#' @param safetytab An object of S3 class \code{safetytab}
#'
#' @param ... Additional parameters passed to [knitr::kable]
#'
#' @importFrom knitr kable
#' @importFrom kableExtra kable_styling add_header_above
#' @export
safety_kable <- function(safetytab, ...) {
  safetytab %>% format() %>% t() %>%
    kable(align='r', ...) %>% kable_styling(position = "left", full_width = FALSE) %>%
    add_header_above(c("Expected counts per toxicity grade"=6, " "=1))
}

#' @export
t.safetytab <- function(x) t(x[,]) # plain-matrix transpose (drops 'safetytab' class)

#' @importFrom escalation num_doses
num_doses.three_plus_three_selector_factory <- function(x, ...) {
  return(x$num_doses)
}

num_doses.dfcrm_selector_factory <- function(x, ...) {
  return(length(x$skeleton))
}

num_doses.boin_selector_factory <- function(x, ...) {
  return(x$num_doses)
}

#' @importFrom escalation dose_indices
#' @export
dose_indices.default <- function(x, ...) {
  n <- num_doses(x)
  if(n > 0) {
    return(1:n)
  } else {
    return(integer(length = 0))
  }
}

