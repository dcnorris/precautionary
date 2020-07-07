## Enhancements to functionality from package 'escalation'

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

