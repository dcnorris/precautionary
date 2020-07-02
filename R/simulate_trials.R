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
#           , def = function(selector_factory, num_sims, true_prob_tox, ...) NULL)

setOldClass("ordtox")
setOldClass("selector_factory")
setOldClass("three_plus_three_selector_factory")

#' @rdname simulate_trials
setMethod("simulate_trials", c(selector_factory="ordtox"),
  function(selector_factory, num_sims, true_prob_tox, ...){
    sim_safety(selector_factory, ...)
  })

#' Syntactic sugar for injecting safety analysis into \code{simulate_trials}
#'
#' @param x A selector_factory, required as first argument to \code{simulate_trials}
#'
#' @return x with class \sQuote{ordtox} prepended
#' @export
#'
#' @examples
#' 
#' # Simulate a 3+3 trial
#' protocol <- get_three_plus_three(num_doses = 6)
#' mu <- 5              # Values estimated in arXiv:2004.12755 [stat, q-bio];
#' sigma <- 1/sqrt(1.3) # see subsection 'Model estimates' and Table 2, pp. 2-3.
#' dose_levels <- c(2, 6, 20, 60, 180, 400)
#' MTDi <- Lognormal$new(meanlog = mu, sdlog = sigma)
#' DLT_probs <- MTDi$cdf(dose_levels)
#' sims <- simulate_trials(protocol, num_sims = 20, true_prob_tox = DLT_probs)
#' 
#' # Do a richer simulation, examining safety qua avoidance of Grade 4 & 5
#' # toxicities, under a default generative model for ordinal toxicities:
#' safety_sims <- simulate_trials(check_safety(protocol))
check_safety <- function(x){
  stopifnot("x must be a selector_factory" = is(x, "selector_factory"))
  class(x) <- c("ordtox", class(x))
  x
}
