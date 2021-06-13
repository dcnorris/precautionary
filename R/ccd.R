## An R6 class for Cumulative-Cohort Designs (CCDs), copied from 'Crm'
## to repurpose its parallelized, recursive Comprehensive Path Enumeration

## TODO: Restructure the class hierarchy, perhaps with CPE as a base class
##       or mixin shared by both 'Ccd' and 'Crm'.

##' @name Ccd-class
##' @title An R6 class encapsulating Cumulative-Cohort Designs
##'
##' @details
##' Syntactically, the method chaining supported by R6 classes makes the
##' invocation of CRM models more transparent. The mutability conferred
##' by reference semantics enables memoization (caching) of results, which
##' can speed up DTP calculations significantly.
##'
##' Presently, this class supports only the 'empiric' (aka 'power') model.
##' But it is hoped that inheritance will assist in rendering other models
##' implemented in package \CRANpkg{dfcrm} clearly, with code reuse.
##' @importFrom R6 R6Class
##' @export
Ccd <- R6Class("Ccd",
               inherit = Cpe,
               public = list(
                 ##' @details
                 ##' Create a new `Ccd` object.
                 ##'
                 ##' @param escalate Escalation boundary
                 ##' @param deescalate Deescalation boundary
                 ##' @param eliminate Elimination boundary
                 ##' @param cohort_max Upper bound on dose-wise enrollment
                 ##' @param enroll_max Upper bound on total enrollment
                 ##' @return A Ccd object.
                 ##'
                 ##' @examples
                 ##' # TODO
                 initialize = function(escalate, deescalate, eliminate, cohort_max, enroll_max) {
                   private$escalate <- escalate
                   private$deescalate <- deescalate
                   ## The 'severe' action of dose elimination may not apply until some
                   ## minimum cohort enrollment (e.g., 3), in which case \CRANpkg{BOIN}
                   ## denotes inapplicable entries by `NA`. But converting these to `Inf`
                   ## renders the 'lifting' of the constraint more elegantly:
                   private$eliminate <- { eliminate[is.na(eliminate)] <- Inf; eliminate }
                   private$cohort_max <- cohort_max
                   private$enroll_max <- enroll_max
                 },
                 ##' @details
                 ##' Return dose recommendation for given tox/no-tox tallies.
                 ##'
                 ##' @param x A dose-wise vector of toxicity counts
                 ##' @param o A dose-wise vector of non-toxicity counts
                 ##' @param last_dose The most recently given dose, as required to implement
                 ##' cumulative-cohort-based escalation decisions.
                 ##' @param max_dose An upper limit on future dose levels
                 ##' @param ... Parameters passed to \code{Crm$esc()}, enabling passthru
                 ##' of required \code{impl} parameter and optional \code{abbrev} flag.
                 ##' @return An object with components:
                 ##' * \code{$stop} - logical value indicating whether stop is indicated
                 ##' * \code{$mtd} - integer value, the recommended dose
                 ##' * \code{$max_dose} - integer value, a dose not to be exceeded henceforth.
                 applied = function(x, o, last_dose, max_dose, ...){
                   tox <- x[last_dose]
                   n <- tox + o[last_dose]
                   rec <- if (tox >= private$eliminate[n])
                            max_dose <- last_dose-1
                          else if (tox >= private$deescalate[n])
                            last_dose-1
                          else if (tox <= private$escalate[n])
                            min(last_dose+1, max_dose)
                          else
                            last_dose

                   stop <- ( rec == 0 ||
                             (x+o)[rec] >= private$cohort_max ||
                             sum(x+o) >= private$enroll_max )

                   return(list(mtd = rec,
                               stop = stop,
                               max_dose = max_dose))
                 } #</applied>
               ), # </public>
               private = list(
                 escalate = NA
               , deescalate = NA
               , eliminate = NA
               , cohort_max = NA
               , enroll_max = NA
               )
               )
