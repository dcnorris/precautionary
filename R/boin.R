##' @name Boin-class
##' @title The BOIN design as a Cumulative-Cohort Design (CCD) subclass
##'
##' @details
##' TODO: Provide references, and note clearly the choice of defaults.
##'       Also provide citations to BOIN paper(s).
##' @importFrom R6 R6Class
##' @importFrom BOIN get.boundary select.mtd
##' @export
Boin <- R6Class("Boin",
               inherit = Ccd,
               public = list(
                 ##' @details
                 ##' Create a new `Boin` object.
                 ##'
                 ##' @param target Target toxicity rate
                 ##' @param cohort_max Upper bound on dose-wise enrollment
                 ##' @param enroll_max Upper bound on total enrollment
                 ##' @return A Boin object.
                 ##'
                 ##' @examples
                 ##' # TODO
                 initialize = function(target, cohort_max, enroll_max) {
                   private$target <- target
                   bdy <- get.boundary(target = 0.25
                                      ,ncohort = enroll_max
                                      ,cohortsize = 1 # TODO: Allow larger cohorts?
                                      ,n.earlystop = cohort_max
                                      ,extrasafe = FALSE
                                      ,offset = 0.05
                                       )$boundary_tab
                   super$initialize(escalate = bdy[2,]
                                   ,deescalate = bdy[3,]
                                   ,eliminate = bdy[4,]
                                   ,cohort_max = cohort_max
                                   ,enroll_max = enroll_max
                                    )
                 } #</initialize>
                 ## TODO: Obtain the *final* BOIN dose recommendations
                 ##       via BOIN:select.mtd(). Note that this seems to
                 ##       require a general design making a distinction
                 ##       between dose-escalation decisions OT1H, and
                 ##       FINAL dose recommendations OTOH.
               ), # </public>
               private = list(
                 target = NA
               )
               )
