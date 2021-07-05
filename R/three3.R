#' @name Cpe3_3-class
#' @title An R6 class encapsulating pre-calculated CPE for 3+3 design
#'
#' @details
#' TODO: Explain the hierarchy of model classes, including connections
#'       with the executable specifications set forth in exec/prolog/ccd.pl.
#' @references
#' 1. Korn EL, Midthune D, Chen TT, Rubinstein LV, Christian MC, Simon RM.
#'    A comparison of two phase I trial designs. *Stat Med*. 1994;13(18):1799-1806.
#'    doi:10.1002/sim.4780131802
#' @importFrom R6 R6Class
#' @export
Cpe3_3 <- R6Class("Cpe3_3",
               inherit = Cpe,
               public = list(
                 #' @details
                 #' Create a new `Cpe3_3` object.
                 #'
                 #' @param D number of pre-specified doses
                 #' @return A `Cpe3_3` object.
                 #'
                 #' @examples
                 #' # TODO
                 initialize = function(D) {
                   stopifnot("3+3 designs precalculated only for 2--8 doses" = D %in% 2:8)
                   private$.max_dose <- D
                   private$T <- precautionary:::T[[D]]
                   private$b <- precautionary:::b[[D]]
                   private$U <- precautionary:::U[[D]]
                 },
                 #' @details
                 #' Query number of doses
                 #' @note This specializes the generic method of Cpe superclass,
                 #' removing the possibility of updating the state.
                 #'
                 #' #param D A positive integer, the highest permissible dose.
                 #' @return Self (invisibly), unless `D` is missing,
                 #' in which case the top dose, an integer, is returned.
                 max_dose = function() {
                   #if (missing(D))
                   return(private$.max_dose)
                   #private$.max_dose <- D
                   #invisible(self)
                 },
                 #' @details
                 #' Get the `b` vector and `U` matrix
                 #' @return Named list with components `b` and `U`
                 bU = function() {
                   list(b = private$b, U = private$U)
                 },
                 #' @details
                 #' Get the number `J` of paths
                 #' @return Integer number of paths
                 J = function() {
                   length(private$b)
                 },
                 #' @details
                 #' No-op specialization of superclass method
                 #'
                 #' Since `Cpe3_3` has cached paths precomputed by Prolog code,
                 #' it does not need to support this method.
                 #'
                 #' @param ... Ignored
                 #' @return Self, invisibly
                 trace_paths = function(...){
                   message("No-op `Cpe3_3$trace_paths` invoked unnecessarily")
                   invisible(self)
                 },
                 #' @details
                 #' Refuse superclass method
                 #'
                 #' @return An error
                 path_matrix = function() {
                   stop("`Cpe3_3$path_matrix()` unimplemented")
                 },
                 #' @details
                 #' Refuse superclass method
                 #'
                 #' @return An error
                 path_rx = function() {
                   stop("`Cpe3_3$path_rx()` unimplemented")
                 }
               ) # </public>
               )
