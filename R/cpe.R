## An R6 class for Complete Path Enumeration (CPE), to be regarded as an
## *abstract* class notwithstanding R6's lack of support for this idea.

#' @name Cpe-class
#' @title An R6 base class for designs requiring Complete Path Enumeration
#'
#' @details
#' TODO: Lots to be said here! The fundamental CPE concept has to be stated,
#' and implementation strategy discussed as well. If it turns out this class
#' ends up providing hooks for cacheing, this also needs to be detailed.
#' @importFrom R6 R6Class
#' @export
Cpe <- R6Class("Cpe",
               public = list(
                 #' @details
                 #' Set or query upper limit on further dosing
                 #' @note This state may in general decrease along a trial path,
                 #' as for example with BOIN's dose elimination. (Designs which
                 #' do not implement this safety feature may clarify this fact
                 #' by overriding (TODO: or deleting?) this method.
                 #'
                 #' @param D A positive integer, the highest permissible dose.
                 #' @return Self (invisibly), unless `D` is missing,
                 #' in which case the top dose, an integer, is returned.
                 max_dose = function(D) {
                   if (missing(D))
                     return(private$.max_dose)
                   private$.max_dose <- D
                   invisible(self)
                 },
                 #' @details
                 #' Get the `b` vector and `U` matrix
                 #' @return Named list with components `b` and `U`
                 bU = function() {
                   if (is.null(private$b) || is.null(private$U))
                     self$path_array()
                   list(b = private$b, U = private$U)
                 },
                 #' @details
                 #' Get the number `J` of paths
                 #' @return Integer number of paths
                 J = function() {
                   length(private$path_list)
                 },
                 #' @details
                 #' Return dose recommendation for given tox/no-tox tallies.
                 #' @note Concrete subclasses must implement this abstract method.
                 #' @param x A dose-wise vector of toxicity counts
                 #' @param o A dose-wise vector of non-toxicity counts
                 #' @param last_dose The most recently given dose, as required to implement
                 #' cumulative-cohort-based escalation decisions.
                 #' @param max_dose An upper limit on future dose levels
                 #' @return An object with components:
                 #' * `$stop` - logical value indicating whether stop is indicated
                 #' * `$mtd` - integer value, the recommended dose
                 #' * `$max_dose` - integer value, a dose not to be exceeded henceforth.
                 applied = function(x, o, last_dose, max_dose){
                   stop("Must override $applied() method of abstract class 'Cpe'")
                 }, #</applied>
                 #' @details
                 #' Hook for concrete subclasses to implement for performance reporting
                 #' from method `Cpe$trace_paths`
                 #' @param ... Optional columns to add to report
                 #' @return NULL
                 report = function(...) {
                   return(NULL)
                 },
                 #' @field performance A vector used for vetting performance
                 performance = NULL,
                 #' @details
                 #' Compute trial paths forward from current tally
                 #'
                 #' The computed paths are saved in a private field, from which variously
                 #' formatted results may be obtained by accessor functions.
                 #'
                 #' @param root_dose The starting dose for tree of paths
                 #' @param cohort_sizes Integer vector giving sizes of future cohorts,
                 #' its length being the maximum number of cohorts to look ahead.
                 #' @param ... Parameters passed ultimately to `mclapply`, presently
                 #' an unexported, specially adapted version of `parallel::mclapply`
                 #' that implements progress reporting.
                 #' @param prog A function of a single integer, the current cumulative
                 #' path count, to be used for progress reporting
                 #' @param unroll Integer; how deep to unroll path tree for parallelism
                 #' @return Self, invisibly
                 #' @seealso `path_matrix`, `path_table`, `path_array`.
                 #' @note If the `parallel` package were to incorporate the necessary
                 #' changes to `mclapply`, I could restore the following import!
                 trace_paths = function(root_dose, cohort_sizes, ..., prog = NULL, unroll = 4){
                   stopifnot("Only constant cohorts_sizes are supported currently" =
                               length(unique(cohort_sizes))==1)
                   ## NB: The reason cohort_sizes must be a constant vector is that
                   ##     $path_array() uses a constant n = max(T) in choose(n,T).
                   private$T <- private$b <- private$U <- NULL # force recalc by $path_array()
                   stopifnot(unroll > 0) # TODO: Handle unroll=0 case gracefully!
                   paths. <- function(n, x, coh, path_m, cohort_sizes, par_t0 = NA){
                     t1. <- Sys.time()
                     path_hash <- new.env(hash = TRUE, size = 100000L) # to collect paths
                     paths_ <- function(n, x, coh, path_m, max_dose, cohort_sizes){
                       ## This recursive accessory function manages PATH BRANCHING at the
                       ## crucial point of INDETERMINACY, where toxicity counts are observed.
                       ## It is called for its SIDE-EFFECT on statically scoped path_hash.
                       d <- path_m["D",coh]
                       ## Handle the terminal case upon entry, simplifying the code to follow.
                       if (coh > length(cohort_sizes) || is.na(d) || d <= 0L) {
                         tox_c <- path_m["T",]
                         key <- paste(tox_c[!is.na(tox_c)], collapse='.')
                         assign(key, as.vector(path_m), envir = path_hash)
                         return()
                       }
                       n[d] <- n[d] + cohort_sizes[coh]
                       x_d <- x[d]
                       for (ntox in 0L:cohort_sizes[coh]) {
                         path_m["T", coh] <- ntox
                         x[d] <- x_d + ntox
                         rec <- self$applied(x=x, o=n-x, last_dose=d, max_dose=max_dose)
                         ## we recommend dose <= 0 to signal stopped path:
                         path_m["D", coh+1] <- ifelse(rec$stop, -1L, 1L) * rec$mtd
                         if (rec$stop)
                           paths_(n, x, coh+1, path_m, rec$max_dose, cohort_sizes[1:coh])
                         else
                           paths_(n, x, coh+1, path_m, rec$max_dose, cohort_sizes)
                       }
                     } #</paths_>
                     paths_(n, x, coh, path_m, self$max_dose(), cohort_sizes)
                     ## Regarding performant nature of the following as.list(env), see
                     ## https://stackoverflow.com/a/29482211/3338147 by Gabor Csardi.
                     path_list <- as.list(path_hash, sorted = FALSE)
                     t2. <- Sys.time()
                     ## Now render t1. and t2. relative to a single event (parallelization)
                     t1 <- 1000*difftime(t1., par_t0, units = "secs") # ~ started
                     attributes(t1)$units = "ms"
                     t2 <- 1000*difftime(t2., par_t0, units = "secs") # ~ finished
                     attributes(t2)$units = "ms"
                     attr(path_list,'performance') <- self$report(J = length(path_list)
                                                                , t1 = round(t1)
                                                                , t2 = round(t2)
                                                                  ## NB: U+0394 is LETTER Delta,
                                                                  ## but U+2206 is needed SYMBOL.
                                                                , "\u2206t" = round(t2 - t1)
                                                                  )
                     return(path_list)
                   } #</paths.>
                   ## 'Unroll' the first few levels of the tree recursion..
                   path_m <- matrix(NA_integer_, nrow=2, ncol=1+length(cohort_sizes),
                                    dimnames=list(c("D","T")))
                   n <- x <- integer(self$max_dose())
                   path_m["D",1] <- as.integer(root_dose)
                   ppe <- paths.(n, x, 1, path_m, cohort_sizes[1:unroll])
                   ppe <- ppe[order(names(ppe))]
                   ## ..set aside any stopped paths..
                   stopped <- sapply(ppe, function(ppe_) {
                     path_m <- matrix(ppe_, nrow=2, dimnames=list(c("D","T")))
                     any(path_m["D",] <= 0, na.rm=TRUE)
                   })
                   cpe_stopped <- ppe[stopped]
                   ppe <- unname(ppe[!stopped]) # NB: unname avoids prefix dup by mclapply
                   ## ..and parallelize over the pending partial paths:
                   par_t0 <- Sys.time() # to report thread-wise times since parallelization
                   FUN <- function(ppe_) {
                     path_m <- matrix(ppe_, nrow=2, dimnames=list(c("D","T")))
                     level <- factor(path_m["D",1:unroll], levels=1:self$max_dose())
                     tox <- path_m["T",1:unroll]
                     enr <- cohort_sizes[1:unroll]
                     n <- as.vector(xtabs(enr ~ level))
                     x <- as.vector(xtabs(tox ~ level))
                     paths.(n, x, unroll+1, path_m, cohort_sizes, par_t0 = par_t0)
                   }
                   if (!is.null(prog))
                     cpe_parts <- mclapply(ppe, FUN, proginit = sum(stopped), progreport = prog, ...)
                   else
                     cpe_parts <- mclapply(ppe, FUN, ...)
                   cpe <- c(cpe_stopped, do.call(c, cpe_parts))
                   private$path_list <- cpe[order(names(cpe))]
                   self$performance <- rbindlist(lapply(cpe_parts, attr, which='performance')
                                               , idcol = 'job')
                   invisible(self)
                 }, # </trace_paths>
                 #' @details
                 #' Return computed trial paths in matrix form
                 #'
                 #' @return An integer matrix with the same column layout as the
                 #' DTP tables of \CRANpkg{dtpcrm}. That is, there is a D0 column
                 #' followed by paired Tc, Dc columns giving the toxicity count
                 #' for cohort c and the resulting dose recommendation *yielded by*
                 #' cohort c -- which is generally the recommendation *for* cohort
                 #' c+1.
                 #' @seealso `trace_paths`, which must already have been invoked
                 #' if this method is to return a meaningful result.
                 path_matrix = function() {
                   stopifnot(length(private$path_list) > 0)
                   Ncol <- max(sapply(private$path_list, length)) - 1L
                   pmx <- matrix(NA_integer_,
                                 nrow=length(private$path_list),
                                 ncol=Ncol)
                   for (j in 1L:nrow(pmx))
                     pmx[j,] <- private$path_list[[j]][1L:ncol(pmx)]
                   Cmax <- (Ncol - 1L)/2L
                   colnames(pmx) <- paste0(c("T","D"), rep(0:Cmax, each=2))[-1]
                   pmx
                 }, #</path_matrix>
                 #' @details
                 #' Return computed trial paths in a 3D array
                 #'
                 #' @param condense Logical value; if FALSE, the returned array has its
                 #' cohorts indexed trial-wise instead of dose-wise. This inflates the
                 #' array more than needed for the matrix computations it must support
                 #' (observe that in Norris2020c Eq. (4), the `c` index is eliminated
                 #' already by summation), but enables the sequence of events along a path
                 #' to be read off directly if this is required e.g. for visualization or
                 #' debugging. Default is TRUE.
                 #' @return For the `j`th path, the C*D matrix `T[j,,]` gives
                 #' the number of toxicities `T[j,c,d]` occurring in the `c`th
                 #' cohort for dose d. In case `condense=FALSE`, see above.
                 #'
                 #' @references
                 #' Norris DC. What Were They Thinking? Pharmacologic priors implicit in
                 #' a choice of 3+3 dose-escalation design. arXiv:2012.05301 \[stat.ME\].
                 #' December 2020. <https://arxiv.org/abs/2012.05301>
                 #' @seealso `trace_paths`, which must already have been invoked
                 #' if this method is to return a result.
                 path_array = function(condense=TRUE) {
                   t0 <- Sys.time()
                   if (is.null(private$T)) {
                     pmx <- self$path_matrix()
                     C <- (ncol(pmx) - 1L)/2L # max number of cohorts enrolled
                     D <- self$max_dose() # dose levels range 1:D
                     ## (Note how the stopping-rec dose<=0 idiom works smoothly here.)
                     I <- outer(pmx[,paste0("D",0:(C-1))], 1:D, FUN = "==")
                     I[I] <- 1L
                     I[!I] <- NA_integer_
                     ## Now I[,,] is a J*C*D indicator array that we may use
                     ## in the manner of a bitmask, to select the toxicities
                     ## into their proper positions within T[j,c,d]:
                     T <- I * outer(pmx[,paste0("T",1:C)], rep(1,D))
                     dimnames(T)[[2]] <- paste0("C",1:C) # cohort number
                     dimnames(T)[[3]] <- paste0("D",1:dim(T)[3]) # dose levels
                     ## Save {T, b, Y, U} to avoid costly recalculation
                     private$T <- T
                     n <- max(T, na.rm=TRUE) # TODO: Handle non-constant cohorts size
                     private$b <- apply(log(choose(n, T)), MARGIN = 1, FUN = sum, na.rm = TRUE)
                     Y <- apply(T, MARGIN = c(1,3), FUN = sum, na.rm = TRUE)
                     Z <- apply(n-T, MARGIN = c(1,3), FUN = sum, na.rm = TRUE)
                     private$U <- cbind(Y, Z)
                   }
                   t1 <- Sys.time()
                   Tdims <- dim(private$T)
                   cat(sprintf("T (%d\u00D7%d\u00D7%d) constructed in %5.3f sec\n",
                               Tdims[1], Tdims[2], Tdims[3], t1 - t0))
                   if (!condense)
                     return(private$T)
                   ## T is at this point larger (sparser) than necessary for the
                   ## matrix computations which it is intended to support.
                   ## Specifically, the cohorts are now indexed within the whole trial,
                   ## whereas they really need only to be indexed on a per-dose basis.
                   ## We next identify the largest number of non-NA entries in any
                   ## dose column of the array; this will be the new, smaller
                   ## C-dimension. Note that for the VIOLA trial example, this
                   ## condenses the C-dimension from 7 to 4. So a massive size
                   ## reduction is not necessarily achieved here, although with
                   ## the absolute size of these T arrays being potentially quite
                   ## large, a 50% reduction could be worthwhile.
                   c_dim <- max(apply(apply(private$T, MARGIN=c(1,3),
                                            FUN=function(x) sum(!is.na(x))),
                                      MARGIN=2, max))
                   dim_ <- dim(private$T); dim_[2] <- c_dim
                   T_ <- array(NA_integer_, dim = dim_)
                   for (j in 1:dim(private$T)[1]) {
                     for (d in 1:dim(private$T)[3]) {
                       x <- private$T[j,,d]
                       x <- x[!is.na(x)]
                       if (length(x) > 0)
                         T_[j,1:length(x),d] <- x
                     }
                   }
                   dimnames(T_)[[3]] <- dimnames(private$T)[[3]]
                   dimnames(T_)[[2]] <- paste0("c", 1:dim(T_)[2]) # NB: little-c
                   t2 <- Sys.time()
                   cat(sprintf("T condensed in %5.3f sec\n", t2 - t1))
                   T_
                 }, # </path_array>
                 #' @details
                 #' Path probabilities for given dose-wise DLT probabilities
                 #'
                 #' The design's paths must already have been completely enumerated by
                 #' `trace_paths`.
                 #' @param probs.DLT Numeric vector of DLT probabilities for the design's
                 #' prespecified doses.
                 #' @return A vector of probabilities for the enumerated paths
                 path_probs = function(probs.DLT) {
                   if (is.null(private$b) || is.null(private$U))
                     self$path_array()
                   log_p <- log(probs.DLT)
                   log_q <- log(1 - probs.DLT)
                   log_pi <- private$b + private$U %*% c(log_p, log_q)
                   as.vector(exp(log_pi))
                 }, # </path_probs>
                 #' @details
                 #' Vector of path-wise final dose recommendations
                 #'
                 #' The design's paths must already have been completely enumerated by
                 #' `trace_paths`. This method is a temporizing measure, bridging
                 #' to a new format for the return value of `path_matrix`, possibly
                 #' a `data.table` where the dose recs would be a column.
                 #' @return Integer vector of final dose recs on the enumerated paths
                 path_rx = function() {
                   pmx <- self$path_matrix()
                   rxcol <- max.col(!is.na(pmx), ties.method="last")
                   rxdose <- abs(pmx[cbind(seq.int(nrow(pmx)),rxcol)])
                 } # </path_recs>
               ), # </public>
               private = list(
                 .max_dose = NA
               , path_list = list()
                 ## See Norris2020c for definitions of T, b, U
               , T = NULL # J*C*D array
               , b = NULL # J-vector
               , U = NULL # J*2D matrix
               )
               )
