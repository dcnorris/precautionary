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
                 ##' Set or query maximum dose
                 ##' @note This private state may decrease along a trial path,
                 ##' due to dose elimination.
                 ##'
                 ##' @param D A positive integer, the highest permissible dose.
                 ##' @return Self (invisibly), unless \code{D} is missing,
                 ##' in which case the top dose, an integer, is returned.
                 max_dose = function(D) {
                   if (missing(D))
                     return(private$.max_dose)
                   private$.max_dose <- D
                   invisible(self)
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
                 }, #</applied>
                 ##' @details
                 ##' Do-nothing implementation to preserve Crm$trace_paths
                 ##' @return NULL
                 report = function() {
                   return(NULL)
                 },
                 ##' @field performance A vector used for vetting performance
                 performance = NULL,
                 ##' @details
                 ##' Compute trial paths forward from current tally
                 ##'
                 ##' The computed paths are saved in a private field, from which variously
                 ##' formatted results may be obtained by accessor functions.
                 ##'
                 ##' @param root_dose The starting dose for tree of paths
                 ##' @param cohort_sizes Integer vector giving sizes of future cohorts,
                 ##' its length being the maximum number of cohorts to look ahead.
                 ##' @param ... Parameters passed ultimately to \code{Crm$esc()},
                 ##' enabling passthru of its required \code{impl} parameter.
                 ##' @param unroll Integer; how deep to unroll path tree for parallelism
                 ##' @return Self, invisibly
                 ##' @seealso \code{path_matrix}, \code{path_table}, \code{path_array}.
                 ##' @importFrom parallel mclapply
                 trace_paths = function(root_dose, cohort_sizes, ..., unroll = 4){
                   stopifnot(unroll > 0) # TODO: Handle unroll=0 case gracefully!
                   paths. <- function(n, x, coh, path_m, cohort_sizes){
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
                         rec <- self$applied(x=x, o=n-x, last_dose=d, max_dose=max_dose, ...)
                         ## we recommend dose <= 0 to signal stopped path:
                         path_m["D", coh+1] <- ifelse(rec$stop, -1L, 1L) * rec$mtd
                         if (rec$stop)
                           paths_(n, x, coh+1, path_m, rec$max_dose, cohort_sizes[1:coh])
                         else
                           paths_(n, x, coh+1, path_m, rec$max_dose, cohort_sizes)
                       }
                     } #</paths_>
                     paths_(n, x, coh, path_m, private$.max_dose, cohort_sizes)
                     ## Regarding performant nature of the following as.list(env), see
                     ## https://stackoverflow.com/a/29482211/3338147 by Gabor Csardi.
                     path_list <- as.list(path_hash, sorted = FALSE)
                     attr(path_list,'performance') <- self$report()
                     return(path_list)
                   } #</paths.>
                   ## 'Unroll' the first few levels of the tree recursion..
                   path_m <- matrix(NA_integer_, nrow=2, ncol=1+length(cohort_sizes),
                                    dimnames=list(c("D","T")))
                   max_dose <- private$.max_dose # copy to local var we may thread thru recursion
                   n <- x <- integer(max_dose)
                   path_m["D",1] <- as.integer(root_dose)
                   ppe <- paths.(n, x, 1, path_m, cohort_sizes[1:unroll])
                   ppe <- ppe[order(names(ppe))]
                   ## ..and parallelize over the pending partial paths:
                   cpe_parts <- mclapply(ppe, function(ppe_) {
                     path_m <- matrix(ppe_, nrow=2, dimnames=list(c("D","T")))
                     ## Let's be sure to skip the stopped paths, tho!
                     if (any(path_m["D",] <= 0, na.rm=TRUE))
                       return(list(ppe_))
                     level <- factor(path_m["D",1:unroll], levels=1:max_dose)
                     tox <- path_m["T",1:unroll]
                     enr <- cohort_sizes[1:unroll]
                     n <- as.vector(xtabs(enr ~ level))
                     x <- as.vector(xtabs(tox ~ level))
                     paths.(n, x, unroll+1, path_m, cohort_sizes)
                   }
                   , mc.preschedule = TRUE)
                   cpe <- do.call(c, cpe_parts)
                   ## attr(cpe,'performance') <- do.call(rbind, lapply(cpe_parts, attr,
                   ##                                                  which='performance'))
                   private$path_list <- cpe[order(names(cpe))]
                   self$performance <- rbindlist(lapply(cpe_parts, attr, which='performance'))
                   invisible(self)
                 }, # </trace_paths>
                 ##' @details
                 ##' Return computed trial paths in matrix form
                 ##'
                 ##' @return An integer matrix with the same column layout as the
                 ##' DTP tables of \CRANpkg{dtpcrm}. That is, there is a D0 column
                 ##' followed by paired Tc, Dc columns giving the toxicity count
                 ##' for cohort c and the resulting dose recommendation \emph{yielded by}
                 ##' cohort c -- which is generally the recommendation \emph{for} cohort
                 ##' c+1.
                 ##' @seealso \code{trace_paths}, which must already have been invoked
                 ##' if this method is to return a meaningful result.
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
                 ##' @details
                 ##' Return computed trial paths in a 3D array
                 ##'
                 ##' @param condense Logical value; if FALSE, the returned array has its
                 ##' cohorts indexed trial-wise instead of dose-wise. This inflates the
                 ##' array more than needed for the matrix computations it must support
                 ##' (observe that in Norris2020c Eq. (4), the \code{c} index is eliminated
                 ##' already by summation), but enables the sequence of events along a path
                 ##' to be read off directly if this is required e.g. for visualization or
                 ##' debugging. Default is TRUE.
                 ##' @return For the \code{j}'th path, the C*D matrix \code{T[j,,]} gives
                 ##' the number of toxicities \code{T[j,c,d]} occurring in the \code{c}'th
                 ##' cohort for dose d. In case \code{condense=FALSE}, see above.
                 ##'
                 ##' @references
                 ##' Norris DC. What Were They Thinking? Pharmacologic priors implicit in
                 ##' a choice of 3+3 dose-escalation design. arXiv:2012.05301 \[stat.ME\].
                 ##' December 2020. \url{https://arxiv.org/abs/2012.05301}
                 ##' @seealso \code{trace_paths}, which must already have been invoked
                 ##' if this method is to return a result.
                 path_array = function(condense=TRUE) {
                   pmx <- self$path_matrix()
                   C <- (ncol(pmx) - 1L)/2L # max number of cohorts enrolled
                   D <- private$.max_dose # dose levels range 1:D
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
                   ## We take a moment to save private b, Y and U calculated from T
                   private$b <- apply(log(choose(2, T)), MARGIN = 1, FUN = sum, na.rm = TRUE)
                   private$Y <- apply(T, MARGIN = c(1,3), FUN = sum, na.rm = TRUE)
                   Z <- apply(2-T, MARGIN = c(1,3), FUN = sum, na.rm = TRUE)
                   private$U <- cbind(private$Y, Z)
                   if (!condense)
                     return(T)
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
                   c_dim <- max(apply(apply(T, MARGIN=c(1,3), FUN=function(x) sum(!is.na(x))),
                                      MARGIN=2, max))
                   dim_ <- dim(T); dim_[2] <- c_dim
                   T_ <- array(NA_integer_, dim = dim_)
                   for (j in 1:dim(T)[1]) {
                     for (d in 1:dim(T)[3]) {
                       x <- T[j,,d]
                       x <- x[!is.na(x)]
                       if (length(x) > 0)
                         T_[j,1:length(x),d] <- x
                     }
                   }
                   dimnames(T_)[[3]] <- dimnames(T)[[3]]
                   dimnames(T_)[[2]] <- paste0("c", 1:dim(T_)[2]) # NB: little-c
                   T_
                 }, # </path_array>
                 ##' @details
                 ##' Path probabilities for given dose-wise DLT probs
                 ##'
                 ##' The design's paths must already have been completely enumerated by
                 ##' \code{trace_paths}, and the associated array \code{T}, matrix \code{U}
                 ##' and vector \code{b} computed by a call to \code{path_array}.
                 ##' @param probs.DLT Numeric vector of DLT probabilities for the design's
                 ##' pre-specified doses.
                 ##' @return A vector of probabilities for the enumerated paths
                 path_probs = function(probs.DLT) {
                   if (is.na(private$b) || is.na(private$U))
                     stop("Must invoke path_array() before log_pi()")
                   log_p <- log(probs.DLT)
                   log_q <- log(1 - probs.DLT)
                   log_pi <- private$b + private$U %*% c(log_p, log_q)
                   exp(log_pi)
                 }, # </path_probs>
                 ##' @details
                 ##' Vector of path-wise final dose recommendations
                 ##'
                 ##' The design's paths must already have been completely enumerated by
                 ##' \code{trace_paths}. This method is a temporizing measure, bridging
                 ##' to a new format for the return value of \code{path_matrix}, possibly
                 ##' a \code{data.table} where the dose recs would be a column.
                 ##' @return Integer vector of final dose recs on the enumerated paths
                 path_rx = function() {
                   pmx <- self$path_matrix()
                   rxcol <- max.col(!is.na(pmx), ties.method="last")
                   rxdose <- abs(pmx[cbind(seq.int(nrow(pmx)),rxcol)])
                 } # </path_recs>
               ), # </public>
               private = list(
                 escalate = NA
               , deescalate = NA
               , eliminate = NA
               , cohort_max = NA
               , enroll_max = NA
               , .max_dose = NA
               , user = c(dfcrm = 0.0,
                          rusti = 0.0,
                          ruste = 0.0)
               , path_list = list()
                 ## See Norris2020c for definitions of T, b, U, Y
               , T = NA # J*C*D array
               , b = NA # J-vector
               , U = NA # J*2D matrix
               , Y = NA # left half of U
               )
               )
