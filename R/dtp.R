## Exploring speedups for DTP computations

##' Execute CRM
##'
##' Run a CRM trial under given design options for dose-escalation decisions.
##' This is identical to \code{\link[dtpcrm]{applied_crm}}, except that it uses
##' \code{precautionary}'s more performant version of \code{crm}.
##'
##' @param prior A vector of prior estimates of toxicity probabilties
##' for the dose levels
##' @param scale The sigma parameter of a lognormal prior (denoted \deqn{g(a)}
##' in the original 1990 O'Quigley et al. paper) over the single parameter
##' of the empiric or logistic CRM models.
##' @param target The target DLT rate
##' @param tox A patient-wise vector of 0/1 toxicity indicators
##' @param level A patient-wise vector of dose-level assignments
##' @param no_skip_esc If FALSE, the method will not enforce no skipping
##' of doses in escalation. Default is TRUE.
##' @param no_skip_deesc If FALSE, the method will not enforce no skipping
##' of doses in de-escalation. Default is TRUE.
##' @param global_coherent_esc If FALSE, the method will not enforce global
##' coherent escalation, that is, escalation if the overall rate of toxicity
##' seen at the current dose level is above the target rate. Default is TRUE.
##' @param stop_func An optional argument to provide a function which will
##' utilised alongside the CRM to determine if the trial should be stopped.
##' @param ... Additional params addressed to the CRM numerical routines.
##' Presently used to select optional implementations via \code{impl}.
##' @return An object of class \code{'mtd'}, as per \CRANpkg{dfcrm}
##' @export
applied_crm <- function (prior, scale, target, tox, level,
                         no_skip_esc = TRUE, no_skip_deesc = TRUE,
                         global_coherent_esc = TRUE, stop_func = NULL, ...)
{
  x <- Crm$new(skeleton = prior, scale = scale, target = target)$
    observe(level = level, tox = tox)$
    est(...) # ... contains impl=___
  if (no_skip_esc & x$mtd > (max(level) + 1)) {
    x$mtd <- max(level) + 1
  }
  if (no_skip_deesc & x$mtd < (min(level) - 1)) {
    x$mtd <- min(level) - 1
  }
  if (global_coherent_esc) {
    last_dose <- utils::tail(level, 1)
    tox_rate_last_dose <- sum(tox[level == last_dose])/sum(level == last_dose)
    if (tox_rate_last_dose > target) {
      x$mtd <- min(x$mtd, last_dose)
    }
  }
  if (!is.null(stop_func)) {
    x = stop_func(x)
  }
  return(x)
}

## Adapted from package:dtpcrm.
## TODO: Refactor this as a method of class 'Crm'. (It is no longer feasible
##       to achieve high performance within the interface set forth in 'dtpcrm'.)
.conduct_dose_finding_cohorts <- function (next_dose, tox_counts, cohort_sizes,
                                           dose_func = applied_crm, ...)
{
  num_cohorts <- length(tox_counts)
  dose_recs <- integer(num_cohorts+1) # used as if a 0-based array
  dose_recs[1] <- next_dose
  D <- length(environment(dose_func)$private$ln_skel) # by the skin of my teeth!
  x <- o <- integer(D)
  for (coh in 1:num_cohorts) {
    ## Work out the recommendation for (coh+1)th cohort ...
    tox <- tox_counts[1:coh]
    level <- factor(dose_recs[1:coh], levels=1:D)
    d <- dose_recs[coh]
    x[d] <- x[d] + tox_counts[coh]
    o[d] <- o[d] + (cohort_sizes[coh] - tox_counts[coh])

    rec <- dose_func(x, o, last_dose = dose_recs[coh], ...)
    dose_recs[coh+1] <- rec$mtd
    if ("stop" %in% names(rec)) {
      if (rec[["stop"]]) {
        dose_recs[1 + coh:num_cohorts] = NA
        break
      }
    }
  }
  return(dose_recs[-1]) # accommodate the quirky aliasing of original dtpcrm code
}

##' Faster Dose Transition Pathways (DTP) calculation
##'
##' This function reimplements \code{dtpcrm::calculate_dtps} with algorithmic
##' improvements that achieve signficant speedups, rendering it suitable for
##' comprehensive enumeration of all paths of a trial.
##'
##' @details
##' In accord with the chosen DTP representation of \pkg{dtpcrm}, we pursue a tabular
##' listing of all paths, inclusive of the degeneracies engendered by early stopping.
##' During the DTP computation, however, these degeneracies are recognized and
##' handled efficiently here, avoiding costly reduplication of CRM model runs.
##' The algorithm employs a top-to-bottom scan of the path table, in the
##' lexicographic order of increasing toxicity counts which naturally arises from
##' \code{expand.grid} used here and in the original \code{dtpcrm::calculate_dtps}.
##' The early termination of any path in the course of this scan may then be carried
##' forward through all subsequent paths which share the same dose assignments and
##' toxicity counts. (In the case of early termination specifically for excessive
##' toxicity, this principle could also be applied 'one level up', possibly capturing
##' some additional efficiencies. But implementing this safely would require greater
##' control over early-termination semantics than \pkg{dtpcrm} --- or indeed perhaps
##' any R package --- seems capable of achieving. So it is not implemented here.)
##'
##' The substantial efficiencies (3.5x speedup of VIOLA trial DTP) obtained in the
##' above-described scan can readily be preserved under coarse-grained parallelization,
##' on account of the DTP table's recursive structure. If the first cohort has size \deqn{n},
##' then it gives rise to \deqn{n+1} independent chunks within which the scan logic applies,
##' one chunk for each possible toxicity count in \deqn{\{0,1,...,n\}}.
##' This principle applies also at the level of the next cohort, and so on. With each
##' such refinement to the chunking of the DTP calculation, some re-work is introduced
##' where the same early termination even is discovered independently in several chunks.
##' But given that few reasonable trial designs will have a substantial number of paths
##' terminating before 6 patients, this re-work should be irrelevant to parallelization
##' of DTP computation over the 4--16 threads of a single modern, multicore processor.
##' (Consider that chunking just the first 2 cohorts of a trial with cohort size 3
##' yields $(3+1)^2=16$ chunks.)
##'
##' @param next_dose The root dose of the trial-path tree to be computed
##' @param cohort_sizes An integer vector of cohort sizes, with length equal to
##' the maximum number of cohorts that may be enrolled. (Thus, \code{sum(cohort_sizes)}
##' is the trial's maximum enrollment.)
##' @param dose_func An adapter function
##' @param ... Ultimately passed to \code{.conduct_dose_finding_cohorts}, whence
##' it ends up informing the \code{dose_func}. Thanks to the method chaining of
##' R6 class \code{Crm}, the role of this argument has been sharply curtailed
##' compared with the original \code{dtpcrm::calculate_dtps}. Here, we use it
##' only for selecting optional CRM numerical implementations via \code{impl}.
##' @param mc.cores Number of logical cores available for parallelizing DTP computation.
##' Setting this to 1 prevents chunking of the computation.
##' @return A \code{data.table} listing trial pathways
##' @examples
##' ## VIOLA trial
##' prior.DLT <- c(0.03, 0.07, 0.12, 0.20, 0.30, 0.40, 0.52)
##' prior.var <- 0.75
##'
##' stop_func <- function(x) {
##'   y <- stop_for_excess_toxicity_empiric(x,
##'                                         tox_lim = 0.3,
##'                                         prob_cert = 0.72,
##'                                         dose = 1)
##'   if(y$stop){
##'     x <- y
##'   } else {
##'     x <- stop_for_consensus_reached(x, req_at_mtd = 12)
##'   }
##' }
##'
##' crm <- Crm$new(skeleton = prior.DLT,
##'                scale = sqrt(prior.var),
##'                target = 0.2)$
##'   stop_func(stop_func)$
##'   no_skip_esc(TRUE)$
##'   no_skip_deesc(FALSE)$
##'   global_coherent_esc(TRUE)
##'
##' if (interactive()) { # don't overtax CRAN servers
##' restore <- options(mc.cores = 1) # you can reasonably set as high as detectCores()
##' dtps <- calculate_dtps(next_dose = 3,
##'                        cohort_sizes = rep(3, 7),
##'                        dose_func = crm$applied,
##'                        impl = "rusti")
##' options(mc.cores = restore)
##' }
##' @author Adapted by David C. Norris from original \code{dtpcrm::calculate_dtps}
##' @export
calculate_dtps <- function (next_dose, cohort_sizes,
                            dose_func = applied_crm, ...,
                            mc.cores = getOption("mc.cores", 2L))
{
  num_cohorts <- length(cohort_sizes)
  feasible_tox_counts <- lapply(cohort_sizes, function(x) 0:x)
  paths <- expand.grid(feasible_tox_counts)
  paths <- paths[do.call(order, as.data.frame(paths)), ]
  row.names(paths) <- 1:nrow(paths)
  ## Blindly applying .conduct_dose_finding_cohorts() to all rows
  ## of the 'paths' data frame duplicates effort, due to degeneracy
  ## from early-stopping. This degeneracy is 3.5x in the VIOLA trial,
  ## so substantial speedup may be attained by skipping over duplicate
  ## rows as they are encountered. Since the stopping criteria may
  ## depend on model outputs, this skipping must in general be decided
  ## on-the-fly. Thus, I will begin by converting this 'apply' to a
  ## for-loop...
  ## dtps <- apply(paths, 1, function(x)
  ##   .make_dtp_row(next_dose,
  ##                 x, dtpcrm:::.conduct_dose_finding_cohorts(next_dose, x, cohort_sizes,
  ##                                                           prev_tox = prev_tox,
  ##                                                           prev_dose = prev_dose,
  ##                                                           dose_func = dose_func,
  ##                                                           ...)))
  ## dtps <- t(dtps)
  scan_dtps <- function(paths, ...) {
    ## NB: Storing the paths in columns improves locality of reference
    dtps <- matrix(NA_integer_, ncol = nrow(paths),
                   nrow = 1+2*ncol(paths)) # preallocate matrix work area
    skipped <- 0 # to keep track of efficiencies
    i <- 0
    while (i < nrow(paths)) {
      i <- i + 1
      x <- as.integer(paths[i,])
      dose_recs <- .conduct_dose_finding_cohorts(next_dose, x, cohort_sizes,
                                                 dose_func = dose_func,
                                                 ...)
      dtps[,i] <- c(next_dose, c(rbind(x, dose_recs))) # NB: c(rbind(a,b)) interleaves a,b
      ## Here is the point where I might recognize an early-termination case,
      ## and propagate it forward (LOCF-like) through the current degeneracy.
      ## The early termination is recognizable from the index of the first NA
      ## value in dose_recs.
      etcoh <- match(NA, dose_recs)
      if (is.na(etcoh) || etcoh==length(dose_recs))
        next
      ## Upon reaching this point, we have detected an early termination (ET).
      degen <- min(prod(1 + cohort_sizes[-(1:etcoh)]) - 1,
                   ncol(dtps) - i) # avoid overrunning bounds
      dtps[,i + (1:degen)] <- dtps[,i]
      i <- i + degen
      skipped <- skipped + degen
    }
    dtps <- data.table(t(dtps))
    attr(dtps,'performance') <- environment(dose_func)$self$report()
    return(dtps)
  }
  chunks <- if (mc.cores == 1 || ncol(paths)<3)
              list(paths) # singleton 'chunk'
            else # TODO: Split into about mc.cores^2 chunks
              split(paths, paths[,3:1], drop=TRUE) # note reversed order of factor columns
  chunked_dtps <- parallel::mclapply(chunks, scan_dtps, ..., mc.cores = mc.cores)
  dtps <- rbindlist(chunked_dtps)
  colnames(dtps) <- c("D0", as.vector(rbind(paste0("T", 1:num_cohorts),
                                            paste0("D", 1:num_cohorts))))
  dtps[t(apply(is.na(dtps), 1, cumsum)) > 0] <- NA
  ## TODO: Make a nice data frame, possibly with PIDs as rownames.
  ##       Alternatively, consider a table or more sophisticated summary
  ##       that explores (e.g.) workload balance across forked processes.
  attr(dtps,'performance') <- do.call(rbind, lapply(chunked_dtps, attr, which='performance'))
  return(dtps)
}

##' A supremely faster version of a function from 'dtpcrm' v0.1.1
##'
##' Originally, the sampling in stats::rnorm() (see inline comments in code)
##' consumed 53% of run-time in a 6-cohort VIOLA DTP. After this change, it
##' doesn't even show up! More importantly, the consumption is now dominated
##' by (at 75%) by the objective function 'f' in integrate().
##' @param x A object of class \code{mtd}
##' @param tox_lim Scalar upper threshold on estimated toxicity rate
##' @param prob_cert Confidence level for threshold exceedance
##' @param dose Integer scalar, the dose being considered
##' @param suppress_dose Logical; if TRUE the MTD is set to \code{NA} when
##' trial stop is recommended.
##' @return
##' The \code{mtd} object x, with stop decision annotated
##' @author Adapted by David C. Norris from original dtpcrm::stop_for_excess_toxicity_empiric
##' @importFrom stats pnorm
##' @export
stop_for_excess_toxicity_empiric <- function (x, tox_lim, prob_cert, dose = 1,
                                              suppress_dose = TRUE) {
  post_beta_mean = x$estimate
  post_beta_var = x$post.var
  ## The following was massively wasteful, with the stats::rnorm()
  ## sampling itself consuming 53% of time in dtp(6)!
  ## > post_beta_samp = stats::rnorm(n = nsamps, mean = post_beta_mean,
  ## >                               sd = sqrt(post_beta_var))
  ## > post_prob_tox_samp = x$prior[dose]^exp(post_beta_samp)
  ## > prob_too_toxic = mean(post_prob_tox_samp > tox_lim)
  ## Yet a bit of math suffices to compute this with ZERO MCSE, via pnorm():
  prob_too_toxic <- pnorm(log(log(tox_lim)/log(x$prior[dose])),
                          mean = post_beta_mean, sd = sqrt(post_beta_var))
  stop_decision = prob_too_toxic > prob_cert
  x$stop = stop_decision
  if (stop_decision) {
    x$stop_reason = paste0("Prob(Prob(Tox[", dose, "]) > ",
                           tox_lim, ") = ", round(prob_too_toxic, 3), " > ",
                           prob_cert)
    if (suppress_dose)
      x$mtd = NA
  }
  return(x)
}

## summaryRprof()$by.total:
##                                          total.time total.pct self.time
## "calculate_dtps"                              55.88    100.00      0.20
## "dtp"                                         55.88    100.00      0.00
## "dtpcrm:::.conduct_dose_finding_cohorts"      53.44     95.63      0.58
## "dose_func"                                   52.52     93.99      0.48
## "dfcrm::crm"                                  49.46     88.51      1.20
## "integrate"                                   47.12     84.32      1.26
## ".External"                                   44.58     79.78      1.10
## "<Anonymous>"                                 43.64     78.10      1.44
## "f"                                           42.18     75.48     42.08
## "["                                            1.94      3.47      0.02
## "[.data.frame"                                 1.92      3.44      0.48
## "stop_func"                                    1.30      2.33      0.22
## "order"                                        0.98      1.75      0.34
## "[["                                           0.96      1.72      0.22
## "stop_for_excess_toxicity_empiric"             0.80      1.43      0.34

## TESTS

## The C=4 default yields 4^4 = 256 paths, with 1.54x degeneracy.
## This offers up plenty of speedup opportunity, while completing
## in just 10s even with the original dtpcrm version.
##' @importFrom dtpcrm stop_for_consensus_reached
##' @importFrom utils Rprof summaryRprof
dtp <- function(C=4, mc.cores=1, profile=FALSE, impl='rusti') {
  ## VIOLA trial set-up
  number.doses <- 7
  start.dose.level <- 3
  max.sample.size <- 21
  target.DLT <- 0.2
  cohort.size <- 3

  prior.DLT <- c(0.03, 0.07, 0.12, 0.20, 0.30, 0.40, 0.52)
  prior.var <- 0.75

  stop_func <- function(x) {
    y <- stop_for_excess_toxicity_empiric(x,
                                          tox_lim = target.DLT + 0.1,
                                          prob_cert = 0.72,
                                          dose = 1)
    if(y$stop){
      x <- y
    } else {
      x <- stop_for_consensus_reached(x, req_at_mtd = 12)
    }
  }

  crm <- Crm$new(skeleton = prior.DLT,
                 scale = sqrt(prior.var),
                 target = target.DLT)$
    stop_func(stop_func)$
    no_skip_esc(TRUE)$
    no_skip_deesc(FALSE)$
    global_coherent_esc(TRUE)

  ## Get timing
  if (profile)
    Rprof()
  else
    t0 <- proc.time()
  viola_dtp <- calculate_dtps(next_dose = start.dose.level,
                              cohort_sizes = rep(cohort.size, C),
                              dose_func = crm$applied,
                              mc.cores = mc.cores,
                              impl = impl
                              )
  if (profile) {
    Rprof(NULL)
    return(summaryRprof())
  } else {
    print(proc.time() - t0)
  }
  invisible(viola_dtp)
}
