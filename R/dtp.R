## Exploring speedups for DTP computations

##' Faster Dose Transition Pathways (DTP) calculation
##'
##' .. content for \details{} ..
##' @title
##' @param next_dose
##' @param cohort_sizes
##' @param prev_tox
##' @param prev_dose
##' @param dose_func
##' @param ...
##' @return
##' @author Adapted by David C. Norris from original dtpcrm::calculate_dtps
calculate_dtps <- function (next_dose, cohort_sizes,
                            prev_tox = c(), prev_dose = c(),
                            dose_func = applied_crm, ...)
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
  dtps <- matrix(99, nrow = nrow(paths),
                 ncol = 1+2*ncol(paths)) # preallocate the matrix
  skipped <- 0 # to keep track of efficiencies
  i <- 0
  while (i < nrow(paths)) {
    i <- i + 1
    x <- as.integer(paths[i,])
    dose_recs <- dtpcrm:::.conduct_dose_finding_cohorts(next_dose, x, cohort_sizes,
                                                        prev_tox = prev_tox,
                                                        prev_dose = prev_dose,
                                                        dose_func = dose_func,
                                                        ...)
    dtps[i,] <- c(next_dose, c(rbind(x, dose_recs))) # NB: interleaves x, dose_recs
    ## Here is the point where I might recognize an early-termination case,
    ## and propagate it forward (LOCF-like) through the current degeneracy.
    ## The early termination is recognizable from the index of the first NA
    ## value in dose_recs.
    etcoh <- match(NA, dose_recs)
    if (is.na(etcoh))
      next
    ## Upon reaching this point, we have detected an early termination (ET).
    while (i < nrow(paths) &&
           all(as.integer(paths[i+1, 1:etcoh]) == x[1:etcoh])) {
             dtps[i+1,] <- dtps[i,]
             i <- i + 1
             skipped <- skipped + 1
           }
  }
  cat("skipped =", skipped, paste0("(", round(100*skipped/nrow(paths)), "%)"), "\n")
  dtps <- data.frame(dtps)
  colnames(dtps) <- c("D0", as.vector(rbind(paste0("T", 1:num_cohorts),
                                            paste0("D", 1:num_cohorts))))
  dtps[t(apply(is.na(dtps), 1, cumsum)) > 0] <- NA
  return(dtps)
}

##' An infinitely faster version of a function from 'dtpcrm' v0.1.1
##'
##' Originally, the sampling in stats::rnorm() (see inline comments in code)
##' consumed 53% of run-time in a 6-cohort VIOLA DTP. After this change, it
##' doesn't even show up! More importantly, the consumption is now dominated
##' by (at 75%) by the objective function 'f' in integrate().
##' @title
##' @param x
##' @param tox_lim
##' @param prob_cert
##' @param dose
##' @param nsamps
##' @param suppress_dose
##' @return
##' @author Adapted by David C. Norris from original dtpcrm::stop_for_excess_toxicity_empiric
stop_for_excess_toxicity_empiric <- function (x, tox_lim, prob_cert, dose = 1,
                                              nsamps = 10^6, suppress_dose = TRUE) {
  post_beta_mean = x$estimate
  post_beta_var = x$post.var
  ## The following was massively wasteful, with the stats::rnorm()
  ## sampling itself consuming 53% of time in dtp(6)!
  ## > post_beta_samp = stats::rnorm(n = nsamps, mean = post_beta_mean,
  ## >                               sd = sqrt(post_beta_var))
  ## > post_prob_tox_samp = x$prior[dose]^exp(post_beta_samp)
  ## > prob_too_toxic = mean(post_prob_tox_samp > tox_lim)
  ## Yet a bit of math suffices to compute this exactly via pnorm():
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
dtp <- function(C=4) {
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
                                          dose = 1,
                                          nsamps = 100000)
    if(y$stop){
      x <- y
    } else {
      x <- stop_for_consensus_reached(x, req_at_mtd = 12)
    }
  }

  ## Get timing
  t0 <- proc.time()
  viola_dtp <- calculate_dtps(next_dose = start.dose.level,
                              cohort_sizes = rep(cohort.size, C),
                              dose_func = applied_crm,
                              prior = prior.DLT,
                              target = target.DLT,
                              stop_func = stop_func,
                              scale = sqrt(prior.var),
                              no_skip_esc = TRUE,
                              no_skip_deesc = FALSE,
                              global_coherent_esc = TRUE
                              )
  print(proc.time() - t0)
  invisible(viola_dtp)
}
