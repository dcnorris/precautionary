library(dfcrm)
library(dtpcrm)
library(precautionary)

## To avoid a NOTE on package check, a dummy '...' formal argument is written.
## But this should not be taken as an invitation to provide actual arguments!
bt_crms <- function(...) {
  if (length(list(...))) stop("Don't invoke bt_crms with arguments; it's a trap!")
  benchtab(crm(prior=c(0.05, 0.10, 0.20, 0.35, 0.50, 0.70),
               target=0.2,
               tox=c(0L, 0L, 1L, 0L, 0L, 1L, 1L, 0L, 0L, 0L),
               level=c(3, 4, 4, 3, 3, 4, 3, 2, 2, 2),
               ...),
           impl=c("dfcrm","rusti"))
}

## This set-up is verbatim from the example in dfcrm::crm
comp_crms <- function(prior = c(0.05, 0.10, 0.20, 0.35, 0.50, 0.70),
                      target = 0.2,
                      level = c(3, 4, 4, 3, 3, 4, 3, 2, 2, 2),
                      y = c(0, 0, 1, 0, 0, 1, 1, 0, 0, 0),
                      impl = c("dfcrm","rusti")) {
  ## TODO: Get this working:
  ##comp_impl(quote(crm(prior, target, y, level, impl)), impl=impl)
  ## Until the above works, do it 'by hand':
  perf <- CJ(impl = impl, sorted = FALSE)
  for (i in seq(nrow(perf))) {
    crm_call <- substitute(crm(prior, target, y, level, impl=impl), perf[i,])
    perf$median.ms[i] <- median(microbenchmark::microbenchmark({
      eval(crm_call, perf[i,])
    })$time/1000000) # convert ns -> ms
    perf[i,'call'] <- deparse(crm_call)
  }
  perf
}

comp_icrm <- function(x=c(1,2,3,2,3,3)*0.1,
                      y=c(0L,0L,1L,0L,0L,1L),
                      w=rep(1,length(y)),
                      s=500) {
  mb_old <- microbenchmark::microbenchmark(
                              old0 <- integrate(dfcrm::crmh, -Inf, Inf, x, y, w, s)[[1]])
  old1 <- integrate(dfcrm::crmht, -Inf, Inf, x, y, w, s)[[1]]
  old2 <- integrate(dfcrm::crmht2, -Inf, Inf, x, y, w, s)[[1]]
  mb_rust <- microbenchmark::microbenchmark(rust0 <- icrm(x, y, w, s, 0))
  rust1 <- icrm(x, y, w, s, 1)
  rust2 <- icrm(x, y, w, s, 2)

  speedup <- mean(mb_old$time) / mean(mb_rust$time)
  cat(paste0("speedup: ", format(speedup, digits=2), "x\n"))

  ## TODO: Undertake a more comprehensive examination of convergence
  data.frame(old0 = old0, rust0 = rust0,
             old1 = old1, rust1 = rust1,
             old2 = old2, rust2 = rust2)
}

## TODO: Test in case where w not identically 1.
comp_crmh <- compiler::cmpfun(function(a=seq(-0.5, 0.5, 0.05),
                      x=c(1,2,3,2,3,3)*0.1,
                      y=c(0L,0L,1L,0L,0L,1L),
                      w=rep(1,length(y)), s=500) {
  old <- dfcrm::crmh(a,x,y,w,s)
  rust <- crmh(a,x,y,w,s)
  deltar <- max(abs(rust - old))
  if (deltar > .Machine$double.eps)
    cat("crmh |rust - old| =", deltar, "\n")
  ##
  old <- dfcrm::crmht(a,x,y,w,s)
  rust <- crmht(a,x,y,w,s)
  tdeltar <- max(abs(rust - old))
  if (tdeltar > .Machine$double.eps)
    cat("crmht |rust - old| =", tdeltar, "\n")
  ##
  old <- dfcrm::crmht(a,x,y,w,s)
  rust <- crmht(a,x,y,w,s)
  t2deltar <- max(abs(rust - old))
  if (t2deltar > .Machine$double.eps)
    cat("crmht2 |rust - old| =", t2deltar, "\n")

  t_old <- microbenchmark::microbenchmark(integrate(dfcrm::crmh, -Inf, Inf, x, y, w, s))
  t_rust <- microbenchmark::microbenchmark(integrate(crmh, -Inf, Inf, x, y, w, s))
  print(t_old)
  print(t_rust)
  speedup_message(t_old, t_rust)
  invisible(list(old=t_old, rust=t_rust))
})

## The defaults select the old routine, but skip.degenerate=TRUE engages other options
prof_viola_dtp <- function(skip.degenerate=FALSE, impl="dfcrm", mc.cores=1) {

  prior.DLT <- c(0.03, 0.07, 0.12, 0.20, 0.30, 0.40, 0.52)
  prior.var <- 0.75

  stop_func <- function(x) {
    y <- stop_for_excess_toxicity_empiric(x,
                                          tox_lim = 0.3,
                                          prob_cert = 0.72,
                                          dose = 1)
    if(y$stop){
      x <- y
    } else {
      x <- stop_for_consensus_reached(x, req_at_mtd = 12)
    }
  }

  system.time(
    if (!skip.degenerate) { # original dfcrm, but with my pnorm fix
      dtpcrm::calculate_dtps(next_dose = 3,
                             cohort_sizes = rep(3, 7),
                             prior = prior.DLT,
                             target = 0.2,
                             stop_func = stop_func,
                             scale = sqrt(prior.var),
                             no_skip_esc = TRUE,
                             no_skip_deesc = FALSE,
                             global_coherent_esc = TRUE)
    } else {
      calculate_dtps(next_dose = 3,
                     cohort_sizes = rep(3, 7),
                     prior = prior.DLT,
                     target = 0.2,
                     stop_func = stop_func,
                     scale = sqrt(prior.var),
                     no_skip_esc = TRUE,
                     no_skip_deesc = FALSE,
                     global_coherent_esc = TRUE,
                     impl = impl,
                     mc.cores = mc.cores)
    }
  )

}

viola_speedup_report <- function() {
  elapsed <- numeric(7)
  names(elapsed) <- c("pnorm", "skipt", "rusti", "core2", "core4", "core6", "core8")

  cat("Running original dtpcrm::calculate_dtps, but with rnorm->pnorm fix...\n")
  elapsed['pnorm'] <- prof_viola_dtp(skip=FALSE)['elapsed']
  cat("Elapsed: ", elapsed['pnorm'], "\n")

  cat("Running precautionary::calculate_dtps(impl = 'dfcrm') ...\n")
  elapsed['skipt'] <- prof_viola_dtp(skip=TRUE)['elapsed']
  cat("Elapsed: ", elapsed['skipt'], "\n")

  cat("Running precautionary::calculate_dtps(impl = 'rusti') ...\n")
  elapsed['rusti'] <- prof_viola_dtp(skip=TRUE, impl="rusti")['elapsed']
  cat("Elapsed: ", elapsed['rusti'], "\n")

  cat("Running precautionary::calculate_dtps(impl = 'rusti', mc.cores = 2) ...\n")
  elapsed['core2'] <- prof_viola_dtp(skip=TRUE, impl="rusti", mc.cores=2)['elapsed']
  cat("Elapsed: ", elapsed['core2'], "\n")

  cat("Running precautionary::calculate_dtps(impl = 'rusti', mc.cores = 4) ...\n")
  elapsed['core4'] <- prof_viola_dtp(skip=TRUE, impl="rusti", mc.cores=4)['elapsed']
  cat("Elapsed: ", elapsed['core4'], "\n")

  cat("Running precautionary::calculate_dtps(impl = 'rusti', mc.cores = 6) ...\n")
  elapsed['core6'] <- prof_viola_dtp(skip=TRUE, impl="rusti", mc.cores=6)['elapsed']
  cat("Elapsed: ", elapsed['core6'], "\n")

  cat("Running precautionary::calculate_dtps(impl = 'rusti', mc.cores = 8) ...\n")
  elapsed['core8'] <- prof_viola_dtp(skip=TRUE, impl="rusti", mc.cores=8)['elapsed']
  cat("Elapsed: ", elapsed['core8'], "\n")

  elapsed
}

## Very interesting results!
##   pnorm   skipt   rusti   core2   core4   core6   core8
## 110.639  82.744  40.744  14.284  10.842   8.531   8.468
##
## The skipping of a 71% degeneracy apparently yields only 110/83 = 1.3x speedup,
## not the 1/0.29 = 3.5x speedup expected from avoiding 71% of the effort!
## The transition from mc.cores=1 to mc.cores=2 however provides a strong hint
## as to the cause of this. The 40./14.3 = 2.85x improvement here (where at most
## 2x would be anticipated) shows that the chunking (suppressed in case mc.cores=1)
## itself accomplishes a great deal. This suggests that my earlier 'rolling' while
## loop may have been the more efficient way to traverse the degeneracy, after all!
## ...
## More learning...
## A more efficient loop (with en bloc degeneracy skip) achieves this:
##   pnorm   skipt   rusti   core2   core4   core6   core8
## 107.239  59.965  21.064  13.445  10.223   8.081   7.945
##
## Now the improvement is 107/60 = 1.8x, and deeper thinking about the algorithm
## reveals why this should not be 3.5x. The more highly degenerate paths have
## fewer crm() calls, so the savings are not proportional to path degeneracy.
## ...
## Now splitting more finely (on first 3 columns) yields this:
##   pnorm   skipt   rusti   core2   core4   core6   core8
## 107.874  60.159  21.056  13.164   9.558   7.805   7.753
##
## Exactly as one might expect, the finer division helps squeeze extra value
## out of more cores.
##
## After switching to a scheme that saves N-D transcendental ops (on the screen,
## that is -- not sure what optimizing compiler might figure out), I obtain
## essentially unchanged results. It will be interesting to see what happens
## as I move the bookkeeping overhead into the caller, and then do the same
## with the D X[d].ln()'s.
##   pnorm   skipt   rusti   core2   core4   core6   core8
## 106.571  57.845  20.872  12.853   9.161   7.789   7.947
##
## Implementing a D+1 txop routine (w==1 case) with pre-bookkeeping, I obtain:
##   pnorm   skipt   rusti   core2   core4   core6   core8
## 105.905  58.509  15.275  10.068   6.656   5.266   5.623


## Calculate potential savings (%) in the cascade
## of function calls during VIOLA DTP
find_dtp_savings <- function(impl="rusti") {
  prof <- prof_viola_dtp(impl)$by.total
  ##prof <- rusti$by.total
  nuisance <- rownames(prof) %in% c('"find_dtp_savings"',
                                    '"prof_viola_dtp"')
  prof <- prof[!nuisance, ]
  savings <- as.data.table(-diff(as.matrix(prof)), keep.rownames="function")
  ## Drop irrelevant columns
  savings[, total.time := NULL]
  savings[, self.time := NULL]
  savings[, self.pct := NULL]
  savings <- savings[1:10,]
  savings
}

## Maybe I'd be better off summarizing self.time across impls?
compare_dtp_effort <- function() {
  dfcrm <- as.data.table(dfcrm$by.self, keep.rownames="routine")[,c('routine','self.time')]
  rusti <- as.data.table(rusti$by.self, keep.rownames="routine")[,c('routine','self.time')]
  setnames(dfcrm, "self.time", "dfcrm")
  setnames(rusti, "self.time", "rusti")
  merged <- merge(dfcrm, rusti, by="routine", all=TRUE)
  merged$routine <- gsub("\"", "", merged$routine, fixed=TRUE)
  ## Calculate a margin row at top
  totals <- colSums(merged[,-1], na.rm=TRUE)
  merged <- rbind(merged[1,], merged)
  merged[1, routine := "Total"]
  merged[1, names(totals)] <- as.list(totals)
  merged[, worst := pmax(dfcrm, rusti, na.rm=TRUE)]
  merged <- merged[order(-worst),]
  merged[, worst := NULL]
  merged
}

## From this effort comparison, I gather a few conclusions:
## (a) Conversion of objective function 'f' in integrate(f...)
##     accounts for the vast majority of the speedup.
## (b) Together f+.External+integrate=3.00 in rusti, nearly
##     cancelling out its .Call advantage of 11.12-7.38=3.74
##     relative to rustq.
## (c) rustq's further 1.04 advantage on <Anonymous> then
##     puts it in the lead, which it never loses to rusti's
##     miscellaneous accumulated R-related costs.
## (d) Given how the costs are spread out, no remaining routine
##     will reward effort invested in performance tuning.
## (e) Given how close rustq is to rusti, there seems to be
##     little reason to export it at this time.
##
##                              routine dfcrm rusti rustq
##  1:                            Total 56.68 20.16 17.38
##  2:                                f 43.64  0.56    NA
##  3:                            .Call    NA  7.38 11.12
##  4:                        .External  1.78  1.16    NA
##  5:                        integrate  1.20  1.28    NA
##  6:                      <Anonymous>  1.20  1.12  0.08
##  7:                              crm  0.20  1.10  0.70
##  8:                       dfcrm::crm  0.78    NA    NA
##  9:                        dose_func  0.42  0.52  0.58
## 10: stop_for_excess_toxicity_empiric  0.44  0.30  0.18
## 11:                               [[  0.42  0.36  0.16
## 12:                     [.data.frame  0.40  0.36  0.24
## 13:                       match.call  0.36  0.38    NA
## 14:                        stopifnot  0.38  0.32    NA
## 15:       stop_for_consensus_reached  0.12  0.26  0.36
## 16:                                $  0.34  0.18  0.24
## 17:    .conduct_dose_finding_cohorts  0.22  0.26  0.34
## 18:                            order  0.30  0.22  0.14
## 19:                    [[.data.frame  0.26  0.26  0.22
## 20:                             %in%  0.12  0.24  0.18
## 21:                              all  0.14  0.24  0.18
## 22:                   calculate_dtps  0.24  0.14  0.16
## 23:                         sys.call  0.18  0.24  0.08
## 24:                      utils::tail  0.10  0.20  0.08
## 25:                           length  0.18  0.02  0.08
## 26:                        match.arg  0.08  0.18  0.08
## 27:                     sys.function  0.12  0.18  0.04
## 28:                            pnorm  0.10  0.04  0.16
## 29:                        stop_func  0.14  0.16  0.04
## 30:                           vapply  0.16  0.14  0.04

## Let's profile the DTP computations to find out why degeneracy skipping
## doesn't yield speedup of 1/degeneracy. (I suspect it's the while loop.)
## So let's compare Rprof output for the two routines
prof_skipping <- function() {
  ## VIOLA trial set-up
  start.dose.level <- 3
  target.DLT <- 0.2

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

  skip_prof <- list()

  restore <- options(mc.cores = 1)

  Rprof()
  dtpcrm::calculate_dtps(
    next_dose = start.dose.level,
    cohort_sizes = rep(3, 7),
    prior = prior.DLT,
    target = target.DLT,
    stop_func = stop_func,
    scale = sqrt(prior.var),
    no_skip_esc = TRUE,
    no_skip_deesc = FALSE,
    global_coherent_esc = TRUE)
  Rprof(NULL)
  skip_prof$all <- summaryRprof()

  Rprof()
  calculate_dtps(
    next_dose = start.dose.level,
    cohort_sizes = rep(3, 7),
    prior = prior.DLT,
    target = target.DLT,
    stop_func = stop_func,
    scale = sqrt(prior.var),
    no_skip_esc = TRUE,
    no_skip_deesc = FALSE,
    global_coherent_esc = TRUE,
    impl = 'dfcrm') # avoid advantage here vs dtpcrm original
  Rprof(NULL)
  skip_prof$skip <- summaryRprof()

  options(mc.cores = restore)

  skip_prof
}
