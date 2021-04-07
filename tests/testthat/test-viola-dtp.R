library(dtpcrm)

test_that(".conduct_dose_finding_cohorts agrees with dtpcrm's version", {
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

  crm <- Crm$new(skeleton = prior.DLT,
                 scale = sqrt(prior.var),
                 target = target.DLT)
  crm$stop_func(stop_func)
  crm$no_skip_esc(TRUE)$no_skip_deesc(FALSE)$global_coherent_esc(TRUE)

  ##path <- c(0, 0, 0, 1, 3, 2, 0)
  ## Let's look for a path with early stopping ...
  path <- c(0, 0, 0, 1, 3, 2, 0)
  cohort_sizes = rep(3, length(path))
  old_recs <- dtpcrm:::.conduct_dose_finding_cohorts(next_dose = 3,
                                                     tox_counts = path,
                                                     cohort_sizes = cohort_sizes,
                                                     prior = prior.DLT,
                                                     target = target.DLT,
                                                     scale = sqrt(prior.var),
                                                     no_skip_esc = TRUE,
                                                     no_skip_deesc = FALSE,
                                                     global_coherent_esc = TRUE
                                                     )
  new_recs <- .conduct_dose_finding_cohorts(next_dose = 3,
                                            tox_counts = path,
                                            cohort_sizes = cohort_sizes,
                                            dose_func = crm$applied,
                                            impl = "rusti"
                                            )

  expect_equal(new_recs, old_recs)
})

test_that("calculate_dtps() yields same VIOLA result as dtpcrm's version", {
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

  crm <- Crm$new(skeleton = prior.DLT,
                 scale = sqrt(prior.var),
                 target = target.DLT)$
    stop_func(stop_func)$
    no_skip_esc(TRUE)$
    no_skip_deesc(FALSE)$
    global_coherent_esc(TRUE)

  new <- calculate_dtps(
    next_dose = start.dose.level,
    cohort_sizes = rep(3, 7),
    dose_func = crm$applied,
    impl = 'rusti')

  data(viola_dtp) # saved for comparison

  attr(new,'performance') <- NULL # drop incomparable attribute
  expect_equal(new, as.data.table(viola_dtp))
})
