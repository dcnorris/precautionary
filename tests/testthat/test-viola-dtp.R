library(dtpcrm)

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

  timings <- list(
    dtpcrm = system.time(
      old <- dtpcrm::calculate_dtps(
                       next_dose = start.dose.level,
                       cohort_sizes = rep(3, 7),
                       prior = prior.DLT,
                       target = target.DLT,
                       stop_func = stop_func,
                       scale = sqrt(prior.var),
                       no_skip_esc = TRUE,
                       no_skip_deesc = FALSE,
                       global_coherent_esc = TRUE)
    )
  , newdtp = system.time(
      new <- calculate_dtps(
        next_dose = start.dose.level,
        cohort_sizes = rep(3, 7),
        dose_func = applied_crm, # i.e., precautionary::applied_crm
        prior = prior.DLT,
        target = target.DLT,
        stop_func = stop_func,
        scale = sqrt(prior.var),
        no_skip_esc = TRUE,
        no_skip_deesc = FALSE,
        global_coherent_esc = TRUE,
        impl = 'rusti')
    )
  )

  with(timings, {
    speedup_message(newdtp, dtpcrm)
  })

  rownames(new) <- rownames(old) # don't compare rownames
  expect_equal(old, new)
})
