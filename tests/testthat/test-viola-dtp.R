library(dtpcrm)

test_that("calculate_dtps() yields same VIOLA result as dtpcrm's version", {
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

  timings <- list()
  dtps <- list()
  for (impl in c('dfcrm','rusti','rustq')) {
    timings[[impl]] <- system.time(
      dtps[[impl]] <- calculate_dtps(next_dose = start.dose.level,
                                     cohort_sizes = rep(cohort.size,
                                                        max.sample.size/cohort.size),
                                     dose_func = applied_crm,
                                     prior = prior.DLT,
                                     target = target.DLT,
                                     stop_func = stop_func,
                                     scale = sqrt(prior.var),
                                     no_skip_esc = TRUE,
                                     no_skip_deesc = FALSE,
                                     global_coherent_esc = TRUE,
                                     impl = impl)
    )
  }

  with(timings, {
    speedup_message(rusti, dfcrm)
    speedup_message(rustq, dfcrm)
  })

  with(dtps, {
    expect_equal(dfcrm, rusti)
    expect_equal(rusti, rustq)
  })
})
