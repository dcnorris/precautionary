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

test_that("Crm$trace_paths() yields same VIOLA result as dtpcrm's version", {
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

  pmx <- crm$trace_paths(root_dose = start.dose.level,
                         cohort_sizes = rep(3, 7),
                         impl = 'rusti')$path_matrix()

  data(viola_dtp) # saved for comparison

  ## Of note, the new Crm$paths() method retains more dose-recommendation info
  ## than the 'dtpcrm' code under its default settings which write 'NA' dose recs
  ## when the trial stops early.
  ## Consequently, to effect a test, we copy NA's over from viola_dtp to pmx:
  viola_paths <- unique(viola_dtp)
  expect_equal(dim(pmx), dim(viola_paths))
  if (all(dim(pmx) == dim(viola_paths))) {
    pmx[is.na(viola_paths)] <- NA_integer_
    expect_equal(as.integer(pmx), as.integer(as.matrix(viola_paths)))
  }
})

test_that("Crm-class path_matrix can be recovered from path_array", {
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

  crm$trace_paths(root_dose = start.dose.level,
                  cohort_sizes = rep(3, 7),
                  impl = 'rusti')

  pmx <- crm$path_matrix()
  Tv <- crm$path_array(condense=FALSE)

  doses <- apply(Tv, MARGIN=1:2, FUN=function(x) which(!is.na(x))[1])
  colnames(doses) <- paste0("D",0:6) # force pmx column names

  ## Since pmx may contain final recommendations that were never dosed
  ## and consequently have no corresponding toxicity counts, which we
  ## can hardly expect to recover from T, we have to mask these cases
  ## in order to effect our comparison:
  M <- !is.na(doses)
  M[M] <- 1L
  M[!M] <- NA_integer_
  expect_equal(M*doses, M*pmx[,colnames(doses)])

})
