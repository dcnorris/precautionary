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
      x <- dtpcrm::stop_for_consensus_reached(x, req_at_mtd = 12)
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

test_that("Crm$trace_paths() answer invariant to unroll depth", {
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
      x <- dtpcrm::stop_for_consensus_reached(x, req_at_mtd = 12)
    }
  }

  crm <- Crm$new(skeleton = prior.DLT,
                 scale = sqrt(prior.var),
                 target = target.DLT)$
    stop_func(stop_func)$
    no_skip_esc(TRUE)$
    no_skip_deesc(FALSE)$
    global_coherent_esc(TRUE)

  pmx1 <- crm$trace_paths(root_dose = start.dose.level,
                          cohort_sizes = rep(3, 7),
                          impl = 'rusti', unroll = 1)$path_matrix()

  pmx2 <- crm$trace_paths(root_dose = start.dose.level,
                          cohort_sizes = rep(3, 7),
                          impl = 'rusti', unroll = 2)$path_matrix()

  pmx3 <- crm$trace_paths(root_dose = start.dose.level,
                          cohort_sizes = rep(3, 7),
                          impl = 'rusti', unroll = 3)$path_matrix()

  pmx4 <- crm$trace_paths(root_dose = start.dose.level,
                          cohort_sizes = rep(3, 7),
                          impl = 'rusti', unroll = 4)$path_matrix()

  expect_equal(dim(pmx2), dim(pmx1))
  expect_equal(dim(pmx3), dim(pmx1))
  expect_equal(dim(pmx4), dim(pmx1))
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
      x <- dtpcrm::stop_for_consensus_reached(x, req_at_mtd = 12)
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
