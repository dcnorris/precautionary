## Manual tests go here, on the way to testthat/ dir

viola_stop_func <- function(x) {
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

viola_prior <- c(0.03, 0.07, 0.12, 0.20, 0.30, 0.40, 0.52)

test_viola_dtp <- function(scale = sqrt(0.75)) {
  timings <- list()
  outputs <- list()
  for (impl in c('dfcrm','rusti')) {
    timings[[impl]] <- system.time(
      outputs[[impl]] <- calculate_dtps(next_dose = 3,
                                        cohort_sizes = rep(3, 7),
                                        dose_func = applied_crm,
                                        prior = viola_prior,
                                        target = 0.2,
                                        stop_func = viola_stop_func,
                                        scale = scale,
                                        no_skip_esc = TRUE,
                                        no_skip_deesc = FALSE,
                                        global_coherent_esc = TRUE,
                                        impl = impl)
    )
  }
  lookhere <<- list(timings = timings, outputs = outputs)
  return(timings)
}


## Utility function
moments <- function(levels, numtox, s = 500, prior = viola_prior) {
  x <- rep(prior[levels], each=3)
  y <- as.integer(t(t(matrix(1:3, nrow=3, ncol=length(levels))) <= numtox))
  w <- rep(1, length(y))

  integrals <- list(
    dfcrm = c(integrate(dfcrm::crmh, -Inf, Inf, x, y, w, s)$value
             ,integrate(dfcrm::crmht, -Inf, Inf, x, y, w, s)$value
             ,integrate(dfcrm::crmht2, -Inf, Inf, x, y, w, s)$value
              )
  , Ri = c(integrate(precautionary:::crmh, -Inf, Inf, x, y, w, s)$value
          ,integrate(precautionary:::crmht, -Inf, Inf, x, y, w, s)$value
          ,integrate(precautionary:::crmht2, -Inf, Inf, x, y, w, s)$value
           )
  , rusti = c(integrate(rcrmh, -Inf, Inf, x, y, w, s)$value
             ,integrate(rcrmht, -Inf, Inf, x, y, w, s)$value
             ,integrate(rcrmht2, -Inf, Inf, x, y, w, s)$value
              )
  )
  integrals <- do.call(rbind, integrals)
  colnames(integrals) <- paste0("m", 0:2)
  quad <- as.data.table(integrals, keep.rownames="impl")
  ## Let's compute some of the crm outputs ...
  quad[,`:=`(
    est = m1/m0
   ,e2 = m2/m0
  )]
  quad[, post.var := e2 - est^2]
  quad
}

## There's a difference for D7 at index 95

test_viola_95 <- function() {
  moments(levels = c(3, 4, 5, 6, 6, 6, 4),
          numtox = c(0, 0, 0, 1, 1, 3, 2))
}

