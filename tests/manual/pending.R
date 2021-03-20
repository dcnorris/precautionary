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
  for (impl in c('dfcrm','Ri','rusti','rustq')) {
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
  viola.dtp <- list(
    dfcrm = calculate_dtps(next_dose = 3,
                           cohort_sizes = rep(3, 7),
                           dose_func = applied_crm,
                           prior = viola_prior,
                           target = 0.2,
                           stop_func = viola_stop_func,
                           scale = scale,
                           no_skip_esc = TRUE,
                           no_skip_deesc = FALSE,
                           global_coherent_esc = TRUE,
                           impl = "dfcrm"
                           )
  , Ri    = calculate_dtps(next_dose = 3,
                           cohort_sizes = rep(3, 7),
                           dose_func = applied_crm,
                           prior = viola_prior,
                           target = 0.2,
                           stop_func = viola_stop_func,
                           scale = scale,
                           no_skip_esc = TRUE,
                           no_skip_deesc = FALSE,
                           global_coherent_esc = TRUE,
                           impl = "Ri"
                           )
  , rusti = calculate_dtps(next_dose = 3,
                           cohort_sizes = rep(3, 7),
                           dose_func = applied_crm,
                           prior = viola_prior,
                           target = 0.2,
                           stop_func = viola_stop_func,
                           scale = scale,
                           no_skip_esc = TRUE,
                           no_skip_deesc = FALSE,
                           global_coherent_esc = TRUE,
                           impl = "rusti"
                           )
  , rustq = calculate_dtps(next_dose = 3,
                           cohort_sizes = rep(3, 7),
                           dose_func = applied_crm,
                           prior = viola_prior,
                           target = 0.2,
                           stop_func = viola_stop_func,
                           scale = scale,
                           no_skip_esc = TRUE,
                           no_skip_deesc = FALSE,
                           global_coherent_esc = TRUE,
                           impl = "rustq"
                           )
    )
  ## Compare vs cached (17-minute) computation. Because the cached original used
  ## dtpcrm's stochastic stop_for_excess_toxicity_empiric(), some irregularities
  ## emerge in a range of indices where high toxicities make a late appearance.
  ## I exclude these indices from the comparison, since they represent errors in
  ## the cached version!
  rust_viola_dtp <<- viola.dtp
  ##testthat::expect_identical(viola.dtp[-(4810:4870),], precautionary::viola_dtp[-(4810:4870),])

}


## A couple of cases where rustq != rusti (at least w.r.t. stopping)
## (See below for output & comment ...)
test_qvi <- function(s = 500) {
  ## > rust_viola_dtp$rusti[c(624,1200),]
  ##      D0 T1 D1 T2 D2 T3 D3 T4 D4 T5 D5 T6 D6 T7 D7
  ## 624   3  0  4  0  5  2  4  1  4  2  2  3  1  3 NA
  ## 1200  3  0  4  1  4  0  5  2  4  2  2  3  1  3 NA
  ## > rust_viola_dtp$rustq[c(624,1200),]
  ##      D0 T1 D1 T2 D2 T3 D3 T4 D4 T5 D5 T6 D6 T7 D7
  ## 624   3  0  4  0  5  2  4  1  4  2  2  3  1  3  1
  ## 1200  3  0  4  1  4  0  5  2  4  2  2  3  1  3  1

  list("path-624"  = moments(levels = c(3, 4, 5, 4, 4, 2, 1),
                             numtox = c(0, 0, 2, 1, 2, 3, 3))
      ,"path-1200" = moments(levels = c(3, 4, 4, 5, 4, 2, 1),
                             numtox = c(0, 1, 0, 2, 2, 3, 3))
       )
}

## > test_qvi(500)
## $`path-624`
##     impl           m0            m1           m2       est       e2  post.var
## 1: dfcrm 3.577549e-08 -4.512669e-08 6.076539e-08 -1.261386 1.698520 0.1074267
## 2:    Ri 3.577549e-08 -4.512669e-08 6.076539e-08 -1.261386 1.698520 0.1074267
## 3: rusti 3.577549e-08 -4.512669e-08 6.076539e-08 -1.261386 1.698520 0.1074267
## 4: rustq 3.577542e-08 -4.502761e-08 6.076712e-08 -1.258619 1.698572 0.1144507

## $`path-1200`
##     impl           m0            m1           m2       est       e2  post.var
## 1: dfcrm 3.577549e-08 -4.512669e-08 6.076539e-08 -1.261386 1.698520 0.1074267
## 2:    Ri 3.577549e-08 -4.512669e-08 6.076539e-08 -1.261386 1.698520 0.1074267
## 3: rusti 3.577549e-08 -4.512669e-08 6.076539e-08 -1.261386 1.698520 0.1074267
## 4: rustq 3.577542e-08 -4.502761e-08 6.076712e-08 -1.258619 1.698572 0.1144507

## ** COMMENT **
## The above shows how even a 0.2% (1 part in 450) error in m1 may translate
## to a large proportional error in post.var := e2 - est^2.
## Firstly, I am surprised that this is even the right calculation of post.var!
## Can we really be certain that it never goes negative?

## ** AHA! **
## As it turns out, increasing the order of Gauss-Kronrod quadrature to 31 points
## suffices to obtain very close accuracy:
## > test_qvi(500) # with 31-point G-K quadrature now
## $`path-624`
##     impl           m0            m1           m2       est      e2  post.var
## 1: dfcrm 3.577549e-08 -4.512669e-08 6.076539e-08 -1.261386 1.69852 0.1074267
## 2:    Ri 3.577549e-08 -4.512669e-08 6.076539e-08 -1.261386 1.69852 0.1074267
## 3: rusti 3.577549e-08 -4.512669e-08 6.076539e-08 -1.261386 1.69852 0.1074267
## 4: rustq 3.577547e-08 -4.512673e-08 6.076569e-08 -1.261388 1.69853 0.1074307

## $`path-1200`
##     impl           m0            m1           m2       est      e2  post.var
## 1: dfcrm 3.577549e-08 -4.512669e-08 6.076539e-08 -1.261386 1.69852 0.1074267
## 2:    Ri 3.577549e-08 -4.512669e-08 6.076539e-08 -1.261386 1.69852 0.1074267
## 3: rusti 3.577549e-08 -4.512669e-08 6.076539e-08 -1.261386 1.69852 0.1074267
## 4: rustq 3.577547e-08 -4.512673e-08 6.076569e-08 -1.261388 1.69853 0.1074307


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
  , rustq = c(icrm(x, y, w, s, 0)
             ,icrm(x, y, w, s, 1)
             ,icrm(x, y, w, s, 2)
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


## Test all moments for all impls

test_quad <- function(s = 500) {
  x <- c(1,2,3,2,3,3)*0.1
  y <- c(0L,0L,1L,0L,0L,1L)
  w <- rep(1,length(y))

  integrals <<- list(
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
  , rustq = c(icrm(x, y, w, s, 0)
             ,icrm(x, y, w, s, 1)
             ,icrm(x, y, w, s, 2)
              )
  )

  expect_equal(integrals$dfcrm, integrals$Ri)
  expect_equal(integrals$dfcrm, integrals$rusti)
  expect_equal(integrals$dfcrm, integrals$rustq, tolerance = 1e-7)
}
