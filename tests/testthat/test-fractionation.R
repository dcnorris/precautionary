test_that("CPE-cognizant tox fractionation matches old nested sims", {

  ## Revisit the 3+3 example from Intro vignette
  design <- get_three_plus_three(num_doses = 5
                                ,allow_deescalate = TRUE)

  mtdi_gen <- hyper_mtdi_lognormal(CV = 1
                                  ,median_mtd = 5
                                  ,median_sdlog = 0.5
                                  ,units="mg/kg"
                                   )
  options(dose_levels = c(0.5, 1, 2, 4, 6)) # specify actual dosing

  set.seed(2021)
  design %>% simulate_trials(
               num_sims = 2000
             , true_prob_tox = mtdi_gen # pull tox probs from a MODEL, not thin air
             ) -> HYPERSIMS

  tox_threshold_scaling <- function(MTDi, r0) {
    MTDi * r0 ^ c(Gr1=-2, Gr2=-1, Gr3=0, Gr4=1, Gr5=2)
  }

  r0 <- 2
  safety <- summary(HYPERSIMS
                   ,ordinalizer = tox_threshold_scaling
                   ,r0 = r0)$safety

  ## ~~ Let's compare results obtained via CPE ~~

  MTDi_gen <- HyperMTDi_lognormal$new(CV = 1
                                     ,median_mtd = 5
                                     ,median_sdlog = 0.5
                                     ,units = "mg/kg"
                                     ,n = 1000)$doses(getOption("dose_levels"))

  bUF <- MTDi_gen$fractionate(Cpe3_3$new(5)
                             ,kappa = log(r0)
                              )

  ## TODO: A more sophisticated view would not put the two MCSE's on equal footing!
  err_se <- (bUF[1,] - safety[1,]) / sqrt(bUF['MCSE',]^2 + safety['MCSE',]^2)
  expect_lt(max(abs(err_se)), 3) # expect all differences under 3 sigma

})
