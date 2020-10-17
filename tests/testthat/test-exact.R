test_that("exact 3+3 ('more common variant') safety matches simulated", {
  design <- get_three_plus_three(
    num_doses = 5,
    allow_deescalate = TRUE # the 'more common' variant per Korn &al Stat Med (1994)
  )
  mtdi_dist <- mtdi_lognormal(CV = 2
                              ,median = 5
                              ,units = "mg/kg"
  )
  options(dose_levels = c(0.5, 1, 2, 4, 6))
  
  set.seed(2020) # avoid gut-wrenching random failures on CRAN
  design %>% simulate_trials(
    num_sims = 200
    , true_prob_tox = mtdi_dist
  ) -> SIMS
  
  options(ordinalizer = function(MTDi, r0 = 1.5) {
    MTDi * r0 ^ c(Gr1=-2, Gr2=-1, Gr3=0, Gr4=1, Gr5=2)
  })
  
  exact(design) %>% simulate_trials(true_prob_tox = mtdi_dist) -> EXACT
  
  SIMS_safety <- summary(SIMS)$safety
  abserrs <- abs(SIMS_safety[1,] - summary(EXACT)$safety)
  relerrs <- abserrs / SIMS_safety['MCSE',]
  
  expect(max(relerrs) < 3,
         "entries in simulated 3+3 safety table are exact ± 3 MCSEs")
  
})
