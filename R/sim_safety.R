# Wait a minute! I see that true_prob_tox does not get attached to the
# protocol, but rather to the escalation::simulations object returned
# by simulate_trials.
# Thus, my own function might look more properly like a generalization
# of escalation::simulate_trials, in which 'true_prob_tox' is an entire
# *prior distribution*, rather than just a vector.
sim_safety <- function(protocol
                       , lambda_CV = 3
                       , median_mtd = protocol$num_doses - 1
                       , median_sd = median_mtd/3
                       , r0 = seq(1.5, 3, 0.5)
                       , K = 40
                       , M = 50
                       ){
  stopifnot(is(protocol, "selector_factory"))
  # 1. Build a data.table of 'true_prob_tox' vectors obtained
  #    by sampling the prior
  true_prob_tox <- data.table(CV = rexp(n=K, rate=lambda_CV))
  true_prob_tox[, `:=`(
    sigma = sqrt(log(CV^2 + 1)) # NB: CV>sigma, w/ near-identity below 0.4
  , mu = rnorm(n=K, mean=median_mtd, sd=median_sd)
  )]
  for(d in seq(protocol$num_doses)){
    true_prob_tox[, paste0("P",d) := pnorm(q=d, mean=mu, sd=sigma)]
  }
  setcolorder(true_prob_tox, paste0("P", 1:protocol$num_doses)) # move P1..Pd to front
  ensembles <<- list()
  for(k in 1:K){
    cat("k =", k, "\n")
    sims <- simulate_trials(protocol, num_sims = M,
                            true_prob_tox = as.numeric(true_prob_tox[k, 1:protocol$num_doses]))
    ensemble <- rbindlist(lapply(sims[[1]], function(.) .[[1]]$fit$outcomes)
                          , idcol = "rep")
    ensemble <- merge(data.frame(r0=r0), ensemble) # cartesian product
    ensemble <- as.data.table(ensemble) # restore data.table after cartesian product
    ensemble[, MTDi3 := qnorm(p = u_i
                              , mean = true_prob_tox[k]$mu
                              , sd = true_prob_tox[k]$sigma
                              )]
    ensemble[, `:=`(
      MTDi1 = MTDi3/r0^2
    , MTDi2 = MTDi3/r0
    , MTDi4 = MTDi3*r0
    , MTDi5 = MTDi3*r0^2
    )]
    ensemble[, toxgrade := (dose>MTDi1) + (dose>MTDi2) + (dose>MTDi3) +
                           (dose>MTDi4) + (dose>MTDi5)]
    stopifnot(with(ensemble, all(tox == as.integer(dose >= MTDi3))))
    # TODO: The following summary is wrong. I need to write the ensemble
    #       to .GlobalEnv, and experiment 'by hand' to find the right
    #       summary scheme.
    # N.B.: Any trial realization ('rep') is just as likely as any other.
    #       Therefore, I have to aggregate with reps as the denominator.
    #       Specifically, pooling all the patients would fail to weight
    #       them inversely to trial size.
    #       Another way to see this is, we are asking about the probability
    #       of a fatality in *this* trial. Alternatively, we can ask for
    #       the expected number of each grade of toxicity in *this* trial.
    #       Note that this should partition the expected size of this trial!
    toxicities <- ensemble[, .(prob = .N/nrow(.SD)), by = .(toxgrade,r0)]
    ensembles[[k]] <<- ensemble
  }
  list(tpt = true_prob_tox, tox = toxicities)
}