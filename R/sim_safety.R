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
                       , r0 = seq(0.5, 2.5, 0.5)
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
  toxicities <<- list()
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
    ensemble[, Tox := factor(toxgrade, levels=0:5, labels=paste0("Gr",0:5))]
    ensembles[[k]] <<- ensemble
    # See https://stackoverflow.com/a/16519612/3338147 explaining below syntax
    toxicities[[k]] <<- ensemble[, .(Tox=ordered(levels(Tox)), N=c(table(Tox)))
                                 , by=r0]
    # Importantly, because these protocols run considering only DLT = Gr>=3,
    # they are r0-agnostic. But I ought to allow for the general case where
    # ordinal toxicities (MTDig for g != 3) affect trial conduct.
    #
    # Initially, let me go simply for *counts* and worry secondarily about
    # obtaining (and dividing by) the denominators.
    # Want a tabulation with several r0 values defining rows, and columns for
    # the counts of toxicity grades. A final column may show total enrollments,
    # which as noted above will all be equal unless ordinal toxicities affect
    # escalation or termination decisions.
  }
  # Indeed, I like the idea of abstracting the calculation of high-level
  # summary statistics into a separate function, possibly even a 'summary'
  # method for a suitably defined R3 class. But for now, let me implement
  # these summaries here, as additional components of the returned list.
  toxdt <- rbindlist(toxicities, idcol = "k")
  counts <- toxdt[, .(n=sum(N)), by=.(r0,Tox)]
  toxtab <- ftable(xtabs(n ~ Tox + r0, counts))
  # It now seems to me there are 2 perspectives on the probabilities here.
  # One POV is the per-trial perspective, which might ask for expected
  # numbers of each grade of toxicity.
  # The second POV is that of the enrolling patient, who might ask what
  # is the probability of each grade of toxicity. But this POV ought to
  # be conditioned on the prior results seen in the trial, and so demands
  # much more sophisticated modeling -- indeed, modeling of the kind I
  # employed in the AFM11 paper.
  # Note in fact that, as soon as such a perspective BEGINS to be acknowledged,
  # I have already 'won' the argument about dose individualization!
  list(tpt = true_prob_tox
      ,toxdt = toxdt
      ,toxtab = toxtab
      ,expect = toxtab/(K*M)
      ,K = K
      ,M = M
      )
}