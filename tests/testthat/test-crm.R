library(microbenchmark)

test_that("Faster objective functions integrate same as 'dfcrm' originals", {

  skeleton <- seq(0.1, 0.6, 0.1)

  ## cohort-wise set-up
  level <- c(1L, 2L, 3L, 2L, 3L, 3L)
  n_tox <- c(0L, 0L, 1L, 0L, 0L, 1L)
  obs <- encode_cohorts(enr = tabulate(level, nbins=8)
                       ,tox = xtabs(n_tox ~ factor(level, levels=1:8),
                                    data=data.frame(n_tox = n_tox,
                                                    level = level))
                        )

  ## The 'dfcrm' & 'rusti' implementations take patient-wise vectors x and y,
  ## and patient-wise weights w that 'break exchangeability'.
  x <- skeleton[rep(level, each=3)]
  coh <- function(tox, n=3) c(rep(0L,n-tox), rep(1L,tox))
  y <- as.integer(sapply(n_tox, coh))
  w <- rep(1,length(y))
  w[y == 1] <- 0.0 # encode y in w for rusti's benefit (dfcrm is unaffected)
  s <- 500

  integrals <- list(
    dfcrm = c(integrate(dfcrm::crmh, -Inf, Inf, x, y, w, s)$value
             ,integrate(dfcrm::crmht, -Inf, Inf, x, y, w, s)$value
             ,integrate(dfcrm::crmht2, -Inf, Inf, x, y, w, s)$value
              )
  , rusti = c(integrate(crmh, -Inf, Inf, log(x), w, s)$value
             ,integrate(crmht, -Inf, Inf, log(x), w, s)$value
             ,integrate(crmht2, -Inf, Inf, log(x), w, s)$value
              )
  , ruste = c(integrate(crmh_ev, -Inf, Inf, obs, log(skeleton), s, 0)$value
             ,integrate(crmh_ev, -Inf, Inf, obs, log(skeleton), s, 1)$value
             ,integrate(crmh_ev, -Inf, Inf, obs, log(skeleton), s, 2)$value
              )
  )

  with(integrals, {
    expect_equal(rusti, dfcrm)
    expect_equal(ruste, dfcrm)
  })

  ## Now again, but with weights not all 1.0:
  w <- runif(n = length(y), min=0.8, max=1.0)
  w[y==1] <- 0.0

  nontrivial_weights <- list(
    dfcrm = c(integrate(dfcrm::crmh, -Inf, Inf, x, y, w, s)$value
             ,integrate(dfcrm::crmht, -Inf, Inf, x, y, w, s)$value
             ,integrate(dfcrm::crmht2, -Inf, Inf, x, y, w, s)$value
              )
  , rusti = c(integrate(crmh, -Inf, Inf, log(x), w, s)$value
             ,integrate(crmht, -Inf, Inf, log(x), w, s)$value
             ,integrate(crmht2, -Inf, Inf, log(x), w, s)$value
              )
    # NB: 'ruste' inapplicable when patients not exchangeable
  )

  with(nontrivial_weights,
       expect_equal(rusti, dfcrm))

})

test_that("Crm$est() yields same result as dfcrm::crm, but faster", {
  ## This set-up is verbatim from the example in dfcrm::crm
  prior <- c(0.05, 0.10, 0.20, 0.35, 0.50, 0.70)
  target <- 0.2
  level <- c(3, 4, 4, 3, 3, 4, 3, 2, 2, 2)
  y     <- c(0, 0, 1, 0, 0, 1, 1, 0, 0, 0)
  s <- sqrt(1.34)
  ##  We provide the levels in a new way to the faster code, however ...
  x <- xtabs(y ~ factor(level, levels=seq_along(prior)))
  o <- xtabs(!y ~ factor(level, levels=seq_along(prior)))
  crm_old <- microbenchmark(old <- dfcrm::crm(prior, target, y, level))
  crm_new <- microbenchmark(new <- crm(prior, target, y, level, impl="rusti"))
  ##r6model <- Crm$new(skeleton = prior, target = target)$tally(x, o)
  r6model <- Crm$new(skeleton = prior, target = target)$observe(level, y)
  crm_ri <- microbenchmark(r6i <- r6model$est(impl="rusti", abbrev=FALSE))
  ## TODO: Let 'ruste' method send dosewise x and o vectors directly,
  ##       without asking the Rust routine to reconstruct these.
  crm_re <- microbenchmark(r6e <- r6model$est(impl="ruste", abbrev=FALSE))

  expect_equal(new, old)

  ## Initially, I won't burden class 'Crm' with preserving identity
  ## where dose-wise observations are exchangeable in the likelihood.
  ## So to effect a check, I simply sort 'old' to match my construction
  ## (o first, then x) of the Crm fields.
  exch <- order(paste(old$tox, old$level, sep='.'))
  old$tox <- old$tox[exch]
  old$level <- old$level[exch]
  expect_equal(r6i, old)

  ##expect_equal(r6e, old)

  ## Expect DEFINITE improvement, with new UPPER quartile < old LOWER:
  expect_lt(summary(crm_new, unit="ms")$uq,
            summary(crm_old, unit="ms")$lq)

  speedup_message(crm_new, crm_old)
  speedup_message(crm_ri, crm_old)
  ##speedup_message(crm_re, crm_old)
})

## test_that("titecrm() yields same result as dfcrm::crm, but faster", {
##   ## This set-up is verbatim from the example in dfcrm::titecrm
##   prior <- c(0.05, 0.10, 0.20, 0.35, 0.50, 0.70)
##   target <- 0.2
##   level <- c(3, 3, 3, 4, 4, 3, 2, 2, 2, 3)
##   y <- c(0, 0, 1, 0, 1, 0, 0, 0, 0, 0)
##   u <- c(178, 181, 168, 181, 24, 181, 179, 102, 42, 3)
##   tau <- 180
##   foo <- titecrm(prior, target, y, level, followup=u, obswin=tau)
##   rec <- foo$mtd  # recommend a dose level for next patient

##   ## An example with adaptive weight
##   foo2 <- titecrm(prior, target, y, level, followup=u, obswin=tau, scheme="adaptive")
##   wts <- foo2$weights

##   ## The 'weights' argument makes 'followup' and 'obswin' obsolete
##   foo3 <- titecrm(prior, target, y, level, weights=wts, followup=u, obswin=tau)
##   \dontrun{plot(foo3, ask=T)}

##   ## Patient time information via 'entry' and 'exit' arguments
##                                         # entry time (days since study begins)
##   entry <- c(7, 29, 49, 76, 92, 133, 241, 303, 363, 402)
##                                         # exit time (days since study begins)
##   exit <- c(185, 210, 217, 257, 116, 314, 420, 405, 405, 405)
##   foo4 <- titecrm(prior, target, y, level, exit=exit, entry=entry, obswin=tau)
##   \dontrun{plot(foo4, ask=T)}

##   prior <- c(0.05, 0.10, 0.20, 0.35, 0.50, 0.70)
##   target <- 0.2
##   level <- c(3, 4, 4, 3, 3, 4, 3, 2, 2, 2)
##   y <- c(0, 0, 1, 0, 0, 1, 1, 0, 0, 0)
##   tite_old <- microbenchmark(old <- dfcrm::crm(prior, target, y, level))
##   tite_new <- microbenchmark(new <- crm(prior, target, y, level, impl="rusti"))

##   ##expect_equal(old$ptox, new$ptox) # TODO: Do a fuller check; consider waldo::compare()
##   expect_equal(old, new) # TODO: Do a fuller check; consider waldo::compare()

##   ## Expect DEFINITE improvement, with new UPPER quartile < old LOWER:
##   expect_lt(summary(tite_new, unit="ms")$uq,
##             summary(tite_old, unit="ms")$lq)

##   speedup_message(tite_new, tite_old)
## })
