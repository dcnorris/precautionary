library(microbenchmark)

test_that("Faster objective functions integrate same as 'dfcrm' originals", {

  x <- c(1,2,3,2,3,3)*0.1
  y <- c(0L,0L,1L,0L,0L,1L)
  w <- rep(1,length(y))
  s <- 500

  integrals <- list(
    dfcrm = integrate(dfcrm::crmh, -Inf, Inf, x, y, w, s)
  , Ri = integrate(precautionary:::crmh, -Inf, Inf, x, y, w, s)
  , rusti = integrate(rcrmh, -Inf, Inf, x, y, w, s)
  , rustq = icrm(x, y, w, s, 0)
  )

  expect_equal(integrals$dfcrm[[1]], integrals$Ri[[1]])
  expect_equal(integrals$dfcrm[[1]], integrals$rusti[[1]])
  expect_equal(integrals$dfcrm[[1]], integrals$rustq[[1]])
})

test_that("crm() yields same result as dfcrm::crm, but faster", {
  ## This set-up is verbatim from the example in dfcrm::crm
  prior <- c(0.05, 0.10, 0.20, 0.35, 0.50, 0.70)
  target <- 0.2
  level <- c(3, 4, 4, 3, 3, 4, 3, 2, 2, 2)
  y <- c(0, 0, 1, 0, 0, 1, 1, 0, 0, 0)
  mb_old <- microbenchmark(old <- dfcrm::crm(prior, target, y, level))
  mb_new <- microbenchmark(new <- crm(prior, target, y, level, impl="Ri"))

  expect_equal(old$ptox, new$ptox) # TODO: Do a fuller check; consider waldo::compare()

  ## Expect DEFINITE improvement, with new UPPER quartile < old LOWER:
  expect_lt(summary(mb_new, unit="ms")$uq,
            summary(mb_old, unit="ms")$lq)

  speedup <- mean(mb_old$time) / mean(mb_new$time)
  message(paste0("speedup: ", format(speedup, digits=2), "x"))
})
