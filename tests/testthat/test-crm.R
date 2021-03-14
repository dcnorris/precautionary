library(microbenchmark)

test_that("crm() yields same result as dfcrm::crm, but faster", {
  ## This set-up is verbatim from the example in dfcrm::crm
  prior <- c(0.05, 0.10, 0.20, 0.35, 0.50, 0.70)
  target <- 0.2
  level <- c(3, 4, 4, 3, 3, 4, 3, 2, 2, 2)
  y <- c(0, 0, 1, 0, 0, 1, 1, 0, 0, 0)
  mb_old <- microbenchmark(old <- dfcrm::crm(prior, target, y, level))
  mb_new <- microbenchmark(new <- precautionary::crm(prior, target, y, level))

  expect_equal(old$ptox, new$ptox) # TODO: Do a fuller check; consider waldo::compare()

  ## Expect DEFINITE improvement, with new UPPER quartile < old LOWER:
  expect_lt(summary(mb_new, unit="ms")$uq,
            summary(mb_old, unit="ms")$lq)

  speedup <- mean(mb_old$time) / mean(mb_new$time)
  message(paste0("speedup: ", format(speedup, digits=2), "x"))
})
