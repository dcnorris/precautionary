## Decision boundaries for 'local BOIN' with pi[0,j] = pi[1,j] = pi[2,j],
## defaulting to the recommended phi_1/phi = 0.6, phi_2/phi = 1.4

boin_lambda <- function(phi, delta = c(-0.4, 0.4)) {
  phi_ <- (1 + delta) * phi

  inv_lambda <- 1 - log(1 + delta)/log((1 - phi_)/(1 - phi))
  1/inv_lambda
}

tabulate_boin_lambdas <- function(phi = seq(0.15, 0.3, 0.05)) {
  for (p in phi)
    print(fractions(boin_lambda(p)))
}

## The 5% quantile of Beta(x+1, o+1) is the value of phi we are
## "95% confident" is exceeded at a dose with tally x/(x+o),
## given a Bayes-Laplace prior B(1,1) is used.
predicate_qbeta05 <- function(nmax=14) {
  for (n in 1:nmax) {
    for (x in 1:n) {
      q <- qbeta(p = 0.05, shape1 = x+1, shape2 = n-x+1)
      cat(paste0("qbeta05_alpha_beta(", fractions(q),", ", x+1, ", ", n-x+1, "). % ", q, "\n"))
    }
  }
}
