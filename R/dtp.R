#' A supremely faster version of a function from 'dtpcrm' v0.1.1
#'
#' Originally, the sampling in stats::rnorm() (see inline comments in code)
#' consumed 53% of run-time in a 6-cohort VIOLA DTP. After this change, it
#' doesn't even show up! More importantly, the consumption is now dominated
#' by (at 75%) by the objective function 'f' in integrate().
#' @param x A object of class \code{mtd}
#' @param tox_lim Scalar upper threshold on estimated toxicity rate
#' @param prob_cert Confidence level for threshold exceedance
#' @param dose Integer scalar, the dose being considered
#' @param suppress_dose Logical; if TRUE the MTD is set to \code{NA} when
#' trial stop is recommended.
#' @return
#' The \code{mtd} object x, with stop decision annotated
#' @author Adapted by David C. Norris from original dtpcrm::stop_for_excess_toxicity_empiric
#' @importFrom stats pnorm
#' @export
stop_for_excess_toxicity_empiric <- function (x, tox_lim, prob_cert, dose = 1,
                                              suppress_dose = TRUE) {
  post_beta_mean = x$estimate
  post_beta_var = x$post.var
  ## The following was massively wasteful, with the stats::rnorm()
  ## sampling itself consuming 53% of time in dtp(6)!
  ## > post_beta_samp = stats::rnorm(n = nsamps, mean = post_beta_mean,
  ## >                               sd = sqrt(post_beta_var))
  ## > post_prob_tox_samp = x$prior[dose]^exp(post_beta_samp)
  ## > prob_too_toxic = mean(post_prob_tox_samp > tox_lim)
  ## Yet a bit of math suffices to compute this with ZERO MCSE, via pnorm():
  prob_too_toxic <- pnorm(log(log(tox_lim)/log(x$prior[dose])),
                          mean = post_beta_mean, sd = sqrt(post_beta_var))
  stop_decision = prob_too_toxic > prob_cert
  x$stop = stop_decision
  if (stop_decision) {
    x$stop_reason = paste0("Prob(Prob(Tox[", dose, "]) > ",
                           tox_lim, ") = ", round(prob_too_toxic, 3), " > ",
                           prob_cert)
    if (suppress_dose)
      x$mtd = NA
  }
  return(x)
}
