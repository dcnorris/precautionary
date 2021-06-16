## BOIN v2.7.2 get.oc() function, reformatted with added commentary
## focused on treatment of 'n.earlystop' parameter documented as follows:

## n.earlystop: the early stopping parameter. If the number of patients
##           treated at the current dose reaches ‘n.earlystop’, stop the
##           trial and select the MTD based on the observed data. The
##           default value ‘n.earlystop=100’ essentially turns off this
##           type of early stopping.

get.oc <- function (target, p.true, ncohort, cohortsize, n.earlystop = 100, 
                    startdose = 1, titration = FALSE, p.saf = 0.6 * target,
                    p.tox = 1.4 * target, cutoff.eli = 0.95, extrasafe = FALSE,
                    offset = 0.05, boundMTD = FALSE, ntrial = 1000, seed = 6) 
{
  if (target < 0.05) {
    stop("the target is too low!")
  }
  if (target > 0.6) {
    stop("the target is too high!")
  }
  if ((target - p.saf) < (0.1 * target)) {
    stop("the probability deemed safe cannot be higher than or too close to the target!")
  }
  if ((p.tox - target) < (0.1 * target)) {
    stop("the probability deemed toxic cannot be lower than or too close to the target!")
  }
  if (offset >= 0.5) {
    stop("the offset is too large!")
  }
  if (n.earlystop <= 6) {
    warning("the value of n.earlystop is too low to ensure good operating characteristics. Recommend n.earlystop = 9 to 18.")
  }
  set.seed(seed)
  if (cohortsize == 1) 
    titration = FALSE
  lambda_e = log((1 - p.saf)/(1 - target))/log(target * (1 - p.saf)/(p.saf * (1 - target)))
  lambda_d = log((1 - target)/(1 - p.tox))/log(p.tox * (1 - target)/(target * (1 - p.tox)))
  ndose = length(p.true)
  npts = ncohort * cohortsize
  Y = matrix(rep(0, ndose * ntrial), ncol = ndose)
  N = matrix(rep(0, ndose * ntrial), ncol = ndose)
  dselect = rep(0, ntrial)
  if (cohortsize > 1) {
    temp = get.boundary(target, ncohort, cohortsize,
                        n.earlystop = ncohort * cohortsize,
                        p.saf, p.tox, cutoff.eli, extrasafe)$full_boundary_tab
  }
  else {
    temp = get.boundary(target, ncohort, cohortsize,
                        n.earlystop = ncohort * cohortsize,
                        p.saf, p.tox, cutoff.eli, extrasafe)$boundary_tab
  }
  b.e = temp[2, ]
  b.d = temp[3, ]
  b.elim = temp[4, ]
  for (trial in 1:ntrial) {
    y <- rep(0, ndose)
    n <- rep(0, ndose)
    earlystop = 0
    d = startdose
    elimi = rep(0, ndose)
    ft = TRUE
    if (titration) {
      z <- (runif(ndose) < p.true)
      if (sum(z) == 0) {
        d = ndose
        n[1:ndose] = 1
      }
      else {
        d = which(z == 1)[1]
        n[1:d] = 1
        y[d] = 1
      }
    }
    for (i in 1:ncohort) {
      if (titration && n[d] < cohortsize && ft) {
        ft = FALSE
        y[d] = y[d] + sum(runif(cohortsize - 1) < p.true[d])
        n[d] = n[d] + cohortsize - 1
      }
      else {
        newcohort = runif(cohortsize) < p.true[d]
        if ((sum(n) + cohortsize) >= npts) {
          nremain = npts - sum(n)
          y[d] = y[d] + sum(newcohort[1:nremain])
          n[d] = n[d] + nremain
          break
        }
        else {
          y[d] = y[d] + sum(newcohort)
          n[d] = n[d] + cohortsize
        }
      }
      if (!is.na(b.elim[n[d]])) {
        if (y[d] >= b.elim[n[d]]) {
          elimi[d:ndose] = 1
          if (d == 1) {
            earlystop = 1
            break
          }
        }
        if (extrasafe) {
          if (d == 1 && n[1] >= 3) {
            if (1 - pbeta(target, y[1] + 1, n[1] - y[1] + 1) > cutoff.eli - offset) {
              earlystop = 1
              break
            }
          }
        }
      }
      ## if (n[d] >= n.earlystop && ((y[d] > b.e[n[d]] && 
      ##                              y[d] < b.d[n[d]]) || (d == 1 && y[d] >= b.d[n[d]]) || 
      ##                             ((d == ndose || elimi[d + 1] == 1) && y[d] <= 
      ##                              b.e[n[d]]))) 
      ##   break
      ##
      ## KEY
      ## d    : most-recently-given dose
      ## n[d] : cumulative enrollment at dose 'd'
      ## y[d] : number of toxicities at dose 'd' (thus n[d]-y[d] is number of non-tox)
      ## b.d  : deescalation boundary (integer vector of denominator-wise toxicity counts)
      ## b.e  : escalation boundary (integer vector of denominator-wise toxicity counts)
      ## ndose : number of pre-specified doses
      ## elimi : 0/1 vector indicating which doses have been eliminated
      if (n[d] >= n.earlystop)
        if ((b.e[n[d]] < y[d] && y[d] < b.d[n[d]]) # ntox strictly between esc & deesc bdys
            || (d == 1 && y[d] >= b.d[n[d]]) # OR deesc indicated but d already at lowest dose
            || ((d == ndose || elimi[d + 1] == 1) && y[d] <= b.e[n[d]]) # OR..
            ## ..escalation is indicated, but dose cannot go any higher
            )
          break
      if (y[d] <= b.e[n[d]] && d != ndose) {
        if (elimi[d + 1] == 0) 
          d = d + 1
      }
      else if (y[d] >= b.d[n[d]] && d != 1) {
        d = d - 1
      }
      else {
        d = d
      }
    }
    Y[trial, ] = y
    N[trial, ] = n
    if (earlystop == 1) {
      dselect[trial] = 99
    }
    else {
      dselect[trial] = select.mtd(target, n, y, cutoff.eli, 
                                  extrasafe, offset, boundMTD = boundMTD, p.tox = p.tox)$MTD
    }
  }
  selpercent = rep(0, ndose)
  nptsdose = apply(N, 2, mean)
  ntoxdose = apply(Y, 2, mean)
  for (i in 1:ndose) {
    selpercent[i] = sum(dselect == i)/ntrial * 100
  }
  if (length(which(p.true == target)) > 0) {
    if (which(p.true == target) == ndose - 1) {
      overdosing60 = mean(N[, p.true > target] > 0.6 * 
                          npts) * 100
      overdosing80 = mean(N[, p.true > target] > 0.8 * 
                          npts) * 100
    }
    else {
      overdosing60 = mean(rowSums(N[, p.true > target]) > 
                          0.6 * npts) * 100
      overdosing80 = mean(rowSums(N[, p.true > target]) > 
                          0.8 * npts) * 100
    }
    out = list(selpercent = selpercent, npatients = nptsdose, 
               ntox = ntoxdose, totaltox = sum(Y)/ntrial, totaln = sum(N)/ntrial, 
               percentstop = sum(dselect == 99)/ntrial * 100, overdose60 = overdosing60, 
               overdose80 = overdosing80, simu.setup = data.frame(target = target, 
                                                                  p.true = p.true, ncohort = ncohort, cohortsize = cohortsize, 
                                                                  startdose = startdose, p.saf = p.saf, p.tox = p.tox, 
                                                                  cutoff.eli = cutoff.eli, extrasafe = extrasafe, 
                                                                  offset = offset, ntrial = ntrial, dose = 1:ndose), 
               flowchart = TRUE, lambda_e = lambda_e, lambda_d = lambda_d)
  }
  else {
    out = list(selpercent = selpercent, npatients = nptsdose, 
               ntox = ntoxdose, totaltox = sum(Y)/ntrial, totaln = sum(N)/ntrial, 
               percentstop = sum(dselect == 99)/ntrial * 100, simu.setup = data.frame(target = target, 
                                                                                      p.true = p.true, ncohort = ncohort, cohortsize = cohortsize, 
                                                                                      startdose = startdose, p.saf = p.saf, p.tox = p.tox, 
                                                                                      cutoff.eli = cutoff.eli, extrasafe = extrasafe, 
                                                                                      offset = offset, ntrial = ntrial, dose = 1:ndose), 
               flowchart = TRUE, lambda_e = lambda_e, lambda_d = lambda_d)
  }
  class(out) <- "boin"
  return(out)
}
