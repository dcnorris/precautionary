## Explorations of quadrature degeneracy in DTP computation

## I'm realizing now the importance of an efficient (and human-readable!) encoding
## for the dose-wise enrollment and toxicity counts. If I limit myself initially
## to no more than 10 separate doses (quite reasonable in practice!) and 9 or fewer
## cohorts (sufficient for VIOLA, with 7), then an R integer can encode each of
## these halves of the quadrature index.

## TODO: In general, it might be possible to ship 'precautionary' with pre-built
##       hashes for each of many types of scenario, or these might be built locally
##       if CRAN complains.

## In the interest of human-readability, I would like to read these integers with
## the LSB on the left. This actually makes decimal places more useful for enrollment!
## But really, all I need from the doubles is a sortable representation.
## The actual display is a separate concern.
## The numeric representation actually allows for 15 digits of precision. Could I
## possibly fit BOTH enrollment and toxicity into this space?
## What are the extreme limits of this possibility?
## * say 8 distinct doses, at most
## * n=3 case first:
##   - say max 5 cohorts enrolled per dose
##   - say max 15 tox at any given dose (only a Therac-25 trial could violate this)
##  => then 8*log2((5+1)*(15+1)) = 52.68 bits needed
##  $$ using the sign bit, or exponent, easily allows this.
## * n=1 case:
##   - enroll at most 15 per dose ==> 8*4 = 32 bits
##   - suppose up to 8 tox may occur at each dose
##   - then 8*log2(8+1) = 27 bits needed
##  $$ again, 32+27=59 just fits if I use 9 bits of the exponent

## Convert a floating point number to mantissa & exponent bits
double_bits <- function(x) {
  man <- x * 2^-ceiling(log2(x))
  exp <- ceiling(log2(x))

  stopifnot(abs(man) < 1.0)
  stopifnot(x == man*2^exp)
  list(man=man, exp=exp)
}

## Initially, let's assume n=3 cohorts. (But it would be ideal to support all of n=1:3
## in the same encoding! I probably have just enough bits to manage that!)
encode_cohorts <- function(tox, enr) {
  ## 'tox' is expected to be an integer vector (length < 9) of dose-wise
  ##       toxicity counts, ordered from lowest to highest doses.
  stopifnot(is.integer(tox))
  stopifnot(length(tox) <= 8)
  stopifnot(all(tox %in% 0:15))
  ## 'enr' is an integer vector of cohort counts, ordered likewise
  stopifnot(is.integer(enr))
  stopifnot(length(enr) <= 8)
  stopifnot(all(enr >= 0))
  stopifnot(all(enr <= 5))

  ## Extend tox & enr to length-8 vectors
  enr8 <- tox8 <- integer(8)
  enr8[seq_along(enr)] <- enr
  tox8[seq_along(tox)] <- tox

  ## Let the toxicity be represented in the MSB, and enrollment in LSBs.
  x <- sum(tox8 * 16^(0:7)) + sum(enr8 * 6^-(1:8))

  check <- decode_cohorts(x)
  if (all(check$tox == tox8))
    if (all(check$enr == enr8))
      return(x)
    else
      message("enr = ", enr8, " != check$enr = ", check$enr)
  else
    message("tox = ", tox8, " != check$tox = ", check$tox)

  return(NULL)

}

decode_cohorts <- function(x) {
  tox <- sapply(1:8, function(d) (x %% 16^d) %/% 16^(d-1)) # base-16 decoding
  .x <- x %% 1 # fractional part
  enr <- sapply(-(0:7), function(d) (.x %% 6^d) %/% 6^(d-1)) # base-6 decoding

  list(tox = tox,
       enr = enr)
}

## Let me initially expect an 'xtabs' object.
## I map to a double-precision float.
## This gives me a guaranteed 15 digits of precision, enabling encoding of up to
## 15 cohorts.
encode_enr <- function(xt) {
  if (!is(xt,'xtabs')) stop("xt should be an 'xtabs' object")
  if (length(xt) > 15) stop("encoding for more than 15 cohorts not yet implemented")

  enr <- as.integer(xt)
  sum(enr * 10^-seq.int(1, length(xt)))
}

encode_tox <- function(xt) {
  if (!is(xt,'xtabs')) stop("xt should be an 'xtabs' object")
  if (length(xt) > 9) stop("encoding for more than 9 cohorts not yet implemented")

}

## Compute (ex post) the quadrature degeneracy of a DTP table.
## For a given skeleton, the quadratures required to run CRM are unique only up to
## the counts of dose assignments.
## Each path gives rise to a sequence of dose tabulations, one for each head slice.
## So my computation of the quadrature degeneracy of a given DTP table
## begins by expanding each path (with R dose recommendations) into R
## rows of a new table.
##
## Q: What is the complete set of inputs that indexes a quadrature?
## A: Numbers of O/X at each of the doses.
##
## Thus, for a D-dose trial (D=7 for VIOLA), a tuple of 2*D integers comprises
## the quadrature index. I might write this as a D-tuple of pairs T/N.
## Thus, there is a dose-recommendation function from (T/N)^D --> 0:D [where 0=NA].
##
## The best way to explore this would be through a function that generates
## all possible quadrature indices for a given row of a DTP table.
quadrature_hash <- function(dtp) {
  vp <- as.matrix(dtp)
  cat("dim(vp) = ", deparse(dim(vp)), "\n")
  C <- (ncol(dtp)-1)/2
  cat(C, " cohorts\n")
  P <- nrow(dtp)
  cat(P, "paths\n")
  vd <- vp[,paste0("D",0:(C-1))]

  ve <- t(apply(vd, MARGIN=1, tabulate, nbins=C))
  colnames(ve) <- paste0("E",1:C) # 'E' for Enrolled
  vex <- rbindlist(rep(list(as.data.table(ve)), C))

  ## Before proceeding, I need to expand (explode!) vp vertically
  ## into a 'vx' that contains partial heads of all sequences.
  ####vx <- as.data.table(vp)
  vx <- rbindlist(rep(list(as.data.table(vp[,-ncol(vp)])), C))
  end <- nrow(vx)
  for (c in 2:C) { # TODO: This loop should zero the 've' copies also
    start <- end - (c-1)*P + 1
    col <- paste0("T",c)
    cat("Blanking column ", col, " from ", start, " to ", end, "\n")
    vx[[col]][start:end] <- NA
    vex[[paste0("E",c)]][start:end] <- 0
  }

  ## Now the trick is to tabulate a SUM of the T_ counts,
  ## on a PER DOSE basis. This requires using *both* sets
  ## of integers.
  dt <- array(as.matrix(vx), dim=c(nrow(vx),2,C))
  dimnames(dt) <- list(path=seq_len(nrow(vx)), c("D","T"), cohort=1:C)
  dt <- aperm(dt, c(3,2,1))

  ## dim(dt) is c(2,7,16384)
  ## Each plane dt(,,p) now should be subject to a computation
  ## that converts it to a dose-wise (D-)vector of toxicity counts.
  ## The essential computation is xtabs.
  ## TODO: Consider parallelizing or otherwise speeding this up!
  tc <- apply(dt, MARGIN=3, function(m) xtabs(T ~ factor(D, levels=1:7), data=m))
  tc <- t(tc)
  dimnames(tc) <- list(path=seq(nrow(tc)), dose=1:7)

  ## At this point, the index is formed from the pair cbind(rep(ve), tc)
  qh <- cbind(vex, tc)
  cat("There are ", nrow(unique(qh)), " unique quadrature indices.\n")

  qh
}


## What about this case?
## I can begin with the number of cohorts enrolled at each dose.
## Because cohort size is fixed in my current DTP application,
## this seems just fine.
## > vd[100,]
## D0 D1 D2 D3 D4 D5 D6
##  3  4  5  6  6  5  5
## > vt[100,]
## T1 T2 T3 T4 T5 T6 T7
##  0  0  0  1  2  0  3
## >
##
## If I want the toxicity counts summed by dose, then perhaps I
## ought to build a contingency matrix, and marginalize.
## One thing that is enormously helpful is converting doses to *factors*
## with levels 1:D.
##
## What about using a 3-dimensional array?
