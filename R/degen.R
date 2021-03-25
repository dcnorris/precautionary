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
## in the same encoding! I may have just enough bits to manage that!)
encode_cohorts <- function(enr, tox) {
  ## 'enr' is a vector of cohort counts, ordered dose-wise
  stopifnot(all(enr %% 1 == 0))
  stopifnot(length(enr) <= 8)
  stopifnot(all(enr >= 0))
  stopifnot(all(enr <= 5))
  ## 'tox' is a integer vector of dose-wise toxicity counts
  stopifnot(all(tox %% 1 == 0))
  stopifnot(length(tox) <= 8)
  stopifnot(all(tox %in% 0:15))

  ## Extend tox & enr to length-8 vectors
  enr8 <- tox8 <- integer(8)
  enr8[seq_along(enr)] <- enr
  tox8[seq_along(tox)] <- tox

  ## Since the base-16 encoding of tox will not spawn continued binary fractions,
  ## we put tox in the fractional part. This should enable recovery of more decimal
  ## digits of representations from the decoder. (A cosmetic issue, I believe.)
  x <- sum(enr8 * 6^(0:7)) + sum(tox8 * 16^-(1:8))

  check <- decode_cohorts(x)
  if (all(check$tox == tox8))
    if (all(check$enr == enr8))
      return(x)
    else
      message("enr = ", enr8, " != check$enr = ", check$enr)
  else
    message("tox = ", tox8, " != check$tox = ", check$tox)

  stop("INCONCEIVABLE!") # This encoding 'cannot fail'
  return(NULL)

}

decode_cohorts <- function(x) {
  enr <- sapply(1:8, function(d) (x %% 6^d) %/% 6^(d-1)) # base-6 decoding
  .x <- x %% 1 # fractional part
  tox <- sapply(-(0:7), function(d) (.x %% 16^d) %/% 16^(d-1)) # base-16 decoding

  list(enr = enr,
       tox = tox)
}

## Encode CRM obs expressed as a vector of interleaved dose assignments and tox counts.
encode_dose_ntox <- function(dt, interleaved=TRUE) {
  ## dt is an integer vector of (dose,ntox) pairs
  stopifnot(all(dt %% 1 == 0, na.rm=TRUE))
  dt <- as.integer(dt)
  stopifnot(length(dt) %% 2 == 0) # must have even number of elements
  stopifnot(length(dt) %/% 2 <= 8)

  if (interleaved) {
    dt <- t(matrix(dt, nrow=2))
  } else {
    dt <- matrix(dt, ncol=2)
  }
  dimnames(dt) <- list(cohort = 1:nrow(dt), c("D","T"))
  df <- as.data.frame(dt)          # TODO: See whether sticking with matrix
  df$D <- factor(df$D, levels=1:8) #       speeds this up when applied.

  enr <- tabulate(df$D, nbins=8)
  tox <- as.integer(xtabs(T ~ D, data=df))
  ## list(enr = enr,
  ##      tox = tox,
  ##      dt  = t(dt), # 2 x 8 matrix makes for easier viewing
  ##      enc = encode_cohorts(enr, tox))

  encode_cohorts(enr, tox)
}

quadcodes <- function(dtp) {
  qc <- matrix(nrow = nrow(dtp),
               ncol = ncol(dtp) %/% 2)
  for (c in seq.int(ncol(qc))) {
    qc[,c] <- apply(dtp[,1:(2*c)], MARGIN=1, encode_dose_ntox)
  }
  qc
}

## I find that I get 3041 unique codes from all steps of all rows of viola_dtp.
## Furthermore, at the console I find:
## > path_length <- apply(viola_dtp[,paste0("D",0:6)], MARGIN=1, function(.) sum(!is.na(.)))
## > xtabs(~path_length)
## path_length
##    2    3    4    5    6    7
## 4096 2560 1792 2272 1540 4124
##
## The transcendental burden of quadratures along a path of length l is (l+1)*l/2.
## > n <- xtabs(~path_length)
## > l <- 2:7
## > (l+1)*l/2 -> txops
## > txops
## [1]  3  6 10 15 21 28
## > sum(n * txops)
## [1] 227460
##
## Even if the cost of every one of the 3041 unique quadratures were 28 (the max!),
## I would still get a 2.67 speedup!
## > 227460 / (28*3041)
## [1] 2.671349

## TODO: Note that the quadcode can also be used to figure out the txcost!
##       Thus, I can do this calculation in a quite refined manner.
txcost <- function(enc) {
  coh <- decode_cohorts(enc)
  ## Let's revisit my former tally of txops!
  ##
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
  tc <<- tc
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
