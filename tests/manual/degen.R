## Explorations of quadrature degeneracy in DTP computation

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

  ## Before proceeding, I need to expand (explode!) vp vertically
  ## into a 'vx' that contains partial heads of all sequences.
  ####vx <- as.data.table(vp)
  vx <- rbindlist(rep(list(as.data.table(vp[,-ncol(vp)])), C))
  end <- nrow(vx)
  for (c in 2:C) {
    start <- end - (c-1)*P + 1
    col <- paste0("T",c)
    cat("Blanking column ", col, " from ", start, " to ", end, "\n")
    vx[[col]][start:end] <- NA
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
  vex <- rbindlist(rep(list(as.data.table(ve)), C))
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
