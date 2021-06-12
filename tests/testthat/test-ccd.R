## Independently corroborate BOIN CPE done in Prolog ...

library(BOIN)

test_that("CPE for BOIN/CCD matches Prolog's BOIN25-<D>-6-24.tab",
{
  D <- 4 # TODO: Get 4 to work, also!
  bdy <- get.boundary(target = 0.25
                     ,ncohort = 24
                     ,cohortsize = 1
                     ,n.earlystop = 12
                     ,extrasafe = FALSE
                     ,offset = 0.05
                      )$boundary_tab

  boin <- Ccd$new(escalate = bdy[2,]
                 ,deescalate = bdy[3,]
                 ,eliminate = bdy[4,]
                 ,cohort_max = 6 # TODO: Try 9 if fast enough
                 ,enroll_max = 24
                  )$max_dose(D)

  paths <- boin$trace_paths(
                  root_dose = 1,
                  cohort_sizes = rep(1, 24),
                  unroll=1
                )$path_matrix()

  ## For comparison with outputs of boin.pl, we transform
  ## the dtpcrm-shaped table 'paths' to tables of REC (T N){D}.

  J <- nrow(paths)
  pathvectors <- matrix(nrow = J, ncol = 1+2*D)
  colnames(pathvectors) <- c("rec", paste0(c("T","N"), rep(1:D, each=2)))

  for (j in 1:J) {
    dtpvec <- paths[j,]
    DTmtx <- matrix(dtpvec[-length(dtpvec)], nrow=2)
    doses <- DTmtx[1,]
    pathvectors[j,'rec'] <- -doses[doses <= 0][1]
    DTmtx <- DTmtx[,!is.na(colSums(DTmtx)), drop=FALSE]
    ## Now DTmtx pairs doses (row 1) with tox counts (row 2)
    N <- tabulate(DTmtx[1,], nbins=D)
    ## Sum the toxicities indexed by dose, including zero counts for non-tried doses:
    dose <- factor(DTmtx[1,], levels=1:D)
    ntox <- DTmtx[2,]
    T <- as.vector(xtabs(ntox ~ dose))
    pathvectors[j,-1] <- as.vector(rbind(T,N))
  }

  ## Since BOIN obtains its dose recommendation from isotonic regression,
  ## we ignore the recommendations yielded by the dose-escalation process.
  ## We focus our checks instead solely on the final cumulative tallies:
  r <- as.data.frame(pathvectors[,-1])

  prolog <- read.table(system.file("testdata", sprintf("BOIN25-%d-6-24.tab", D)
                                    , package = "precautionary", mustWork = TRUE))[,-1]
  colnames(prolog) <- colnames(r)

  ## TODO: Generalize these constructions to arbitrary D
  rsols <- with(r, paste(paste(T1, N1, sep="/")
                        ,paste(T2, N2, sep="/")
                        ,paste(T3, N3, sep="/")
                        ,paste(T4, N4, sep="/")
                         ))
  psols <- with(prolog, paste(paste(T1, N1, sep="/")
                             ,paste(T2, N2, sep="/")
                             ,paste(T3, N3, sep="/")
                             ,paste(T4, N4, sep="/")
                               ))

  expect_identical(sort(rsols), sort(psols))

})
