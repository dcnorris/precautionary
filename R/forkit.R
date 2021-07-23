## Experiments with IPC for forked processes

library(parallel)

## Let's have a test using the target computation!
crm_b2020 <- function() {
  d1_maxn <- 5
  cum_maxn <- 10
  calmod <- Crm$new(skeleton = c(0.03, 0.11, 0.25, 0.42, 0.58, 0.71),
                  scale = 0.85, # aka 'sigma'
                  target = 0.25)$
    no_skip_esc(TRUE)$    # compare Braun's 'restrict = T'
    no_skip_deesc(FALSE)$
      stop_func(function(x) {
        enrolled <- tabulate(x$level, nbins = length(x$prior))
        x$stop <- enrolled[1] >= d1_maxn || max(enrolled) >= cum_maxn
        x
      })
}

testprog <- function(C = 13, ...) {
  mod <- crm_b2020()
  mod$trace_paths(1, rep(2, C), unroll = 4, ...)
  cat("Done.\n")
  mod$performance
}


## This design begins to look promising. But I need a clear nomenclature!
## Let me say that a body of 'work' is parallelized into 'tasks' that are
## 'parceled' into 'jobs', and that these jobs are 'assigned' to 'workers'.
## We can also say the tasks are 'parceled out' to the workers.

## TODO: Should I reserve 'job' for work that is underway by a worker?

## This is represented by a 'task matrix' with workers (identified by pid)
## assigned (as rownames) to rows that represent the parcels as sequences
## of 'task indexes' (possibly NA) into the 'task list'.

## These matrices, although 'ragged' and generally quite likely to contain
## NA indexes (because work can't be evenly divided), do give handy means
## to represent and compute with the parcellation.

#' @importFrom stats na.exclude
#' @importFrom parallel mcparallel mccollect
#' @importFrom utils getFromNamespace
proglapply <- function(X, FUN
                     , parcellator = NULL # use default defined in body
                     , metric = length # progress metric on task results
                     , init = 0L
                     , prog = function(p) cat(sprintf("progress so far: %d\n", p))
                     , workers = parallelly::availableCores(omit = 2)
                       ) {

  sendMaster <- getFromNamespace('sendMaster', 'parallel')
  selectChildren <- getFromNamespace('selectChildren', 'parallel')
  readChild <- getFromNamespace('readChild', 'parallel')

  ## Parcel out the work to the workers
  N <- ((length(X) + workers - 1) %/% workers) * workers
  m <- matrix(1:N, nrow = workers)
  m[m > length(X)] <- NA # task matrix is generally 'ragged'

  ## Get workers started on their jobs, with per-task progress messaging
  jobs <- lapply(1:workers, function(w)
    mcparallel(lapply(X[na.exclude(m[w,])], function(x) {
      result <- FUN(x)
      sendMaster(metric(result), raw.asis=FALSE)
      result
    })))
  rownames(m) <- sapply(jobs, function(.).$pid)

  cat("Work was parcelled out as follows:\n")
  print(m)

  ## We track progress reports due per worker, to avoid `readChild`-ing
  ## the final results, which would steal them from `mccollect`:
  progressReportsDue <- rowSums(!is.na(m)) # count non-NA tasks for each worker
  cat("Expecting the following progress reports:\n")
  print(progressReportsDue)

  cat("Watching progress notifications...\n")
  progress <- init
  while (sum(progressReportsDue) > 0) {
    dueToReport = names(progressReportsDue)[progressReportsDue > 0]
    ## TODO: The following is not restricting polling to the desired values!
    pids <- selectChildren(as.integer(dueToReport))
    if (!is.integer(pids))
      next
    msg <- readChild(pids[1])
    if (is.raw(msg)) {
      pid <- as.character(attr(msg,'pid'))
      if (!(pid %in% dueToReport)) {
        cat("ERROR!\n")
        print(dueToReport)
      }
      progressReportsDue[pid] <- progressReportsDue[pid] - 1
      progress <- progress + unserialize(msg)
      prog(progress)
    } else { # this branch should never happen now
      if (is.integer(msg))
        stop("Received notice from exiting pid:", msg, "\n")
    }
  }

  do.call(c, unname(mccollect(jobs))) # unname() to strips pids
}
