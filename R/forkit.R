## Experiments with IPC for forked processes

library(parallel)

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

proglapply <- function(X, FUN
                     , parcellator = NULL # use default defined in body
                     , metric = length # progress metric on task results
                     , handler = NULL # use body's default progress handler
                     , ... # arguments to pass to progress reporting?
                     , workers = parallelly::availableCores(omit = 2)
                       ) {
  ## Parcel out the work to the workers
  N <- ((length(X) + workers - 1) %/% workers) * workers
  m <- matrix(1:N, nrow = workers)
  m[m > length(X)] <- NA # task matrix is generally 'ragged'

  ## Get workers started on their jobs, with per-task progress messaging
  jobs <- lapply(1:workers, function(w)
    mcparallel(lapply(X[na.exclude(m[w,])], function(x) {
      result <- FUN(x)
      parallel:::sendMaster(metric(result), raw.asis=FALSE)
      result
    }), name = w))
  rownames(m) <- sapply(jobs, function(.).$pid)

  cat("Work was parcelled out as follows:\n")
  print(m)

  ## We track progress reports due per worker, to avoid `readChild`-ing
  ## the final results, which would steal them from `mccollect`:
  progressReportsDue <- rowSums(!is.na(m)) # count non-NA tasks for each worker
  cat("Expecting the following progress reports:\n")
  print(progressReportsDue)

  cat("Watching progress notifications...\n")
  while (sum(progressReportsDue) > 0) {
    dueToReport = as.integer(names(progressReportsDue)[progressReportsDue > 0])
    ## TODO: The following is not restricting polling to the desired values!
    pids <- parallel:::selectChildren(dueToReport)
    if (!is.integer(pids))
      next
    msg <- parallel:::readChild(pids[1])
    if (is.raw(msg)) {
      pid <- attr(msg,'pid')
      if (!(pid %in% dueToReport)) {
        cat("ERROR!\n")
        print(dueToReport)
      }
      pid <- as.character(attr(msg,'pid'))
      cat("Received message from pid", pid, "\n")
      cat("Message says:", unserialize(msg), "\n")
      ##cat(sprintf("Sum: %d\n", sum)) # Here's where I'd do Shiny UI update!
      progressReportsDue[pid] <- progressReportsDue[pid] - 1
      print(progressReportsDue)
    } else { # this branch should never happen now
      if (is.integer(msg))
        stop("Received notice from exiting pid:", msg, "\n")
    }
  }

  mccollect(jobs)
}

## Can I now do same as above, only with fewer, composite jobs?
## The key to this may be sendMaster(), which replaces the fifo
## and allows more direct communication of progress.
sumComp <- function(n=10, mc.cores = 5) {
  ## Here's the 'work' to be done in parallel:
  work <- function(x) {
    pid <- Sys.getpid()
    cat(sprintf("pid %d received work list of length %d\n", pid, length(x)))
    x2 <- numeric(length(x))
    for (i in seq(length(x))) {
      ##cat(sprintf("%d sleeping for %d seconds...\n", pid, 1))
      Sys.sleep(1)
      ##cat(sprintf("%d waking up...\n", pid))
      ## Contribute to the answer...
      x2[i] <- x[[i]]^2
      ##cat(sprintf("%d sending {%d} to master.\n", pid, x2[i]))
      parallel:::sendMaster(x2[i], raw.asis=FALSE)
    }
    x2 # Does child process cease to exist after returning this value?
  }

  ## Assign work to `mc.cores` jobs. (In applications, this could be
  ## a point where good heuristics can capture non-negligible gains.)
  m <- matrix(1:n, nrow = mc.cores)
  print(m)

  ## Use the 'naive parallel lapply' from `?mccollect`
  tasks <- as.list(1:n)
  jobs <- lapply(1:mc.cores, function(p) mcparallel(work(tasks[m[p,]]), name = p))

  rownames(m) <- sapply(jobs, function(job) job$pid)
  cat("The following jobs were created:")
  print(m)

  progressReportsDue <- rowSums(!is.na(m))
  cat("Expecting the following progress reports:\n")
  print(progressReportsDue)

  cat("Watching results come in...\n")

  ## Aha! At this point, I need to know who got what work,
  ## and stop polling workers once they've signalled progress
  ## on each assigned item. This way, workers are left 'holding
  ## the bag' until mccollect stage!

  tasksDone <- 0L
  while (sum(progressReportsDue) > 0) {
    dueToReport = as.integer(names(progressReportsDue)[progressReportsDue > 0])
    ## TODO: The following is not restricting polling to the desired values!
    pids <- parallel:::selectChildren(dueToReport)
    if (!is.integer(pids))
      next
    msg <- parallel:::readChild(pids[1])
    if (is.raw(msg)) {
      pid <- attr(msg,'pid')
      if (!(pid %in% dueToReport)) {
        cat("ERROR!\n")
        print(dueToReport)
      }
      pid <- as.character(attr(msg,'pid'))
      cat("Received message from pid", pid, "\n")
      cat("Message says:", unserialize(msg), "\n")
      ## TODO: Figure out how to read an integer msg!
      ##cat(sprintf("Sum: %d\n", sum)) # Here's where I'd do Shiny UI update!
      tasksDone <- tasksDone + 1
      progressReportsDue[pid] <- progressReportsDue[pid] - 1
    print(progressReportsDue)
      cat(sprintf("Tasks complete so far: %d\n", tasksDone))
    } else { # this branch should never happen now
      if (is.integer(msg))
        cat("Received notice from exiting pid:", msg, "\n")
    }
  }

  Sys.sleep(0.1)
  cat("Collecting results...\n")
  mccollect(jobs)
}
