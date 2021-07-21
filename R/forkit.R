## Experiments with IPC for forked processes

## See https://stackoverflow.com/a/27729791/3338147 by @fotNelton

library(parallel)

## TODO: I'd like to do something like sumSquaresTo below,
##       except that the polling loop should be in the main
##       process, and the main process should manage the UI
##       progress reporting, then (when work done) receive
##       the results of mclapply().
##
## This may very well be too much to ask! Perhaps I must do
## the work of mclapply myself, using mcparallel + mccollect.

sumSqool <- function(n) {
  f <- fifo(tempfile(), open="w+b", blocking=TRUE)
  on.exit(close(f))

  ## Here's the 'work' to be done in parallel:
  work <- function(x) {
    Sys.sleep(x/2)
    # Contribute to the answer...
    x2 <- x^2
    writeBin(x2, f)
    x2
  }

  ## Use the 'naive parallel lapply' from `?mccollect`
  jobs <- lapply(1:n, function(x) mcparallel(work(x), name = x))

  cat("Watching results come in...\n")
  sum <- 0
  jobsDone <- 0L
  while (jobsDone < length(jobs) && !isIncomplete(f)) {
    msg <- readBin(f, "double")
    sum <- sum + as.numeric(msg)
    cat(sprintf("Sum: %d\n", sum)) # Here's where I'd do Shiny UI update!
    jobsDone <- jobsDone + 1
  }

  cat("Collecting results...\n")
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


## Add the numbers (1:n)^2 by posting them to a FIFO

sumSquaresTo <- function(n) {
  f <- fifo(tempfile(), open="w+b", blocking=TRUE)

  ## Anticipating that the parent process will suspend while
  ## waiting for the mclapply call below to return, we spawn
  ## a child dedicated to collecting the results from those
  ## parallel workers:
  if (inherits(parallel:::mcfork(), "masterProcess")) { # 'spawning' idiom
    ## Child 'collector' process
    sum <- 0
    answer <- sum((1:n)^2) # TODO: How can I avoid needing to know this?
    while (sum < answer && !isIncomplete(f)) {
      msg <- readBin(f, "double")
      sum <- sum + as.numeric(msg)
      cat(sprintf("Sum: %d\n", sum)) # Q: Can a child process update Shiny UI?
      ## I rather suspect that only the main event-loop can update the
      ## Shiny UI for me. So perhaps I need to handle this collecting
      ## task in the main function, and delegate rather the mclapply()
      ## invocation to a child?
    }
    ## Q: What is the role of !isIncomplete(f) above? Presumably,
    ##    the while loop should exit in case isIncomplete(f); why?
    ##    The documentation states that isIncomplete() returns
    ##    "whether the last read attempt was blocked," which here
    ##    would presumably be the readBin(f) from the previous
    ##    iteration? What would it mean that this read blocked?
    ##    Would that mean the fifo has been closed?
    ##    Or does this possibly serve to catch exceptional conditions,
    ##    related to the 'adjustments' to blocking described here?
    ##
    ## > Blocking:
    ## >   ...
    ## >   In blocking mode, functions using the connection do not return to
    ## >   the R evaluator until the read/write is complete.  In non-blocking
    ## >   mode, operations return as soon as possible, so on input they will
    ## >   return with whatever input is available (possibly none) and for
    ## >   output they will return whether or not the write succeeded.
    ## >   ...
    ## >   Even when a connection is in blocking mode, attempts are made to
    ## >   ensure that it does not block the event loop and hence the
    ## >   operation of GUI parts of R.  These do not always succeed, and the
    ## >   whole R process will be blocked during a DNS lookup on Unix, for
    ## >   example.
    parallel:::mcexit()
  }

  ## Now the main work begins:
  numJobs <- n
  result <- mclapply(1:numJobs, function(x) {
    Sys.sleep(1)
    # Contribute to the answer...
    x2 <- x^2
    writeBin(x2, f)
    x2
  }
  , mc.cores = 5)
  close(f)
  result
}


## TODO: Why not use non-blocking fifos, per the default behavior?
##       Given that my message sizes are 'atomic' (just 1 integer!),
##       don't I have a pretty good shot at making this work?
