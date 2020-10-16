## Exact simulation for the 3+3 design

setOldClass(c("exact","three_plus_three_selector_factory","tox_selector_factory","selector_factory"))

#' A wrapper function supporting exact simulation of 3+3 design
#' 
#' @param selector_factory An object of type
#'  \code{\link[escalation:get_three_plus_three]{three_plus_three_selector_factory}},
#'  with \code{allow_deescalation = TRUE}.
#'
#' @export
exact <- function(selector_factory) {
  stopifnot("Class 'exact' applies only to a three_plus_three_selector_factory"
            = is(selector_factory,"three_plus_three_selector_factory"))
  stopifnot("Exact sim is currently implemented only for 3+3 with 'allow_deescalate=TRUE'"
            = selector_factory$allow_deescalate)
  prependClass("exact", selector_factory)
}

# Calculate the G matrix from an mtdi_dist
# NB: The ... may be used to pass r0 thru to ordinalizer
G <- function(mtdi_dist, ordinalizer = getOption("ordinalizer"), ...) {
  # 1. Get a (D x 5) matrix of doses rescaled by r0^(3-g), g in 1:5
  # 2. Get the CDF values on this matrix
  r0_spread <- 1/ordinalizer(1, ...)
  G5 <- mtdi_dist@dist$cdf(getOption("dose_levels") %o% r0_spread)
  # 3. Attach trivial columns of 1's and 0's; label for clarity
  G7 <- cbind(None=rep(1,nrow(G5)), G5, rep(0,nrow(G5)))
  dimnames(G7)[1] <- list(D = paste(getOption('dose_levels'), mtdi_dist@units))
  # The (now) 7 columns give the probabilities of toxicity Gr=0+, Gr=1+, ..., Gr=6+.
  # The trivial leftmost and rightmost columns assist with difference-taking below.
  # 4. By taking differences, we obtain probabilities for the grades *categorically*
  G <- G7[,1:6] - G7[,2:7]
  # 5. Form the (2D x 6) blocked matrix (note, this remains a very small matrix!)
  Z <- matrix(0, nrow = nrow(G), ncol = 3)
  dimnames(Z) <- dimnames(G[,1:3])
  G <- rbind(cbind(Z, G[,4:6]),
             cbind(G[,1:3], Z))
  # 6. Normalize the rows
  G <- G / rowSums(G)
  return(G)
}

# Calculate the length-2D vector c(log(p), log(1-p))
log_pq <- function(mtdi_dist) {
  p <- mtdi_dist@dist$cdf(getOption("dose_levels"))
  q <- 1 - p
  log_pq <- c(log(p), log(q))
  names(log_pq) <- rep(paste(getOption('dose_levels'), mtdi_dist@units), 2)
  log_pq
}

# Exact calculation of safety table
exact_safety <- function(mtdi_dist, ordinalizer = getOption("ordinalizer"), ...) {
  D <- length(getOption("dose_levels"))
  log_pi <- b[[D]] + U[[D]] %*% log_pq(mtdi_dist)
  safety <- t(exp(log_pi)) %*% U[[D]] %*% G(mtdi_dist, ...)
}

#' @importFrom escalation prob_recommend
prob_recommend.exact <- function(x, ...) {
  prob_recs <- NextMethod()
  names(prob_recs)[-1] <- paste0(x$dose_levels, x$dose_units)
  prob_recs
}

#' @importFrom escalation prob_administer
prob_administer.exact <- function(x, ...) {
  prob_admin <- NextMethod()
  names(prob_admin) <- paste0(x$dose_levels, x$dose_units)
  prob_admin
}

#' Specialize print method for objects of class \code{\link[escalation]{simulations}}
#'
#' @param x An object of class  c("exact","precautionary","simulations")
#'
#' @param ... Additional arguments; ignored
#'
#' @importFrom escalation num_patients num_tox trial_duration
print.exact <- function(x, ...) {
  
  cat('Number of iterations:', length(x$fits), '\n')
  cat('\n')
  
  cat('Number of doses:', length(x$dose_levels), '\n')
  cat('\n')
  
  cat('True probability of toxicity:\n')
  print(x$true_prob_tox, digits = 3)
  cat('\n')
  
  cat('Probability of recommendation:\n')
  print(prob_recommend(x), digits = 3)
  cat('\n')
  
  cat('Probability of administration:\n')
  print(prob_administer(x), digits = 3)
  cat('\n')
  
  cat('Sample size:\n')
  print(summary(num_patients(x)))
  cat('\n')
  
  cat('Total toxicities:\n')
  print(summary(num_tox(x)))
  cat('\n')
  
  cat('Trial duration:\n')
  print(summary(trial_duration(x)))
  cat('\n')
}


#' Specialize print method defined for class \code{\link[escalation]{simulations}}
#'
#' @param x An object of class  c("hyper","precautionary","simulations")
#'
#' @param ... Additional arguments; ignored
#'
#' @importFrom escalation num_patients num_tox trial_duration
print.hyper <- function(x, ...) {
  
  cat('Number of iterations:', length(x$fits), '\n')
  cat('\n')
  
  cat('Number of doses:', length(x$dose_levels), '\n')
  cat('\n')
  
  cat('Average probability of toxicity:\n')
  print(x$avg_prob_tox, digits = 3)
  cat('\n')
  
  cat('Probability of recommendation:\n')
  print(prob_recommend(x), digits = 3)
  cat('\n')
  
  cat('Probability of administration:\n')
  print(prob_administer(x), digits = 3)
  cat('\n')
  
  cat('Sample size:\n')
  print(summary(num_patients(x)))
  cat('\n')
  
  cat('Total toxicities:\n')
  print(summary(num_tox(x)))
  cat('\n')
  
  cat('Trial duration:\n')
  print(summary(trial_duration(x)))
  cat('\n')
}


#' Convert an object of class c('exact','precautionary','simulations') to a data.table
#'
#' TODO: Actually implement this!
#' @param x An object of class c('precautionary','simulations')
#'
#' @param keep.rownames Unused; retained for S3 generic/method consistency
#' @param ordinalizer If not NULL, this is a function mapping the threshold
#'  dose ('MTDi') at which an individual experiences a binary toxicity (as
#'  recognized by the dose-escalation design) to a named vector giving dose
#'  thresholds for multiple grades of toxicity. The names of this vector will
#'  be taken as designations of the toxicity grades.
#' @param ... Additional parameters passed to the \code{ordinalizer}
#'
#' @export
as.data.table.exact <- function(x, keep.rownames = FALSE
                                        , ordinalizer = getOption('ordinalizer')
                                        , ...) {
  extractor <- ifelse(is(x$fits[[1]][[1]]$fit, "derived_dose_selector")
                     ,function(.) .[[1]]$fit$parent$outcomes
                     ,function(.) .[[1]]$fit$outcomes
                     )
  ensemble <- rbindlist(lapply(x$fits, extractor), idcol = "rep")
  if( is.null(ordinalizer) )
    return(ensemble)
  # TODO: Do add these columns to 'ensemble', so that the whole table
  #       may later be inspected by user to improve understanding.
  MTDig <- t(sapply(ensemble$MTDi, ordinalizer, ...))
  if( is.null(colnames(MTDig)) ) {
    warning("Ordinalizer returns unnamed vector; using default names for toxicity grades.")
    colnames(MTDig) <- paste("Grade", 1:ncol(MTDig))
  }
  tox_grades <- colnames(MTDig)
  # Compare with actual dose to obtain toxicity grade indicator matrix
  tox_ind <- ( x$dose_levels[ensemble$dose] > MTDig )
  # Tally the thresholds crossed to obtain integer toxgrade
  ensemble$toxgrade <- rowSums(tox_ind)
  # Convert toxgrade to an ordered factor Tox
  ensemble$Tox <- ordered(ensemble$toxgrade+1
                          , levels=seq(1+length(tox_grades))
                          , labels=c('None', tox_grades))
  ensemble
}

#' Summarize an exact treatment of a dose-escalation design
#' 
#' Algorithmic (or 'rule-based') dose-escalation designs admit exact computation
#' of their outcomes. This method summarizes such an exact treatment in a manner
#' roughly parallel to that of \code{summary.precautionary}.
#' 
#' @param object An object of class 'exact' 
#'
#' @param ordinalizer An ordinalizer function
#' @param ... Additional parameters passed to the ordinalizer
#'
#' @importFrom dplyr mutate rename_with select everything
#' @importFrom stats xtabs addmargins sd var
#' @importFrom rlang .data
#' @export
summary.exact <- function(object, ordinalizer = getOption('ordinalizer'), ...) {
  # NB: We are severing the call chain, with a special implementation here
  ##summary <- NextMethod()
  # TODO: Consider implementing the 'escalation' & 'toxTab' components.
  #       (I believe both are obtainable in the exact formulation.)
  # I need to deal with the 'hyper' case separately.
  summary <- list(escalation = "not yet implemented",
                  safety = NULL, # to be calculated ...
                  toxTab = "not yet implemented")
  if (is(object, 'hyper')) {
    toxTab <- object$hyper$safety %>%
      addmargins(margin = 2, FUN = list(Total=sum))
    summary$safety <- rbind(
      "Expected participants" = colMeans(toxTab)
      ,"MCSE" = apply(toxTab, MARGIN = 2, FUN = sd) / sqrt(nrow(toxTab))
    )
  } else {
    summary$safety <- object$safety %>%
      addmargins(margin = 2, FUN = list(Total=sum))
  }
  summary
}

