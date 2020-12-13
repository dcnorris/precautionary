## Exact simulation for the 3+3 design

setOldClass(c("exact","three_plus_three_selector_factory","tox_selector_factory","selector_factory"))

#' A wrapper function supporting exact simulation of dose-escalation trials.
#' 
#' Implemented currently only for the (most?) common variant of the 3+3 design,
#' which requires that at least 6 patients be treated at a dose before it may be
#' declared as \sQuote{the} MTD.
#' 
#' @param selector_factory An object of type
#'  \code{\link[escalation:get_three_plus_three]{three_plus_three_selector_factory}},
#'  with \code{allow_deescalation = TRUE}.
#'
#' @details 
#' In any given realization of a 3+3 design, each of the \eqn{D} prespecified doses
#' will enroll 0, 1 or 2 cohorts, each with 3 patients. Each cohort will result in
#' a tally of 0--3 dose-limiting toxicities (DLTs), and these may be recorded in a
#' \eqn{2 \times D}{2 x D} matrix. Moreover, the 3+3 dose-escalation rules allow for
#' only one path through any such matrix. For example, the matrix
#' \preformatted{
#'    d
#' c   D1 D2 D3 D4
#'   1  0  1  2 NA
#'   2 NA  0 NA NA
#' }
#' represents the path in a 4-dose 3+3 trial, where the following events occur:
#' \enumerate{
#' \item Initial cohort at \eqn{d=1} results 0/3
#' \item Escalation to \eqn{d=2} results 1/3
#' \item Additional cohort at \eqn{d=2} results 0/3 for net 1/6 at this dose
#' \item Escalation to \eqn{d=3} results 2/3; MTD declared as \eqn{d=1}.
#' }
#' (Indeed, as you may verify at the R prompt, the above matrix is the 262nd of 442
#' such paths enumerated comprehensively in the \eqn{2 \times 4 \times 442}{2 x 4 x 442}
#' array \code{precautionary:::T[[4]]}.)
#' 
#' As detailed in Norris 2020c (below), these matrices may be used to construct simple
#' matrix computations that altogether eliminate the need for discrete-event simulation
#' of the 3+3 design. For each \eqn{D = 3,...,8}, the \code{precautionary} package has
#' pretabulated a \eqn{J \times 2D}{J x 2D} matrix \code{precautionary:::U[[D]]} and
#' \eqn{J}-vector \code{precautionary:::b[[D]]} such that the eqn{J}-vector \eqn{\pi}
#' of path probabilities is given by:
#' \deqn{
#' log(\pi) = b + U [log(p), log(q)]',
#' }
#' where \eqn{p} is the \eqn{D}-vector of DLT probabilities at the prespecified
#' doses, and \eqn{q \equiv 1-p}{q := 1-p} is its complement. See Eq. (4) of
#' Norris (2020c).
#' 
#' For details on the enumeration itself, please see the Prolog code in directory
#' \code{exec/} of the installed package.
#' @references 
#' Norris DC. What Were They Thinking? Pharmacologic priors implicit in a choice
#' of 3+3 dose-escalation design. arXiv:2012.05301 \[stat.ME\]. December 2020.
#' \url{http://arxiv.org/abs/2012.05301}
#' @examples 
#' # Run an exact version of the simulation from FDA-proactive vignette
#' design <- get_three_plus_three(
#'   num_doses = 6
#' , allow_deescalate = TRUE)
#' old <- options(
#'   dose_levels = c(2, 6, 20, 60, 180, 400)
#' , ordinalizer = function(MTDi, r0 = 1.5)
#'     MTDi * r0 ^ c(Gr1=-2, Gr2=-1, Gr3=0, Gr4=1, Gr5=2)
#' )
#' mtdi_gen <- hyper_mtdi_lognormal(CV = 0.5
#'                                 ,median_mtd = 180
#'                                 ,median_sdlog = 0.6
#'                                 ,units="ng/kg/week")
#' exact(design) %>% simulate_trials(
#'   num_sims = 1000
#' , true_prob_tox = mtdi_gen) -> EXACT
#' summary(EXACT)$safety
#' if (interactive()) { # runs too long for CRAN servers
#'   # Compare with discrete-event-simulation trials
#'   design %>% simulate_trials(
#'     num_sims = 1000
#'   , true_prob_tox = mtdi_gen) -> DISCRETE
#'   summary(DISCRETE)$safety[,]
#'   # Note the larger MCSEs in this latter simulation, reflecting combined noise
#'   # from hyperprior sampling and the nested discrete-event trial simulations.
#'   # The MCSE of the former simulation is purely from the hyperprior sampling,
#'   # since the nested trial simulation is carried out by an exact computation.
#' }
#' options(old)
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
  # Note that the following assumes ordinalizers multiply by a 5-vector
  # of the form c(r_2, r_1, 1, r1, r2), where r_2 < r_1 < 1 < r1 < r2.
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
  # If we do not bound CV away from zero, then it is possible to draw an mtdi_dist
  # that e.g. makes non-DLTs *impossible* at some high dose -- i.e., all toxicities
  # are of grade 3 or higher. In such extreme cases, the normalizing factor rowSums(G)
  # will have zero in components corresponding to such doses. (The converse may also
  # occur, in which low doses *cannot* cause a DLT.) This results in NaN's in the G
  # matrix, which cascade into results despite multiplication by zero probabilities
  # (corresponding to impossible events) in the U matrix. To avoid this, we set these
  # elements of the matrix to 0.
  G[is.nan(G)] <- 0
  return(G)
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

# TODO: Incorporate this logic to map 3+3 sim outcomes indices j in A[[D]][,,j]
#' @importFrom stats aggregate
haystack <- function(sims) {
  fit <- sims$fits[[1]][[1]]$fit
  m <- function(fit) {
    oc <- fit$outcomes
    ag <- aggregate(oc$tox, by=oc[,c("cohort","dose")], FUN=sum)
    ag$cohort <- NULL # drop this column
    # For any doses appearing just once, add an NA dose
    jo <- which(tabulate(ag$dose) == 1) # jo = 'just once'
    if (length(jo))
      ag <- rbind(ag, data.frame(dose = jo, x = NA))
    # For doses that never appeared, add 2 NAs
    doses <- 1:fit$num_doses
    nd <- doses[doses > max(ag$dose)]
    if (length(nd))
      ag <- rbind(ag, data.frame(dose = rep(nd, 2), x = NA))
    m <- matrix(ag[order(ag$dose),"x"], nrow = 2)
    dimnames(m) <- list(c = 1:2, d = paste0("D",doses))
    m
  }
  j <- function(m) {
    D <- ncol(m)
    J <- dim(A[[D]])[3]
    which(sapply(1:J, function(j)
      identical(m, A[[D]][,,j])))
  }
  j(m(fit))
}
