## Test environments
* local OS X install, R 4.0.1
* win-builder (devel)

## R CMD check --as-cran results
There were no ERRORs or WARNINGs

There was 1 NOTE:

* checking dependencies in R code ... NOTE
  Unexported object imported by a ':::' call: ‘escalation:::spruce_outcomes_df’
    See the note in ?`:::` about the use of this operator.

  This unexported function from 'escalation' v0.1-3 is defined as follows:
  
  # Make sure columns in an outcomes data-frame are integers.
  spruce_outcomes_df <- function(df) {
    df$dose <- as.integer(df$dose)
    df$tox <- as.integer(df$tox)
    if('cohort' %in% colnames(df)) df$cohort <- as.integer(df$cohort)
    if('patient' %in% colnames(df)) df$patient <- as.integer(df$patient)
    df
  }
  
  It would be easy enough to provide a duplicated version in 'precautionary',
  but I think it more sensible to accept this NOTE as a reminder of a 'wart'
  to be resolved through collaboration with the author of 'escalation'. 
