#' Complete dose transition pathways (DTP) table for the VIOLA trial.
#'
#' This data set caches a long computation (17 minutes on a 2.6GHz i7)
#' needed to build one package vignette.
#'
#' @format A data frame with 4^7 = 16384 rows and 15 columns, each row
#' representing one possible path the trial could take:
#' \describe{
#'   \item{D0}{Initial dose level}
#'   \item{T1}{Number of toxicities observed in first cohort}
#'   \item{D1}{Dose recommendation after the first cohort}
#'   \item{T2}{Number of toxicities observed in second cohort}
#'   \item{D2}{Dose recommendation after the second cohort}
#'   \item{T3}{Number of toxicities observed in third cohort}
#'   \item{D3}{Dose recommendation after the third cohort}
#'   \item{T4}{Number of toxicities observed in fourth cohort}
#'   \item{D4}{Dose recommendation after the fourth cohort}
#'   \item{T5}{Number of toxicities observed in fifth cohort}
#'   \item{D5}{Dose recommendation after the fifth cohort}
#'   \item{T6}{Number of toxicities observed in sixth cohort}
#'   \item{D6}{Dose recommendation after the sixth cohort}
#'   \item{T7}{Number of toxicities observed in seventh cohort}
#'   \item{D7}{Dose recommendation after the seventh cohort}
#' }
#'
#' @seealso
#' Documentation of the \code{\link[dtpcrm]{calculate_dtps}} function.
#'
"viola_dtp"
