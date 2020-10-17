#' @name precautionary-package
#' @title Safety Diagnostics for Dose-Escalation Trial Designs
#' @description Enhances various R packages that support the design and simulation
#' of phase 1 dose-escalation trials, adding diagnostics to examine the safety
#' characteristics of these designs in light of expected inter-individual variation
#' in pharmacokinetics and pharmacodynamics.
#' @aliases precautionary-package precautionary
#' @docType package
#' @author David C. Norris (<https://orcid.org/0000-0001-9593-6343>)
#' 
#' @references
#' 1. Norris DC. Dose Titration Algorithm Tuning (DTAT) should supersede
#'    \sQuote{the} Maximum Tolerated Dose (MTD) in oncology dose-finding trials.
#'    \emph{F1000Research}. 2017;6:112. \doi{10.12688/f1000research.10624.3}.
#'    \url{https://f1000research.com/articles/6-112/v3}
#' 
#' 2. Norris DC. Costing \sQuote{the} MTD. \emph{bioRxiv}. August 2017:150821.
#'    \doi{10.1101/150821}.
#'    \url{https://www.biorxiv.org/content/10.1101/150821v3}
#' 
#' 3. Norris DC. Precautionary Coherence Unravels Dose Escalation Designs.
#'    \emph{bioRxiv}. December 2017:240846. \doi{10.1101/240846}.
#'    \url{https://www.biorxiv.org/content/10.1101/240846v1}
#' 
#' 4. Norris DC. One-size-fits-all dosing in oncology wastes money, innovation
#'    and lives. \emph{Drug Discov Today}. 2018;23(1):4-6.
#'    \doi{10.1016/j.drudis.2017.11.008}.
#'    \url{https://www.sciencedirect.com/science/article/pii/S1359644617303586}
#' 
#' 5. Norris DC. Costing \sQuote{the} MTD ... in 2-D. \emph{bioRxiv}. July 2018:370817.
#'    \doi{10.1101/370817}.
#'    \url{https://www.biorxiv.org/content/10.1101/370817v1}
#' 
#' 6. Norris DC. Ethical Review and Methodologic Innovation in Phase 1 Cancer Trials.
#'    \emph{JAMA Pediatrics}. 2019;173(6):609
#'    \doi{10.1001/jamapediatrics.2019.0811}.
#'    
#' 7. Norris DC. Comment on Wages et al., Coherence principles in interval-based
#'    dose finding. Pharmaceutical Statistics 2019, DOI: 10.1002/pst.1974.
#'    \emph{Pharmaceutical Statistics}. March 2020.
#'    \doi{10.1002/pst.2016}.
#'    
#' 8. Norris DC. Retrospective analysis of a fatal dose-finding trial.
#'    arXiv:2004.12755 \[stat.ME\]. April 2020.
#'    \url{https://arxiv.org/abs/2004.12755}
#' 
#' @import methods
#' @import magrittr
#' @import data.table
NULL

#' Ordinalizers
#' 
#' TODO: Explain ordinalization and ordinalizers here.
#' 
#' @section Limitations:
#' TODO: Be sure to note the highly restrictive assumptions made about
#' ordinalizer functions, especially as noted e.g. in connection with
#' the internal function G (defined in 'exact.R').
#' 
#' @section Validity:
#' The validity of an ordinalizer might well be programmatically testable.
#' If so, all this discussion might well be carried out in documentation
#' for a validity-testing function.
#' 
#' @name ordinalizer
#' @aliases ordinalizers ordinalization
NULL
