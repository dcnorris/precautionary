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
#'    *F1000Research*. 2017;6:112. \doi{10.12688/f1000research.10624.3}.
#'    <https://f1000research.com/articles/6-112/v3>
#'
#' 2. Norris DC. Costing \sQuote{the} MTD. *bioRxiv*. August 2017:150821.
#'    \doi{10.1101/150821}.
#'    <https://www.biorxiv.org/content/10.1101/150821v3>
#'
#' 3. Norris DC. Precautionary Coherence Unravels Dose Escalation Designs.
#'    *bioRxiv*. December 2017:240846. \doi{10.1101/240846}.
#'    <https://www.biorxiv.org/content/10.1101/240846v1>
#'
#' 4. Norris DC. One-size-fits-all dosing in oncology wastes money, innovation
#'    and lives. *Drug Discov Today*. 2018;23(1):4-6.
#'    \doi{10.1016/j.drudis.2017.11.008}.
#'    <https://www.sciencedirect.com/science/article/pii/S1359644617303586>
#'
#' 5. Norris DC. Costing \sQuote{the} MTD ... in 2-D. *bioRxiv*. July 2018:370817.
#'    \doi{10.1101/370817}.
#'    <https://www.biorxiv.org/content/10.1101/370817v1>
#'
#' 6. Norris DC. Ethical Review and Methodologic Innovation in Phase 1 Cancer Trials.
#'    *JAMA Pediatrics*. 2019;173(6):609
#'    \doi{10.1001/jamapediatrics.2019.0811}.
#'
#' 7. Norris DC. Comment on Wages et al., Coherence principles in interval-based
#'    dose finding. Pharmaceutical Statistics 2019, DOI: 10.1002/pst.1974.
#'    *Pharmaceutical Statistics*. March 2020.
#'    \doi{10.1002/pst.2016}.
#'
#' 8. Norris DC. Retrospective analysis of a fatal dose-finding trial.
#'    arXiv:2004.12755 \[stat.ME\]. April 2020.
#'    <https://arxiv.org/abs/2004.12755>
#'
#' 9. Norris DC. What Were They Thinking? Pharmacologic priors implicit in a choice
#'    of 3+3 dose-escalation design. arXiv:2012.05301 \[stat.ME\]. December 2020.
#'    <https://arxiv.org/abs/2012.05301>
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

#' Plan
#'
#' @section Version 0.2.6-1:
#' * Improve feedback during long calculations
#'   - Set J = '...' pending first progress report
#'   - Analyze where `Cpe$path_array()` spends time
#'   - Give suitable feedback on `path_array` progress
#'
#' @section Version 0.2.7:
#' * Demote \pkg{escalation} to a 'Suggests'
#' * Thoroughly rewrite 'Intro' and 'FDA Clinical Hold' vignettes
#'
#' @section Version 0.2.8:
#' * BOIN recs via isotonic regression
#' * Testthat CPE matches BOIN::get.oc() sim
#' * Implement CRM logistic model
#'
#' @section Version 0.2.9:
#' * Implement TITE CRM?
#'   - Is TITE even amenable to path enumeration?
#'   - If not, there may be little sense in implementing it here!
#'   - If TITE is truly off the table, can existing numerics be
#'     sped up further?
#'   - OTOH, implementing TITE would help place the numerics
#'     for enumerable designs in context.
#'
#' @section Version 0.3.0:
#' * Stop exposing the `impl` parameter
#'
#' @section Version 0.3.1:
#' * Native Rust CPE
#'
#' @section Dependencies:
#' * Eliminate dependence on/adherence to dtpcrm layout
#'   - Package dtpcrm makes many design decisions unsuited to
#'     comprehensive enumeration of whole trials.
#'   - Allowing the Crm class to evolve along now-'natural' lines
#'     probably makes a suitable 'plan' for now.
#'   - Expunging superseded code/dependencies ASAP will facilitate
#'     this evolution.
#'   - Even the term 'DTP' seems not quite right anymore for a complete
#'     path enumeration (CPE?), and should be abandoned.
#' * Eliminate dependence on 'escalation' and 'dfcrm'
#'   - With the special emphasis on *speed* created by CPE,
#'     I now must implement all underlying trial designs in Rust.
#' * What may be retained in each case are 'Suggests:'-type relations,
#'   wherever regression tests are helpful -- esp. wrt dfcrm.
#'
#' @section Fast CRM:
#' * Benchmark; try mul_add()'s
#'
#' @section Document:
#' * Add examples to the documented Rust functions
#' * Add a vignette applying DTP to mTPI and BOIN
#' * Expose some visual numerics checks, via vignette or documented function
#'
#' @section Refactor:
#' * Remove the `$safety` component of exact trials?
#'   - Perhaps this ought to be calculated 'on the fly' by the summary method.
#'   - On-the-fly calculation would postpone use of ordinalizer, in keeping with
#'     the pattern established already for (non-exact) simulations.
#' * Should summary(EXACT)$safety bear class 'safetytab'?
#' * Example for 'as.data.table.exact'
#' * What is role of G function in exact.R?
#'
#' @section Extend:
#' * Implement exact 3+3 variant with `allow_deescalation=FALSE`
#' * Implement rolling 6
#' * Allow an accelerated titration phase
#'   - Note that this requires access to graded toxicities at simulation time,
#'     and therefore constitutes a substantial challenge to the generalizability
#'     of this software design.
#' * Index `sims$fits` to exact outcomes in `A[[D]]` where appropriate
#'   - See the `haystack` function in `exact.R`
#'
#' @section Robustify:
#' * Tests comparing results from multiple CRAN packages
#'
#' @name plan
#' @aliases todo
NULL
