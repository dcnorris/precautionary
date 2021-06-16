# precautionary <img src="man/figures/logo.svg" align="right" alt="LOGO" width="120" />

<!-- badges: start -->

[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![R-CMD-check](https://github.com/dcnorris/precautionary/workflows/R-CMD-check/badge.svg)](https://github.com/dcnorris/precautionary/actions)
<!-- badges: end -->

`precautionary` implements new layers of patient-centered safety analysis
for phase 1 dose-escalation trials, adding diagnostics to examine
the safety characteristics of these designs in light of expected
inter-individual variation in pharmacokinetics and pharmacodynamics.
See Norris (2020b), "Retrospective analysis of a fatal dose-finding trial"
[arXiv:2004.12755](https://arxiv.org/abs/2004.12755) and (2020c)
"What Were They Thinking? Pharmacologic priors implicit in a choice of 3+3
dose-escalation design" [arXiv:2012.05301](https://arxiv.org/abs/2012.05301).

## Installation

Releases starting with 0.2.3 incorporate fast numerics implemented in
[Rust](https://www.rust-lang.org), a modern programming language that
emphasizes performance and reliability---attributes crucial to
applications such as the analysis of clinical trials.

Until CRAN maintainers have worked out suitable practices for deploying
R packages containing Rust code, the newest features of `precautionary`
will be available only here on GitHub.

``` r
# Install release version from GitHub
remotes::install_github("dcnorris/precautionary")

# Install obsolete version from CRAN (where review of new Rust library remains pending)
install.package("precautionary")
```

To date, those features of `precautionary` which depend on the Prolog code
in `exec/prolog/` have been pre-built into the package, for example as the
arrays `T[,,,]` written into `R/sysdata.rda` by `exec/make_sysdata_TUb.R`.
Methodologists who wish to examine, recompute and verify these arrays are
advised to install [Scryer Prolog](https://github.com/mthom/scryer-prolog).

It is a near-term goal for `precautionary` to reveal more transparently
Prolog's special contributions to its analysis of dose-escalation designs.

## Usage

Please see the vignettes under the [Articles](#) tab above.

## References

The `precautionary` package is the pointy end of the spear in a larger
[DTAT research programme](https::/precisionmethods.guru), of which the
following are key outputs. Several of these citations have accompanying
online resources such as web applications. For the key references,
lay explanations are available.

<ol>
<li>Norris DC. Dose Titration Algorithm Tuning (DTAT) should supersede &lsquo;the&rsquo; Maximum Tolerated Dose (MTD) in oncology dose-finding trials. <i>F1000Research.</i> 2017;6:112. doi:<a href="https://f1000research.com/articles/6-112/v3">10.12688/f1000research.10624.3</a>. [<a href="https://precisionmethods.guru/2019/04/16/a-new-concept-may-help-us-at-last-abandon-one-size-fits-all-dosing-of-cancer-treatment-drugs/">lay explanation</a>]
</li>
<li>
&ndash;&ndash;&ndash;&ndash;&ndash;. Dose Titration Algorithm Tuning (DTAT) should supplant &lsquo;the&rsquo; MTD. May 2017. [podium presentation] Symposium on Dose Selection for Cancer Treatment Drugs, Stanford Center for Innovative Study Design (CISD) May 12, 2017. doi:<a href="https://f1000research.com/slides/6-854">10.7490/f1000research.1114209.1</a>.
</li>
<li>
&ndash;&ndash;&ndash;&ndash;&ndash;. Costing &lsquo;the&rsquo; MTD. <i>bioRxiv.</i> August 2017:150821. doi:<a href="https://www.biorxiv.org/content/early/2017/08/22/150821">10.1101/150821</a>. [<a href="https://precisionmethods.guru/2019/04/16/one-size-fits-all-dosing-of-cancer-treatment-drugs-how-much-does-it-cost-society/">lay explanation</a>]
</li>
<li>
&ndash;&ndash;&ndash;&ndash;&ndash;. Costing &lsquo;the&rsquo; MTD: What Is the Economic and Human Cost of 1-Size-Fits-All Dose Finding in Oncology? [poster] Presented at 8th American Conference on Pharmacometrics (ACoP8), October 16, 2017. doi:<a href="https://f1000research.com/posters/6-1861">10.7490/f1000research.1114988.1</a>.
</li>
<li>
&ndash;&ndash;&ndash;&ndash;&ndash;. One-size-fits-all dosing in oncology wastes money, innovation and lives. <i>Drug Discovery Today.</i> 2018;23(1):4-6. doi:<a href="Norris (2018) One-size-fits-all dosing in oncology wastes money, innovation and lives.pdf">10.1016/j.drudis.2017.11.008</a>. [<a href="https://precision-methodologies.shinyapps.io/thecost/">Shiny app</a>]
</li>
<li>
&ndash;&ndash;&ndash;&ndash;&ndash;. Precautionary Coherence Unravels Dose Escalation Designs. <i>bioRxiv.</i> December 2017:240846. doi:<a href="https://www.biorxiv.org/content/early/2017/12/29/240846">10.1101/240846</a>. [<a href="https://precisionmethods.guru/2019/04/14/the-conduct-of-most-first-in-human-oncology-drug-trials-is-conceptually-incoherent-and-unethical/">lay explanation</a>] [<a href="../3+3/PC/">D3 app</a>]
</li>
<li>
&ndash;&ndash;&ndash;&ndash;&ndash;. Costing &lsquo;the&rsquo; MTD ... in 2-D. <i>bioRxiv.</i> July 2018:370817. doi:<a
href="https://www.biorxiv.org/content/early/2018/07/17/370817">10.1101/370817</a> [<a href="https://precisionmethods.guru/2019/04/16/clinicians-must-regain-control-over-phase-1-cancer-combination-therapy-trials/">lay explanation</a>]
</li>
<li>
&ndash;&ndash;&ndash;&ndash;&ndash;. Ethical Review and Methodologic Innovation in Phase 1 Cancer Trials. <i>JAMA Pediatrics.</i> April 2019. doi:<a href="https://dx.doi.org/10.1001/jamapediatrics.2019.0811">10.1001/jamapediatrics.2019.0811</a> [<a href="https://precisionmethods.guru/2019/04/25/precautionary-coherence-for-irbs/">2-minute video</a>]
</li>
<li>
&ndash;&ndash;&ndash;&ndash;&ndash;. Impeachment of One-Size-Fits-All Dosing for Obstruction of Synergism. [working paper]
December 4, 2019. doi:<a href="https://osf.io/3hcdb/">10.17605/OSF.IO/S7XDU</a>. [<a href="https://precisionmethods.guru/2020/01/13/therapeutic-synergism-and-the-statistician/">2-minute video</a>]
</li>
<li>
&ndash;&ndash;&ndash;&ndash;&ndash;. Comment on Wages et al., Coherence principles in interval-based dose finding. Pharmaceutical Statistics 2019, DOI: 10.1002/pst.1974. <i>Pharmaceutical Statistics.</i> March 2020. doi:<a href="https://onlinelibrary.wiley.com/doi/full/10.1002/pst.2016">10.1002/pst.2016</a>
[<a href="https://precisionmethods.guru/2019/12/02/comment-on-wages-et-al-coherence-principles-in-interval-based-dose-finding/">additional background</a>]
</li>
<li>
&ndash;&ndash;&ndash;&ndash;&ndash;. Retrospective analysis of a fatal dose-finding trial. <a href="https://arxiv.org/abs/2004.12755">arXiv:2004.12755 [stat, q-bio]</a>. April 2020. [<a href="https://threadreaderapp.com/thread/1255095770627428352.html">Tweetorial</a>]
</li>
<li>
Norris DC, Sen S, Groisberg R, Subbiah V. Patient-Centered, Physician-Investigator Friendly Pragmatic Phase I/II Trial Designs&mdash;The 4P Model. <i>Mayo Clinic Proceedings.</i> 2020;95(11):2566-2568. doi:<a href="https://www.mayoclinicproceedings.org/article/S0025-6196(20)31039-9/fulltext">10.1016/j.mayocp.2020.09.009</a>
</li>
<li>
Norris DC. What Were They Thinking? Pharmacologic priors implicit in a choice of 3+3 dose-escalation design. <a href="https://arxiv.org/abs/2012.05301">arXiv:2012.05301 [stat, q-bio]</a>. December 9, 2020. [<a href="https://threadreaderapp.com/thread/1339219770730799106.html">Tweetorial</a>]
</li>
</ol>
