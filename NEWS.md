# precautionary 0.2.6

## Changes

* EscRisk app displays CPE progress via 'odometer'
* EscRisk app lets user set max enrollment in range 20:30
* Restrict to OS_type: unix, as feasible CPE demands fork-based parallelism

# precautionary 0.2.5

## Changes

* Overhauled EscRisk shiny app to use CPE

# precautionary 0.2.4

## Changes

* Refined version numbering scheme to 0.<tranche>.<feature>-<patch>
* Implemented CPE as R6 superclass 'Cpe'
* Implemented Cumulative-Cohort Designs (CCDs) in 'Cpe' subclass 'Ccd'
* Implemented BOIN design as 'Boin' subclass of 'Ccd'
* Executable specification for CCD and BOIN in exec/prolog/
* Demoted package dtpcrm from a 'Depends' to a 'Suggests'

# precautionary 0.2-3

## Changes

* Parallelized Crm$trace_paths method
* Added 'MCSE-Free' calibration vignette

# precautionary 0.2-2

## Changes

* Implements R6 class to wrap CRM functions, with caching to speed DTP
* Implements faster versions of certain objective functions from package 'dtpcrm'
* Explicitly integrate() over c(-Inf,Inf) in 'crm' function (as per documentation)
* Implements faster version of dtpcrm::stop_for_excess_toxicity_empiric
* Corrected overloaded use of 'C' variable in text of DTP vignette

# precautionary 0.2-1

## Changes

* Added vignette 'Generalized dose-escalation safety schematics via DTP'
* Corrected a few typos

# precautionary 0.2

## Changes

* Added a Prolog program to exhaustively enumerate 3+3 trial outcomes

# precautionary 0.1-5

## Changes

* Accommodate new 'true_prob_eff' arg in escalation::simulate_trials

# precautionary 0.1-4

## Changes

* Added vignette for a regulatory application

# precautionary 0.1-3

## Changes

* Added the EscRisk shiny app, runnable by demo(EscRisk)

# precautionary 0.1-2

## Changes

* Simulations are now extensible, via extend.(precautionary|hyper) methods
* No longer using or depending on distr6, eliminating object creation overhead
* Hyper-simulations no longer have a superfluous 'K' (within-hyperprior-draw) dimension
* Sim summary()'s $safety component has prepended class 'safetytab'
* The format.safetytab() method shows significant digits in accordance with MCSEs

# precautionary 0.1-1

* First public release
