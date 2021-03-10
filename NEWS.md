# precautionary 0.2-2

## Changes

* Corrected overloaded use of 'C' variable in DTP vignette

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
