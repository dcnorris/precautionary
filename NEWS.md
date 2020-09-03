# precautionary 0.1-3

TODO: Ping Piotr Juszczak <piotr.juszczak@roche.com> upon release.

# precautionary 0.1-2

## Changes

* Simulations are now extensible, via extend.(precautionary|hyper) methods
* No longer using or depending on distr6, eliminating object creation overhead
* Hyper-simulations no longer have a superfluous 'K' (within-hyperprior-draw) dimension
* Sim summary()'s $safety component has prepended class 'safetytab'
* The format.safetytab() method shows significant digits in accordance with MCSEs

# precautionary 0.1-1

* First public release
