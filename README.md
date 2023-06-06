# Backreaction Module for Dustpy

Includes the dynamical dust backreaction into [DustPy](https://github.com/stammler/dustpy) (Stammler and Birnstiel, in prep.), following the Garate et al.([2019](https://ui.adsabs.harvard.edu/abs/2019ApJ...871...53G%2F/abstract), [2020](https://ui.adsabs.harvard.edu/abs/2020A%26A...635A.149G/abstract)) implementation.


To setup the backreaction module, add the following lines after initialization.

`from setup_backreaction import setup_backreaction`

`setup_backreaction(sim, vertical_setup = False)`



See the `run_backreaction.py` code for an example.
If you use this module, please cite [Garate et al.(2020)](https://ui.adsabs.harvard.edu/abs/2020A%26A...635A.149G/abstract)


This module is currently being implemented into the [Dustpylib](https://github.com/stammler/dustpylib) under the `dustpylib.dynamical.backreaction` module, with additional documentation and examples.
