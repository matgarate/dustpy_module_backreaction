import numpy as np
import sys
import os


import dustpy
from dustpy import Simulation
from dustpy import constants as c


sim = Simulation()

################################
# DISK SETUP
################################

# Relevant Dust parameters for backreaction experiments
sim.ini.gas.alpha = 1.e-3
sim.ini.dust.d2gRatio = 0.01
sim.ini.dust.vfrag = 1000.0

################################
# GRID SETUP
################################

# Radial Grid Parameters
sim.ini.grid.Nr = 200
sim.ini.grid.rmin = 5 * c.au
sim.ini.grid.rmax = 500 * c.au

sim.initialize()

################################
# BACKREACTION SETUP
################################

# Add the next two lines to your script after "initialize()"" to setup the backreaction updaters
# To account for the effect of vertical dust settling use "vertical_setup = True"
from setup_backreaction import setup_backreaction
setup_backreaction(sim, vertical_setup = False)


################################
# RUN SIMULATION
################################
sim.writer.datadir = "./Simulation/"
sim.t.snapshots = np.linspace(0.5, 5.0, 10) * 1.e5 * c.year

sim.writer.dumping = True
sim.writer.overwrite = True
sim.verbosity = 2

sim.run()
