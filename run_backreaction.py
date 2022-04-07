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
# Star and disk parameters
sim.ini.star.M = 1.0 * c.M_sun
sim.ini.gas.Mdisk = 0.05 * sim.ini.star.M
sim.ini.gas.SigmaRc = 60 * c.au
sim.ini.gas.gamma = 1.0


# Relevant Dust parameters
sim.ini.gas.alpha = 1.e-3
sim.ini.dust.d2gRatio = 0.01
sim.ini.dust.vfrag = 1000.0

################################
# GRID SETUP
################################

# Mass Grid Parameters
sim.ini.grid.Nmbpd = 7
sim.ini.grid.mmin = 1.e-12
sim.ini.grid.mmax = 1.e5


# Radial Grid Parameters
sim.ini.grid.Nr = 250
sim.ini.grid.rmin = 5 * c.au
sim.ini.grid.rmax = 500 * c.au

sim.initialize()


################################
# BACKREACTION SETUP
################################

from functions_backreaction import BackreactionCoefficients
from functions_backreaction import Backreaction_A, Backreaction_B
from functions_backreaction import dustDiffusivity_Backreaction


# Set the backreaction coefficients
sim.dust.backreaction.addfield("AB", np.array([np.ones_like(sim.grid.r), np.zeros_like(sim.grid.r)]),  description = "Backreaction Coefficients (joint - internal)")

sim.dust.backreaction.updater = ["AB","A", "B"]
sim.dust.backreaction.AB.updater = BackreactionCoefficients
sim.dust.backreaction.A.updater = Backreaction_A
sim.dust.backreaction.B.updater = Backreaction_B

# Update the dust diffusivity to account for high dust-to-gas ratios
sim.dust.D.updater = dustDiffusivity_Backreaction

# Update all
sim.update()



################################
# UPDATE VELOCITIES (?)
################################
'''
from functions_backreaction import update_RadialVelocities
# In the implicit integration scheme the dust.v.rad and gas.v are updated in the finalizer,
# But I think that in the explicit integration scheme these need to be updated simultaneously in the dust.v.disastole
sim.dust.v.diastole = update_RadialVelocities
'''


################################
# RUN SIMULATION
################################
print("Running Simulation")

outputDir = "./Simulation/"
sim.writer.datadir = outputDir
sim.t.snapshots = np.linspace(0.5, 5.0, 11) * 1.e5 * c.year
sim.writer.overwrite = True
sim.verbosity = 2

sim.run()
