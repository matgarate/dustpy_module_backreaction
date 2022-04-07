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
# ADVANCED BACKREACTION SETUP
# (this setup considers the vertical structure)
# (each dust species have their independent pair of backreaction coefficients)
################################

from functions_backreaction import BackReactionCoefficients_VerticalStructure
from functions_backreaction import Backreaction_A, Backreaction_B, Backreaction_A_VerticalStructure, Backreaction_B_VerticalStructure
from functions_backreaction import vrad_dust_BackreactionVerticalStructure
from functions_backreaction import dustDiffusivity_Backreaction


# Set the backreaction coefficients
sim.dust.backreaction.addfield("AB", np.ones((2 * (sim.grid.Nm[0] + 1), sim.grid.Nr[0])),  description = "Backreaction Coefficients (joint - internal)")

# Additional back-reaction coefficients for the dust
sim.dust.backreaction.addfield("A_vertical", np.ones_like(sim.dust.a) ,  description = "Backreaction Coefficient A, considering dust vertical settling")
sim.dust.backreaction.addfield("B_vertical", np.zeros_like(sim.dust.a),  description = "Backreaction Coefficient B, considering dust vertical settling")


sim.dust.backreaction.updater = ["AB","A", "B", "A_vertical", "B_vertical"]
sim.dust.backreaction.AB.updater = BackReactionCoefficients_VerticalStructure

# The standard backreaction coefficients are used for the gas dynamics
sim.dust.backreaction.A.updater = Backreaction_A
sim.dust.backreaction.B.updater = Backreaction_B
# The backreaction coefficients A_vertical and B_vertical are used for the dust dynamics
sim.dust.backreaction.A_vertical.updater = Backreaction_A_VerticalStructure
sim.dust.backreaction.B_vertical.updater = Backreaction_B_VerticalStructure


# Redefine the radial dust velocity to consider one pair of backreaction coefficients per dust species
sim.dust.v.rad.updater = vrad_dust_BackreactionVerticalStructure

# Update the dust diffusivity to account for high dust-to-gas ratios
sim.dust.D.updater = dustDiffusivity_Backreaction

# Update all
sim.update()



################################
# UPDATE VELOCITIES (?)
################################
'''
# In the implicit integration scheme the dust.v.rad and gas.v are updated in the finalizer,
# But I think that in the explicit integration scheme these need to be updated simultaneously in the dust.v.disastole
from functions_backreaction import update_RadialVelocities
sim.dust.v.diastole = update_RadialVelocities
sim.update()
'''


################################
# RUN SIMULATION
################################
print("Running Simulation")

outputDir = "./Simulation_VerticalStructure/"
sim.writer.datadir = outputDir
sim.t.snapshots = np.linspace(0.5, 5.0, 11) * 1.e5 * c.year
sim.writer.overwrite = True
sim.verbosity = 2

sim.run()
