import numpy as np
import sys
import os


import dustpy
from dustpy import Simulation
from simframe.frame import Field
from dustpy import constants as c


################################
# BENCHMARK
################################
'''
There is a power-law steady-state solution for a disk with back-reaction under certain conditions.
This disk behaves as standard power-law viscous accretion, with a reduced equivalent alpha viscosity.
The description of the can be found in Garate et al.(2020) Appendix A.
'''

def get_Simulation(build_up = False, SigmaDust = None):
    sim = Simulation()

    ################################
    # DISK SETUP
    ################################
    # Star and disk parameters
    sim.ini.star.M = 1.0 * c.M_sun
    sim.ini.gas.gamma = 1.0


    # Relevant Dust parameters
    sim.ini.gas.alpha = 1.e-2   # # High alpha to make the disk gas evolution noticeable
    sim.ini.dust.d2gRatio = 0.25 # High dust-to-gas ratio
    sim.ini.dust.vfrag = 1000 # Extremely low fragmentation velocity

    ################################
    # GRID SETUP
    ################################

    # Mass Grid Parameters
    sim.ini.grid.Nmbpd = 7
    sim.ini.grid.mmin = 1.e-12
    sim.ini.grid.mmax = 1.e0

    # Radial Grid Parameters
    sim.ini.grid.Nr = 200
    sim.ini.grid.rmin = 5 * c.au
    sim.ini.grid.rmax = 500 * c.au

    # Make the grids in advace to set the gas surface density using the power law distribution
    sim.makegrids()


    ################################
    # INITIAL SURFACE DENSITY PROFILE
    ################################

    def Sigma_PowerLaw(r, Sigma0, R0, exp):
        return Sigma0 * (r / R0)**exp

    SigmaGas = np.array(Sigma_PowerLaw(sim.grid.r, 30., 10 * c.au, -1))
    sim.gas.Sigma = Field(sim, SigmaGas, description="Surface density [g/cm²]")

    if not build_up:
        sim.dust.Sigma = Field(sim, SigmaDust, description="Surface density per mass bin [g/cm²]")

    sim.initialize()




    ########################################################################
    # ADJUST Deltas
    ########################################################################
    # Settling and turbulence can match the gas viscosity, but there should not be any radial diffusivity
    sim.dust.delta.turb = sim.gas.alpha
    sim.dust.delta.vert = sim.gas.alpha
    sim.dust.delta.rad = 1.e-8
    sim.update()

    ########################################################################
    # ADJUST FRAGMENTATION VELOCITY
    ########################################################################

    # The fragmentation limit needs to be constant in order to reach the steady-state solution for dust backreaction
    # We set the Stokes number that we wish to obtain, and then compute the fragmentation velocity accordingly
    def fragmentationVelocity_ConstantStokes(sim):
        St_target = 2.5e-3
        v_frag = np.sqrt(3 * St_target * sim.dust.delta.turb) * sim.gas.cs
        return v_frag

    sim.dust.v.frag.updater = fragmentationVelocity_ConstantStokes
    sim.update()


    ################################
    # BACKREACTION SETUP
    ################################
    # Set up back-reaction (assuming a vertically uniform dust-to-gas ratio)
    sys.path.append("../")
    from setup_backreaction import setup_backreaction
    setup_backreaction(sim, vertical_setup = False, velocity_update = False)


    ################################################################
    # TURN-OFF GAS EVOLUTION AND DUST ADVECTION DURING BUILDUP
    ################################################################
    if build_up:
        # Set the dust advection to zero
        sim.dust.v.rad = 0
        sim.dust.D = 0.
        sim.dust.v.rad.updater = None
        sim.dust.D.updater = None

        # Delete gas evolution during build-up
        del(sim.integrator.instructions[1])
    else:
        # Adjust the outer boundary condition here...


        # Update the velocities so that they are written correctly in the first snapshot
        sim.update()
        sim.gas.v.update()
        sim.dust.v.update()
        sim.dust.v.rad.update()



    sim.update()
    return sim

################################
# RUN SIMULATION
################################

# Make the build-up simulation to get the dust distribution
sim_build_up = get_Simulation(build_up = True)
sim_build_up.writer.datadir = "./Simulation_BuildUp/"
sim_build_up.t.snapshots = np.linspace(0.5, 2.5, 5) * 1.e4 * c.year
sim_build_up.writer.overwrite = True
sim_build_up.run()


# Let the simulation evolve normally with gas and dust advection, starting from the fully grown dust distribution
sim = get_Simulation(build_up = False, SigmaDust = sim_build_up.dust.Sigma)
sim.writer.datadir = "./Simulation_PowerLaw/"
sim.t.snapshots = np.linspace(0.1, 1.0, 10) * 1.e5 * c.year
sim.writer.overwrite = True
sim.run()