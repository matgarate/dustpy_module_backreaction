import dustpy
import numpy as np

from functions_backreaction import BackreactionCoefficients, BackreactionCoefficients_VerticalStructure
from functions_backreaction import Backreaction_A, Backreaction_B, Backreaction_A_VerticalStructure, Backreaction_B_VerticalStructure
from functions_backreaction import vrad_dust_BackreactionVerticalStructure
from functions_backreaction import dustDiffusivity_Backreaction

################################
# Helper routine to add backreaction to your Simulation object in one line.
################################
def setup_backreaction(sim, vertical_setup = False, velocity_update = False):
    '''
    Add the backreaction setup to your simulation object.
    Call the backreaction setup function after the initialization and then run, as follows:

    sim.initialize()
    setup_backreaction(sim)
    sim.run()
    '''


    if vertical_setup:
        # Set the backreaction coefficients
        sim.dust.backreaction.addfield("AB", np.ones((2 * (sim.grid.Nm[0] + 1), sim.grid.Nr[0])),  description = "Backreaction Coefficients (joint - internal)")

        # Additional back-reaction coefficients for the dust
        sim.dust.backreaction.addfield("A_vertical", np.ones_like(sim.dust.a) ,  description = "Backreaction Coefficient A, considering dust vertical settling")
        sim.dust.backreaction.addfield("B_vertical", np.zeros_like(sim.dust.a),  description = "Backreaction Coefficient B, considering dust vertical settling")


        sim.dust.backreaction.updater = ["AB","A", "B", "A_vertical", "B_vertical"]
        sim.dust.backreaction.AB.updater = BackreactionCoefficients_VerticalStructure

        # The standard backreaction coefficients are used for the gas dynamics
        sim.dust.backreaction.A.updater = Backreaction_A
        sim.dust.backreaction.B.updater = Backreaction_B
        # The backreaction coefficients A_vertical and B_vertical are used for the dust dynamics
        sim.dust.backreaction.A_vertical.updater = Backreaction_A_VerticalStructure
        sim.dust.backreaction.B_vertical.updater = Backreaction_B_VerticalStructure

        # Redefine the radial dust velocity to consider one pair of backreaction coefficients per dust species
        sim.dust.v.rad.updater = vrad_dust_BackreactionVerticalStructure

    else:
        # Set the backreaction coefficients
        sim.dust.backreaction.addfield("AB", np.array([np.ones_like(sim.grid.r), np.zeros_like(sim.grid.r)]),  description = "Backreaction Coefficients (joint - internal)")

        sim.dust.backreaction.updater = ["AB","A", "B"]
        sim.dust.backreaction.AB.updater = BackreactionCoefficients
        sim.dust.backreaction.A.updater = Backreaction_A
        sim.dust.backreaction.B.updater = Backreaction_B



    # Update the dust diffusivity to account for high dust-to-gas ratios
    sim.dust.D.updater = dustDiffusivity_Backreaction


    # Update the gas and dust velocities in dust.v.diastole
    # Not necessary for implicit simulations?
    if velocity_update:
        from functions_backreaction import update_RadialVelocities
        sim.dust.v.diastole = update_RadialVelocities

    # Update all
    sim.update()
