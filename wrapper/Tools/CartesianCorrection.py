# Script to calculate the analytical correction for releasing Cartesian
# restraints - this includes the standard state correction.
# @author: Finlay Clark
# Based very closely on StandardState.py by Stefano Bosisio and Julien Michel

import numpy as np
import os
from math import pi, sin, log
from scipy.special import dawsn

from Sire.Tools import Parameter, resolveParameters
#from Sire.Units import *
from Sire.Tools.OpenMMMD import *

# Constants
v0 = 1660.53907 # A^3, the standard state volume
R = 0.00198720425864083 # kcal mol-1, the molar gas constant

#temperature = Parameter("temperature", 25 * celsius, """Simulation temperature""")

@resolveParameters
def run():
    try:
        host = os.environ['HOSTNAME']
    except KeyError:
        host = "unknown"
    print("### Running Standard state correction calculation on %s ###" % host)

    if verbose.val:
        print("###====================Utilised Parameters=====================###")
        print(temperature)
        print(cartesian_restraints_dict)
        print ("###===========================================================###\n")

    # Get Cartesian restraint dict in dict form
    cart_dict = dict(cartesian_restraints_dict.val)

    # Params
    T = temperature.val.value() # K

    #force_constants = list(boresch_dict["force_constants"].values()) # kcal mol-1 A-2 or rad-2
    #prod_force_constants = np.prod(force_constants)

    prefactor = 8*(pi**2)*v0 # Divide this to account for force constants of 0
    force_constants = []
    k_theta = 0
    k_theta = cart_dict["force_constants"]["k_theta"]

    # Loop through and correct for angle force constants of zero,
    # which break the analytical correction
    for k, val in cart_dict["force_constants"].items():
        if val == 0:
            if k[3] == "r":
                print("Error: Positional restraints must not be zero")
                sys.exit(-1)
            if k in ["k_phi", "k_psi"]:
                prefactor /= 2*pi
            if k == "k_theta":
                prefactor /= 2
        else:
            force_constants.append(val)

    n_nonzero_k = len(force_constants)
    prod_force_constants = np.prod(force_constants)

    # Calculation
    numerator = np.sqrt(prod_force_constants)
    denominator_1 = ((2*pi*R*T)**(n_nonzero_k/2))/pi**((bool(k_theta))/2) # If one of the force consts 
                                                                          # is k_theta, divide by sqrt(pi)
    if bool(k_theta):
        denominator_2 = dawsn(((R*T)**(0.5))/((2*k_theta)**0.5))
    else:
        denominator_2 = 1

    dg = -R*T*log((prefactor*numerator)/(denominator_1 * denominator_2))

    print(f"Analytical correction for releasing Cartesian restraints = {dg:.2f} kcal mol-1")