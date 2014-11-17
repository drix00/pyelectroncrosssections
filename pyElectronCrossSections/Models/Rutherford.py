#!/usr/bin/env python
"""
.. py:currentmodule:: Models.Rutherford
.. moduleauthor:: Hendrix Demers <hendrix.demers@mail.mcgill.ca>

Rutherford electron elastic cross section models.
"""

# Script information for the file.
__author__ = "Hendrix Demers (hendrix.demers@mail.mcgill.ca)"
__version__ = ""
__date__ = ""
__copyright__ = "Copyright (c) 2012 Hendrix Demers"
__license__ = ""

# Standard library modules.

# Third party modules.
import numpy as np
from scipy.constants import e, epsilon_0
# Local modules.

# Project modules

# Globals and constants variables.

def partialWilliamsCarter(atomicNumber, energy_eV, theta_rad):
    """
    """
    Z = atomicNumber
    E0 = energy_eV
    e0 = epsilon_0

    factor = e*e*e*e*Z*Z/(16.0*(4.0*np.pi*e0*E0)**2)
    denominator = np.power(np.sin(theta_rad/2.0), 4.0)

    dThetadOmega_m2_sr = factor/denominator

    return dThetadOmega_m2_sr

def partialScreenedWilliamsCarter(atomicNumber, energy_eV, theta_rad):
    """
    """
    Z = atomicNumber
    E0 = energy_eV

    theta0 = 0.117*np.power(Z, 1.0/3.0)/np.sqrt(E0)
    lambdaRelativistic = 0.0
    a0 = 0.0

    factor = Z*Z*np.power(lambdaRelativistic, 4.0)/(64.0*np.power(np.pi, 4.0)*a0*a0)
    denominator = np.power(np.power(np.sin(theta_rad/2.0), 2.0) + theta0*theta0/4.0, 2.0)

    dThetadOmega_m2_sr = factor/denominator

    return dThetadOmega_m2_sr

def totalWilliamsCarter(atomicNumber, energy_eV, theta_rad):
    """
    """
    Z = atomicNumber
    E0 = energy_eV

    sigma = 1.62e-24 * np.power(Z/E0, 2.0) /np.power(np.tan(theta_rad/2.0), 2.0)

    return sigma

if __name__ == '__main__': #pragma: no cover
    import pyHendrixDemersTools.Runner as Runner
    Runner.Runner().run(runFunction=None)
