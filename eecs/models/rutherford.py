#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. py:currentmodule:: eecs.models.rutherford
.. moduleauthor:: Hendrix Demers <hendrix.demers@mail.mcgill.ca>

Rutherford electron elastic cross section models.
"""

###############################################################################
# Copyright 2021 Hendrix Demers
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
###############################################################################

# Standard library modules.

# Third party modules.
import numpy as np
from scipy.constants import e, epsilon_0

# Local modules.

# Project modules.

# Globals and constants variables.


def partial_williams_carter(atomic_number, energy_eV, theta_rad):
    """
    """
    Z = atomic_number
    E0 = energy_eV
    e0 = epsilon_0

    factor = e*e*e*e*Z*Z/(16.0*(4.0*np.pi*e0*E0)**2)
    denominator = np.power(np.sin(theta_rad/2.0), 4.0)

    dtheta_domega_m2_sr = factor / denominator

    return dtheta_domega_m2_sr


def partial_screened_williams_carter(atomic_number, energy_eV, theta_rad):
    """
    """
    Z = atomic_number
    E0 = energy_eV

    theta0 = 0.117*np.power(Z, 1.0/3.0)/np.sqrt(E0)
    lambda_relativistic = 0.0
    a0 = 0.0

    factor = Z*Z*np.power(lambda_relativistic, 4.0)/(64.0*np.power(np.pi, 4.0)*a0*a0)
    denominator = np.power(np.power(np.sin(theta_rad/2.0), 2.0) + theta0*theta0/4.0, 2.0)

    dtheta_domega_m2_sr = factor/denominator

    return dtheta_domega_m2_sr


def total_williams_carter(atomic_number, energy_eV, theta_rad):
    """
    """
    Z = atomic_number
    E0 = energy_eV

    sigma = 1.62e-24 * np.power(Z/E0, 2.0) / np.power(np.tan(theta_rad/2.0), 2.0)

    return sigma
