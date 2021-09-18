#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. py:currentmodule:: tests.test_element_properties
.. moduleauthor:: Hendrix Demers <hendrix.demers@mail.mcgill.ca>

Tests for the :py:mod:`eecs.element_properties` module.
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
from pytest import approx

# Local modules.

# Project modules.
from eecs.element_properties import *

# Globals and constants variables.


def test_is_discovered():
    """
    Test used to validate the file is included in the tests
    by the test framework.
    """
    # assert False
    assert True


def test_get_mass_density_g_cm3():
    assert 7.1900 == getMassDensity_g_cm3(24)


def test_get_atomic_mass_g_mol():
    assert 51.996000 == getAtomicMass_g_mol(24)


def test_get_fermi_energy_eV():
    assert 1.000 == getFermiEnergy_eV(24)

    assert 4.700 == getFermiEnergy_eV(3)


def test_get_k_fermi_eV():
    assert 7.00E7 == getKFermi_eV(24)

    assert 1.10E8 == getKFermi_eV(3)


def test_get_plasmon_energy_eV():
    assert 24.9 == getPlasmonEnergy_eV(24)


def test_compute_atomic_density_atom_cm3():
    mass_density_g_cm3 = getMassDensity_g_cm3(13)
    atomic_mass_g_mol = getAtomicMass_g_mol(13)

    value = computeAtomicDensity_atom_cm3(mass_density_g_cm3, atomic_mass_g_mol)
    value *= 1.0E-22

    assert 6.02617011482666 == approx(value, 4)


def test_get_k_ratio_correction():
    assert 0.80707 == approx(getKRatioCorrection(13), 4)


def test_get_mean_ionization_energy_eV():
    assert 149.5 == getMeanIonizationEnergy_eV(13)


def test_get_symbol():
    assert 'Al' == getSymbol(13)


def test_get_name():
    assert 'Aluminium' == getName(13)
