#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. py:currentmodule:: tests.test_numeric_conversion
.. moduleauthor:: Hendrix Demers <hendrix.demers@mail.mcgill.ca>

Tests for the :py:mod:`eecs.numeric_conversion` module.
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
from eecs.numeric_conversion import *

# Globals and constants variables.


def test_is_discovered():
    """
    Test used to validate the file is included in the tests
    by the test framework.
    """
    # assert False
    assert True


def test_cm2_to_nm2():
    value_cm2 = 8.90454E-16
    value_ref_nm2 = 0.0890454
    value_nm2 = cm2_to_nm2(value_cm2)
    assert value_ref_nm2 == approx(value_nm2)


def test_cm2_to_um2():
    value_ref_cm2 = 1.0
    value_ref_um2 = 100000000

    value_um2 = cm2_to_um2(value_ref_cm2)
    assert value_ref_um2 == approx(value_um2)

    value_cm2 = um2_to_cm2(value_ref_um2)
    assert value_ref_cm2 == approx(value_cm2)


def test_cm3_to_um3():
    value_ref_cm3 = 1.0
    value_ref_um3 = 1.0e12

    value_um3 = cm3_to_um3(value_ref_cm3)
    assert value_ref_um3 == approx(value_um3)

    value_cm3 = um3_to_cm3(value_ref_um3)
    assert value_ref_cm3 == approx(value_cm3)


def test_1_cm3_to_1_um3():
    value_ref_cm3 = 1.0/1.0
    value_ref_um3 = 1.0/1.0e12

    value_um3 = 1.0/cm3_to_um3(1.0/value_ref_cm3)
    assert value_ref_um3 == approx(value_um3)

    value_cm3 = 1.0/um3_to_cm3(1.0/value_ref_um3)
    assert value_ref_cm3 == approx(value_cm3)


def test_MeV_to_eV():
    value_ref_MeV = 2.3
    value_ref_eV = 2.3e6

    value_eV = MeV_to_eV(value_ref_MeV)
    assert value_ref_eV == approx(value_eV)

    value_MeV = eV_to_MeV(value_ref_eV)
    assert value_ref_MeV == approx(value_MeV)


def test_keV_to_eV():
    value_ref_keV = 2.3
    value_ref_eV = 2.3e3

    value_eV = keV_to_eV(value_ref_keV)
    assert value_ref_eV == approx(value_eV)

    value_MeV = eV_to_keV(value_ref_eV)
    assert value_ref_keV == approx(value_MeV)


def test_barn_to_cm2():
    value_ref_barn = 2.34
    value_ref_cm2 = 2.34e-24

    value_cm2 = barn_to_cm2(value_ref_barn)
    assert value_ref_cm2 == approx(value_cm2)

    value_barn = cm2_to_barn(value_ref_cm2)
    assert value_ref_barn == approx(value_barn)


def test_barn_to_nm2():
    value_ref_barn = 2.34
    value_ref_nm2 = 2.34e-10

    value_nm2 = barn_to_nm2(value_ref_barn)
    assert value_ref_nm2 == approx(value_nm2)

    value_barn = nm2_to_barn(value_ref_nm2)
    assert value_ref_barn == approx(value_barn)


def test_m_to_fm():
    value_ref_m = 12.5
    value_ref_fm = 12.5e15

    value_cm = m_to_fm(value_ref_m)
    assert value_ref_fm == approx(value_cm)

    value_m = fm_to_m(value_ref_fm)
    assert value_ref_m == approx(value_m)


def test_m_to_cm():
    value_ref_m = 10.0
    value_ref_cm = 1000.0

    value_cm = m_to_cm(value_ref_m)
    assert value_ref_cm == approx(value_cm)

    value_m = cm_to_m(value_ref_cm)
    assert value_ref_m == approx(value_m)


def test_m_to_mm():
    value_ref_m = 10.0
    value_ref_mm = 10000.0

    value_mm = m_to_mm(value_ref_m)
    assert value_ref_mm == approx(value_mm)

    value_m = mm_to_m(value_ref_mm)
    assert value_ref_m == approx(value_m)


def test_m_to_nm():
    value_ref_m = 12.5
    value_ref_nm = 12.5e9

    value_cm = m_to_nm(value_ref_m)
    assert value_ref_nm == approx(value_cm)

    value_m = nm_to_m(value_ref_nm)
    assert value_ref_m == approx(value_m)


def test_m_to_um():
    value_ref_m = 12.5
    value_ref_um = 12.5e6

    value_um = m_to_um(value_ref_m)
    assert value_ref_um == approx(value_um)

    value_m = um_to_m(value_ref_um)
    assert value_ref_m == approx(value_m)


def test_m_to_pm():
    value_ref_m = 12.5
    value_ref_pm = 12.5e12

    value_cm = m_to_pm(value_ref_m)
    assert value_ref_pm == approx(value_cm)

    value_m = pm_to_m(value_ref_pm)
    assert value_ref_m == approx(value_m)


def test_cm_to_mm():
    value_ref_cm = 12.5
    value_ref_mm = 125.0

    value_mm = cm_to_mm(value_ref_cm)
    assert value_ref_mm == approx(value_mm)

    value_cm = mm_to_cm(value_ref_mm)
    assert value_ref_cm == approx(value_cm)


def test_cm_to_um():
    value_ref_cm = 12.5
    value_ref_um = 12.5e4

    value_um = cm_to_um(value_ref_cm)
    assert value_ref_um == approx(value_um)

    value_cm = um_to_cm(value_ref_um)
    assert value_ref_cm == approx(value_cm)


def test_nm_to_um():
    value_ref_nm = 12.5e3
    value_ref_um = 12.5

    value_um = nm_to_um(value_ref_nm)
    assert value_ref_um == approx(value_um)

    value_nm = um_to_nm(value_ref_um)
    assert value_ref_nm == approx(value_nm)


def test_nm_to_cm():
    value_ref_nm = 12.5e7
    value_ref_cm = 12.5

    value_um = nm_to_cm(value_ref_nm)
    assert value_ref_cm == approx(value_um)

    value_nm = cm_to_nm(value_ref_cm)
    assert value_ref_nm == approx(value_nm)


def test_rad_to_mrad():
    value_ref_rad = 0.1
    value_ref_mrad = 100.0

    value_mrad = rad_to_mrad(value_ref_rad)
    assert value_ref_mrad == approx(value_mrad)

    value_rad = mrad_to_rad(value_ref_mrad)
    assert value_ref_rad == approx(value_rad)
