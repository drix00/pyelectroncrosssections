#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. py:currentmodule:: eecs.numeric_conversion
.. moduleauthor:: Hendrix Demers <hendrix.demers@mail.mcgill.ca>

Numeric conversion functions.
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

# Local modules.

# Project modules.

# Globals and constants variables.


def cm2_to_nm2(value_cm2):
    value_nm2 = value_cm2 * 1.0e-4 * 1.0e18
    return value_nm2


def cm2_to_um2(value_cm2):
    value_um2 = value_cm2 * 1.0e8
    return value_um2


def um2_to_cm2(value_um2):
    value_cm2 = value_um2 * 1.0e-8
    return value_cm2


def cm3_to_um3(value_cm3):
    value_um3 = value_cm3 * 1.0e12
    return value_um3


def um3_to_cm3(value_um3):
    value_cm3 = value_um3 * 1.0e-12
    return value_cm3


def MeV_to_eV(value_MeV):
    value_eV = value_MeV * 1.0e6
    return value_eV


def eV_to_MeV(value_eV):
    value_MeV = value_eV * 1.0e-6
    return value_MeV


def keV_to_eV(value_MeV):
    value_eV = value_MeV * 1.0e3
    return value_eV


def eV_to_keV(value_eV):
    value_MeV = value_eV * 1.0e-3
    return value_MeV


def barn_to_cm2(value_barn):
    value_cm2 = value_barn * 1.0e-24
    return value_cm2


def cm2_to_barn(value_cm2):
    value_barn = value_cm2 * 1.0e24
    return value_barn


def barn_to_nm2(value_barn):
    value_nm2 = value_barn * 1.0e-10
    return value_nm2


def nm2_to_barn(value_nm2):
    value_barn = value_nm2 * 1.0e10
    return value_barn


def m_to_fm(value_m):
    value_fm = value_m * 1.0e15
    return value_fm


def fm_to_m(value_fm):
    value_m = value_fm * 1.0e-15
    return value_m


def m_to_cm(value_m):
    value_cm = value_m * 1.0e2
    return value_cm


def cm_to_m(value_cm):
    value_m = value_cm * 1.0e-2
    return value_m


def m_to_mm(value_m):
    value_mm = value_m * 1.0e3
    return value_mm


def mm_to_m(value_mm):
    value_m = value_mm * 1.0e-3
    return value_m


def m_to_pm(value_m):
    value_pm = value_m * 1.0e12
    return value_pm


def pm_to_m(value_pm):
    value_m = value_pm * 1.0e-12
    return value_m


def m_to_nm(value_m):
    value_nm = value_m * 1.0e9
    return value_nm


def nm_to_m(value_nm):
    value_m = value_nm * 1.0e-9
    return value_m


def m_to_um(value_m):
    value_um = value_m * 1.0e6
    return value_um


def um_to_m(value_um):
    value_m = value_um * 1.0e-6
    return value_m


def cm_to_mm(value_cm):
    value_m = cm_to_m(value_cm)
    value_mm = m_to_mm(value_m)
    return value_mm


def mm_to_cm(value_mm):
    value_m = mm_to_m(value_mm)
    value_cm = m_to_cm(value_m)
    return value_cm


def cm_to_um(value_cm):
    value_m = cm_to_m(value_cm)
    value_um = m_to_um(value_m)
    return value_um


def um_to_cm(value_um):
    value_m = um_to_m(value_um)
    value_cm = m_to_cm(value_m)
    return value_cm


def nm_to_um(value_nm):
    value_m = nm_to_m(value_nm)
    value_um = m_to_um(value_m)
    return value_um


def um_to_nm(value_um):
    value_m = um_to_m(value_um)
    value_nm = m_to_nm(value_m)
    return value_nm


def nm_to_cm(value_nm):
    value_m = nm_to_m(value_nm)
    value_cm = m_to_cm(value_m)
    return value_cm


def cm_to_nm(value_cm):
    value_m = cm_to_m(value_cm)
    value_nm = m_to_nm(value_m)
    return value_nm


def rad_to_mrad(value_rad):
    value_mrad = value_rad * 1.0e3
    return value_mrad


def mrad_to_rad(value_mrad):
    value_rad = value_mrad * 1.0e-3
    return value_rad
