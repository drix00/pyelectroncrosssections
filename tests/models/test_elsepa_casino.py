#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. py:currentmodule:: tests.models.test_elsepa_casino
.. moduleauthor:: Hendrix Demers <hendrix.demers@mail.mcgill.ca>

Tests for the :py:mod:`eecs.models.elsepa_casino` module.
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
import os.path

# Third party modules.
import pytest
from pytest import approx

# Local modules.

# Project modules.
import eecs.models.elsepa_casino
from eecs import get_current_module_path

# Globals and constants variables.


def test_is_discovered():
    """
    Test used to validate the file is included in the tests
    by the test framework.
    """
    # assert False
    assert True


@pytest.fixture
def zip_file_path():
    file_path = "../../test_data/ELSEPA_AbsorptionCorrection_LinearInterpolationTabulation_0.1.zip"
    file_path = get_current_module_path(__file__, file_path)
    if not os.path.isfile(file_path):
        pytest.skip("No file: {}".format(file_path))
    return file_path


@pytest.fixture
def casino_cross_section(zip_file_path):
    cross_section = eecs.models.elsepa_casino.ElsepaCasino(zip_file_path)
    return cross_section


def test_zip_filepath():
    """
    First test to check if the testcase is working with the testing framework.
    """

    file_name = "../../test_data/ELSEPA_AbsorptionCorrection_LinearInterpolationTabulation_0.1.zip"
    zip_filepath = get_current_module_path(__file__, file_name)
    cross_section = eecs.models.elsepa_casino.ElsepaCasino(zip_filepath)
    assert zip_filepath == cross_section.zip_filepath

    assert os.path.isfile(cross_section.zip_filepath)


def test_total_nm2(casino_cross_section):
    """
    Tests for method :py:meth:`total_nm2`.
    """

    total_ref_nm2 = 0.00120146

    total_nm2 = casino_cross_section.total_nm2(6, 5.0e3)

    assert total_ref_nm2 == approx(total_nm2)


def test_angle_deg(casino_cross_section):
    """
    Tests for method :py:meth:`angle_deg`.
    """

    angle_ref_deg = 4.345130434680232

    angle_deg = casino_cross_section.angle_deg(6, 5.0e3, 0.5)

    assert angle_ref_deg == approx(angle_deg)


def test_partial_nm2_sr(casino_cross_section):
    """
    Tests for method :py:meth:`partial_nm2_sr`.
    """

    partial_ref_nm2_sr = 0.017259500000000004

    partial_nm2_sr = casino_cross_section.partial_nm2_sr(6, 5.0e3, 4.35)

    assert partial_ref_nm2_sr == approx(partial_nm2_sr)


def test_read_total_data(casino_cross_section):
    """
    Tests for method :py:meth:`read_total_data`.
    """

    atomic_number = 6
    assert atomic_number not in casino_cross_section.total_functions
    casino_cross_section.read_total_data(atomic_number)
    assert atomic_number in casino_cross_section.total_functions


def test_read_angle_data(casino_cross_section):
    """
    Tests for method :py:meth:`read_angle_data`.
    """

    atomic_number = 6
    assert atomic_number not in casino_cross_section.angle_functions
    casino_cross_section.read_angle_data(atomic_number)
    assert atomic_number in casino_cross_section.angle_functions


def test_read_partial_data(casino_cross_section):
    """
    Tests for method :py:meth:`read_partial_data`.
    """

    atomic_number = 6
    assert atomic_number not in casino_cross_section.partial_functions
    casino_cross_section.read_partial_data(atomic_number)
    assert atomic_number in casino_cross_section.partial_functions
