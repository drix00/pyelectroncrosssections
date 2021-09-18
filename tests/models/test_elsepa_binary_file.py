#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. py:currentmodule:: tests.models.test_elsepa_binary_file
.. moduleauthor:: Hendrix Demers <hendrix.demers@mail.mcgill.ca>

Tests for the :py:mod:`eecs.models.elsepa_binary_file` module.
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

# Globals and constants variables.


def test_is_discovered():
    """
    Test used to validate the file is included in the tests
    by the test framework.
    """
    # assert False
    assert True


def test_init(el29_file_path, el29_file):
    assert el29_file_path == el29_file._filepath


def test_read_file(el29_file_path, el29_file):
    file = open(el29_file_path, 'rb')
    el29_file.read_file(file)

    assert 48 == len(el29_file._els_cs_info_list)

    els_cs_info = el29_file._els_cs_info_list[0]
    assert 30105020 == els_cs_info._file_version
    assert 29 == els_cs_info._atomic_number
    assert 1.000000e-01 == approx(els_cs_info._energy_keV)
    assert 2.870967e-02 == approx(els_cs_info._total_cs_nm2)

    assert 606 == els_cs_info._number_points

    assert 0.0 == approx(els_cs_info._ratios[0])
    assert 0.0 == approx(els_cs_info._theta_rad[0])

    assert 1.0 == approx(els_cs_info._ratios[-1])
    assert 3.1415926539 == approx(els_cs_info._theta_rad[-1])

    els_cs_info = el29_file._els_cs_info_list[-1]
    assert 30105020 == els_cs_info._file_version
    assert 29 == els_cs_info._atomic_number
    assert 500.0 == approx(els_cs_info._energy_keV)
    assert 0.000247381469262322 == approx(els_cs_info._total_cs_nm2)

    assert 606 == els_cs_info._number_points

    assert 0.0 == approx(els_cs_info._ratios[0])
    assert 0.0 == approx(els_cs_info._theta_rad[0])

    assert 1.0 == approx(els_cs_info._ratios[-1])
    assert 3.1415926539 == approx(els_cs_info._theta_rad[-1])

    assert 0.999999999788251 == approx(els_cs_info._ratios[-1])
    assert 3.13286600763917 == approx(els_cs_info._theta_rad[-2])
