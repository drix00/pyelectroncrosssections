#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. py:currentmodule:: tests.models.test_elsepa_cross_section_info
.. moduleauthor:: Hendrix Demers <hendrix.demers@mail.mcgill.ca>

Tests for the :py:mod:`eecs.models.elsepa_cross_section_info` module.
"""

###############################################################################
# Copyright 2019 Hendrix Demers
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
from eecs.models.elsepa_cross_section_info import ElsepaCrossSectionInfo

# Globals and constants variables.


def test_is_discovered():
    """
    Test used to validate the file is included in the tests
    by the test framework.
    """
    # assert False
    assert True


def test_read_file(el29_file_path):
    file = open(el29_file_path, 'rb')
    els_cs_info  = ElsepaCrossSectionInfo()
    els_cs_info .read_file(file)

    assert 30105020 == els_cs_info ._file_version
    assert 29 == els_cs_info ._atomic_number
    assert 1.000000e-01 == approx(els_cs_info ._energy_keV)
    assert 2.870967e-02 == approx(els_cs_info ._total_cs_nm2)

    assert 606 == els_cs_info ._number_points

    assert 0.0 == approx(els_cs_info ._ratios[0])
    assert 0.0 == approx(els_cs_info ._theta_rad[0])

    assert 1.0 == approx(els_cs_info ._ratios[-1])
    assert 3.1415926539 == approx(els_cs_info ._theta_rad[-1])
