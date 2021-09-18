#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. py:currentmodule:: tests.models.test_browning
.. moduleauthor:: Hendrix Demers <hendrix.demers@mail.mcgill.ca>

Tests for the module :py:mod:`eecs.models.browning`.
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

# Project modules
from eecs.models.browning import total_elastic_cross_section_browning1994_cm2, polar_angle_rad

# Globals and constants variables.


def test_is_discovered():
    """
    Test used to validate the file is included in the tests
    by the test framework.
    """
    # assert False
    assert True


def test_total_elastic_cross_section_browning1994_cm2_atomic_number():
    """
    Tests for method :py:meth:`total_elastic_cross_section_browning1994_cm2`.
    """

    total_css_A2_ref = {(6, 1.0): 0.55816651729734323,
                        (6, 10.0): 0.061015054708437634,
                        (6, 20.0): 0.030813100110956096,
                        (6, 100.0): 0.0062434494508425758,
                        (6, 1000.0): 0.0006288343892704956,
                        (13, 1.0): 1.5556811447079389,
                        (13, 10.0): 0.20830407903592069,
                        (13, 20.0): 0.10785138063688208,
                        (26, 1.0): 2.7799108935122673,
                        (26, 10.0): 0.53846911168966925,
                        (26, 20.0): 0.29584540228492651,
                        (47, 1.0): 3.4646680196969259,
                        (47, 10.0): 0.97140022086807613,
                        (47, 20.0): 0.58144059494789129,
                        (79, 1.0): 3.6626344033951284,
                        (79, 10.0): 1.3288532834221803,
                        (79, 20.0): 0.86137481241591683,
                        }

    for atomic_number, energy_keV in total_css_A2_ref:
        total_cs_A2_ref = total_css_A2_ref[(atomic_number, energy_keV)]
        total_cs_cm2 = total_elastic_cross_section_browning1994_cm2(atomic_number, energy_keV)
        total_cs_A2 = total_cs_cm2*1.0e16
        assert total_cs_A2_ref == approx(total_cs_A2, 6)


def test_total_elastic_cross_section_browning1994_cm2():
    """
    Tests for method :py:meth:`total_elastic_cross_section_browning1994_cm2`.
    """

    energy_keV = 100.0

    total_css_A2_ref = {6: 0.006243,
                        13: 0.02260,
                        26: 0.0676705532253224,
                        47: 0.15470732661851455,
                        79: 0.2734892890031806}

    for atomic_number in total_css_A2_ref:
        total_cs_A2_ref = total_css_A2_ref[atomic_number]
        total_cs_cm2 = total_elastic_cross_section_browning1994_cm2(atomic_number, energy_keV)
        total_cs_A2 = total_cs_cm2*1.0e16
        assert total_cs_A2_ref == approx(total_cs_A2, 6)


def test_angle_rad():
    angles_rad_ref = {(6, 1.0, 0.0, 0.0): 0.0,
                      (6, 1.0, 0.0, 0.1): 0.0,
                      (6, 1.0, 0.0, 0.5): 0.0,
                      (6, 1.0, 0.0, 0.9): 0.0,
                      (6, 1.0, 0.0, 1.0): 0.0,
                      (6, 1.0, 0.1, 0.0): 0.055568829753649025,
                      (6, 1.0, 0.1, 0.1): 0.055568829753649025,
                      (6, 1.0, 0.1, 0.5): 0.055568829753649025,
                      (6, 1.0, 0.1, 0.9): 0.055568829753649025,
                      (6, 1.0, 0.1, 1.0): 0.64350110879328426,
                      (6, 1.0, 0.5, 0.5): 0.16636462643999517,
                      (6, 1.0, 0.9, 0.5): 0.49019091996002279,
                      (6, 1.0, 1.0, 0.5): 3.1415926535897931,
                      (6, 20.0, 0.0, 0.5): 0.0,
                      (6, 20.0, 0.1, 0.5): 0.012469847640606605,
                      (6, 20.0, 0.5, 0.5): 0.037405665625122229,
                      (6, 20.0, 0.9, 0.5): 0.1121125047652253,
                      (6, 20.0, 1.0, 0.5): 3.1415917435203582,
                      (6, 30.0, 0.0, 0.5): 0.0,
                      (6, 30.0, 0.1, 0.5): 0.010182225703635034,
                      (6, 30.0, 0.5, 0.5): 0.030544566044577722,
                      (6, 30.0, 0.9, 0.5): 0.091576770150234882,
                      (6, 30.0, 1.0, 0.5): 3.1415926535897931,
                      (6, 50.0, 0.0, 0.5): 0.0,
                      (6, 50.0, 0.1, 0.5): 0.0078875133756613411,
                      (6, 50.0, 0.5, 0.5): 0.023661558793672852,
                      (6, 50.0, 0.9, 0.5): 0.070958200143830599,
                      (6, 50.0, 1.0, 0.5): 3.1415917435203582,
                      (6, 100.0, 0.0, 0.5): 0.0,
                      (6, 100.0, 0.1, 0.5): 0.0055775238406077154,
                      (6, 100.0, 0.5, 0.5): 0.016732224515452667,
                      (6, 100.0, 0.9, 0.5): 0.050187307871050975,
                      (6, 100.0, 1.0, 0.5): 3.1415917435203582,
                      }

    for key in angles_rad_ref:
        atomic_number, energy_keV, random_number1, random_number2 = key
        angle_rad = polar_angle_rad(atomic_number, energy_keV, random_number1, random_number2)
        assert angles_rad_ref[key] == approx(angle_rad, 6)
