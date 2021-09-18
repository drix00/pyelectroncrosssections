#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. py:currentmodule:: eecs.models.elsepa_cross_section_info
.. moduleauthor:: Hendrix Demers <hendrix.demers@mail.mcgill.ca>

Description
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
from casinotools.file_format.file_reader_writer_tools import read_int, read_double

# Project modules.

# Globals and constants variables.


class ElsepaCrossSectionInfo:
    def __init__(self):
        self._file_version = None
        self._atomic_number = None
        self._energy_keV = None
        self._total_cs_nm2 = None

        self._number_points = None
        self._ratios = None
        self._theta_rad = None

    def read_file(self, file):
        self._file_version = read_int(file)
        self._atomic_number = int(read_double(file))
        self._energy_keV = read_double(file)
        self._total_cs_nm2 = read_double(file)

        self._number_points = read_int(file)
        self._ratios = []
        self._theta_rad = []
        for dummy in range(self._number_points):
            ratio = read_double(file)
            theta_rad = read_double(file)

            self._ratios.append(ratio)
            self._theta_rad.append(theta_rad)

        assert len(self._ratios) == self._number_points
        assert len(self._theta_rad) == self._number_points
