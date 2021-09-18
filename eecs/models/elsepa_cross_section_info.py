#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. py:currentmodule:: pyElectronCrossSection.models.ElsepaCrossSectionInfo

.. moduleauthor:: Hendrix Demers <hendrix.demers@mail.mcgill.ca>


"""

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

# Standard library modules.

# Third party modules.

# Local modules.
from casinotools.file_format.file_reader_writer_tools import read_int, read_double

# Project modules.

# Globals and constants variables.


class ElsepaCrossSectionInfo:
    def read_file(self, file):
        self._fileVersion = read_int(file)
        self._atomicNumber = int(read_double(file))
        self._energy_keV = read_double(file)
        self._totalCS_nm2 = read_double(file)

        self._numberPoints = read_int(file)
        self._ratios = []
        self._theta_rad = []
        for dummy in range(self._numberPoints):
            ratio = read_double(file)
            theta_rad = read_double(file)

            self._ratios.append(ratio)
            self._theta_rad.append(theta_rad)

        assert len(self._ratios) == self._numberPoints
        assert len(self._theta_rad) == self._numberPoints
