#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. py:currentmodule:: eecs.models.elsepa_casino
.. moduleauthor:: Hendrix Demers <hendrix.demers@mail.mcgill.ca>

Read ELSEPA CASINO files.
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
from zipfile import ZipFile

# Third party modules.
from scipy.interpolate import interp1d

# Local modules.

# Project modules.
from eecs.element_properties import get_symbol

# Globals and constants variables.
PREFIX_ANGLE = "A_"
PREFIX_PARTIAL = "P_"
PREFIX_TOTAL = "T_"
SUFFIX = ".txt"


class AngleElement:
    def __init__(self, energies_keV, random_numbers, angles_deg):
        self.angle_functions = {}
        for row_id, energy_keV in enumerate(energies_keV):
            self.angle_functions[energy_keV] = interp1d(random_numbers[row_id], angles_deg, kind='linear')

    def __call__(self, energy_keV, random_number):
        angle_deg = self.angle_functions[energy_keV](random_number)
        return angle_deg


class PartialElement:
    def __init__(self, energies_keV, angles_deg, partials_nm2_sr_2D):
        self.partial_functions = {}
        for row_id, energy_keV in enumerate(energies_keV):
            self.partial_functions[energy_keV] = interp1d(angles_deg, partials_nm2_sr_2D[row_id], kind='linear')

    def __call__(self, energy_eV, angle_deg):
        angle_deg = self.partial_functions[energy_eV](angle_deg)
        return angle_deg


class ElsepaCasino:
    def __init__(self, zip_filepath):
        self.zip_filepath = zip_filepath

        self.total_functions = {}
        self.angle_functions = {}
        self.partial_functions = {}

    def read_total_data(self, atomic_number):
        zip_file = ZipFile(self.zip_filepath, mode='r')

        symbol = get_symbol(atomic_number).capitalize()
        name = PREFIX_TOTAL + symbol + SUFFIX

        lines = zip_file.read(name).decode('ascii').split('\r\n')

        energies_keV = []
        totals_nm2 = []

        for line in lines[1:]:
            items = line.split('\t')
            if len(items) > 0:
                try:
                    energy_keV = float(items[0])
                    total_nm2 = float(items[1])
                    energies_keV.append(energy_keV)
                    totals_nm2.append(total_nm2)
                except ValueError:
                    pass

        self.total_functions[atomic_number] = interp1d(energies_keV, totals_nm2, kind='linear')

    def read_angle_data(self, atomic_number):
        zip_file = ZipFile(self.zip_filepath, mode='r')

        symbol = get_symbol(atomic_number).capitalize()
        name = PREFIX_ANGLE + symbol + SUFFIX

        lines = zip_file.read(name).decode('ascii').split('\r\n')

        energies_keV = []
        random_numbers_2D = []

        line = lines[0]
        angles_deg = [float(angle_deg) for angle_deg in line.split('\t')[1:]]

        for line in lines[1:]:
            items = line.strip().split('\t')
            if len(items) > 0:
                try:
                    energy_keV = float(items[0])
                    random_numbers = []
                    for item in items[1:]:
                        random_number = float(item)
                        random_numbers.append(random_number)

                    energies_keV.append(energy_keV)
                    random_numbers_2D.append(random_numbers)
                except ValueError:
                    pass

        assert len(energies_keV) == len(random_numbers_2D)
        assert len(angles_deg) == len(random_numbers_2D[0])
        assert len(angles_deg) == len(random_numbers_2D[-1])

        self.angle_functions[atomic_number] = AngleElement(energies_keV, random_numbers_2D, angles_deg)

    def read_partial_data(self, atomic_number):
        zip_file = ZipFile(self.zip_filepath, mode='r')

        symbol = get_symbol(atomic_number).capitalize()
        name = PREFIX_PARTIAL + symbol + SUFFIX

        lines = zip_file.read(name).decode('ascii').split('\r\n')

        energies_keV = []
        partials_nm2_sr_2D = []

        line = lines[0]
        angles_deg = [float(angle_deg) for angle_deg in line.split('\t')[1:]]

        for line in lines[1:]:
            items = line.strip().split('\t')
            if len(items) > 0:
                try:
                    energy_keV = float(items[0])
                    partials_nm2_sr = []
                    for item in items[1:]:
                        partial_nm2_sr = float(item)
                        partials_nm2_sr.append(partial_nm2_sr)

                    energies_keV.append(energy_keV)
                    partials_nm2_sr_2D.append(partials_nm2_sr)
                except ValueError:
                    pass

        assert len(energies_keV) == len(partials_nm2_sr_2D)
        assert len(angles_deg) == len(partials_nm2_sr_2D[0])
        assert len(angles_deg) == len(partials_nm2_sr_2D[-1])

        self.partial_functions[atomic_number] = PartialElement(energies_keV, angles_deg, partials_nm2_sr_2D)

    def total_nm2(self, atomic_number, energy_eV):
        if atomic_number not in self.total_functions:
            self.read_total_data(atomic_number)

        total_nm2 = self.total_functions[atomic_number](energy_eV)

        return total_nm2

    def angle_deg(self, atomic_number, energy_eV, random_number):
        if atomic_number not in self.angle_functions:
            self.read_angle_data(atomic_number)

        angle_deg = self.angle_functions[atomic_number](energy_eV, random_number)
        return angle_deg

    def partial_nm2_sr(self, atomic_number, energy_eV, angle_deg):
        if atomic_number not in self.partial_functions:
            self.read_partial_data(atomic_number)

        partial_nm2_sr = self.partial_functions[atomic_number](energy_eV, angle_deg)
        return partial_nm2_sr
