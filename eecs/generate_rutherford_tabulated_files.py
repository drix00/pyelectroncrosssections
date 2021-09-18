#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. py:currentmodule:: eecs.generate_rutherford_tabulated_files
.. moduleauthor:: Hendrix Demers <hendrix.demers@mail.mcgill.ca>

Generate tabulated file for Rutherford model.
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
import logging
import os

# Third party modules.

# Local modules.

# Project modules.
from eecs.generate_interpolation_points import RunnerGrid
from eecs.models.rutherford_reimer_tem import total_relativistic_screened_elastic_cross_section_henoc_maurice_nm2

# Globals and constants variables.


class RutherfordRunnerGrid(RunnerGrid):
    def __init__(self, error_percentage, initial_grid, atomic_number):
        super().__init__(initial_grid)

        self._start = 10.0
        self._end = 500.0e3
        self._error_percentage = error_percentage
        self._atomic_number = atomic_number

    def run(self):
        self._run()

        return self._g_int_points.get_points()

    def total_nm2(self, energy_eV):
        return total_relativistic_screened_elastic_cross_section_henoc_maurice_nm2(self._atomic_number, energy_eV)


class GenerateRutherfordTabulatedFiles:
    def __init__(self, atomic_number):
        logging.info("GenerateRutherfordTabulatedFiles for %i", atomic_number)
        self._atomic_number = atomic_number
        self._error_percentage = "0.1"

        self._input_path = None
        self._output_path = None

    def set_input_path(self, path):
        self._input_path = path

    def set_output_path(self, path):
        self._output_path = path

    def run(self):
        logging.info("run")
        self._generate_interpolation_energy_grid()

#        self._generate_total_file()
#        self._generate_partial_file()
#        self._generate_partial_angles_file()
#        self._generate_binary_files()

    def _generate_interpolation_energy_grid(self):
        self._initialEnergiesGrid_eV = self._create_initial_energy_grid_eV()
        runner_grid = RutherfordRunnerGrid(float(self._error_percentage), self._initialEnergiesGrid_eV, self._atomic_number)

        self._energiesGrid_eV, self._totals_nm2 = runner_grid.run()

        logging.info("Original number of points: %i", len(self._initialEnergiesGrid_eV))
        logging.info("Nre grid number of points: %i", len(self._energiesGrid_eV))

    @staticmethod
    def _create_initial_energy_grid_eV():
        energy_grid_eV = []
        energy_grid_eV.extend(range(10, 100, 10))
        energy_grid_eV.extend(range(100, 1000, 100))
        energy_grid_eV.extend(range(1000, 10000, 1000))
        energy_grid_eV.extend(range(10000, 100000, 10000))
        energy_grid_eV.extend(range(100000, 500000+1, 20000))

        logging.debug("number of Initial Energy Grid: %i", len(energy_grid_eV))
        return energy_grid_eV

    def get_energies_grid(self):
        return self._energiesGrid_eV


def run_carbon():
    atomic_number = 6
    print(_run_element(atomic_number))


def run_all_elements():
    energies_grid_list = {}

    for atomic_number in range(1, 99+1):
        energies_grid = _run_element(atomic_number)
        energies_grid_list[atomic_number] = energies_grid

    output_path = _get_output_path()

    filepath = os.path.join(output_path, "RutherfordEnergiesGridList.txt")
    file = open(filepath, 'w')

    for atomic_number in sorted(energies_grid_list.keys()):
        line = "%i" % atomic_number
        for energy in energies_grid_list[atomic_number]:
            line += f"\t{energy:g}"
        file.write(line + "\n")

    file.close()


def _run_element(atomic_number):
    output_path = _get_output_path()
    tabulated_files = GenerateRutherfordTabulatedFiles(atomic_number)
    tabulated_files.set_output_path(output_path)
    tabulated_files.run()

    return tabulated_files.get_energies_grid()


def _get_output_path():
    # configuration_filepath = get_current_module_path(__file__, "eecs.cfg")
    output_path = "calculations/Rutherford"

    if not os.path.isdir(output_path):
        os.makedirs(output_path)

    return output_path
