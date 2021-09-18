#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. py:currentmodule:: eecs.generate_tabulated_files
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
import os.path
import logging
import csv
import math

# Third party modules.
import scipy.integrate as integrate
import matplotlib.pyplot as plt

# Local modules.
from casinotools.file_format.casino3.models.cross_section_file import generate_raw_binary_files

# Project modules.
from eecs.element_properties import get_symbol

from eecs import csv2txt

# Globals and constants variables.


def generate_filename_with_symbol(prefix, extension, atomic_numbers):
    symbol = get_symbol(atomic_numbers)
    filename = prefix + f"_{symbol}" + extension
    return filename


def generate_filename(prefix, extension, atomic_numbers):
    filename = prefix + f"_{atomic_numbers:02d}" + extension
    return filename


class GenerateTabulatedFiles:
    def __init__(self):
        self._outputPath = None
        self._atomic_number = None

    def _generate_total_file(self, pathname, energies_grid_eV, totals_nm2):
        filename = self._generate_total_filename()

        cvs_file = self._create_csv_file(pathname, filename)

        row = ["Energy (eV)", "total (nm2)"]
        cvs_file.writerow(row)

        for energy_eV in energies_grid_eV:
            total_nm2 = totals_nm2[energy_eV]
            row = [energy_eV, total_nm2]
            cvs_file.writerow(row)

        del cvs_file

        path = os.path.join(self._outputPath, pathname)
        filepath = os.path.join(path, filename)
        csv2txt(filepath)

        plt.figure()
        x = energies_grid_eV
        y = [totals_nm2[energy_eV] for energy_eV in energies_grid_eV]
        plt.loglog(x, y)
        plt.xlabel(r"$E_{0}$ eV")
        plt.ylabel(r"$\sigma_{T}$ (nm$^{2}$)")

        figure_path = filepath[:-4] + ".pdf"
        plt.savefig(figure_path)

    def _create_csv_file(self, pathname, filename):
        path = os.path.join(self._outputPath, pathname)
        if not os.path.isdir(path):
            os.makedirs(path)

        filepath = os.path.join(path, filename)
        logging.info(filepath)

        csv_file = csv.writer(open(filepath, 'wb'))
        return csv_file

    def _generate_partial_file(self, pathname, polar_angle_grid_deg, energies_grid_eV, partials_nm2_sr):
        filename = self._generate_partial_filename()

        csv_file = self._create_csv_file(pathname, filename)

        row = ["Energy (eV)"]
        for angle_deg in polar_angle_grid_deg:
            row.append(angle_deg)
        csv_file.writerow(row)

        for energy_eV in energies_grid_eV:
            row = [energy_eV]
            for partial_nm2_sr in partials_nm2_sr[energy_eV]:
                row.append(partial_nm2_sr)

            csv_file.writerow(row)

        path = os.path.join(self._outputPath, pathname)
        filepath = os.path.join(path, filename)
        csv2txt(filepath)

    def _generate_partial_angles_file(self, pathname, polar_angles_grid_deg,
                                      energies_grid_eV, totals_nm2, partials_nm2_sr):
        filename = self._generate_partial_angles_filename()

        csv_file = self._create_csv_file(pathname, filename)

        row = ["Energy (eV)"]

        for angle_deg in polar_angles_grid_deg:
            row.append(angle_deg)

        csv_file.writerow(row)

        for energy_eV in energies_grid_eV:
            row = [energy_eV]
            partial_sin_thetas_nm2_sr = []
            polar_angle_grid_rad = []
            for angle_deg, partial_nm2_sr in zip(polar_angles_grid_deg, partials_nm2_sr[energy_eV]):
                angle_rad = math.radians(angle_deg)
                polar_angle_grid_rad.append(angle_rad)
                partial_sin_thetas_nm2_sr.append(partial_nm2_sr*math.sin(angle_rad)*2.0*math.pi)

            # logging.debug("%0.1f", energy_eV)
            # total_nm2 = totals_nm2[energy_eV]
            # logging.debug("%0.3e", total_nm2)
            # logging.debug("%0.3e", integrate.trapz(partial_sin_thetas_nm2_sr, polar_angle_grid_rad))
            # logging.debug("%0.3e", integrate.simps(partial_sin_thetas_nm2_sr, polar_angle_grid_rad))

            computed_total_nm2 = integrate.trapz(partial_sin_thetas_nm2_sr, polar_angle_grid_rad)
            for index in range(1, len(partial_sin_thetas_nm2_sr)+1):
                x = polar_angle_grid_rad[:index]
                y = partial_sin_thetas_nm2_sr[:index]
                ratio = integrate.trapz(y, x)/computed_total_nm2
                row.append(ratio)

            csv_file.writerow(row)

        path = os.path.join(self._outputPath, pathname)
        filepath = os.path.join(path, filename)
        csv2txt(filepath)

    def _generate_binary_files(self, path_name, atomic_number, energies_grid_eV,
                               totals_nm2, polar_angles_grid_deg, partial_nm2_sr):
        filepath = self._generate_binary_filepath(path_name)
        totals_list_nm2 = [totals_nm2[energy_eV] for energy_eV in energies_grid_eV]
        generate_raw_binary_files(filepath, atomic_number, energies_grid_eV, totals_list_nm2,
                                  polar_angles_grid_deg, partial_nm2_sr)

    def _generate_binary_filepath(self, path_name):
        filename = self._generate_binary_filename()
        path = os.path.join(self._outputPath, path_name)
        filepath = os.path.join(path, filename)
        return filepath

    def _generate_total_filename(self):
        return generate_filename_with_symbol("T", ".csv", self._atomic_number)

    def _generate_partial_filename(self):
        return generate_filename_with_symbol("P", ".csv", self._atomic_number)

    def _generate_partial_angles_filename(self):
        return generate_filename_with_symbol("A", ".csv", self._atomic_number)

    def _generate_binary_filename(self):
        return generate_filename_with_symbol("CS", ".bin", self._atomic_number)
