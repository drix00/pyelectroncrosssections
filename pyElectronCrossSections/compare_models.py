#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. py:currentmodule:: pyElectronCrossSection.compare_modles

.. moduleauthor:: Hendrix Demers <hendrix.demers@mail.mcgill.ca>

Compare elastic cross section models.
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
import os.path

# Third party modules.
import numpy as np
import matplotlib.pyplot as plt

# Local modules.

# Project modules.
from pyElectronCrossSections.Models.Browning import total_elastic_cross_section_browning1991a_cm2, total_elastic_cross_section_rutherford_cm2
from pyElectronCrossSections.Models.ElsepaBinaryFile import ElsepaBinaryFile
from pyElectronCrossSections import current_module_path

# Globals and constants variables.


def compare_total():
    path = r"C:\hdemers\bin\casino3\latest\Tabulated\Elsepa"

    energies_keV = np.arange(0.1, 100.0, 0.1)
    element = 'C'
    atomic_number = 6

    plt.figure()
    plt.title(element)

    Z = atomic_number

    crossSections_cm2 = [total_elastic_cross_section_rutherford_cm2(Z, E_keV) for E_keV in energies_keV]
    crossSections_nm2 = np.array(crossSections_cm2)*1.0e14
    plt.plot(energies_keV, crossSections_nm2, '--', label="Rutherford")

    crossSections_cm2 = [total_elastic_cross_section_browning1991a_cm2(Z, E_keV) for E_keV in energies_keV]
    crossSections_nm2 = np.array(crossSections_cm2)*1.0e14
    plt.plot(energies_keV, crossSections_nm2, label="Browning")

    filename = "CS_{}.bin".format(element)
    filepath = os.path.join(path, filename)
    elsepa_file = ElsepaBinaryFile(filepath)

    casino_energies_keV = []
    totals_nm2 = []
    for data in elsepa_file._elsCSInfoList:
        energy_keV = data._energy_keV
        total_nm2 = data._totalCS_nm2

        casino_energies_keV.append(energy_keV)
        totals_nm2.append(total_nm2)

    crossSections_nm2 = np.array(totals_nm2)
    plt.plot(casino_energies_keV, crossSections_nm2, label="ELSEPA CASINO")

    plt.xlabel("Electron Energy (keV)")
    plt.ylabel(r"Total Elastic Cross Section $nm$^{2}$)")
    plt.legend()

    plt.xlim((5, 100))
    plt.ylim((3.0e-5, 2.0e-3))
    plt.yscale('log')

    plt.savefig("compare_total_C.png")


if __name__ == '__main__':  # pragma: no cover
    compare_total()

    plt.show()
