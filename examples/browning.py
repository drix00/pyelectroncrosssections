#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. py:currentmodule:: browning
.. moduleauthor:: Hendrix Demers <hendrix.demers@mail.mcgill.ca>

Plot Browning cross section results.
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
import matplotlib.pyplot as plt
import numpy as np

# Local modules.
import pyHendrixDemersTools.Colors as Colors
from eecs.numeric_conversion import cm2_to_nm2

# Project modules.
from eecs.models.browning import total_elastic_cross_section_rutherford_cm2, average_scattering_angle_rutherford_deg, \
    total_elastic_cross_section_browning1991a_cm2, average_scattering_angle_rutherford_decreased_screening_deg, \
    differential_cross_section_browning1991_cm2_sr, ratio_browning1994, ratio_browning1994_mcx_ray, \
    compute_polar_angle_one_random_number_rad, compute_polar_angle_two_random_numbers_rad

# Globals and constants variables.


def plot_figure2_browning1991a():
    energies_keV = np.arange(0.1, 100.0, 0.1)
    elements = ['U', 'Mo', 'Al', 'C', 'H']
    atomic_numbers = {'U': 92, 'Mo': 42, 'Al': 13, 'C': 6, 'H': 1}

    colors = Colors.BaseColors(len(elements))
    plt.figure()

    for element in elements:
        Z = atomic_numbers[element]
        color = colors.next()
        cross_sections_cm2 = [total_elastic_cross_section_rutherford_cm2(Z, E_keV) for E_keV in energies_keV]
        cross_sections_pm2 = np.array(cross_sections_cm2)*1.0e16
        plt.loglog(energies_keV, cross_sections_pm2, '--', color=color)
        cross_sections_cm2 = [total_elastic_cross_section_browning1991a_cm2(Z, E_keV) for E_keV in energies_keV]
        cross_sections_pm2 = np.array(cross_sections_cm2)*1.0e16
        plt.loglog(energies_keV, cross_sections_pm2, color=color, label=element)

    plt.xlabel("Electron Energy (keV)")
    plt.ylabel(r"Total Elastic Cross Section ($10^{-16}$cm$^{2}$)")
    plt.legend(loc='best')

    plt.ylim((0.001, 10.0))


def plot_figure3_browning1991():
    energies_keV = np.arange(0.1, 100.0, 0.1)
    elements = ['Au']
    atomic_numbers = {'Au': 79}

    colors = Colors.BaseColors(len(elements))
    plt.figure()

    for element in elements:
        Z = atomic_numbers[element]
        color = colors.next()
        average_angles_deg = [average_scattering_angle_rutherford_deg(Z, E_keV) for E_keV in energies_keV]
        plt.plot(energies_keV, average_angles_deg, color=color, label=element)

        color = colors.next()
        average_angles_deg = [average_scattering_angle_rutherford_decreased_screening_deg(Z, E_keV)
                              for E_keV in energies_keV]
        plt.plot(energies_keV, average_angles_deg, color=color)

    plt.xlabel("Energy (keV)")
    plt.ylabel(r"Average Scattering Angle (Degrees)")
    plt.legend(loc='best')

    plt.ylim((0.0, 40.0))


def plot_figure4_browning1991():
    energy_keV = 1.0
    angles_rad = np.arange(0.0, np.pi, 0.01)
    angles_deg = np.degrees(angles_rad)
    elements = ['Au']
    atomic_numbers = {'Au': 79}

    colors = Colors.BaseColors(len(elements))
    plt.figure()

    for element in elements:
        Z = atomic_numbers[element]
        color = colors.next()
        differential_cross_sections = [differential_cross_section_browning1991_cm2_sr(Z, energy_keV, angle_rad)
                                       for angle_rad in angles_rad]
        plt.plot(angles_deg, differential_cross_sections, color=color, label=element)

    plt.xlabel("Angle (Degrees)")
    plt.ylabel(r"Cross Section (A$^{2}$)")
    plt.legend(loc='best')

    # plt.ylim((1.0e-4, 1.0e1))


def plot_comparison_ratio_mcxray_error():
    energies_keV = np.arange(0.1, 100.0, 0.1)
    elements = ['Au', 'Cu', 'C']
    atomic_numbers = {'Au': 79, 'Cu': 29, 'C': 6}

    colors = Colors.BaseColors(len(elements))
    plt.figure()

    for element in elements:
        Z = atomic_numbers[element]
        color = colors.next()
        ratios = np.array([ratio_browning1994(Z, E_keV) for E_keV in energies_keV])
        ratios = ratios / (1.0 + ratios)
        label = "{} Browning".format(element)
        plt.semilogx(energies_keV, ratios, color=color, label=label)

        color = colors.next()
        ratios = np.array([ratio_browning1994_mcx_ray(Z, E_keV) for E_keV in energies_keV])
        ratios = ratios / (1.0 + ratios)
        label = "{} MCXRay".format(element)
        plt.semilogx(energies_keV, ratios, color=color, label=label)

    plt.xlabel("Energy (keV)")
    plt.ylabel(r"Ratio / (1 + Ratio)")
    plt.legend(loc='best')

    # plt.ylim((0.0, 40.0))


def plot_comparison_polar_angle_random_number():
    number_samples = 10000000
    atomic_number = 79
    energy_keV = 1.0

    ratio = ratio_browning1994(atomic_number, energy_keV)
    print("Ratio = {:.4f}".format(ratio))
    print("Ratio/(1 + Ratio) = %.4f" % (ratio/(1.0 + ratio)))

    polar_angles_one_random_number_rad = [compute_polar_angle_one_random_number_rad(atomic_number, energy_keV)
                                          for _i in range(number_samples)]

    polar_angles_two_random_numbers_rad = [compute_polar_angle_two_random_numbers_rad(atomic_number, energy_keV)
                                           for _i in range(number_samples)]

    plt.figure()

    plt.hist(polar_angles_two_random_numbers_rad, bins=100, normed=True, label='2 RNs', histtype='step')
    plt.hist(polar_angles_one_random_number_rad, bins=100, normed=True, label='1 RN', histtype='step')

    plt.xlabel(r"Scattering Angle (rad)")
    plt.ylabel(r"Probabilities")

    plt.legend(loc='best')
    filename = "ComparaisonPolarAngleRandomNumber.pdf"
    filepath = os.path.join(r"J:\hdemers\work\documents\labbooks\e-labbook\graphics\BrowningCrossSection", filename)
    plt.savefig(filepath)


def plot_compute_differential_cross_section():
    number_bins = 50
    number_samples = 1000000
    atomic_number = 79
    energy_keV = 1.0

    total_cm2 = total_elastic_cross_section_browning1991a_cm2(atomic_number, energy_keV)
    total_nm2 = cm2_to_nm2(total_cm2)
    total_A2 = total_nm2 * 1.0e2

    ratio = ratio_browning1994(atomic_number, energy_keV)
    print("Ratio = {:.4f}" .format(ratio))
    print("Ratio/(1 + Ratio) = %.4f" % (ratio/(1.0 + ratio)))

    polar_angles_two_random_numbers_rad = [compute_polar_angle_two_random_numbers_rad(atomic_number, energy_keV)
                                           for _i in range(number_samples)]

    histogram, bin_edges = np.histogram(polar_angles_two_random_numbers_rad, numbins=number_bins)
    bin_size = bin_edges[1] - bin_edges[0]
    histogram /= number_samples*bin_size
    histogram *= total_A2

    max_range = bin_edges + bin_size*number_bins
    thetas_rad = np.arange(bin_edges, max_range, bin_size)

    plt.figure()

    # plt.hist(polar_angles_two_random_numbers_rad, bins=100, normed=True, histtype='step')
    plt.semilogy(thetas_rad, histogram)

    plt.xlabel(r"Scattering Angle (rad)")
    plt.ylabel(r"Probabilities (A2/sr)")
    plt.xlim(xmin=1.0*bin_size)

    # plt.legend(loc='best')
    filename = "ComputeDifferentialCrossSection.pdf"
    filepath = os.path.join(r"J:\hdemers\work\documents\labbooks\e-labbook\graphics\BrowningCrossSection", filename)
    plt.savefig(filepath)



def run():
    plot_figure2_browning1991a()

    # plot_figure3_browning1991()
    # plot_figure4_browning1991()

    plot_comparison_ratio_mcxray_error()

    # plot_comparison_polar_angle_random_number()
    # plot_compute_differential_cross_section()

    plt.show()


if __name__ == '__main__':  # pragma: no cover
    run()
