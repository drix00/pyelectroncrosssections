#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. py:currentmodule:: Models.Browning

.. moduleauthor:: Hendrix Demers <hendrix.demers@mail.mcgill.ca>

Cross section models from Browning.
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
import math
import os.path
import csv
from math import acos

# Third party modules.
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import scipy.integrate

# Local modules.
import pyHendrixDemersTools.Colors as Colors
from pyHendrixDemersTools.NumericConversion import cm2Tonm2

# Project modules

# Globals and constants variables.


def total_elastic_cross_section_rutherford_cm2(atomic_number, electron_energy_keV):
    """
    From browning1991a ref 1: J
    D. C. Joy, Proceedings of the 9th European Congress of Electron Microscopy,
    edited by P. J. Goodhew and H. G. Dickinson (Institute of Physics,
    Bristol, UK, 1988), p. 22.

    Without relativistic effects correction.
    Screened Rutherford cross section in cm2

    """

    factor = 5.21e-21
    term_a = atomic_number * atomic_number / (electron_energy_keV * electron_energy_keV)

    alpha = compute_screening_parameter(atomic_number, electron_energy_keV)
    term_b = math.pi/(alpha * (1.0 + alpha))

    cross_section_cm2 = factor * term_a * term_b

    return cross_section_cm2


def compute_screening_parameter(atomic_number, energy_keV):
    factor = 3.4e-3
    term_a = math.pow(atomic_number, 0.67) / energy_keV

    alpha = factor * term_a
    return alpha


def average_scattering_angle_rutherford_deg(atomic_number, energy_keV):
    alpha = compute_screening_parameter(atomic_number, energy_keV)

    average_angle_rad = math.pi * math.sqrt(alpha) * math.sqrt(1.0 + alpha) - math.pi * alpha
    average_angle_deg = math.degrees(average_angle_rad)
    return average_angle_deg


def average_scattering_angle_rutherford_decreased_screening_deg(atomic_number, energy_keV):
    alpha = compute_decreased_screening_parameter(atomic_number, energy_keV)

    average_angle_rad = math.pi * math.sqrt(alpha) * math.sqrt(1.0 + alpha) - math.pi * alpha
    average_angle_deg = math.degrees(average_angle_rad)
    return average_angle_deg


def compute_decreased_screening_parameter(atomic_number, energy_keV):
    alpha_rutherford = compute_screening_parameter(atomic_number, energy_keV)
    alpha = (0.6 - 0.0035 * energy_keV) * alpha_rutherford
    return alpha


def total_elastic_cross_section_browning1991a_cm2(atomic_number, energy_keV):
    """
    From browning1991a
    Valid in the range 1 to 100 keV for all elements.
    """
    Z = atomic_number
    E = energy_keV
    factor = 4.7e-18
    nominator = math.pow(Z, 1.33) + 0.032*Z*Z
    denominator = E + 0.0155 * math.pow(Z, 1.33) * math.pow(E, 0.5)
    term_a = nominator/denominator

    u = compute_factor_u(atomic_number, energy_keV)
    denominator = 1.0 - 0.02 * math.pow(Z, 0.5) * math.exp(-u*u)
    term_b = 1.0/denominator

    cross_section_cm2 = factor * term_a * term_b

    return cross_section_cm2


def compute_factor_u(atomic_number, energy_keV):
    u = math.log10(8.0) * energy_keV * math.pow(atomic_number, -1.33)
    return u


def total_elastic_cross_section_browning1994_cm2(atomic_number, energy_keV):
    """
    From browning1994
    Valid in the range 100 eV to 30 keV for elements 1 to 92.
    """
    Z = atomic_number
    E = energy_keV
    factor = 3.0e-18
    power_z = math.pow(Z, 1.7)
    power_e = math.pow(E, 0.5)
    nominator = factor*power_z
    denominator = E + 0.005 * power_z * power_e + 0.0007 * Z * Z / power_e
    cross_section_cm2 = nominator/denominator

    return cross_section_cm2


def differential_cross_section_browning1991_cm2_sr(atomic_number, energy_keV, theta_rad):
    factor = 5.21e-21
    Z = atomic_number
    E = energy_keV
    term_a = Z * Z / (E * E)

    alpha = compute_screening_parameter_browning1991(atomic_number, energy_keV)
    denominator = 1.0 - math.cos(theta_rad) - alpha
    term_b = 1.0 / denominator

    nominator = alpha * (alpha + 1.0)
    denominator = 4.2 * math.pow(E, 1.1)
    term_c = nominator / denominator

    differential_cross_section_cm2_sr = factor * term_a * (term_b + term_c)

    return differential_cross_section_cm2_sr


def compute_screening_parameter_browning1991(atomic_number, energy_keV):
    alpha = 5.5e-4 * math.pow(atomic_number, 0.67) / energy_keV
    return alpha


def ratio_browning1994(atomic_number, energy_keV):
    Z = atomic_number
    E = energy_keV
    term_a = 300.0 * math.pow(E, 1.0 - Z/2000.0) / Z

    term_b = math.pow(Z, 3.0) / (3.0e5 * E)

    ratio = term_a + term_b
    return ratio


def ratio_browning1994_mcx_ray(atomic_number, energy_keV):
    Z = atomic_number
    E = energy_keV
    term_a = 300.0 * math.pow(E, 1.0 - Z/2000.0) / Z

    term_b = math.pow(Z, 3.0) / 3.0e5 * E

    ratio = term_a + term_b
    return ratio


def compute_polar_angle_two_random_numbers_rad(atomic_number, energy_keV):
    random_number1 = np.random.random()

    return polar_angle_rad(atomic_number, energy_keV, random_number1, random_number1)


def compute_polar_angle_one_random_number_rad(atomic_number, energy_keV):
    random_number1 = np.random.random()

    return polar_angle_rad(atomic_number, energy_keV, random_number1, random_number1)


def polar_angle_rad(atomic_number, energy_keV, random_number1, random_number2):
    ratio = ratio_browning1994(atomic_number, energy_keV)
    if random_number2 <= ratio / (1.0 + ratio):
        alpha = 7.0e-3 / energy_keV
        cos_theta = 1.0 - 2.0 * alpha * random_number1 / (1.0 + alpha - random_number1)
    else:
        cos_theta = 1.0 - 2.0 * random_number1

    if cos_theta > 1.0:
        return 0.0
    elif cos_theta < -1.0:
        return math.pi

    theta_rad = acos(cos_theta)
    return theta_rad


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
    total_nm2 = cm2Tonm2(total_cm2)
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


def read_data(filepath):
    energies_eV = []
    totals_nm2 = {}
    mean_thetas_rad = {}

    reader = csv.reader(open(filepath, 'rb'))

    reader.next()

    for row in reader:
        energy_eV = float(row[0])
        meanTheta_rad = float(row[1])
        total_nm2 = float(row[2])

        energies_eV.append(energy_eV)
        mean_thetas_rad[energy_eV] = meanTheta_rad
        totals_nm2[energy_eV] = total_nm2

    del reader

    return mean_thetas_rad, totals_nm2


def save_data(filepath, energies_eV, meanThetas_rad, totals_nm2):
    writer = csv.writer(open(filepath, 'wb'))

    row_header = ["Energy (eV)", "Mean Theta (rad)", "Total (nm2)"]
    writer.writerow(row_header)

    for energy_eV in energies_eV:
        row = [energy_eV, meanThetas_rad[energy_eV], totals_nm2[energy_eV]]
        writer.writerow(row)


def compute_mean_theta_total_browning(atomic_number, energy_eV):
    number_bins = 50
    number_samples = 1000000

    energy_keV = energy_eV*1.0e-3
    total_cm2 = total_elastic_cross_section_browning1991a_cm2(atomic_number, energy_keV)
    total_nm2 = cm2Tonm2(total_cm2)

    polar_angles_two_random_numbers_rad = [compute_polar_angle_two_random_numbers_rad(atomic_number, energy_keV)
                                           for _i in range(number_samples)]

    histogram, bin_edges = np.histogram(polar_angles_two_random_numbers_rad, numbins=number_bins)
    bin_size = bin_edges[1] - bin_edges[0]
    histogram /= number_samples*bin_size
    histogram *= total_nm2

    max_range = bin_edges + bin_size*number_bins
    thetas_rad = np.arange(bin_edges, max_range, bin_size)

    # partial_sin_thetas_nm2_sr = [2.0*np.pi*x1*np.sin(x2) for x1,x2 in zip(histogram, thetas_rad)]
    partial_sin_thetas_nm2_sr = [x1 for x1, x2 in zip(histogram, thetas_rad)]
    total_calculated_nm2 = scipy.integrate.trapz(partial_sin_thetas_nm2_sr, thetas_rad)

    # partial_sin_thetas_nm2_sr = [2.0*np.pi*x1*x2*np.sin(x2) for x1,x2 in zip(histogram, thetas_rad)]
    partial_sin_thetas_nm2_sr = [x1*x2 for x1, x2 in zip(histogram, thetas_rad)]
    mean_theta_rad = scipy.integrate.trapz(partial_sin_thetas_nm2_sr, thetas_rad)
    mean_theta_rad = mean_theta_rad/total_calculated_nm2

    return mean_theta_rad, total_calculated_nm2


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
