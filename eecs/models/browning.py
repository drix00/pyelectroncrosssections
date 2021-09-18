#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. py:currentmodule:: eecs.models.Browning
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
import csv
from math import acos

# Third party modules.
import numpy as np
import scipy.stats
import scipy.integrate

# Local modules.
from eecs.numeric_conversion import cm2_to_nm2

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


def read_data(filepath):
    energies_eV = []
    totals_nm2 = {}
    mean_thetas_rad = {}

    reader = csv.reader(open(filepath, 'rb'))

    next(reader)

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
    total_nm2 = cm2_to_nm2(total_cm2)

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
