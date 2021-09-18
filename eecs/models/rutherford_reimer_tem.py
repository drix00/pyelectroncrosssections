#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. py:currentmodule:: eecs.models.rutherford_reimer_tem
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
import matplotlib.pyplot as plt
import scipy.constants as constants
import numpy as np


# Local modules.

# Project modules.

# Globals and constants variables.


# TODO: Born Wentzel does not work, dimension miss match.
def total_elastic_cross_section_born_wentzel(atomic_number, energy_eV):
    r"""
    From Reimer TEM book's page 152.

    .. math:

    \sigma_{el} &= \frac{Z^{2}R^{2}\lambda^{2}\left(1 + \frac{E}{E_{0}}\right)^{2}}{\pi a^{2}_{H}} \\
    &= \frac{h^{2} Z^{4/3}}{\pi E_{0}^2 \beta^{2}}

    """
    pi = constants.pi
    # h = constants.Planck
    h_eVs = 4.13566733e-15

    E0_eV = 511.0e3
    Z = atomic_number
    beta2 = _compute_beta2(energy_eV)

    nominator = h_eVs ** 2 * np.power(Z, 4.0 / 3.0)
    denominator = pi * E0_eV ** 2 * beta2

    total_nm2 = nominator / denominator
    return total_nm2


def total_elastic_cross_section_quantum_approximation_nm2(atomic_number, energy_eV):
    r"""
    From Reimer TEM book's page 153.

    .. math:

    \sigma_{el} &= \frac{1.5\time 10^{-6}}{\beta^{2}} Z^{3/2} \left(1 - 0.23 \frac{Z}{137 \beta}\right).

    for .. math:`\frac{Z}{137 \beta} < 1.2`.

    """
    Z = atomic_number
    beta = _compute_beta(energy_eV)

    limit_factor = Z / (137.0 * beta)
    # if limit_factor > 1.2:
    #     logging.warning("Model used out of his validity range: %i, %f", atomicNumber, energy_eV)

    factor_a = 1.5e-6 * np.power(Z, 3.0 / 2.0) / beta ** 2
    factor_b = 1.0 - 0.23 * limit_factor

    total_nm2 = factor_a * factor_b
    return total_nm2


def _compute_beta(energy_eV):
    beta2 = _compute_beta2(energy_eV)
    beta = np.sqrt(beta2)
    return beta


def _compute_beta2(energy_eV):
    E0_eV = 511.0e3
    temp = (1.0 / (1.0 + energy_eV / E0_eV))
    beta2 = 1.0 - temp * temp
    return beta2


def differential_rutherford_small_angle_nm2_str(atomicNumber, energy_eV, angle_rad):
    r"""
    From Reimer TEM book's page 151.

    .. math:

    \frac{d\sigma_{el}}{d\Omega} = \frac{4 Z^{2} R^{4} \left(\left[1 +
    \frac{E}{E_{0}}\right)^{2} }{a_{H}^{2}} \frac{1}{1 + \left(\frac{\theta}{\theta_{0}}\right)^{2}\right]^{2}}

    with .. math:`\theta_{0} = \frac{\lambda}{2 \pi R}` and .. math:`R = a_{H}Z^{-1/3}`.
    """
    Z = atomicNumber
    R_nm = _compute_r_nm(atomicNumber)
    aH_nm = get_bohr_radius_nm()
    E0_eV = 511.0e3
    theta0 = _compute_characteristic_angle(atomicNumber, energy_eV)
    E_eV = energy_eV

    nominator = 4.0 * Z ** 2 * R_nm ** 4 * (1.0 + E_eV / E0_eV) ** 2
    denominator = aH_nm ** 2
    factor_a = nominator / denominator

    denominator = (1.0 + (angle_rad / theta0) ** 2) ** 2
    factor_b = 1.0 / denominator

    differential_nm2_str = factor_a * factor_b

    return differential_nm2_str


def _compute_r_nm(atomic_number):
    aH_nm = get_bohr_radius_nm()
    R_nm = aH_nm * np.power(atomic_number, -1.0 / 3.0)
    return R_nm


def _compute_characteristic_angle(atomic_number, energy_eV):
    R_nm = _compute_r_nm(atomic_number)
    lambda_nm = _compute_lambda_nm(energy_eV)

    theta0_rad = lambda_nm / (2.0 * np.pi * R_nm)
    return theta0_rad


def _compute_lambda_nm(energy_eV):
    h_eVs = 4.13566733e-15
    c_m_s = constants.c
    E0_eV = 511.0e3
    E_eV = energy_eV

    nominator = h_eVs * c_m_s
    sqrt_arg = 2.0 * E_eV * E0_eV + E_eV ** 2
    denominator = np.sqrt(sqrt_arg)

    lambda_m = nominator / denominator

    lambda_nm = lambda_m * 1.0e9
    return lambda_nm


def get_bohr_radius_nm():
    return 0.0529


def total_relativistic_screened_elastic_cross_section_henoc_maurice_nm2(atomic_number, energy_eV):
    r"""
    From Joy et al. "Principles of Analytical Electron Microscopy" book's page 5.

    .. math:

    \sigma_{el} &= \frac{Z^{2}\lambda_{R}^{4}}{16 \pi^{3}a_{0}^{2}} \frac{1}{\delta\left(\delta + 1\right)}
    &= \frac{7.20\times 10^{13} Z^{2}\lambda_{R}^{4}}{\delta\left(\delta + 1\right)}

    where .. math:`\delta = \frac{\theta_{0}^{2}}{4}`.

    """
    factor_1_cm2 = _compute_total_factor_1_cm2()
    Z2 = atomic_number ** 2
    lambda_1_cm = _compute_lambda_cm(energy_eV)
    delta = _compute_delta(atomic_number, energy_eV)

    nominator = factor_1_cm2 * Z2 * lambda_1_cm ** 4
    denominator = delta * (delta + 1.0)

    total_cm2 = nominator / denominator
    total_nm2 = total_cm2 * 1.0e14

    return total_nm2


def _compute_total_factor_1_cm2():
    r"""
    From Joy et al. "Principles of Analytical Electron Microscopy" book's page 4.

    .. math:

    a_{0} &= \frac{\epsilon_{0} h^{2}}{\pi m_{0} e^{2}}
    &= 5.29\times10^{-9} (cm).

    """
    a0_cm = 5.29e-9
    factor = 1.0 / (16.0 * np.pi ** 3 * a0_cm ** 2)
    return factor


def _compute_lambda_cm(energy_eV):
    energy_keV = energy_eV * 1.0e-3

    nominator = 3.87e-9
    denominator = np.sqrt(energy_keV) * np.sqrt(1.0 + 9.79e-4 * energy_keV)

    lambda_cm = nominator / denominator

    return lambda_cm


def _compute_delta(atomic_number, energy_eV):
    energy_keV = energy_eV * 1.0e-3
    delta = 3.4e-3 * np.power(atomic_number, 2.0 / 3.0) / energy_keV

    return delta


def total_cross_section_greater_than_angle_nm2(atomic_number, energy_eV, angle_rad):
    r"""
    From Reimer TEM book's page 151.

    .. math:

    \frac{\sigma_{el}}{\alpha_{0}} = \frac{Z^{2} R^{2} \lambda^{2} \left(1 +
    \frac{E}{E_{0}}\right)^{2} }{\pi a_{H}^{2}} \right)^{2} \frac{1}{1 + \left(\frac{\alpha_{0}}{\theta_{0}}\right)^{2}}

    with .. math:`\theta_{0} = \frac{\lambda}{2 \pi R}` and .. math:`R = a_{H}Z^{-1/3}`.
    """
    Z = atomic_number
    R_nm = _compute_r_nm(atomic_number)
    aH_nm = get_bohr_radius_nm()
    E0_eV = 511.0e3
    theta0 = _compute_characteristic_angle(atomic_number, energy_eV)
    E_eV = energy_eV
    lambda_nm = _compute_lambda_nm(energy_eV)

    nominator = Z ** 2 * R_nm ** 2 * lambda_nm ** 2 * (1.0 + E_eV / E0_eV) ** 2
    denominator = np.pi * aH_nm ** 2
    factor_a = nominator / denominator

    denominator = (1.0 + (angle_rad / theta0) ** 2)
    factor_b = 1.0 / denominator

    total_nm2 = factor_a * factor_b

    return total_nm2


def figure_total():
    plt.figure()

    energy_eV = 200.0e3
    atomic_numbers = range(1, 100)

    eecs_models = {"Quantum approximation": total_elastic_cross_section_quantum_approximation_nm2,
                   "Born Wentzel": total_elastic_cross_section_born_wentzel}

    for label in eecs_models:
        eecs_model = eecs_models[label]
        totals_nm2 = [eecs_model(Z, energy_eV) for Z in atomic_numbers]
        plt.loglog(atomic_numbers, totals_nm2, label=label)

    plt.xlabel("Atomic number")
    plt.ylabel(r"$\sigma_{el}$ (nm$^{2}$)")
    plt.subplots_adjust(left=0.16)
    plt.legend(loc='best')

    plt.show()


def _figure_differential():
    plt.figure()
    atomic_number = 18
    energy_eV = 25.0e3

    angles_rad = np.logspace(np.log10(1.0e-4), np.log10(np.pi), num=1000)
    diff_str = []
    aH2 = get_bohr_radius_nm() ** 2

    for angle_rad in angles_rad:
        diff_str.append(differential_rutherford_small_angle_nm2_str(atomic_number, energy_eV, angle_rad) / aH2)

    plt.loglog(angles_rad, diff_str)
    plt.show()


def _figure_total():
    atomic_numbers = [6, 29, 79]
    models = {'QuantumApproximation': total_elastic_cross_section_quantum_approximation_nm2,
              'RSRutherford': total_relativistic_screened_elastic_cross_section_henoc_maurice_nm2}
    # models['BornWentzel'] = total_elastic_cross_section_born_wentzel

    start = np.log10(1.0)
    stop = np.log10(5.0e6)
    number_point = 10000
    energies_eV = np.logspace(start, stop, number_point)

    for atomicNumber in atomic_numbers:
        plt.figure()
        plt.title(atomicNumber)

        for modelName in models:
            model = models[modelName]
            total_nm2 = model(atomicNumber, energies_eV)
            plt.loglog(energies_eV, total_nm2, label=modelName)

        plt.legend(loc='best')
        plt.xlabel("Energy (eV)")
        plt.ylabel(r"$\sigma_{el}$ (nm$^{2}$)")

    plt.show()


def run():
    # _figure_differential()
    _figure_total()


if __name__ == '__main__':  # pragma: no cover
    run()
