#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. py:currentmodule:: module_name
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
from pathlib import Path

# Third party modules.
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import interpolate
from numpy.random import default_rng
from scipy import integrate

# Local modules.

# Project modules.

# Globals and constants variables.


def plot_total(path):
    file_path = path / "T_C.txt"

    data = np.loadtxt(file_path, skiprows=1)

    energies_eV = data[:, 0]
    totals_nm = data[:, 1]

    f_total = interpolate.interp1d(energies_eV, totals_nm)

    plt.figure()

    plt.plot(energies_eV, totals_nm, '.', label="data")

    energies_new_eV = np.linspace(energies_eV[0], energies_eV[-1], 100000)
    plt.plot(energies_new_eV, f_total(energies_new_eV), 'o', label="int")

    plt.xscale("log")
    plt.yscale("log")

    plt.xlabel("Energy (eV)")
    plt.ylabel("Total (nm)")


def skip_str(item):
    try:
        value = float(item)
    except ValueError:
        value = 0.0
    return value


def plot_partial(path):
    file_path = path / "P_C.txt"

    data = np.loadtxt(file_path, converters={0: skip_str}, delimiter="\t")

    angles_deg = data[0, 1:]
    energies_eV = data[1:, 0]
    partial_nm_sr = data[1:, 1:]

    print(data.shape)
    print(energies_eV.shape)
    print(angles_deg.shape)
    print(partial_nm_sr.shape)

    energy_id = 75

    plt.figure()
    title = f"Energy: {energies_eV[energy_id]:.1f} eV"
    plt.title(title)

    plt.plot(angles_deg, partial_nm_sr[energy_id, :], '.', label="data")

    # plt.xscale("log")
    plt.yscale("log")

    plt.xlabel("Angle (deg)")
    plt.ylabel("Partial (nm/sr)")


def plot_angle(path):
    file_path = path / "A_C.txt"

    data = np.loadtxt(file_path, converters={0: skip_str}, delimiter="\t")

    angles_deg = data[0, 1:]
    energies_eV = data[1:, 0]
    r_values = data[1:, 1:]

    print(data.shape)
    print(energies_eV.shape)
    print(angles_deg.shape)
    print(r_values.shape)

    energy_id = 75

    plt.figure()
    title = f"Energy: {energies_eV[energy_id]:.1f} eV"
    plt.title(title)

    plt.plot(angles_deg, r_values[energy_id, :], '.', label="data")

    # plt.xscale("log")
    # plt.yscale("log")

    plt.xlabel("Angle (deg)")
    plt.ylabel("R")

    X, Y = np.meshgrid(angles_deg, energies_eV)
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    # ax.plot_wireframe(X, Y, r_values)
    ax.plot_surface(X, Y, r_values, cmap=cm.coolwarm)
    ax.set_ylim(bottom=0.0, top=30.0e3)


def plot_angle_random(path):
    file_path = path / "P_C.txt"

    data_partial = np.loadtxt(file_path, converters={0: skip_str}, delimiter="\t")

    angles_deg = data_partial[0, 1:]
    energies_eV = data_partial[1:, 0]
    partial_nm_sr = data_partial[1:, 1:]

    print(data_partial.shape)
    print(energies_eV.shape)
    print(angles_deg.shape)
    print(partial_nm_sr.shape)

    file_path = path / "A_C.txt"
    data = np.loadtxt(file_path, converters={0: skip_str}, delimiter="\t")
    r_values = data[1:, 1:]

    energy_id = 75

    f_angle = interpolate.interp1d(r_values[energy_id, :], angles_deg)

    rng = default_rng()
    number_points = 1000000
    randoms = rng.random(number_points)
    values = f_angle(randoms)
    print(values.shape)

    plt.figure()
    title = f"Energy: {energies_eV[energy_id]:.1f} eV"
    plt.title(title)

    area = integrate.simpson(partial_nm_sr[energy_id, :], angles_deg)
    print(area)

    plt.plot(angles_deg, partial_nm_sr[energy_id, :], '-', label="data")
    plt.yscale("log")
    plt.hist(values, bins=len(angles_deg), label="random")

    # plt.xscale("log")
    plt.yscale("log")

    plt.xlabel("Angle (deg)")
    plt.ylabel("Partial (nm/sr)")
    plt.legend()


def main():
    path = Path(r"../test_data")

    # plot_total(path)
    # plot_partial(path)
    plot_angle(path)
    # plot_angle_random(path)


if __name__ == '__main__':
    main()

    plt.show()
