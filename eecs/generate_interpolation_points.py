#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. py:currentmodule:: eecs.generate_interpolation_points
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
import logging
import math
import os.path

# Third party modules.
import numpy as np
from scipy import interpolate
from scipy import integrate
import matplotlib.pyplot as plt

# Local modules.

# Project modules.
from eecs.models.rutherford_reimer_tem import total_relativistic_screened_elastic_cross_section_henoc_maurice_nm2

# Globals and constants variables.
SCALE_LINEAR = "linear"
SCALE_LOG10 = "log10"


class GenerateInterpolationPoints:
    def __init__(self, initial_grid=None):
        self._scale = SCALE_LINEAR
        self._x_minimum = None
        self._x_maximum = None
        self._function = None
        self._interpolation_model = None
        self._interpolation_function = None

        self._number_initial_points = 20
        self._initial_grid = initial_grid
        self._number_integration_calls = None

        self._x_grid = None
        self._y_grid = None

    def set_scale(self, scale):
        self._scale = scale

    def set_x_range(self, x_minimum, x_maximum):
        self._x_minimum = x_minimum
        self._x_maximum = x_maximum

    def set_function(self, function):
        self._function = function

    def set_interpolation_model(self, model):
        self._interpolation_model = model
        self._interpolation_function = None

    def generate(self, error_percentage=5.0, number_integration_calls=None):
        logging.info("generate interpolation points")
        self._number_integration_calls = number_integration_calls

        error_fraction = error_percentage / 100.0

        x = self._generate_initial_grid()
        stop_generation = False
        iteration = 0
        max_iteration = 500

        while not stop_generation:
            y_true = self._function(x)
            self._interpolation_function = self._interpolation_model(x, y_true)

            errors = self._compute_errors(x)

            old_x = x
            x = self._generate_new_point(errors, error_fraction, x)

            iteration += 1

            logging.info("Iteration: %i -> %f", iteration, max(errors))
            if max(errors) < error_fraction or iteration >= max_iteration:
                stop_generation = True

            if np.all(old_x == x):
                stop_generation = True

        y_true = self._function(x)
        assert len(x) == len(y_true)
        self._x_grid = x
        self._y_grid = y_true

    @staticmethod
    def _generate_new_points(errors, error_fraction, x):
        new_x = []
        index = 0
        for xx, error in zip(x[:-1], errors):
            # logging.info("%10.1e \t %g", xx, error)
            if error > error_fraction:
                new_xx = xx + (x[index + 1] - xx) / 2.0
                if int(round(new_xx)) != int(round(xx)) and int(round(new_xx)) != int(round(x[index + 1])):
                    new_xx = round(new_xx)
                    new_x.append(new_xx)
            index += 1

        if len(new_x) > 0:
            x = x.tolist()
            x.extend(new_x)
            x = list(set(x))
            x.sort()
            x = np.array(x)

        return x

    def _generate_new_point(self, errors, error_fraction, x):
        max_error = max(errors)

        if max_error > error_fraction:
            index = np.argmax(errors)
            if self._scale == SCALE_LOG10:
                new_log_x = math.log(x[index]) + (math.log(x[index + 1]) - math.log(x[index])) / 2.0
                new_x = math.exp(new_log_x)
            else:
                new_x = x[index] + (x[index + 1] - x[index]) / 2.0
            if int(round(new_x)) != int(round(x[index])) and int(round(new_x)) != int(round(x[index + 1])):
                new_x = round(new_x)
                x = np.unique(np.sort(np.insert(x, index + 1, new_x)))
            logging.info("Add point: %f", new_x)

        return x

    def _generate_initial_grid(self):
        logging.info("generateInitialGrid")
        if self._initial_grid is None:
            if self._scale == SCALE_LINEAR:
                start = self._x_minimum
                stop = self._x_maximum
                number = self._number_initial_points
                x = np.linspace(start, stop, number)

                return x
            elif self._scale == SCALE_LOG10:
                start = math.log10(self._x_minimum)
                stop = math.log10(self._x_maximum)
                number = self._number_initial_points
                x = np.logspace(start, stop, number)

                return x
            else:
                raise NotImplementedError
        else:
            return self._initial_grid

    def _compute_errors(self, x_array):
        errors = []
        for index in range(len(x_array[:-1])):
            xi = x_array[index]
            xi1 = x_array[index + 1]

            error = self._compute_error(xi, xi1)
            # x = xi + (xi1 - xi)/2.0
            value = self._compute_value(xi, xi1)
            error_relative = error / value
            errors.append(error_relative)

        return errors

    def _compute_error(self, a, b):
        def func(x):
            return np.abs(self._function(x) - self._interpolation_function(x))

        if self._number_integration_calls is None:
            error, integration_error, info_dict, message, explain = integrate.quad(func, a, b)
        else:
            error, integration_error = integrate.fixed_quad(func, a, b, n=self._number_integration_calls)
        logging.debug(integration_error)

        return error

    def _compute_value(self, a, b):
        def func(x):
            return self._function(x)

        if self._number_integration_calls is None:
            error, integration_error, info_dict, message, explain = integrate.quad(func, a, b)
        else:
            error, integration_error = integrate.fixed_quad(func, a, b, n=self._number_integration_calls)
        logging.debug(integration_error)

        return error

    def get_points(self):
        return self._x_grid, self._y_grid

    def get_errors(self):
        x = []
        for index in range(len(self._x_grid[:-1])):
            xx = self._x_grid[index] + (self._x_grid[index + 1] - self._x_grid[index]) / 2.0
            x.append(xx)

        errors = self._compute_errors(self._x_grid)
        return x, errors


class RunnerGrid:
    def __init__(self, initial_grid=None):
        self._start = 10.0
        self._end = 5.0e6
        self._atomic_number = 79
        self._number_function_calls = 0
        self._number_integration_calls = 2
        self._error_percentage = 0.1
        self._initial_grid = initial_grid

    def total_nm2(self, energy_eV):
        return total_relativistic_screened_elastic_cross_section_henoc_maurice_nm2(self._atomic_number, energy_eV)

    def _total_nm2(self, energy_eV):
        self._number_function_calls += 1
        return self.total_nm2(energy_eV)

    def run(self):
        self._run()

        start = math.log10(self._start)
        stop = math.log10(self._end)
        x = np.logspace(start, stop, 10000)
        y_true = function(x)

        x_grid, y_grid = self._g_int_points.get_points()
        int_function = self._interpolationModel(x_grid, y_grid)
        y_int = int_function(x)
        logging.info("Number grid points: %i", len(x_grid))

        self._save_grid(x_grid, y_grid)

        _graphic_grid(x, y_true, y_int, self._g_int_points)

        _graphic_error(x, y_true, y_int, self._g_int_points, self._error_percentage)

        _graphic_value_variation(x_grid, y_grid)

        plt.show()

    def _run(self):
        self._g_int_points = GenerateInterpolationPoints(self._initial_grid)
        self._g_int_points.set_scale(SCALE_LOG10)

        self._g_int_points.set_x_range(self._start, self._end)

        vfunc = vectorize1(self._total_nm2, args=(), vec_func=False)
        function = vfunc
        self._g_int_points.set_function(function)

        self._interpolationModel = interpolate.interp1d
        # interpolationModel = interpolate.InterpolatedUnivariateSpline

        self._g_int_points.set_interpolation_model(self._interpolationModel)

        self._g_int_points.generate(error_percentage=self._error_percentage,
                                    number_integration_calls=self._number_integration_calls)

        logging.info("Number function calls: %i", self._number_function_calls)

    def _save_grid(self, x_grid, y_grid):
        if self._number_integration_calls is None:
            filename = f"FullIntegration_{self._error_percentage:04.2f}.txt" % ()
        else:
            filename = "Integration_N%02i_%04.2f.txt" % (self._number_integration_calls, self._error_percentage)

        filepath = os.path.join(r"C:\hdemers\tmp", filename)
        file = open(filepath, 'w')

        for x, y in zip(x_grid, y_grid):
            line = "%f\t%f\n" % (x, y)
            file.write(line)

        file.close()


def vectorize1(func, args=(), vec_func=False):
    if vec_func:
        def vfunc(x):
            return func(x, *args)
    else:
        def vfunc(x):
            if np.isscalar(x):
                return func(x, *args)
            x = np.asarray(x)
            # call with first point to get output type
            y0 = func(x[0], *args)
            n = len(x)
            if hasattr(y0, 'dtype'):
                output = np.empty((n,), dtype=y0.dtype)
            else:
                output = np.empty((n,), dtype=type(y0))
            output[0] = y0
            for i in range(1, n):
                output[i] = func(x[i], *args)
            return output
    return vfunc


def _graphic_grid(x, y_true, y_int, g_int_points):
    x_grid, y_grid = g_int_points.get_points()
    plt.figure()
    plt.loglog(x, y_true, '.', label="True")
    plt.loglog(x_grid, y_grid, 'o', label='Grid')
    plt.loglog(x, y_int, label="Interpolation")

    plt.legend(loc='best')
    plt.xlabel(r"$x$")
    plt.ylabel(r"$p(x)$")


def _graphic_error(x, y_true, y_int, g_int_points, error_percentage):
    x_errors, y_errors = g_int_points.get_errors()
    plt.figure()
    plt.semilogx(x, (y_true - y_int) / y_true, label="error")
    plt.semilogx(x_errors, y_errors, label="y_errors")

    error_fraction = error_percentage / 100.0
    plt.axhline(-error_fraction)
    plt.axhline(error_fraction)

    plt.xlabel(r"$x$")
    plt.ylabel(r"error $\epsilon(x)$")
    plt.legend(loc='best')


def _graphic_value_variation(x, y):
    plt.figure()
    new_x = x[:-1]

    new_y = []
    delta_x = []
    for index in range(len(y[:-1])):
        difference = (y[index] - y[index + 1]) / y[index]
        new_y.append(difference)

        difference_x = (x[index] - x[index + 1]) / x[index]
        delta_x.append(difference_x)

    plt.semilogx(new_x, new_y, label=r"$\frac{\Delta y}{y}$")
    plt.semilogx(new_x, delta_x, label=r"$\frac{\Delta x}{x}$")

    error_fraction = 0.05
    plt.axhline(-error_fraction)
    plt.axhline(error_fraction)

    plt.xlabel(r"$x$")
    plt.ylabel(r"Value variation")
    plt.legend(loc='best')


def run_salvat_analytic_function():
    g_int_points = GenerateInterpolationPoints()
    g_int_points.set_x_range(0.0, 5)

    def function(local_x):
        value = 7.0 * x * np.exp(-4.0 * local_x) + 0.6 * np.exp(-12.5 * (local_x - 3.5) ** 2)
        return value

    g_int_points.set_function(function)

    interpolation_model = interpolate.interpolate.interp1d
    g_int_points.set_interpolation_model(interpolation_model)

    error_percentage = 0.01
    g_int_points.generate(error_percentage=error_percentage)

    x = np.linspace(0.0, 5.0, 10000)
    y_true = function(x)

    x_grid, y_grid = g_int_points.get_points()
    int_function = interpolation_model(x_grid, y_grid)
    y_int = int_function(x)

    plt.figure()
    plt.plot(x, y_true, label="True")
    plt.plot(x_grid, y_grid, 'o', label='Grid')
    plt.plot(x, y_int, label="Interpolation")

    plt.legend(loc='best')
    plt.xlabel(r"$x$")
    plt.ylabel(r"$p(x)$")

    plt.figure()
    plt.plot(x, (y_true - y_int) / y_true)

    plt.xlabel(r"$x$")
    plt.ylabel(r"error $\epsilon(x)$")
    error_fraction = error_percentage / 100.0
    plt.axhline(-error_fraction)
    plt.axhline(error_fraction)

    plt.show()


def run():
    RunnerGrid().run()
