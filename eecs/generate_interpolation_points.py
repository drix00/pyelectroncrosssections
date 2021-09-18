#!/usr/bin/env python
""" """

# Script information for the file.
__author__ = "Hendrix Demers (hendrix.demers@mail.mcgill.ca)"
__version__ = ""
__date__ = ""
__copyright__ = "Copyright (c) 2009 Hendrix Demers"
__license__ = ""

# Subversion informations for the file.
__svnRevision__ = "$Revision$"
__svnDate__ = "$Date$"
__svnId__ = "$Id$"

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
from eecs.models.rutherford_reimer_tem import totalRelativisticScreenedElasticCrossSectionHenocMaurice_nm2

# Globals and constants variables.
SCALE_LINEAR = "linear"
SCALE_LOG10 = "log10"

class InterpolationModel(object):
    def __init__(self, x, y):
        pass

class GenerateInterpolationPoints(object):
    def __init__(self, initialGrid=None):
        self._scale = SCALE_LINEAR
        self._xMinimum = None
        self._xMaximum = None
        self._function = None
        self._interpolationModel = None
        self._interpolationFunction = None

        self._numberInitialPoints = 20
        self._initialGrid = initialGrid

    def setScale(self, scale):
        self._scale = scale

    def setXRange(self, xMinimum, xMaximum):
        self._xMinimum = xMinimum
        self._xMaximum = xMaximum

    def setFunction(self, function):
        self._function = function

    def setInterpolationModel(self, model):
        self._interpolationModel = model
        self._interpolationFunction = None

    def generate(self, errorPercentage=5.0, numberIntegrationCalls=None):
        logging.info("generate interpolation points")
        self._numberIntegrationCalls = numberIntegrationCalls

        errorFraction = errorPercentage/100.0

        x = self._generateInitialGrid()
        stopGeneration = False
        iteration = 0
        maxIteration = 500

        while not stopGeneration:
            yTrue = self._function(x)
            self._interpolationFunction = self._interpolationModel(x, yTrue)

            errors = self._computeErrors(x)

            oldX = x
            x = self._generateNewPoint(errors, errorFraction, x)

            iteration += 1

            logging.info("Iteration: %i -> %f", iteration, max(errors))
            if max(errors) < errorFraction or iteration >= maxIteration:
                stopGeneration = True

            if np.all(oldX == x):
                stopGeneration = True

        yTrue = self._function(x)
        assert len(x) == len(yTrue)
        self._xGrid = x
        self._yGrid = yTrue

    def _generateNewPoints(self, errors, errorFraction, x):
        newX = []
        index = 0
        for xx, error in zip(x[:-1], errors):
            #logging.info("%10.1e \t %g", xx, error)
            if error > errorFraction:
                newXX = xx + (x[index+1] - xx)/2.0
                if int(round(newXX)) != int(round(xx)) and int(round(newXX)) != int(round(x[index+1])):
                    newXX = round(newXX)
                    newX.append(newXX)
            index += 1

        if len(newX) > 0:
            x = x.tolist()
            x.extend(newX)
            x = list(set(x))
            x.sort()
            x = np.array(x)

        return x

    def _generateNewPoint(self, errors, errorFraction, x):
        maxError = max(errors)

        if maxError > errorFraction:
            index = np.argmax(errors)
            if self._scale == SCALE_LOG10:
                newLogX = math.log(x[index]) + (math.log(x[index+1]) - math.log(x[index]))/2.0
                newX = math.exp(newLogX)
            else:
                newX = x[index] + (x[index+1] - x[index])/2.0
            if int(round(newX)) != int(round(x[index])) and int(round(newX)) != int(round(x[index+1])):
                newX = round(newX)
                x = np.unique(np.sort(np.insert(x, index+1, newX)))
            logging.info("Add point: %f", newX)

        return x

    def _generateInitialGrid(self):
        logging.info("generateInitialGrid")
        if self._initialGrid == None:
            if self._scale == SCALE_LINEAR:
                start = self._xMinimum
                stop = self._xMaximum
                number = self._numberInitialPoints
                x = np.linspace(start, stop, number)

                return x
            elif self._scale == SCALE_LOG10:
                start = math.log10(self._xMinimum)
                stop = math.log10(self._xMaximum)
                number = self._numberInitialPoints
                x = np.logspace(start, stop, number)

                return x
            else:
                raise NotImplementedError
        else:
            return self._initialGrid

    def _computeErrors(self, xArray):
        errors = []
        for index in range(len(xArray[:-1])):
            xi = xArray[index]
            xi1 = xArray[index+1]

            error = self._computeError(xi, xi1)
            #x = xi + (xi1 - xi)/2.0
            value = self._computeValue(xi, xi1)
            errorRelative = error/value
            errors.append(errorRelative)

        return errors

    def _computeError(self, a, b):
        def func(x):
            return np.abs(self._function(x) - self._interpolationFunction(x))

        if self._numberIntegrationCalls == None:
            error, integrationError, info_dict, message, explain = integrate.quad(func, a, b)
        else:
            error, integrationError = integrate.fixed_quad(func, a, b, n=self._numberIntegrationCalls)
        logging.debug(integrationError)

        return error

    def _computeValue(self, a, b):
        def func(x):
            return self._function(x)

        if self._numberIntegrationCalls == None:
            error, integrationError, info_dict, message, explain = integrate.quad(func, a, b)
        else:
            error, integrationError = integrate.fixed_quad(func, a, b, n=self._numberIntegrationCalls)
        logging.debug(integrationError)

        return error

    def getPoints(self):
        return self._xGrid, self._yGrid

    def getErrors(self):
        x = []
        for index in range(len(self._xGrid[:-1])):
            xx = self._xGrid[index] + (self._xGrid[index+1] - self._xGrid[index])/2.0
            x.append(xx)

        errors = self._computeErrors(self._xGrid)
        return x, errors

class RunnerGrid(object):
    def __init__(self, initialGrid=None):
        self._start = 10.0
        self._end = 5.0e6
        self._atomicNumber = 79
        self._numberFunctionCalls = 0
        self._numberIntegrationCalls = 2
        self._errorPercentage = 0.1
        self._initialGrid = initialGrid

    def total_nm2(self, energy_eV):
        return totalRelativisticScreenedElasticCrossSectionHenocMaurice_nm2(self._atomicNumber, energy_eV)

    def _total_nm2(self, energy_eV):
        self._numberFunctionCalls += 1
        return self.total_nm2(energy_eV)

    def run(self, function):
        self._run()

        start = math.log10(self._start)
        stop = math.log10(self._end)
        x = np.logspace(start, stop, 10000)
        yTrue = function(x)

        xGrid, yGrid = self._gIntPoints.getPoints()
        intFunction = self._interpolationModel(xGrid, yGrid)
        yInt = intFunction(x)
        logging.info("Number grid points: %i", len(xGrid))

        self._saveGrid(xGrid, yGrid)

        _graphicGrid(x, yTrue, yInt, self._gIntPoints)

        _graphicError(x, yTrue, yInt, self._gIntPoints, self._errorPercentage)

        _graphicValueVariation(xGrid, yGrid)

        plt.show()

    def _run(self):
        self._gIntPoints = GenerateInterpolationPoints(self._initialGrid)
        self._gIntPoints.setScale(SCALE_LOG10)

        self._gIntPoints.setXRange(self._start, self._end)

        vfunc = vectorize1(self._total_nm2, args=(), vec_func=False)
        function = vfunc
        self._gIntPoints.setFunction(function)

        self._interpolationModel = interpolate.interp1d
        #interpolationModel = interpolate.InterpolatedUnivariateSpline

        self._gIntPoints.setInterpolationModel(self._interpolationModel)

        self._gIntPoints.generate(errorPercentage=self._errorPercentage, numberIntegrationCalls=self._numberIntegrationCalls)

        logging.info("Number function calls: %i", self._numberFunctionCalls)


    def _saveGrid(self, xGrid, yGrid):
        if self._numberIntegrationCalls == None:
            filename = "FullIntegration_%04.2f.txt" % (self._errorPercentage)
        else:
            filename = "Integration_N%02i_%04.2f.txt" % (self._numberIntegrationCalls, self._errorPercentage)

        filepath = os.path.join(r"C:\hdemers\tmp", filename)
        file = open(filepath, 'wb')

        for x, y in zip(xGrid, yGrid):
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

def _graphicGrid(x, yTrue, yInt, gIntPoints):
    xGrid, yGrid = gIntPoints.getPoints()
    plt.figure()
    plt.loglog(x, yTrue, '.', label="True")
    plt.loglog(xGrid, yGrid, 'o', label='Grid')
    plt.loglog(x, yInt, label="Interpolation")

    plt.legend(loc='best')
    plt.xlabel(r"$x$")
    plt.ylabel(r"$p(x)$")

def _graphicError(x, yTrue, yInt, gIntPoints, errorPercentage):
    xErrors, yErrors = gIntPoints.getErrors()
    plt.figure()
    plt.semilogx(x, (yTrue-yInt)/yTrue, label="error")
    plt.semilogx(xErrors, yErrors, label="yErrors")

    errorFraction = errorPercentage/100.0
    plt.axhline(-errorFraction)
    plt.axhline(errorFraction)

    plt.xlabel(r"$x$")
    plt.ylabel(r"error $\epsilon(x)$")
    plt.legend(loc='best')

def _graphicValueVariation(x, y):
    plt.figure()
    newX = x[:-1]

    newY = []
    deltaX = []
    for index in range(len(y[:-1])):
        difference = (y[index] - y[index+1])/y[index]
        newY.append(difference)

        differenceX = (x[index] - x[index+1])/x[index]
        deltaX.append(differenceX)

    plt.semilogx(newX, newY, label=r"$\frac{\Delta y}{y}$")
    plt.semilogx(newX, deltaX, label=r"$\frac{\Delta x}{x}$")

    errorFraction = 0.05
    plt.axhline(-errorFraction)
    plt.axhline(errorFraction)

    plt.xlabel(r"$x$")
    plt.ylabel(r"Value variation")
    plt.legend(loc='best')

def runSalvatAnalyticFunction():
    gIntPoints = GenerateInterpolationPoints()
    gIntPoints.setXRange(0.0, 5)

    def function(x):
        value = 7.0*x*np.exp(-4.0*x) + 0.6*np.exp(-12.5*(x - 3.5)**2)
        return value

    gIntPoints.setFunction(function)

    interpolationModel = interpolate.interpolate.interp1d
    gIntPoints.setInterpolationModel(interpolationModel)

    errorPercentage = 0.01
    gIntPoints.generate(errorPercentage=errorPercentage)

    x = np.linspace(0.0, 5.0, 10000)
    yTrue = function(x)

    xGrid, yGrid = gIntPoints.getPoints()
    intFunction = interpolationModel(xGrid, yGrid)
    yInt = intFunction(x)

    plt.figure()
    plt.plot(x, yTrue, label="True")
    plt.plot(xGrid, yGrid, 'o', label='Grid')
    plt.plot(x, yInt, label="Interpolation")

    plt.legend(loc='best')
    plt.xlabel(r"$x$")
    plt.ylabel(r"$p(x)$")

    plt.figure()
    plt.plot(x, (yTrue-yInt)/yTrue)

    plt.xlabel(r"$x$")
    plt.ylabel(r"error $\epsilon(x)$")
    errorFraction = errorPercentage/100.0
    plt.axhline(-errorFraction)
    plt.axhline(errorFraction)

    plt.show()

def run():
    #logging.getLogger().setLevel(logging.DEBUG)
    RunnerGrid().run()
