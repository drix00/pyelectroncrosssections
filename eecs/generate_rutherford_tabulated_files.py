#!/usr/bin/env python
""" """

# Script information for the file.
__author__ = "Hendrix Demers (hendrix.demers@mail.mcgill.ca)"
__version__ = ""
__date__ = ""
__copyright__ = "Copyright (c) 2009 Hendrix Demers"
__license__ = ""

# Subversion informations for the file.
__svnRevision__ = "$Revision: 2292 $"
__svnDate__ = "$Date: 2011-03-21 11:29:50 -0400 (Mon, 21 Mar 2011) $"
__svnId__ = "$Id: generate_rutherford_tabulated_files.py 2292 2011-03-21 15:29:50Z hdemers $"

# Standard library modules.
import logging
import os

# Third party modules.

# Local modules.
import pyHendrixDemersTools.Files as Files
import eecs.generate_interpolation_points as GenerateInterpolationPoints
from eecs.models.rutherford_reimer_tem import totalRelativisticScreenedElasticCrossSectionHenocMaurice_nm2

# Globals and constants variables.

class RutherfordRunnerGrid(GenerateInterpolationPoints.RunnerGrid):
    def __init__(self, errorPercentage, initialGrid, atomicnumber):
        super(RutherfordRunnerGrid, self).__init__(initialGrid)

        self._start = 10.0
        self._end = 500.0e3
        self._errorPercentage = errorPercentage
        self._atomicNumber = atomicnumber

    def run(self):
        self._run()

        return self._gIntPoints.getPoints()

    def total_nm2(self, energy_eV):
        return totalRelativisticScreenedElasticCrossSectionHenocMaurice_nm2(self._atomicNumber, energy_eV)

class GenerateRutherfordTabulatedFiles(object):
    def __init__(self, atomicNumber):
        logging.info("GenerateRutherfordTabulatedFiles for %i", atomicNumber)
        self._atomicNumber = atomicNumber
        self._errorPercentage = "0.1"

    def setInputPath(self, path):
        self._inputPath = path

    def setOutputPath(self, path):
        self._outputPath = path

    def run(self):
        logging.info("run")
        self._generateInterpolationEnergyGrid()

#        self._generateTotalFile()
#        self._generatePartialFile()
#        self._generatePartialAnglesFile()
#        self._generateBinaryFiles()

    def _generateInterpolationEnergyGrid(self):
        self._initialEnergiesGrid_eV = self._createInitialEnergyGrid_eV()
        runnerGrid = RutherfordRunnerGrid(float(self._errorPercentage), self._initialEnergiesGrid_eV, self._atomicNumber)

        self._energiesGrid_eV, self._totals_nm2 = runnerGrid.run()

        logging.info("Original number of points: %i", len(self._initialEnergiesGrid_eV))
        logging.info("Nre grid number of points: %i", len(self._energiesGrid_eV))

    def _createInitialEnergyGrid_eV(self):
        energyGrid_eV = []
        energyGrid_eV.extend(range(10, 100, 10))
        energyGrid_eV.extend(range(100, 1000, 100))
        energyGrid_eV.extend(range(1000, 10000, 1000))
        energyGrid_eV.extend(range(10000, 100000, 10000))
        energyGrid_eV.extend(range(100000, 500000+1, 20000))

        logging.debug("number of Initial Energy Grid: %i", len(energyGrid_eV))
        return energyGrid_eV

    def getEnergiesGrid(self):
        return self._energiesGrid_eV

def runCarbon():
    atomicNumber = 6

    print(_runElement(atomicNumber))

def runAllElements():
    energiesGridList = {}

    for atomicNumber in range(1,99+1):
        energiesGrid = _runElement(atomicNumber)
        energiesGridList[atomicNumber] = energiesGrid

    outputPath = _getOutputPath()

    filepath = os.path.join(outputPath, "RutherfordEnergiesGridList.txt")
    file = open(filepath, 'wb')

    for atomicNumber in sorted(energiesGridList.keys()):
        line = "%i" % atomicNumber
        for energy in energiesGridList[atomicNumber]:
            line += "\t%g" % (energy)
        file.write(line+"\n")

    file.close()

def _runElement(atomicNumber):
    outputPath = _getOutputPath()
    tabulatedFiles = GenerateRutherfordTabulatedFiles(atomicNumber)
    tabulatedFiles.setOutputPath(outputPath)
    tabulatedFiles.run()

    return tabulatedFiles.getEnergiesGrid()

def _getOutputPath():
    configurationFilepath = Files.getCurrentModulePath(__file__, "eecs.cfg")
    outputPath = Files.getResultsSherbrookePath(configurationFilepath, "calculations/Rutherford")

    if not os.path.isdir(outputPath):
        os.makedirs(outputPath)

    return outputPath

if __name__ == '__main__':    #pragma: no cover
    import pyHendrixDemersTools.Runner as Runner
    Runner.Runner().run(runFunction=runAllElements)
