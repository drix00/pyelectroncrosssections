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
import os.path
import logging
import csv
import math

# Third party modules.
import scipy.integrate as integrate
import matplotlib.pyplot as plt

# Local modules.
from casinotools.file_format.casino3.models.cross_section_file import generate_raw_binary_files
from eecs.element_properties import getSymbol

from eecs import csv2txt

# Globals and constants variables.

class GenerateTabulatedFiles(object):
    def _generateTotalFile(self, pathname, energiesGrid_eV, totals_nm2):
        filename = self._generateTotalFilename()

        cvsFile = self._createCsvFile(pathname, filename)

        row = ["Energy (eV)", "total (nm2)"]
        cvsFile.writerow(row)

        for energy_eV in energiesGrid_eV:
            total_nm2 = totals_nm2[energy_eV]
            row = [energy_eV, total_nm2]
            cvsFile.writerow(row)

        del cvsFile

        path = os.path.join(self._outputPath, pathname)
        filepath = os.path.join(path, filename)
        csv2txt(filepath)

        plt.figure()
        x = energiesGrid_eV
        y = [totals_nm2[energy_eV] for energy_eV in energiesGrid_eV]
        plt.loglog(x, y)
        plt.xlabel(r"$E_{0}$ eV")
        plt.ylabel(r"$\sigma_{T}$ (nm$^{2}$)")

        figurePath = filepath[:-4] + ".pdf"
        plt.savefig(figurePath)

    def _createCsvFile(self, pathname, filename):
        path = os.path.join(self._outputPath, pathname)
        if not os.path.isdir(path):
            os.makedirs(path)

        filepath = os.path.join(path, filename)
        logging.info(filepath)

        csvFile = csv.writer(open(filepath, 'wb'))
        return csvFile

    def _generatePartialFile(self, pathname, polarAngleGrid_deg, energiesGrid_eV, partials_nm2_sr):
        filename = self._generatePartialFilename()

        csvFile = self._createCsvFile(pathname, filename)

        row = ["Energy (eV)"]
        for angle_deg in polarAngleGrid_deg:
            row.append(angle_deg)
        csvFile.writerow(row)

        for energy_eV in energiesGrid_eV:
            row = [energy_eV]
            for partial_nm2_sr in partials_nm2_sr[energy_eV]:
                row.append(partial_nm2_sr)

            csvFile.writerow(row)

        path = os.path.join(self._outputPath, pathname)
        filepath = os.path.join(path, filename)
        csv2txt(filepath)

    def _generatePartialAnglesFile(self, pathname, polarAnglesGrid_deg,
                                                                 energiesGrid_eV, totals_nm2, partials_nm2_sr):
        filename = self._generatePartialAnglesFilename()

        csvFile = self._createCsvFile(pathname, filename)

        row = ["Energy (eV)"]

        for angle_deg in polarAnglesGrid_deg:
            row.append(angle_deg)

        csvFile.writerow(row)

        for energy_eV in energiesGrid_eV:
            total_nm2 = totals_nm2[energy_eV]
            row = [energy_eV]
            partialSinThetas_nm2_sr = []
            polarAngleGrid_rad = []
            for angle_deg, partial_nm2_sr in zip(polarAnglesGrid_deg, partials_nm2_sr[energy_eV]):
                angle_rad = math.radians(angle_deg)
                polarAngleGrid_rad.append(angle_rad)
                partialSinThetas_nm2_sr.append(partial_nm2_sr*math.sin(angle_rad)*2.0*math.pi)

            #logging.debug("%0.1f", energy_eV)
            #logging.debug("%0.3e", total_nm2)
            #logging.debug("%0.3e", integrate.trapz(partialSinThetas_nm2_sr, polarAngleGrid_rad))
            #logging.debug("%0.3e", integrate.simps(partialSinThetas_nm2_sr, polarAngleGrid_rad))

            computedTotal_nm2 = integrate.trapz(partialSinThetas_nm2_sr, polarAngleGrid_rad)
            for index in range(1, len(partialSinThetas_nm2_sr)+1):
                x = polarAngleGrid_rad[:index]
                y = partialSinThetas_nm2_sr[:index]
                ratio = integrate.trapz(y, x)/computedTotal_nm2
                row.append(ratio)

            csvFile.writerow(row)

        path = os.path.join(self._outputPath, pathname)
        filepath = os.path.join(path, filename)
        csv2txt(filepath)

    def _generateBinaryFiles(self, pathName, atomicNumber, energiesGrid_eV,
                                                     totals_nm2, polarAnglesGrid_deg, partial_nm2_sr):
        filepath = self._generateBinaryFilepath(pathName)
        totalsList_nm2 = [totals_nm2[energy_eV] for energy_eV in energiesGrid_eV]
        generate_raw_binary_files(filepath, atomicNumber, energiesGrid_eV, totalsList_nm2,
                                                polarAnglesGrid_deg, partial_nm2_sr)

    def _generateBinaryFilepath(self, pathName):
        filename = self._generateBinaryFilename()
        path = os.path.join(self._outputPath, pathName)
        filepath = os.path.join(path, filename)
        return filepath

    def _generateTotalFilename(self):
        return self._generateFilenameWithSymbol("T", ".csv", self._atomicNumber)

    def _generatePartialFilename(self):
        return self._generateFilenameWithSymbol("P", ".csv", self._atomicNumber)

    def _generatePartialAnglesFilename(self):
        return self._generateFilenameWithSymbol("A", ".csv", self._atomicNumber)

    def _generateBinaryFilename(self):
        return self._generateFilenameWithSymbol("CS", ".bin", self._atomicNumber)

    def _generateFilename(self, prefix, extension, atomicNumbers):
        filename =    prefix + "_%02i" % (atomicNumbers) + extension
        return filename

    def _generateFilenameWithSymbol(self, prefix, extension, atomicNumbers):
        symbol = getSymbol(atomicNumbers)
        filename = prefix + "_%s" % (symbol) + extension
        return filename
