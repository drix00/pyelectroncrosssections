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
__svnId__ = "$Id: analyze_czyzewski90.py 2292 2011-03-21 15:29:50Z hdemers $"

# Standard library modules.
import logging
import math

# Third party modules.
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from scipy import interpolate

# Local modules.
import DatabasesTools.Czyzewski90.MottTabulatedFactory as MottTabulatedFactory
import eecs.models.rutherford_reimer_tem as RutherfordReimerTem

# Globals and constants variables.

def _multipagesMultifiguresAtomicNumbers(atomicNumbers, energies_eV, mottTabulated):
    logging.info("_multipagesMultifiguresAtomicNumbers")

    pdf = PdfPages('Czyzewski90_total_energy.pdf')

    for atomicNumber in atomicNumbers[::6]:
        for index in xrange(6):
            atomicNumberIndex = atomicNumber + index
            if atomicNumberIndex <= atomicNumbers[-1]:
                total_nm2 = mottTabulated.getTotals_nm2(atomicNumber=atomicNumberIndex)

                assert len(energies_eV) == len(total_nm2)

                plotNumber = index + 1
                axes = plt.subplot(3, 2, plotNumber)
                plt.loglog(energies_eV, total_nm2, '.')

                plt.text(0.9, 0.9, atomicNumberIndex,
                                 horizontalalignment='center',
                                 verticalalignment='center',
                                 transform = axes.transAxes)

                plt.subplots_adjust(wspace=0.35, hspace=0.35)

                if plotNumber in [5, 6]:
                    plt.xlabel("Energy (eV)")
                if plotNumber in [1, 3, 5]:
                    plt.ylabel(r"$\sigma_{el}$ (nm$^{2}$)")

        pdf.savefig()
        plt.close()

    pdf.close()

def _multipagesMultifiguresEnergies(atomicNumbers, energies_eV, mottTabulated):
    logging.info("_multipagesMultifiguresEnergies")

    pdf = PdfPages('Czyzewski90_total_atomicNumber.pdf')

    totals_nm2 = {}
    for atomicNumber in atomicNumbers:
        totals_nm2[atomicNumber] = mottTabulated.getTotals_nm2(atomicNumber=atomicNumber)

    maxBigIndex = len(energies_eV)
    for indexByStep in xrange(0, maxBigIndex, 6):
        for indexSmall in xrange(6):
            index = indexByStep + indexSmall
            if index < len(energies_eV):
                y = [totals_nm2[Z][index] for Z in atomicNumbers]

                plotNumber = indexSmall + 1
                axes = plt.subplot(3, 2, plotNumber)
                plt.plot(atomicNumbers, y, '.')

                energy_eV = energies_eV[index]
                if energy_eV < 1.0e3:
                    label = "%i eV" % (energy_eV)
                else:
                    label = "%i keV" % (energy_eV/1.0e3)

                plt.text(0.15, 0.9, label,
                                 horizontalalignment='center',
                                 verticalalignment='center',
                                 transform = axes.transAxes)

                plt.subplots_adjust(wspace=0.35, hspace=0.35)

                if plotNumber in [5, 6]:
                    plt.xlabel("Atomic number")
                if plotNumber in [1, 3, 5]:
                    plt.ylabel(r"$\sigma_{el}$ (nm$^{2}$)")

        pdf.savefig()
        plt.close()

    pdf.close()

def _interpolationComparison(atomicNumbers, energies_eV, mottTabulated):
    logging.info("_multipagesMultifiguresEnergies")

    pdf = PdfPages('Interpolation_Czyzewski90_total_atomicNumber.pdf')

    totals_nm2 = {}
    for atomicNumber in atomicNumbers:
        totals_nm2[atomicNumber] = mottTabulated.getTotals_nm2(atomicNumber=atomicNumber)

    interpolationModels = ['linear', 'cubic']

    for dummy_index, atomicNumber in enumerate(atomicNumbers):
        x = np.array(energies_eV)
        y = np.array(totals_nm2[atomicNumber])

        interpolations = {}
        for interpolationModel in interpolationModels:
            interpolations[interpolationModel] = interpolate.interp1d(x, y, kind=interpolationModel)

#        interpolations['univariateSpline'] = interpolate.UnivariateSpline(x, y)
#        logging.info("univariateSpline residual: %f", interpolations['univariateSpline'].get_residual())
#        interpolationModels.append('univariateSpline')

        interpolations['interpolatedUnivariateSpline'] = interpolate.InterpolatedUnivariateSpline(x, y)
        logging.info("interpolatedUnivariateSpline residual: %f", interpolations['interpolatedUnivariateSpline'].get_residual())
        interpolationModels.append('interpolatedUnivariateSpline')

#        interpolations['LSQUnivariateSpline'] = interpolate.LSQUnivariateSpline(x, y, t=x[1:-1:2])
#        logging.info("LSQUnivariateSpline residual: %f", interpolations['LSQUnivariateSpline'].get_residual())
#        interpolationModels.append('LSQUnivariateSpline')

#        tck = interpolate.splrep(x, y, s=0)
#        def interpolateSpline(x):
#            return interpolate.splev(x, tck, der=0)
#        interpolations['spline'] = interpolateSpline
#        interpolationModels.append('spline')

        start = math.log10(energies_eV[0])
        stop = math.log10(energies_eV[-1]-1.0)
        numberPoint = 1000
        newX = np.logspace(start, stop, numberPoint)

        fig = plt.figure()
        axes = fig.add_subplot(211)
        left = fig.subplotpars.left
        bottom = 0.3 + fig.subplotpars.hspace
        bottom = 0.3 + fig.subplotpars.hspace/2.0
        width = fig.subplotpars.right - fig.subplotpars.left
        height = fig.subplotpars.top - bottom
        box = [left, bottom, width, height]
        axes.set_position(box)
        plt.loglog(x, y, 'o', label='Raw data')

        for interpolationModel in interpolationModels:
            newY = interpolations[interpolationModel](newX)
            plt.loglog(newX, newY, label=str(interpolationModel))

        plt.text(0.9, 0.9, atomicNumber,
                         horizontalalignment='center',
                         verticalalignment='center',
                         transform = axes.transAxes)

        plt.legend(loc='best')
        plt.ylabel(r"$\sigma_{el}$ (nm$^{2}$)")

        axes = fig.add_subplot(212, sharex=axes)
        bottom = fig.subplotpars.bottom
        width = fig.subplotpars.right - fig.subplotpars.left
        height = 0.3 - bottom
        box = [left, bottom, width, height]
        axes.set_position(box)

        interpolationModelRef = 'cubic'
        refY = interpolations[interpolationModelRef](x)
        error = (refY - y)/refY
        plt.semilogx(x, error, 'o', label=str(interpolationModelRef))
        logging.info("%s sum error: %f", interpolationModelRef, sum(error))

        refY = interpolations[interpolationModelRef](newX)
        for interpolationModel in interpolationModels:
            error = (refY - interpolations[interpolationModel](newX))/refY
            plt.semilogx(newX, error, label=str(interpolationModel))
            logging.info("%s sum error: %f", interpolationModel, sum(error))

        plt.xlabel("Energy (eV)")
        plt.ylabel("Error relative")
        #plt.ylim(ymin=-1, ymax=1)
        #pdf.savefig()
        #plt.close()

    pdf.close()

def _interpolationComparisonRutherford(atomicNumbers, energies_eV):
    logging.info("_interpolationComparisonRutherford")

    #pdf = PdfPages('Interpolation_Czyzewski90_total_atomicNumber.pdf')

    totalModels = {}
    #totalModels['BornWentzel'] = RutherfordReimerTem.totalElasticCrossSectionBornWentzel
    #totalModels['QuantumApproximation'] = RutherfordReimerTem.totalElasticCrossSectionQuantumApproximation_nm2
    totalModels['RSRutherford'] = RutherfordReimerTem.totalRelativisticScreenedElasticCrossSectionHenocMaurice_nm2

    for totalModelName in totalModels:
        logging.info("Total model: %s", totalModelName)
        totalModel = totalModels[totalModelName]

        for atomicNumber in atomicNumbers:
            logging.info("Atomic number: %s", atomicNumber)
            x = np.array(energies_eV)
            y = totalModel(atomicNumber, x)

            interpolationModels = ['linear', 'cubic']
            #interpolationModels = ['cubic']
            interpolations = {}
            for interpolationModel in interpolationModels:
                interpolations[interpolationModel] = interpolate.interp1d(x, y, kind=interpolationModel)

    #        interpolations['univariateSpline'] = interpolate.UnivariateSpline(x, y)
    #        logging.info("univariateSpline residual: %f", interpolations['univariateSpline'].get_residual())
    #        interpolationModels.append('univariateSpline')

            interpolations['interpolatedUnivariateSpline'] = interpolate.InterpolatedUnivariateSpline(x, y)
            #logging.info("interpolatedUnivariateSpline residual: %f", interpolations['interpolatedUnivariateSpline'].get_residual())
            interpolationModels.append('interpolatedUnivariateSpline')

    #        interpolations['LSQUnivariateSpline'] = interpolate.LSQUnivariateSpline(x, y, t=x[1:-1:2])
    #        logging.info("LSQUnivariateSpline residual: %f", interpolations['LSQUnivariateSpline'].get_residual())
    #        interpolationModels.append('LSQUnivariateSpline')

#            tck = interpolate.splrep(x, y)
#            def interpolateSpline(x):
#                return interpolate.splev(x, tck, der=0)
#            interpolations['spline'] = interpolateSpline
#            interpolationModels.append('spline')

            start = math.log10(energies_eV[0])
            stop = math.log10(energies_eV[-1]-1.0)
            numberPoint = 10000
            newX = np.logspace(start, stop, numberPoint)

            fig = plt.figure()
            fig.suptitle("%s %i" % (totalModelName, atomicNumber))

            axes = fig.add_subplot(211)
            left = fig.subplotpars.left
            bottom = 0.3 + fig.subplotpars.hspace
            bottom = 0.3 + fig.subplotpars.hspace/2.0
            width = fig.subplotpars.right - fig.subplotpars.left
            height = fig.subplotpars.top - bottom
            box = [left, bottom, width, height]
            axes.set_position(box)
            for interpolationModel in interpolationModels:
                newY = interpolations[interpolationModel](newX)
                plt.loglog(newX, newY, label=str(interpolationModel))

            plt.loglog(x, y, 'o', label='Raw data')
            refY = totalModel(atomicNumber, newX)
            plt.loglog(newX, refY, 'k', label='Data')

            plt.text(0.9, 0.9, atomicNumber,
                             horizontalalignment='center',
                             verticalalignment='center',
                             transform = axes.transAxes)

            plt.legend(loc='best')
            plt.ylabel(r"$\sigma_{el}$ (nm$^{2}$)")

            axes = fig.add_subplot(212, sharex=axes)
            bottom = fig.subplotpars.bottom
            width = fig.subplotpars.right - fig.subplotpars.left
            height = 0.3 - bottom
            box = [left, bottom, width, height]
            axes.set_position(box)

            refY = totalModel(atomicNumber, newX)
            for interpolationModel in interpolationModels:
                error = (refY - interpolations[interpolationModel](newX))/refY
                plt.semilogx(newX, error, label=str(interpolationModel))
                chi2 = sum((refY - interpolations[interpolationModel](newX))**2)
                logging.info("%s sum error: %f", interpolationModel, chi2)

            plt.xlabel("Energy (eV)")
            plt.ylabel("Error relative")

            #plt.ylim(ymin=-1, ymax=1)
            #pdf.savefig()
            #plt.close()

    #pdf.close()

def run():
    logging.getLogger().setLevel(logging.DEBUG)

    configurationFilepath = "eecs.cfg"
    mottTabulated = MottTabulatedFactory.getMottTabulated(configurationFilepath=configurationFilepath)

    atomicNumbers = mottTabulated.getAtomicNumbers()
    logging.debug("Atomic numbers read: %s", atomicNumbers)

    atomicNumbers = mottTabulated.getAvailableAtomicNumbers()
    logging.debug("Atomic numbers available: %s", atomicNumbers)

    energies_eV = mottTabulated.getEnergiesGrid_eV()
    logging.debug("Energy available: %s", energies_eV)

    #_multipagesMultifiguresAtomicNumbers(atomicNumbers, energies_eV, mottTabulated)
    #_multipagesMultifiguresEnergies(atomicNumbers, energies_eV, mottTabulated)

    atomicNumbers = [6, 12, 33, 55, 79, 85, 92]
    atomicNumbers = [79]
    _interpolationComparison(atomicNumbers, energies_eV, mottTabulated)
    #_interpolationComparisonRutherford(atomicNumbers, energies_eV)

    plt.show()

if __name__ == '__main__':    #pragma: no cover
    import pyHendrixDemersTools.Runner as Runner
    Runner.Runner().run(runFunction=run)
