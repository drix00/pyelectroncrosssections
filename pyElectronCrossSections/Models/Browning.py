#!/usr/bin/env python
"""
.. py:currentmodule:: Models.Browning
.. moduleauthor:: Hendrix Demers <hendrix.demers@mail.mcgill.ca>

Cross section models from Browning.
"""

# Script information for the file.
__author__ = "Hendrix Demers (hendrix.demers@mail.mcgill.ca)"
__version__ = ""
__date__ = ""
__copyright__ = "Copyright (c) 2012 Hendrix Demers"
__license__ = ""

# Standard library modules.
import math
import os.path
import csv

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

def totalElasticCrossSectionRutherford_cm2(atomicNumber, electronEnergy_keV):
    """
    From browning1991a ref 1: J
    D. C. Joy, Proceedings of the 9th European Congress of Electron Microscopy,
    edited by P. J. Goodhew and H. G. Dickinson (Institute of Physits,
    Bristol, UK, 1988), p. 22.

    Without relativistic effects correction.
    Screened Rutherford cross section in cm2

    """

    factor = 5.21e-21
    termA = atomicNumber*atomicNumber/(electronEnergy_keV*electronEnergy_keV)

    alpha = computeScreeningParameter(atomicNumber, electronEnergy_keV)
    termB = math.pi/(alpha * (1.0 + alpha))

    crossSection_cm2 = factor * termA * termB

    return crossSection_cm2

def computeScreeningParameter(atomicNumber, electronEnergy_keV):
    factor = 3.4e-3
    termA = math.pow(atomicNumber, 0.67) / electronEnergy_keV

    alpha = factor * termA
    return alpha

def averageScatteringAngleRutherford_deg(atomicNumber, electronEnergy_keV):
    alpha = computeScreeningParameter(atomicNumber, electronEnergy_keV)

    averageAngle_rad = math.pi * math.sqrt(alpha) * math.sqrt(1.0 + alpha) - math.pi * alpha
    averageAngle_deg = math.degrees(averageAngle_rad)
    return averageAngle_deg

def averageScatteringAngleRutherfordDecreasedScreening_deg(atomicNumber, electronEnergy_keV):
    alpha = computeDecreasedScreeningParameter(atomicNumber, electronEnergy_keV)

    averageAngle_rad = math.pi * math.sqrt(alpha) * math.sqrt(1.0 + alpha) - math.pi * alpha
    averageAngle_deg = math.degrees(averageAngle_rad)
    return averageAngle_deg

def computeDecreasedScreeningParameter(atomicNumber, electronEnergy_keV):
    alphaRutherford = computeScreeningParameter(atomicNumber, electronEnergy_keV)
    alpha = (0.6 - 0.0035 * electronEnergy_keV) * alphaRutherford
    return alpha

def totalElasticCrossSectionBrowning1991a_cm2(atomicNumber, electronEnergy_keV):
    """
    From browning1991a
    Valid in the range 1 to 100 keV for all elments.
    """
    Z = atomicNumber
    E = electronEnergy_keV
    factor = 4.7e-18
    nominator = math.pow(Z, 1.33) + 0.032*Z*Z
    denominator = E + 0.0155 * math.pow(Z, 1.33) * math.pow(E, 0.5)
    termA = nominator/denominator

    u = computeFactorU(atomicNumber, electronEnergy_keV)
    denominator = 1.0 - 0.02 * math.pow(Z, 0.5) * math.exp(-u*u)
    termB = 1.0/denominator

    crossSection_cm2 = factor * termA * termB

    return crossSection_cm2

def computeFactorU(atomicNumber, electronEnergy_keV):
    u = math.log10(8.0) * electronEnergy_keV * math.pow(atomicNumber, -1.33)
    return u

def differentialCrossSectionBrowning1991_cm2_sr(atomicNumber, electronEnergy_keV, theta_rad):
    factor = 5.21e-21
    Z = atomicNumber
    E = electronEnergy_keV
    termA = Z * Z / (E * E)

    alpha = computeScreeningParameterBrowning1991(atomicNumber, electronEnergy_keV)
    denominator = 1.0 - math.cos(theta_rad) - alpha
    termB = 1.0 / denominator

    nominator = alpha * (alpha + 1.0)
    denominator = 4.2 * math.pow(E, 1.1)
    termC = nominator / denominator

    differentialCrossSection_cm2_sr = factor * termA * (termB + termC)

    return differentialCrossSection_cm2_sr

def computeScreeningParameterBrowning1991(atomicNumber, electronEnergy_keV):
    alpha = 5.5e-4 * math.pow(atomicNumber, 0.67) / electronEnergy_keV
    return alpha

def ratioBrowning1994(atomicNumber, electronEnergy_keV):
    Z = atomicNumber
    E = electronEnergy_keV
    termA = 300.0 * math.pow(E, 1.0 - Z/2000.0) / Z

    termB = math.pow(Z, 3.0) / (3.0e5 * E)

    ratio = termA + termB
    return ratio

def ratioBrowning1994MCXRay(atomicNumber, electronEnergy_keV):
    Z = atomicNumber
    E = electronEnergy_keV
    termA = 300.0 * math.pow(E, 1.0 - Z/2000.0) / Z

    termB = math.pow(Z, 3.0) / 3.0e5 * E

    ratio = termA + termB
    return ratio

def computePolarAngleTwoRandomNumbers_rad(atomicNumber, electronEnergy_keV):
    ratio = ratioBrowning1994(atomicNumber, electronEnergy_keV)

    r1 = np.random.random()
    r2 = np.random.random()

    if r1 <= ratio / (1.0 + ratio):
        alpha = 7.0e-3 / electronEnergy_keV
        cosTheta = 1.0 - 2.0 * alpha * r2 / (1.0 + alpha - r2)
        theta_rad = np.arccos(cosTheta)
    else:
        cosTheta = 1.0 - 2.0 * r2
        theta_rad = np.arccos(cosTheta)

    return theta_rad

def computePolarAngleOneRandomNumber_rad(atomicNumber, electronEnergy_keV):
    ratio = ratioBrowning1994(atomicNumber, electronEnergy_keV)

    r1 = np.random.random()

    if r1 <= ratio / (1.0 + ratio):
        alpha = 7.0e-3 / electronEnergy_keV
        cosTheta = 1.0 - 2.0 * alpha * r1 / (1.0 + alpha - r1)
        theta_rad = np.arccos(cosTheta)
    else:
        cosTheta = 1.0 - 2.0 * r1
        theta_rad = np.arccos(cosTheta)

    return theta_rad

def plotFigure2Browning1991a():
    energies_keV = np.arange(0.1, 100.0, 0.1)
    elements = ['U', 'Mo', 'Al', 'C', 'H']
    atomicNumbers = {'U': 92, 'Mo': 42, 'Al': 13, 'C': 6, 'H': 1}

    colors = Colors.BaseColors(len(elements))
    plt.figure()

    for element in elements:
        Z = atomicNumbers[element]
        color = colors.next()
        crossSections_cm2 = [totalElasticCrossSectionRutherford_cm2(Z, E_keV) for E_keV in energies_keV]
        crossSections_pm2 = np.array(crossSections_cm2)*1.0e16
        plt.loglog(energies_keV, crossSections_pm2, '--', color=color)
        crossSections_cm2 = [totalElasticCrossSectionBrowning1991a_cm2(Z, E_keV) for E_keV in energies_keV]
        crossSections_pm2 = np.array(crossSections_cm2)*1.0e16
        plt.loglog(energies_keV, crossSections_pm2, color=color, label=element)

    plt.xlabel("Electron Energy (keV)")
    plt.ylabel(r"Total Elastic Cross Section ($10^{-16}$cm$^{2}$)")
    plt.legend(loc='best')

    plt.ylim((0.001, 10.0))

def plotFigure3Browning1991():
    energies_keV = np.arange(0.1, 100.0, 0.1)
    elements = ['Au']
    atomicNumbers = {'Au': 79}

    colors = Colors.BaseColors(len(elements))
    plt.figure()

    for element in elements:
        Z = atomicNumbers[element]
        color = colors.next()
        averageAngles_deg = [averageScatteringAngleRutherford_deg(Z, E_keV) for E_keV in energies_keV]
        plt.plot(energies_keV, averageAngles_deg, color=color, label=element)

        color = colors.next()
        averageAngles_deg = [averageScatteringAngleRutherfordDecreasedScreening_deg(Z, E_keV) for E_keV in energies_keV]
        plt.plot(energies_keV, averageAngles_deg, color=color)

    plt.xlabel("Energy (keV)")
    plt.ylabel(r"Average Scattering Angle (Degrees)")
    plt.legend(loc='best')

    plt.ylim((0.0, 40.0))

def plotFigure4Browning1991():
    energy_keV = 1.0
    angles_deg = np.arange(0.0, 180.0, 1.0)
    angles_rad = np.arange(0.0, np.pi, 0.01)
    angles_deg = np.degrees(angles_rad)
    elements = ['Au']
    atomicNumbers = {'Au': 79}

    colors = Colors.BaseColors(len(elements))
    plt.figure()

    for element in elements:
        Z = atomicNumbers[element]
        color = colors.next()
        differentialCrossSections = [differentialCrossSectionBrowning1991_cm2_sr(Z, energy_keV, angle_rad) for angle_rad in angles_rad]
        plt.plot(angles_deg, differentialCrossSections, color=color, label=element)

    plt.xlabel("Agnle (Degrees)")
    plt.ylabel(r"Cross Section (A$^{2}$)")
    plt.legend(loc='best')

    #plt.ylim((1.0e-4, 1.0e1))

def plotComparaisonRatioMCXRayError():
    energies_keV = np.arange(0.1, 100.0, 0.1)
    elements = ['Au', 'Cu', 'C']
    atomicNumbers = {'Au': 79, 'Cu': 29, 'C': 6}

    colors = Colors.BaseColors(len(elements))
    plt.figure()

    for element in elements:
        Z = atomicNumbers[element]
        color = colors.next()
        ratios = np.array([ratioBrowning1994(Z, E_keV) for E_keV in energies_keV])
        ratios = ratios / (1.0 + ratios)
        label = "%s Browning" % (element)
        plt.semilogx(energies_keV, ratios, color=color, label=label)

        color = colors.next()
        ratios = np.array([ratioBrowning1994MCXRay(Z, E_keV) for E_keV in energies_keV])
        ratios = ratios / (1.0 + ratios)
        label = "%s MCXRay" % (element)
        plt.semilogx(energies_keV, ratios, color=color, label=label)

    plt.xlabel("Energy (keV)")
    plt.ylabel(r"Ratio / (1 + Ratio)")
    plt.legend(loc='best')

    #plt.ylim((0.0, 40.0))

def plotComparaisonPolarAngleRandomNumber():
    numberSamples = 10000000
    atomicNumber = 79
    energy_keV = 1.0

    ratio = ratioBrowning1994(atomicNumber, energy_keV)
    print("Ratio = %.4f" % (ratio))
    print("Ratio/(1 + Ratio) = %.4f" % (ratio/(1.0 + ratio)))

    polarAnglesOneRandomNumber_rad = [computePolarAngleOneRandomNumber_rad(atomicNumber, energy_keV) for _i in range(numberSamples)]

    polarAnglesTwoRandomNumbers_rad = [computePolarAngleTwoRandomNumbers_rad(atomicNumber, energy_keV) for _i in range(numberSamples)]

    plt.figure()

    plt.hist(polarAnglesTwoRandomNumbers_rad, bins=100, normed=True, label='2 RNs', histtype='step')
    plt.hist(polarAnglesOneRandomNumber_rad, bins=100, normed=True, label='1 RN', histtype='step')

    plt.xlabel(r"Scattering Angle (rad)")
    plt.ylabel(r"Probabilities")

    plt.legend(loc='best')
    filename = "ComparaisonPolarAngleRandomNumber.pdf"
    filepath = os.path.join(r"J:\hdemers\work\documents\labbooks\e-labbook\graphics\BrowningCrossSection", filename)
    plt.savefig(filepath)

def plotComputeDifferentialCrossSection():
    numberBins = 50
    numberSamples = 1000000
    atomicNumber = 79
    energy_keV = 1.0

    total_cm2 = totalElasticCrossSectionBrowning1991a_cm2(atomicNumber, energy_keV)
    total_nm2 = cm2Tonm2(total_cm2)
    total_A2 = total_nm2 * 1.0e2

    ratio = ratioBrowning1994(atomicNumber, energy_keV)
    print("Ratio = %.4f" % (ratio))
    print("Ratio/(1 + Ratio) = %.4f" % (ratio/(1.0 + ratio)))

    polarAnglesTwoRandomNumbers_rad = [computePolarAngleTwoRandomNumbers_rad(atomicNumber, energy_keV) for _i in range(numberSamples)]

    histogram, low_range, binsize, _extrapoints = scipy.stats.histogram(polarAnglesTwoRandomNumbers_rad, numbins=numberBins)

    histogram /= numberSamples*binsize
    histogram *= total_A2

    max_range = low_range + binsize*numberBins
    thetas_rad = np.arange(low_range, max_range, binsize)

    plt.figure()

    #plt.hist(polarAnglesTwoRandomNumbers_rad, bins=100, normed=True, histtype='step')
    plt.semilogy(thetas_rad, histogram)

    plt.xlabel(r"Scattering Angle (rad)")
    plt.ylabel(r"Probabilities (A2/sr)")
    plt.xlim(xmin=1.0*binsize)

    #plt.legend(loc='best')
    filename = "ComputeDifferentialCrossSection.pdf"
    filepath = os.path.join(r"J:\hdemers\work\documents\labbooks\e-labbook\graphics\BrowningCrossSection", filename)
    plt.savefig(filepath)

def readData(filepath):
    energies_eV = []
    totals_nm2 = {}
    meanThetas_rad = {}

    reader = csv.reader(open(filepath, 'rb'))

    row = reader.next()
    for row in reader:
        energy_eV = float(row[0])
        meanTheta_rad = float(row[1])
        total_nm2 = float(row[2])

        energies_eV.append(energy_eV)
        meanThetas_rad[energy_eV] = meanTheta_rad
        totals_nm2[energy_eV] = total_nm2

    del reader

    return meanThetas_rad, totals_nm2

def saveData(filepath, energies_eV, meanThetas_rad, totals_nm2):
    writer = csv.writer(open(filepath, 'wb'))

    rowHeader = ["Energy (eV)", "Mean Theta (rad)", "Total (nm2)"]
    writer.writerow(rowHeader)

    for energy_eV in energies_eV:
        row = [energy_eV, meanThetas_rad[energy_eV], totals_nm2[energy_eV]]
        writer.writerow(row)

def computeMeanThetaTotalBrowning(atomicNumber, energy_eV):
    numberBins = 50
    numberSamples = 1000000

    energy_keV = energy_eV*1.0e-3
    total_cm2 = totalElasticCrossSectionBrowning1991a_cm2(atomicNumber, energy_keV)
    total_nm2 = cm2Tonm2(total_cm2)

    polarAnglesTwoRandomNumbers_rad = [computePolarAngleTwoRandomNumbers_rad(atomicNumber, energy_keV) for _i in range(numberSamples)]

    histogram, low_range, binsize, _extrapoints = scipy.stats.histogram(polarAnglesTwoRandomNumbers_rad, numbins=numberBins)

    histogram /= numberSamples*binsize
    histogram *= total_nm2

    max_range = low_range + binsize*numberBins
    thetas_rad = np.arange(low_range, max_range, binsize)

    #partialSinThetas_nm2_sr = [2.0*np.pi*x1*np.sin(x2) for x1,x2 in zip(histogram, thetas_rad)]
    partialSinThetas_nm2_sr = [x1 for x1,x2 in zip(histogram, thetas_rad)]
    totalCalculated_nm2 = scipy.integrate.trapz(partialSinThetas_nm2_sr, thetas_rad)

    #partialSinThetas_nm2_sr = [2.0*np.pi*x1*x2*np.sin(x2) for x1,x2 in zip(histogram, thetas_rad)]
    partialSinThetas_nm2_sr = [x1*x2 for x1,x2 in zip(histogram, thetas_rad)]
    meanTheta_rad = scipy.integrate.trapz(partialSinThetas_nm2_sr, thetas_rad)
    meanTheta_rad = meanTheta_rad/totalCalculated_nm2

    return meanTheta_rad, totalCalculated_nm2

def run():
    #plotFigure2Browning1991a()

    #plotFigure3Browning1991()
    #plotFigure4Browning1991()

    #plotComparaisonRatioMCXRayError()

    #plotComparaisonPolarAngleRandomNumber()
    plotComputeDifferentialCrossSection()

    plt.show()

if __name__ == '__main__': #pragma: no cover
    import pyHendrixDemersTools.Runner as Runner
    Runner.Runner().run(runFunction=run)
