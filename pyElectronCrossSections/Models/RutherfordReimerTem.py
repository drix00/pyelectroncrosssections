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

# Third party modules.
import matplotlib.pyplot as plt
import scipy.constants as constants
import numpy as np

# Local modules.

# Globals and constants variables.

# TODO: Born Wentzel does not work, dimension miss match.
def totalElasticCrossSectionBornWentzel(atomicNumber, energy_eV):
    """
    From Reimer TEM book's page 152.

    .. math:

    \sigma_{el} &= \frac{Z^{2}R^{2}\lambda^{2}\left(1 + \frac{E}{E_{0}}\right)^{2}}{\pi a^{2}_{H}} \\
    &= \frac{h^{2} Z^{4/3}}{\pi E_{0}^2 \beta^{2}}

    """
    pi = constants.pi
    #h = constants.Planck
    h_eVs = 4.13566733e-15

    E0_eV = 511.0e3
    Z = atomicNumber
    beta2 = _computeBeta2(energy_eV)

    nominator = h_eVs**2 * np.power(Z, 4.0/3.0)
    denominator = pi * E0_eV**2 * beta2

    total_nm2 = nominator/denominator
    return total_nm2

def totalElasticCrossSectionQuantumApproximation_nm2(atomicNumber, energy_eV):
    """
    From Reimer TEM book's page 153.

    .. math:

    \sigma_{el} &= \frac{1.5\time 10^{-6}}{\beta^{2}} Z^{3/2} \left(1 - 0.23 \frac{Z}{137 \beta}\right).

    for .. math:`\frac{Z}{137 \beta} < 1.2`.

    """
    Z = atomicNumber
    beta = _computeBeta(energy_eV)

    limitFactor = Z/(137.0*beta)
#    if limitFactor > 1.2:
#        logging.warning("Model used out of his validity range: %i, %f", atomicNumber, energy_eV)

    factorA = (1.5e-6) * np.power(Z, 3.0/2.0) / beta**2
    factorB = 1.0 - 0.23*limitFactor

    total_nm2 = factorA*factorB
    return total_nm2

def _computeBeta(energy_eV):
    beta2 = _computeBeta2(energy_eV)
    beta = np.sqrt(beta2)
    return beta

def _computeBeta2(energy_eV):
    E0_eV = 511.0e3
    temp = (1.0 / (1.0 + energy_eV / E0_eV))
    beta2 = 1.0 - temp * temp
    return beta2

def differentialRutherfordSmallAngle_nm2_str(atomicNumber, energy_eV, angle_rad):
    """
    From Reimer TEM book's page 151.

    .. math:

    \frac{d\sigma_{el}}{d\Omega} = \frac{4 Z^{2} R^{4} \left(\left[1 + \frac{E}{E_{0}}\right)^{2} }{a_{H}^{2}} \frac{1}{1 + \left(\frac{\theta}{\theta_{0}}\right)^{2}\right]^{2}}

    with .. math:`\theta_{0} = \frac{\lambda}{2 \pi R}` and .. math:`R = a_{H}Z^{-1/3}`.
    """
    Z = atomicNumber
    R_nm = _computeR_nm(atomicNumber)
    aH_nm = getBohrRadius_nm()
    E0_eV = 511.0e3
    theta0 = _computeCharateristicAngle(atomicNumber, energy_eV)
    E_eV = energy_eV

    nominator = 4.0 * Z**2 * R_nm**4 * (1.0 + E_eV/E0_eV)**2
    denominator = aH_nm**2
    factorA = nominator/denominator

    denominator = (1.0 + (angle_rad/theta0)**2)**2
    factorB = 1.0/denominator

    differential_nm2_str = factorA*factorB

    return differential_nm2_str

def _computeR_nm(atomicNumber):
    aH_nm = getBohrRadius_nm()
    R_nm = aH_nm * np.power(atomicNumber, -1.0/3.0)
    return R_nm

def _computeCharateristicAngle(atomicNumber, energy_eV):
    R_nm = _computeR_nm(atomicNumber)
    lambda_nm = _computeLambda_nm(energy_eV)

    theta0_rad = lambda_nm/(2.0*np.pi*R_nm)
    return theta0_rad

def _computeLambda_nm(energy_eV):
    h_eVs = 4.13566733e-15
    c_m_s = constants.c
    E0_eV = 511.0e3
    E_eV = energy_eV

    nominator = h_eVs * c_m_s
    sqrtArg = 2.0*E_eV*E0_eV + E_eV**2
    denominator = np.sqrt(sqrtArg)

    lambda_m = nominator/denominator

    lambda_nm = lambda_m*1.0e9
    return lambda_nm

def getBohrRadius_nm():
    return 0.0529

def totalRelativisticScreenedElasticCrossSectionHenocMaurice_nm2(atomicNumber, energy_eV):
    """
    From Joy et al. "Principles of Analytical Electron Microscopy" book's page 5.

    .. math:

    \sigma_{el} &= \frac{Z^{2}\lambda_{R}^{4}}{16 \pi^{3}a_{0}^{2}} \frac{1}{\delta\left(\delta + 1\right)}
    &= \frac{7.20\times 10^{13} Z^{2}\lambda_{R}^{4}}{\delta\left(\delta + 1\right)}

    where .. math:`\delta = \frac{\theta_{0}^{2}}{4}`.

    """
    factor_1_cm2 = _computeTotalFactor_1_cm2()
    Z2 = atomicNumber**2
    lambda_1_cm = _computeLambda_cm(energy_eV)
    delta = _computeDelta(atomicNumber, energy_eV)

    nominator = factor_1_cm2 * Z2 * lambda_1_cm**4
    denominator = delta*(delta + 1.0)

    total_cm2 = nominator/denominator
    total_nm2 = total_cm2*1.0e14

    return total_nm2

def _computeTotalFactor_1_cm2():
    """
    From Joy et al. "Principles of Analytical Electron Microscopy" book's page 4.

    .. math:

    a_{0} &= \frac{\epsilon_{0} h^{2}}{\pi m_{0} e^{2}}
    &= 5.29\times10^{-9} (cm).

    """
    a0_cm = 5.29e-9
    factor = 1.0/(16.0 * np.pi**3 * a0_cm**2)
    return factor

def _computeLambda_cm(energy_eV):
    energy_keV = energy_eV*1.0e-3

    nominator = 3.87e-9
    denominator = np.sqrt(energy_keV)*np.sqrt(1.0 + 9.79e-4 * energy_keV)

    lambda_cm = nominator/denominator

    return lambda_cm

def _computeDelta(atomicNumber, energy_eV):
    energy_keV = energy_eV*1.0e-3
    delta = 3.4e-3 * np.power(atomicNumber, 2.0/3.0) / energy_keV

    return delta

def totalCrossSectionGreaterThanAngle_nm2(atomicNumber, energy_eV, angle_rad):
    """
    From Reimer TEM book's page 151.

    .. math:

    \frac{\sigma_{el}}{\alpha_{0}} = \frac{Z^{2} R^{2} \lambda^{2} \left(1 + \frac{E}{E_{0}}\right)^{2} }{\pi a_{H}^{2}} \right)^{2} \frac{1}{1 + \left(\frac{\alpha_{0}}{\theta_{0}}\right)^{2}}

    with .. math:`\theta_{0} = \frac{\lambda}{2 \pi R}` and .. math:`R = a_{H}Z^{-1/3}`.
    """
    Z = atomicNumber
    R_nm = _computeR_nm(atomicNumber)
    aH_nm = getBohrRadius_nm()
    E0_eV = 511.0e3
    theta0 = _computeCharateristicAngle(atomicNumber, energy_eV)
    E_eV = energy_eV
    lambda_nm = _computeLambda_nm(energy_eV)

    nominator = Z**2 * R_nm**2 * lambda_nm**2 * (1.0 + E_eV/E0_eV)**2
    denominator = np.pi * aH_nm**2
    factorA = nominator/denominator

    denominator = (1.0 + (angle_rad/theta0)**2)
    factorB = 1.0/denominator

    total_nm2 = factorA*factorB

    return total_nm2

def figureTotal():
    plt.figure()

    energy_eV = 200.0e3
    atomicNumbers = range(1, 100)

    eccsModels = {"Quantum approximation": totalElasticCrossSectionQuantumApproximation_nm2,
                                "Born Wentzel": totalElasticCrossSectionBornWentzel}

    for label in eccsModels:
        eccs = eccsModels[label]
        totals_nm2 = [eccs(Z, energy_eV) for Z in atomicNumbers]
        plt.loglog(atomicNumbers, totals_nm2, label=label)

    plt.xlabel("Atomic number")
    plt.ylabel(r"$\sigma_{el}$ (nm$^{2}$)")
    plt.subplots_adjust(left=0.16)
    plt.legend(loc='best')

    plt.show()

def _figureDifferential():
    plt.figure()
    atomicNumber = 18
    energy_eV = 25.0e3

    angles_rad = np.logspace(np.log10(1.0e-4), np.log10(np.pi), num=1000)
    diff_str = []
    aH2 = getBohrRadius_nm()**2

    for angle_rad in angles_rad:
        diff_str.append(differentialRutherfordSmallAngle_nm2_str(atomicNumber, energy_eV, angle_rad)/aH2)

    plt.loglog(angles_rad, diff_str)
    #plt.xlim(xmin=1.0e-4)
    plt.show()

def _figureTotal():
    atomicNumbers = [6, 29, 79]
    models = {}
    #models['BornWentzel'] = totalElasticCrossSectionBornWentzel
    models['QuantumApproximation'] = totalElasticCrossSectionQuantumApproximation_nm2
    models['RSRutherford'] = totalRelativisticScreenedElasticCrossSectionHenocMaurice_nm2

    start = np.log10(1.0)
    stop = np.log10(5.0e6)
    numberPoint = 10000
    energies_eV = np.logspace(start, stop, numberPoint)

    for atomicNumber in atomicNumbers:
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
    #_figureDifferential()
    _figureTotal()

if __name__ == '__main__':    #pragma: no cover
    import DrixUtilities.Runner as Runner
    Runner.Runner().run(runFunction=run)
