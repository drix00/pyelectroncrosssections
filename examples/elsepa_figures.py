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

# Third party modules.
import matplotlib.pyplot as plt

# Local modules.
import DatabasesTools.Casino3.ElsepaBinaryFile as ElsepaBinaryFile
from eecs.models.elsepa_casino import ElsepaCasino
from eecs import current_module_path

# Globals and constants variables.

def run():
    path = r"D:\work\data\Casino3\ELDB"
    zipfilepath = current_module_path(__file__, "../testdata/ELSEPA_AbsorptionCorrection_LinearInterpolationTabulation_0.1.zip")

    atomicNumbers = range(1, 99+1)
    totalCS_nm2 = {}
    totalCS2_nm2 = {}

    elsepa_casino = ElsepaCasino(zipfilepath)
    for atomicNumber in atomicNumbers:
        totalCS_nm2[atomicNumber] = {}
        totalCS2_nm2[atomicNumber] = {}

        filename = "EL%i.els" % (atomicNumber)
        filepath = os.path.join(path, filename)
        elsepaFile = ElsepaBinaryFile.ElsepaBinaryFile(filepath)

        for data in elsepaFile._els_cs_info_list:
            energy_keV = data._energy_keV
            total = data._total_cs_nm2

            totalCS_nm2[atomicNumber][energy_keV] = total

            energy_eV = energy_keV * 1.0e3
            total_nm2 = elsepa_casino.total_nm2(atomicNumber, energy_eV)
            totalCS2_nm2[atomicNumber][energy_keV] = total_nm2

    energy_keV = 20.0
    y = []
    y2 = []
    for atomicNumber in atomicNumbers:
        y.append(totalCS_nm2[atomicNumber][energy_keV])
        y2.append(totalCS2_nm2[atomicNumber][energy_keV])

    plt.figure()
    plt.title(r"$E_{0} = %.1f$" % (energy_keV))
    plt.plot(atomicNumbers, y)
    plt.plot(atomicNumbers, y2, '.')

    plt.xlabel('Z')
    plt.ylabel(r"$\sigma_{T}$ (nm$^{2}$)")

    energies_keV = sorted(totalCS_nm2[atomicNumber].keys())
    atomicNumber = 29
    y = []
    y2 = []
    for energy_keV in energies_keV:
        y.append(totalCS_nm2[atomicNumber][energy_keV])
        y2.append(totalCS2_nm2[atomicNumber][energy_keV])

    plt.figure()
    plt.title("Z = %i" % (atomicNumber))
    plt.loglog(energies_keV, y)
    plt.loglog(energies_keV, y2, '.')
    plt.xlabel('energy (keV)')
    plt.ylabel(r"$\sigma_{T}$ (nm$^{2}$)")

    plt.show()


if __name__ == '__main__':  # pragma: no cover
    run()
