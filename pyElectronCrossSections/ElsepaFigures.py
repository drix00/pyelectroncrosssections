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

# Globals and constants variables.

def run():
    path = r"C:\hdemers\data\Casino3\ELDB"

    atomicNumbers = range(1, 99+1)
    totalCS_nm2 = {}
    for atomicNumber in atomicNumbers:
        totalCS_nm2[atomicNumber] = {}
        filename = "EL%i.els" % (atomicNumber)
        filepath = os.path.join(path, filename)
        elsepaFile = ElsepaBinaryFile.ElsepaBinaryFile(filepath)

        for data in elsepaFile._elsCSInfoList:
            energy_keV = data._energy_keV
            total = data._totalCS_nm2

            totalCS_nm2[atomicNumber][energy_keV] = total


    energy_keV = 200.0
    y = []
    for atomicNumber in atomicNumbers:
        y.append(totalCS_nm2[atomicNumber][energy_keV])

    plt.figure()

    plt.plot(atomicNumbers, y)

    plt.show()

if __name__ == '__main__':    #pragma: no cover
    import pyHendrixDemersTools.Runner as Runner
    Runner.Runner().run(runFunction=run)
