#!/usr/bin/env python
"""Regression testing for the project."""

# Script information for the file.
__author__ = "Hendrix Demers (hendrix.demers@mail.mcgill.ca)"
__version__ = ""
__date__ = ""
__copyright__ = "Copyright (c) 2007 Hendrix Demers"
__license__ = ""

# Subversion informations for the file.
__svnRevision__ = "$Revision: 2180 $"
__svnDate__ = "$Date: 2011-01-27 23:14:22 -0500 (Thu, 27 Jan 2011) $"
__svnId__ = "$Id: tests.py 2180 2011-01-28 04:14:22Z hdemers $"

# Standard library modules.

# Third party modules.

# Local modules.
import pyHendrixDemersTools.Testings as Testings

# Globals and constants variables.

if __name__ == "__main__":
    Testings.runTestSuiteWithCoverage(packageName=__file__)
