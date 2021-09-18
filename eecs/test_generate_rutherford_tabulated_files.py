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
__svnId__ = "$Id: test_generate_rutherford_tabulated_files.py 2292 2011-03-21 15:29:50Z hdemers $"

# Standard library modules.
import unittest
import logging

# Third party modules.

# Local modules.
import eecs.generate_rutherford_tabulated_files as GenerateRutherfordTabulatedFiles

# Globals and constants variables.

class TestGenerateRutherfordTabulatedFiles(unittest.TestCase):

    def setUp(self):
        unittest.TestCase.setUp(self)

    def tearDown(self):
        unittest.TestCase.tearDown(self)

    def testSkeleton(self):
        #self.fail("Test if the testcase is working.")
        self.assert_(True)

if __name__ == '__main__':    #pragma: no cover
    logging.getLogger().setLevel(logging.DEBUG)
    unittest.main()
