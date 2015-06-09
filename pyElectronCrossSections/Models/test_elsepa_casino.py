#!/usr/bin/env python
"""
.. py:currentmodule:: pyElectronCrossSections.Models.test_elsepa_casino
   :synopsis: Tests for the module :py:mod:`elsepa_casino`

.. moduleauthor:: Hendrix Demers <hendrix.demers@mail.mcgill.ca>

Tests for the module :py:mod:`pyElectronCrossSections.Models.elsepa_casino`.
"""

# Script information for the file.
__author__ = "Hendrix Demers (hendrix.demers@mail.mcgill.ca)"
__version__ = "0.1"
__date__ = "Jun 1, 2015"
__copyright__ = "Copyright (c) 2015 Hendrix Demers"
__license__ = "GPL 3"

# Standard library modules.
import unittest
import os.path

# Third party modules.

# Local modules.

# Project modules
import pyElectronCrossSections.Models.elsepa_casino
from pyElectronCrossSections import current_module_path

# Globals and constants variables.

class TestElsepaCasino(unittest.TestCase):
    """
    TestCase class for the module `ElsepaCasino`.
    """

    def setUp(self):
        """
        Setup method.
        """

        unittest.TestCase.setUp(self)

        zipfilepath = current_module_path(__file__, "../../testdata/ELSEPA_AbsorptionCorrection_LinearInterpolationTabulation_0.1.zip")
        self.cross_section = pyElectronCrossSections.Models.elsepa_casino.ElsepaCasino(zipfilepath)

    def tearDown(self):
        """
        Teardown method.
        """

        unittest.TestCase.tearDown(self)

    def testSkeleton(self):
        """
        First test to check if the testcase is working with the testing framework.
        """

        #self.fail("Test if the testcase is working.")
        self.assert_(True)

    def test_zipfilepath(self):
        """
        First test to check if the testcase is working with the testing framework.
        """

        zipfilepath = current_module_path(__file__, "../../testdata/ELSEPA_AbsorptionCorrection_LinearInterpolationTabulation_0.1.zip")
        cross_section = pyElectronCrossSections.Models.elsepa_casino.ElsepaCasino(zipfilepath)
        self.assertEqual(zipfilepath, cross_section.zipfilepath)

        self.assertTrue(os.path.isfile(cross_section.zipfilepath))

        #self.fail("Test if the testcase is working.")
        self.assert_(True)

    def test_total_nm2(self):
        """
        Tests for method :py:meth:`total_nm2`.
        """

        totalRef_nm2 = 0.00120146

        total_nm2 = self.cross_section.total_nm2(6, 5.0e3)

        self.assertAlmostEqual(totalRef_nm2, total_nm2)

        #self.fail("Test if the testcase is working.")

    def test_angle_deg(self):
        """
        Tests for method :py:meth:`angle_deg`.
        """

        angleRef_deg = 4.345130434680232

        angle_deg = self.cross_section.angle_deg(6, 5.0e3, 0.5)

        self.assertAlmostEqual(angleRef_deg, angle_deg)

        #self.fail("Test if the testcase is working.")

    def test_partial_nm2_sr(self):
        """
        Tests for method :py:meth:`partial_nm2_sr`.
        """

        partialRef_nm2_sr = 0.017259500000000004

        partial_nm2_sr = self.cross_section.partial_nm2_sr(6, 5.0e3, 4.35)

        self.assertAlmostEqual(partialRef_nm2_sr, partial_nm2_sr)

        #self.fail("Test if the testcase is working.")

    def test_read_total_data(self):
        """
        Tests for method :py:meth:`read_total_data`.
        """

        atomic_number = 6
        self.assertNotIn(atomic_number, self.cross_section.total_functions)
        self.cross_section.read_total_data(atomic_number)
        self.assertIn(atomic_number, self.cross_section.total_functions)

        #self.fail("Test if the testcase is working.")

    def test_read_angle_data(self):
        """
        Tests for method :py:meth:`read_angle_data`.
        """

        atomic_number = 6
        self.assertNotIn(atomic_number, self.cross_section.angle_functions)
        self.cross_section.read_angle_data(atomic_number)
        self.assertIn(atomic_number, self.cross_section.angle_functions)

        #self.fail("Test if the testcase is working.")

    def test_read_partial_data(self):
        """
        Tests for method :py:meth:`read_partial_data`.
        """

        atomic_number = 6
        self.assertNotIn(atomic_number, self.cross_section.partial_functions)
        self.cross_section.read_partial_data(atomic_number)
        self.assertIn(atomic_number, self.cross_section.partial_functions)

        #self.fail("Test if the testcase is working.")

if __name__ == '__main__':  #pragma: no cover
    import nose
    import sys
    argv = sys.argv
    #argv.append("--with-coverage")
    argv.append("--cover-package=pyElectronCrossSections.Models.elsepa_casino")
    nose.runmodule(argv=argv)
