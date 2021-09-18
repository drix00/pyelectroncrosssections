#!/usr/bin/env python
"""
.. py:currentmodule:: eecs.models.test_elsepa_casino
   :synopsis: Tests for the module :py:mod:`elsepa_casino`

.. moduleauthor:: Hendrix Demers <hendrix.demers@mail.mcgill.ca>

Tests for the module :py:mod:`eecs.models.elsepa_casino`.
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
import eecs.models.elsepa_casino
from eecs import current_module_path

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

        file_name = "../../testdata/ELSEPA_AbsorptionCorrection_LinearInterpolationTabulation_0.1.zip"
        zip_filepath = current_module_path(__file__, file_name)
        if not os.path.isfile(zip_filepath):
            self.skipTest("No file: {}".format(zip_filepath))

        self.cross_section = eecs.Models.elsepa_casino.ElsepaCasino(zip_filepath)

    def tearDown(self):
        """
        Teardown method.
        """

        unittest.TestCase.tearDown(self)

    def testSkeleton(self):
        """
        First test to check if the testcase is working with the testing framework.
        """

        # self.fail("Test if the testcase is working.")
        self.assert_(True)

    def test_zip_filepath(self):
        """
        First test to check if the testcase is working with the testing framework.
        """

        file_name = "../../testdata/ELSEPA_AbsorptionCorrection_LinearInterpolationTabulation_0.1.zip"
        zip_filepath = current_module_path(__file__, file_name)
        cross_section = eecs.Models.elsepa_casino.ElsepaCasino(zip_filepath)
        self.assertEqual(zip_filepath, cross_section.zipfilepath)

        self.assertTrue(os.path.isfile(cross_section.zipfilepath))

        # self.fail("Test if the testcase is working.")
        self.assert_(True)

    def test_total_nm2(self):
        """
        Tests for method :py:meth:`total_nm2`.
        """

        total_ref_nm2 = 0.00120146

        total_nm2 = self.cross_section.total_nm2(6, 5.0e3)

        self.assertAlmostEqual(total_ref_nm2, total_nm2)

        # self.fail("Test if the testcase is working.")

    def test_angle_deg(self):
        """
        Tests for method :py:meth:`angle_deg`.
        """

        angle_ref_deg = 4.345130434680232

        angle_deg = self.cross_section.angle_deg(6, 5.0e3, 0.5)

        self.assertAlmostEqual(angle_ref_deg, angle_deg)

        # self.fail("Test if the testcase is working.")

    def test_partial_nm2_sr(self):
        """
        Tests for method :py:meth:`partial_nm2_sr`.
        """

        partial_ref_nm2_sr = 0.017259500000000004

        partial_nm2_sr = self.cross_section.partial_nm2_sr(6, 5.0e3, 4.35)

        self.assertAlmostEqual(partial_ref_nm2_sr, partial_nm2_sr)

        # self.fail("Test if the testcase is working.")

    def test_read_total_data(self):
        """
        Tests for method :py:meth:`read_total_data`.
        """

        atomic_number = 6
        self.assertNotIn(atomic_number, self.cross_section.total_functions)
        self.cross_section.read_total_data(atomic_number)
        self.assertIn(atomic_number, self.cross_section.total_functions)

        # self.fail("Test if the testcase is working.")

    def test_read_angle_data(self):
        """
        Tests for method :py:meth:`read_angle_data`.
        """

        atomic_number = 6
        self.assertNotIn(atomic_number, self.cross_section.angle_functions)
        self.cross_section.read_angle_data(atomic_number)
        self.assertIn(atomic_number, self.cross_section.angle_functions)

        # self.fail("Test if the testcase is working.")

    def test_read_partial_data(self):
        """
        Tests for method :py:meth:`read_partial_data`.
        """

        atomic_number = 6
        self.assertNotIn(atomic_number, self.cross_section.partial_functions)
        self.cross_section.read_partial_data(atomic_number)
        self.assertIn(atomic_number, self.cross_section.partial_functions)

        # self.fail("Test if the testcase is working.")
