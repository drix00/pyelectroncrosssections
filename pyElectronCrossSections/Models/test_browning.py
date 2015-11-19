#!/usr/bin/env python
"""
.. py:currentmodule:: pyElectronCrossSections.Models.test_browning
   :synopsis: Tests for the module :py:mod:`Browning`

.. moduleauthor:: Hendrix Demers <hendrix.demers@mail.mcgill.ca>

Tests for the module :py:mod:`Browning`.
"""

# Script information for the file.
__author__ = "Hendrix Demers (hendrix.demers@mail.mcgill.ca)"
__version__ = "0.1"
__date__ = "Nov 19, 2015"
__copyright__ = "Copyright (c) 2015 Hendrix Demers"
__license__ = "GPL 3"

# Standard library modules.
import unittest

# Third party modules.

# Local modules.

# Project modules
import pyElectronCrossSections.Models.Browning as Browning

# Globals and constants variables.

class TestBrowning(unittest.TestCase):
    """
    TestCase class for the module `Browning`.
    """

    def setUp(self):
        """
        Setup method.
        """

        unittest.TestCase.setUp(self)

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

    def test_totalElasticCrossSectionBrowning1994_cm2(self):
        """
        Tests for method :py:meth:`totalElasticCrossSectionBrowning1994_cm2`.
        """

        energy_keV = 100.0

        total_css_A2_ref = {6: 0.006243,
                            13: 0.02260,
                            26: 0.0676705532253224,
                            47: 0.15470732661851455,
                            79: 0.2734892890031806}

        for atomic_number in total_css_A2_ref:
            total_cs_A2_ref = total_css_A2_ref[atomic_number]
            total_cs_cm2 = Browning.totalElasticCrossSectionBrowning1994_cm2(atomic_number, energy_keV)
            total_cs_A2 = total_cs_cm2*1.0e16
            self.assertAlmostEqual(total_cs_A2_ref, total_cs_A2, 6)

        #self.fail("Test if the testcase is working.")

if __name__ == '__main__':  #pragma: no cover
    import nose
    import sys
    argv = sys.argv
    argv.append("--with-coverage")
    argv.append("--cover-package=pyElectronCrossSections.Models.Browning")
    nose.runmodule(argv=argv)
