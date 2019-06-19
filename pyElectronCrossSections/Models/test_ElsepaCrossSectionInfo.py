#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. py:currentmodule:: pyElectronCrossSection.Models.test_ElsepaCrossSectionInfo

.. moduleauthor:: Hendrix Demers <hendrix.demers@mail.mcgill.ca>


"""

###############################################################################
# Copyright 2019 Hendrix Demers
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
###############################################################################

# Standard library modules.
import unittest
import os.path

# Third party modules.

# Local modules.
from pyElectronCrossSections.Models.ElsepaCrossSectionInfo import ElsepaCrossSectionInfo
from pyElectronCrossSections import current_module_path

# Globals and constants variables.


class TestElsepaCrossSectionInfo(unittest.TestCase):

    def setUp(self):
        unittest.TestCase.setUp(self)
        self.filepath = current_module_path(__file__, "../testData/Casino3/EL29.els")
        if not os.path.isfile(self.filepath):
            self.skipTest("No file: {}".format(self.filepath))

    def tearDown(self):
        unittest.TestCase.tearDown(self)

    def testSkeleton(self):
        # self.fail("Test if the testcase is working.")
        self.assert_(True)

    def test_readFile(self):
        file = open(self.filepath, 'rb')
        self.elsCSInfo = ElsepaCrossSectionInfo()
        self.elsCSInfo.read_file(file)

        self.assertEquals(30105020, self.elsCSInfo._fileVersion)
        self.assertEquals(29, self.elsCSInfo._atomicNumber)
        self.assertAlmostEquals(1.000000e-01, self.elsCSInfo._energy_keV)
        self.assertAlmostEquals(2.870967e-02, self.elsCSInfo._totalCS_nm2)

        self.assertEquals(606, self.elsCSInfo._numberPoints)

        self.assertAlmostEquals(0.0, self.elsCSInfo._ratios[0])
        self.assertAlmostEquals(0.0, self.elsCSInfo._theta_rad[0])

        self.assertAlmostEquals(1.0, self.elsCSInfo._ratios[-1])
        self.assertAlmostEquals(3.1415926539, self.elsCSInfo._theta_rad[-1])

        # self.fail("Test if the testcase is working.")
