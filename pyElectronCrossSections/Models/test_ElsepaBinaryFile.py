#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. py:currentmodule:: pyElectronCrossSection.Models.test_ElsepaBinaryFile

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
from pyElectronCrossSections.Models.ElsepaBinaryFile import ElsepaBinaryFile
from pyElectronCrossSections import current_module_path

# Globals and constants variables.


class TestElsepaBinaryFile(unittest.TestCase):

    def setUp(self):
        unittest.TestCase.setUp(self)
        self.filepath = current_module_path(__file__, "../testData/Casino3/EL29.els")
        if not os.path.isfile(self.filepath):
            self.skipTest("No file: {}".format(self.filepath))
        self.elsFile = ElsepaBinaryFile(self.filepath)

    def tearDown(self):
        unittest.TestCase.tearDown(self)

    def testSkeleton(self):
        # self.fail("Test if the testcase is working.")
        self.assert_(True)

    def test_init(self):
        self.assertEquals(self.filepath, self.elsFile._filepath)

        # self.fail("Test if the testcase is working.")

    def test_readFile(self):
        file = open(self.filepath, 'rb')
        self.elsFile.read_file(file)

        self.assertEquals(48, len(self.elsFile._elsCSInfoList))

        els_cs_info = self.elsFile._elsCSInfoList[0]
        self.assertEquals(30105020, els_cs_info._fileVersion)
        self.assertEquals(29, els_cs_info._atomicNumber)
        self.assertAlmostEquals(1.000000e-01, els_cs_info._energy_keV)
        self.assertAlmostEquals(2.870967e-02, els_cs_info._totalCS_nm2)

        self.assertEquals(606, els_cs_info._numberPoints)

        self.assertAlmostEquals(0.0, els_cs_info._ratios[0])
        self.assertAlmostEquals(0.0, els_cs_info._theta_rad[0])

        self.assertAlmostEquals(1.0, els_cs_info._ratios[-1])
        self.assertAlmostEquals(3.1415926539, els_cs_info._theta_rad[-1])

        els_cs_info = self.elsFile._elsCSInfoList[-1]
        self.assertEquals(30105020, els_cs_info._fileVersion)
        self.assertEquals(29, els_cs_info._atomicNumber)
        self.assertAlmostEquals(500.0, els_cs_info._energy_keV)
        self.assertAlmostEquals(0.000247381469262322, els_cs_info._totalCS_nm2)

        self.assertEquals(606, els_cs_info._numberPoints)

        self.assertAlmostEquals(0.0, els_cs_info._ratios[0])
        self.assertAlmostEquals(0.0, els_cs_info._theta_rad[0])

        self.assertAlmostEquals(1.0, els_cs_info._ratios[-1])
        self.assertAlmostEquals(3.1415926539, els_cs_info._theta_rad[-1])

        self.assertAlmostEquals(0.999999999788251, els_cs_info._ratios[-1])
        self.assertAlmostEquals(3.13286600763917, els_cs_info._theta_rad[-2])

        # self.fail("Test if the testcase is working.")
