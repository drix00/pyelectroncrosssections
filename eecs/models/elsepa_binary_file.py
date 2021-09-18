#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. py:currentmodule:: eecs.models.elsepa_binary_file
.. moduleauthor:: Hendrix Demers <hendrix.demers@mail.mcgill.ca>

Description
"""

###############################################################################
# Copyright 2021 Hendrix Demers
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
import struct
import logging

# Third party modules.

# Local modules.

# Project modules.
from eecs.models.elsepa_cross_section_info import ElsepaCrossSectionInfo

# Globals and constants variables.


class ElsepaBinaryFile:
    def __init__(self, filepath):
        self._filepath = filepath
        self._els_cs_info_list = []

        file = open(self._filepath, 'rb')
        self.read_file(file)

    def read_file(self, file):
        self._els_cs_info_list = []
        while file:
            try:
                elsepa_cs_info = ElsepaCrossSectionInfo()
                elsepa_cs_info.read_file(file)
                self._els_cs_info_list.append(elsepa_cs_info)
            except struct.error as message:
                logging.error(message)
                return
