#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. py:currentmodule:: eecs
.. moduleauthor:: Hendrix Demers <hendrix.demers@mail.mcgill.ca>

pyElectronCrossSections python package for different elastic scattering cross section of electrons models and tools.
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
import os.path

# Third party modules.

# Local modules.

# Globals and constants variables.
__author__ = """Hendrix Demers"""
__email__ = 'hendrix.demers@mail.mcgill.ca'
__version__ = '0.2'
__project_name__ = "pyElectronCrossSections"


def get_current_module_path(module_path, relative_path=""):
    """
    Extract the current module path and combine it with the relative path and return it.

    :param str module_path: Pass the `__file__` python keyword for this parameter
    :param str relative_path: The relative path to combine with the module path
    :return: The path obtained when combine the module path and relative path
    :rtype: str
    """
    base_path = os.path.dirname(module_path)
    file_path = os.path.join(base_path, relative_path)
    file_path = os.path.abspath(file_path)
    file_path = os.path.normpath(file_path)

    return file_path


def create_path(path):
    """
    Create a path from the input string if does not exists.
    Does not try to distinct between file and directory in the input string.
    path = "dir1/filename.ext" => "dir1/filename.ext/"
    where the new directory "filename.ext" is created.
    @param[in] path input string.
    @return the path with the path separator at the end.
    """
    path = os.path.normpath(path)
    if not os.path.exists(path):
        os.makedirs(path)

    if len(path) > 0 and path[-1] != os.sep:
        path += os.sep

    return path
