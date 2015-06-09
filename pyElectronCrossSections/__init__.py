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

import os.path

def current_module_path(modulePath, relativePath=""):
    basepath = os.path.dirname(modulePath)

    filepath = os.path.join(basepath, relativePath)
    filepath = os.path.normpath(filepath)

    return filepath
