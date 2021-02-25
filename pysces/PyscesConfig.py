"""
PySCeS - Python Simulator for Cellular Systems (http://pysces.sourceforge.net)

Copyright (C) 2004-2020 B.G. Olivier, J.M. Rohwer, J.-H.S Hofmeyr all rights reserved,

Brett G. Olivier (bgoli@users.sourceforge.net)
Triple-J Group for Molecular Cell Physiology
Stellenbosch University, South Africa.

Permission to use, modify, and distribute this software is given under the
terms of the PySceS (BSD style) license. See LICENSE.txt that came with
this distribution for specifics.

NO WARRANTY IS EXPRESSED OR IMPLIED.  USE AT YOUR OWN RISK.
Brett G. Olivier
"""
from __future__ import division, print_function
from __future__ import absolute_import
from __future__ import unicode_literals

from pysces.version import __version__

__doc__ = '''
            PyscesConfig
            ------------

            This module contains templates for the default configuration files
            on POSIX and WIN32 systems as well as utility functions for reading
            and writing them.
            '''

try:
    import configparser  # Py 3
except ImportError:
    import ConfigParser as configparser  # Py 2

import string
import os

if os.sys.platform == 'win32':
    __DefaultWin = {
        "install_dir": os.path.join(os.sys.prefix, 'lib', 'site-packages', 'pysces'),
        "model_dir": os.path.join(os.getenv('USERPROFILE'), 'Pysces', 'psc'),
        "output_dir": os.path.join(os.getenv('USERPROFILE'), 'Pysces'),
        "gnuplot_dir": None,
        "pitcon": True,
        "nleq2": True,
        "gnuplot": False,
        "matplotlib": True,
        "matplotlib_backend": 'TKagg',
        "silentstart": False,
        "change_dir_on_start": False,
    }
    __DefaultWinUsr = {
        "model_dir": os.path.join(os.getenv('USERPROFILE'), 'Pysces', 'psc'),
        "output_dir": os.path.join(os.getenv('USERPROFILE'), 'Pysces'),
        "silentstart": False,
        "change_dir_on_start": False,
    }
else:
    if hasattr(os.sys, 'lib'):
        lib = os.sys.lib
    else:
        lib = 'lib'
    __DefaultPosix = {
        "install_dir": os.path.join(
            os.sys.prefix,
            lib,
            "python%d.%d" % tuple(os.sys.version_info[:2]),
            'site-packages',
            'pysces',
        ),
        "model_dir": os.path.join(os.path.expanduser('~'), 'Pysces', 'psc'),
        "output_dir": os.path.join(os.path.expanduser('~'), 'Pysces'),
        "gnuplot_dir": None,
        "pitcon": True,
        "nleq2": True,
        "gnuplot": False,
        "matplotlib": True,
        "matplotlib_backend": 'TKagg',
        "silentstart": False,
        "change_dir_on_start": False,
    }
    __DefaultPosixUsr = {
        "model_dir": os.path.join(os.path.expanduser('~'), 'Pysces', 'psc'),
        "output_dir": os.path.join(os.path.expanduser('~'), 'Pysces'),
        "silentstart": False,
        "change_dir_on_start": False,
    }
    # OSX patch by AF
    if os.sys.platform == 'darwin':
        __DefaultPosix["matplotlib_backend"] = 'MacOSX'


def ReadConfig(file_path, config={}):
    """
    Read a PySCeS configuration file

    - *file_path*	full path to file
    - *config [default={}]* configuration data
    """
    filein = open(file_path, 'r')
    cp = configparser.ConfigParser()
    cp.readfp(filein)
    for sec in cp.sections():
        name = sec.lower()
        for opt in cp.options(sec):
            config[opt.lower()] = cp.get(sec, opt).strip()
    filein.close()
    return config


def WriteConfig(file_path, config={}, section='Pysces'):
    """
    Write a PySCeS configuration file

    - *file_path* full path to file
    - *config* [default={}]: config dictionary
    - *section* [default='Pysces']: default man section name

    """
    cfgfile = open(file_path, 'w')
    cp = configparser.ConfigParser()
    cp.add_section(section)
    for key in config:
        cp.set(section, key, str(config[key]))
    cp.write(cfgfile)
    cfgfile.close()
