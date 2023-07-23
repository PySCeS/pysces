#!/usr/bin/env python

"""
PySCeS - Python Simulator for Cellular Systems (http://pysces.sourceforge.net)

Copyright (C) 2004-2023 B.G. Olivier, J.M. Rohwer, J.-H.S Hofmeyr all rights reserved,

Brett G. Olivier (bgoli@users.sourceforge.net)
Triple-J Group for Molecular Cell Physiology
Stellenbosch University, South Africa

Permission to use, modify, and distribute this software is given under the
terms of the PySceS (BSD style) license.  See LICENSE.txt that came with
this distribution for specifics.

NO WARRANTY IS EXPRESSED OR IMPLIED.  USE AT YOUR OWN RISK.
Brett G. Olivier
"""
from __future__ import division, print_function
from __future__ import absolute_import

__doc__ = "PySCeS: the Python Simulator for Cellular Systems setup file"

# get __version__ from version.txt
with open('pysces/version.txt') as f:
    __version__ = f.read().strip()

# avoid duplication
with open('requirements.txt') as f:
    requirements = f.read().splitlines()

import os
import re
import sysconfig

try:
    import configparser  # Py 3
except ImportError:
    import ConfigParser as configparser  # Py 2

try:
    import setuptools

    print('setuptools is available.')
except Exception as ex:
    print('setuptools not available.')

try:
    from skbuild import setup
except Exception as ex:
    print(ex)
    print("PySCeS requires scikit-build to build.\n")
    os.sys.exit(-1)

try:
    import distutils.command.bdist_conda
    print('bdist_conda is available.')
except Exception as ex:
    print('bdist_conda not available.')


########## User configuration section ##########

# Install extension modules (zero for skip module) - brett 20040310
pitcon = 1

##Special licence libraries see the relevant readme.txt file for details
nleq2 = 1  # pysces/nleq2/readme.txt

use = re.split('\s+', os.getenv('PYSCES_USE', ''))
for e in use:
    if e == 'pitcon':
        pitcon = 1
    elif e == 'nopitcon':
        pitcon = 0
    elif e == 'nleq2':
        nleq2 = 1
    elif e == 'nonleq2':
        nleq2 = 0
'''
# this is now obsolete with nleq2 4.3
    elif e.startswith('nleq2_byteorder='):
        nleq2_byteorder_override = int(e[16:])
# detects and uses IEEE fp big/little endian
'''
### End user configuration section

########## From here on it's up to scikit-build ##########

# get the dir of setup.py
local_path = os.path.dirname(os.path.abspath(__file__))
os.chdir(local_path)

mypackage_data = {}

# add some model files into pscmodels
mypackage_data['pysces.pscmodels'] = ['*.psc']

# Default configurations for the pyscfg.ini files
config = {
    "install_dir": os.path.join(
        sysconfig.get_path('platlib'),
        'pysces',
    ),
    "model_dir": os.path.join(os.path.expanduser('~'), 'Pysces', 'psc'),
    "output_dir": os.path.join(os.path.expanduser('~'), 'Pysces'),
    "gnuplot_dir": None,
    "silentstart": False,
    "change_dir_on_start": False,
    "custom_datatype": None,
}

def writeConfig(local_path, config={}):
    cfgfile = open(os.path.join(local_path, 'pysces', 'pyscfg.ini'), 'w')
    cp = configparser.ConfigParser()
    # PySCeS internal setup
    cp.add_section('Pysces')
    for key in config:
        print(repr(key) + ' :: ' + str(config[key]))
        cp.set('Pysces', key, str(config[key]))
    # add configuration data
    cp.add_section('PyscesConfig')
    cp.set('PyscesConfig', 'matplotlib', 'True')
    # OSX patch thanks to AF
    if os.sys.platform == 'darwin':
        cp.set('PyscesConfig', 'matplotlib_backend', 'MacOSX')
    else:
        cp.set('PyscesConfig', 'matplotlib_backend', 'TkAgg')
    cp.set('PyscesConfig', 'gnuplot', 'False')
    # Built in modules
    cp.add_section('PyscesModules')
    if pitcon:
        cp.set('PyscesModules', 'pitcon', 'True')
    else:
        cp.set('PyscesModules', 'pitcon', 'False')
    # PySCeS external module setup
    cp.add_section('ExternalModules')
    if nleq2:
        cp.set('ExternalModules', 'nleq2', 'True')
    else:
        cp.set('ExternalModules', 'nleq2', 'False')
    cp.write(cfgfile)
    cfgfile.close()

writeConfig(local_path, config)
print('Default configuration file installed')

# my subpackage list
mypackages = [
    'pysces',
    'pysces.tests',
    'pysces.examples',
    'pysces.pscmodels',
    'pysces.docs',
    'pysces.lib',
    'pysces.sandbox',
    'pysces.contrib',
    'pysces.contrib.demo',
    'pysces.core2',
    'pysces.kraken',
    'pysces.kraken.controllers',
]

if pitcon:
    print('\nBuilding pitcon')
    mypackages.append('pysces.pitcon')
else:
    print('\nSkipping pitcon')


if nleq2:
    print('\nBuilding nleq2')
    mypackages.append('pysces.nleq2')
else:
    print('\nSkipping nleq2')
mypackage_data['pysces.nleq2'] = ['nleq2_readme.txt', 'readme.txt']

# Data files to copy
mypackage_data['pysces'] = ['pyscfg.ini', 'version.txt']
mypackage_data['pysces.docs'] = ['userguide.pdf']
mypackage_data['pysces.examples'] = ['*.ipy']

os.chdir(local_path)

setup(
    name="pysces",
    version=__version__,
    description="The Python Simulator for Cellular Systems - simulation and analysis tools for modelling biological systems",
    long_description="""
 PySCeS is developed by the Triple-J Group for Molecular Cell Physiology
 in order to try model and understand the complex processes and systems
 which make up the living cell.

    PySCeS features, amongst other things:
    - A text based model description language.
    - A structural analysis module.
    - Integrators for time simulation
    - Non-linear solvers for steady-state analysis
    - A module for performing Metabolic Control Analysis
    - A bifurcation module for systems which exhibit multiple steady states
    - A variety of extra utilites for parameter scans, data output and plotting.
    - A dynamic module loading framework.
    - SBML import and export capability.
    """,
    author="Brett G. Olivier and Johann M. Rohwer",
    author_email="pysces@googlegroups.com",
    maintainer="Brett G. Olivier and Johann M. Rohwer",
    maintainer_email="pysces@googlegroups.com",
    url="http://pysces.sourceforge.net",
    download_url="https://pypi.org/project/pysces/#files",
    license="New BSD style",
    keywords="computational systems biology, modelling, simulation, systems biology",
    zip_safe=False,
    install_requires=requirements,
    extras_require={
        'parscan': ['ipyparallel'],
        'cvode': ['assimulo'],
        'sbml': ['python_libsbml'],
        'all': ['ipyparallel', 'assimulo', 'python_libsbml'],
    },
    platforms=["Windows", "POSIX", "Max OSX"],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Development Status :: 6 - Mature',
        'Environment :: Console',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Fortran',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
    ],
    packages=mypackages,
    package_data=mypackage_data,
)
