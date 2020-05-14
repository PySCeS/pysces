#!/usr/bin/env python

"""
PySCeS - Python Simulator for Cellular Systems (http://pysces.sourceforge.net)

Copyright (C) 2004-2020 B.G. Olivier, J.M. Rohwer, J.-H.S Hofmeyr all rights reserved,

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
__version__ = '0.9.8post1'

import os, re
import fileinput
import shutil
try:
    import configparser  # Py 3
except ImportError:
    import ConfigParser as configparser   # Py 2

try:
    print('Building an egg? %s.' % FRYING_EGGS)
except:
    FRYING_EGGS = False

try:
    import setuptools
    print('setuptools is available.')
except Exception as ex:
    print('setuptools not available.')

try:
    from numpy.distutils.core import setup, Extension
except Exception as ex:
    print(ex)
    print("PySCeS requires NumPy and SciPy 0.6x+\n")
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
nleq2 = 1                       # pysces/nleq2/readme.txt
# this is now obsolete with nleq2 4.3
##  nleq2_byteorder_override = 0    # allow a user supplied nleq2.f (mach_spec) otherwise PySCeS

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
# this is now obsolete with nleq2 4.3
    ##  elif e.startswith('nleq2_byteorder='):
        ##  nleq2_byteorder_override = int(e[16:])

# detects and uses IEEE fp big/little endian
### End user configuration section

########## From here on it's up to distutils ##########

# get the dir of setup.py
local_path = os.path.dirname(os.path.abspath(os.sys.argv[0]))
os.chdir(local_path)

myscripts = []
mydata_files = []

#add some model files into pscmodels
modfold = os.path.join('pysces', 'pscmodels')
mods = os.listdir(modfold)
alist = []
for x in mods:
    if x[-4:] != '.psc':
        pass
    else:
        alist.append(os.path.join(modfold,x))
mydata_files.append((os.path.join('pysces','pscmodels'), alist))

# Default configurations for the pyscfg.ini files
if os.sys.platform == 'win32':
    if FRYING_EGGS:
        eggdir = 'pysces-%s-py%s.%s-%s.egg' %(__version__, os.sys.version_info[0],\
        os.sys.version_info[1], os.sys.platform)
        installdir = os.path.join(os.sys.prefix,'lib','site-packages',eggdir,'pysces')
    else:
        installdir = os.path.join(os.sys.prefix,'lib','site-packages','pysces')
    config = {
    "install_dir"  : installdir,
    "model_dir"    : os.path.join(os.getenv('USERPROFILE'),'Pysces','psc'),
    "output_dir"   : os.path.join(os.getenv('USERPROFILE'),'Pysces'),
    "gnuplot_dir"  : None,
    "silentstart"  : False
    }
else:
    if hasattr(os.sys, 'lib'):
        lib = os.sys.lib
    else:
        lib = 'lib'
    config = {
    "install_dir"  : os.path.join(os.sys.prefix,lib,"python%d.%d" % tuple(os.sys.version_info[:2]) ,'site-packages','pysces'),
    "model_dir"    : os.path.join(os.path.expanduser('~'),'Pysces','psc'),
    "output_dir"   : os.path.join(os.path.expanduser('~'),'Pysces'),
    "gnuplot_dir"  : None,
    "silentstart"  : False
    }

def writeConfig(local_path, config={}):
    cfgfile = open(os.path.join(local_path,'pysces','pyscfg.ini'),'w')
    cp = configparser.ConfigParser()
    # PySCeS internal setup
    cp.add_section('Pysces')
    for key in config:
        print(repr(key) + ' :: ' + str(config[key]))
        cp.set('Pysces', key, str(config[key]))
    #add configuration data
    cp.add_section('PyscesConfig')
    cp.set('PyscesConfig','matplotlib', 'True')
    # OSX patch thanks to AF
    if os.sys.platform == 'darwin':
        cp.set('PyscesConfig','matplotlib_backend', 'MacOSX')
    else:
        cp.set('PyscesConfig','matplotlib_backend', 'TkAgg')
    cp.set('PyscesConfig','gnuplot', 'False')
    # Built in modules
    cp.add_section('PyscesModules')
    if pitcon:
        cp.set('PyscesModules','pitcon', 'True')
    else:
        cp.set('PyscesModules','pitcon', 'False')
    #PySCeS external module setup
    cp.add_section('ExternalModules')
    if nleq2:
        cp.set('ExternalModules','nleq2', 'True')
        mydata_files.append((os.path.join('pysces','nleq2'), [os.path.join('pysces','nleq2','nleq2_readme.txt')]))
    else:
        cp.set('ExternalModules','nleq2', 'False')
        mydata_files.append((os.path.join('pysces','nleq2'), [os.path.join('pysces','nleq2','readme.txt')]))
    cp.write(cfgfile)
    cfgfile.close()

writeConfig(local_path,config)
print('Default configuration file installed')

# my subpackage list
mypackages= ['pysces','pysces.tests','pysces.lib','pysces.pitcon',\
            'pysces.sandbox', 'pysces.contrib','pysces.contrib.demo', 'pysces.core2',\
            'pysces.kraken','pysces.kraken.controllers']

#PySCeS modules
mymodules = []

if pitcon:
    print('\nBuilding pitcon')
    extpath = os.path.join('pysces', 'pitcon')
    pitcon = Extension('pysces.pitcon.pitcon',[os.path.join(extpath,'pitcon.pyf'),os.path.join(extpath,'pcon61subd.f'),os.path.join(extpath,'dpcon61.f'),os.path.join(extpath,'dpcon61w.f')])
    mymodules.append(pitcon)
    #mydata_files.append((os.path.join('pysces','pitcon'), [os.path.join(local_path, 'pysces', 'pitcon','readme.txt'), os.path.join(local_path, 'pysces', 'pitcon','readme.txt')]))
else:
    print('\nSkipping pitcon')


if nleq2:
    print('\nBuilding nleq2')
    # this is now obsolete with nleq2 4.3 ... i hope !
    ##  print 'System ByteOrder', os.sys.byteorder
    ##  if os.path.exists(os.path.join(local_path, 'pysces', 'nleq2','nleq2.f')) and nleq2_byteorder_override:
        ##  print 'INFO: using user supplied nleq2.f'
    ##  else:
        ##  if os.sys.byteorder == 'little':
            ##  shutil.copyfile(os.path.join(extpath,'nleq2_little.f'), os.path.join(extpath,'nleq2.f'))
        ##  elif os.sys.byteorder == 'big':
            ##  shutil.copyfile(os.path.join(extpath,'nleq2_big.f'), os.path.join(extpath,'nleq2.f'))
    extpath = os.path.join('pysces', 'nleq2')
    nleq2 = Extension('pysces.nleq2.nleq2',[os.path.join(extpath,'nleq2.pyf'),\
            os.path.join(extpath,'nleq2.f'), os.path.join(extpath,'linalg_nleq2.f'),\
            os.path.join(extpath,'zibmon.f'), os.path.join(extpath,'zibsec.f'),\
            os.path.join(extpath,'zibconst.f'), os.path.join(extpath,'wnorm.f')
            ])
    mymodules.append(nleq2)
    mypackages.append('pysces.nleq2')
else:
    print('\n')

if len(mymodules) == 0:
    noext = Extension('None',[],None) # Not ideal but seems safe
    mymodules.append(noext)

# Data files to copy
mydata_files.append((os.path.join('pysces'), [os.path.join('pysces','pyscfg.ini')]))
mydata_files.append(('',[os.path.join('pysces','pysces.pth')]))
mydata_files.append((os.path.join('pysces','docs'), [os.path.join('pysces','docs','userguide.pdf')]))
mydata_files.append((os.path.join('pysces','examples'), [os.path.join('pysces','examples',examplefile) for examplefile in os.listdir(os.path.join(local_path,'pysces','examples'))]))
##not sure if this is necessary anymore, removed to test
if os.sys.platform == 'win32':
    mydata_files.append((os.path.join('pysces','win32'), [os.path.join('pysces','win32','libquadmath-0.dll'), os.path.join('pysces','win32','libgfortran-3.dll'),\
    os.path.join('pysces','win32','libgcc_s_seh-1.dll'), os.path.join('pysces','win32','libwinpthread-1.dll')]))

os.chdir(local_path)
# Install packages and the metatool binaries as "data"
setup(name="pysces",
    version = __version__,
    description = "The Python Simulator for Cellular Systems - simulation and analysis tools for modelling biological systems",
    long_description = """
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
    author = "Brett G. Olivier",
    author_email = "bgoli@users.sourceforge.net",
    maintainer = "Brett G. Olivier",
    maintainer_email = "bgoli@users.sourceforge.net",
    url = "http://pysces.sourceforge.net",
    download_url = "http://pysces.sourceforge.net/download.html",
    license = "New BSD style",
    keywords = "computational systems biology, modelling, simulation, systems biology" ,
    zip_safe = False,
    #requires = ['numpy','scipy','matplotlib'],
    install_requires = ['numpy','scipy','matplotlib','python-libsbml','nose'],
    platforms = ["Windows", "Linux", "macOS"],
    classifiers = [
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
    'Topic :: Scientific/Engineering :: Chemistry'],
    packages = mypackages,
    data_files = mydata_files,
    ext_modules = mymodules
    )

