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
from .version import __version__

__doc__ = """
            PySCeS: the Python Simulator for Cellular Systems
            -------------------------------------------------
            PySCeS is developed by the Triple-J Group for Molecular Cell Physiology
            in order to try model and understand the complex processes and systems
            which make up the living cell. PySCeS features, amongst other things:

            - A text based Model Description Language.
            - A structural analysis module.
            - Integrators for time simulation
            - Non-linear solvers for steady-state analysis
            - A module for performing Metabolic Control Analysis
            - A bifurcation module for systems which exhibit multiple steady states
            - A variety of extra utilites for parameter scans, data output and plotting.
            - A dynamic module loading framework.
            - SBML import and export capability.
            """

import os, time
from pkg_resources import get_build_platform
from . import PyscesConfig
from . import PyscesParse
from . import PyscesLink as link
from . import PyscesSED as SED

from .PyscesUtils import str2bool
from .PyscesModelMap import ModelMap

# TODO get rid unused imports
from .PyscesWeb import PyscesHTML

html = PyscesHTML()

DEBUG = False

inipath = None
lpath = None
install_dir = None
output_dir = None
model_dir = None
pitcon_switch = False
nleq2_switch = False
__USE_MATPLOTLIB__ = True
__MATPLOTLIB_BACKEND__ = 'TkAgg'
__USE_GNUPLOT__ = False
__SILENT_START__ = False

extra_dll_dir = os.path.join(os.path.dirname(__file__), '.libs')
if os.sys.platform == 'win32' and os.path.isdir(extra_dll_dir):
    os.environ["PATH"] += os.pathsep + extra_dll_dir
    if hasattr(os, 'add_dll_directory'):
        os.add_dll_directory(extra_dll_dir)

if os.sys.platform == 'win32':
    __PyscesConfigDefault = PyscesConfig.__DefaultWin
else:
    __PyscesConfigDefault = PyscesConfig.__DefaultPosix

if DEBUG:
    print(time.strftime('1-%H:%M:%S'))

eggdir = 'pysces-%s-py%s.%s-%s.egg' % (
    __version__,
    os.sys.version_info[0],
    os.sys.version_info[1],
    get_build_platform(),
)
for path in os.sys.path:
    chkPath = path.split(os.path.sep)[-1]
    if chkPath == 'pysces' and path != os.getcwd():
        if os.path.isdir(os.path.join(path, 'pysces')):
            # for in-place development with setup.py develop
            install_dir = os.path.join(path, 'pysces')
        else:
            install_dir = path
        inipath = os.path.join(install_dir, 'pyscfg.ini')
        break
    elif chkPath == eggdir:
        install_dir = os.path.join(path, 'pysces')
        inipath = os.path.join(install_dir, 'pyscfg.ini')
        break
if inipath == None:
    for k in os.sys.path_importer_cache:
        if k.split(os.path.sep)[-1] == 'pysces':
            install_dir = k
            inipath = os.path.join(install_dir, 'pyscfg.ini')
            break
        elif k.split(os.path.sep)[-1] == eggdir:
            install_dir = os.path.join(k, 'pysces')
            inipath = os.path.join(install_dir, 'pyscfg.ini')
            break
del eggdir
if DEBUG:
    print(time.strftime('2-%H:%M:%S'))

try:
    __config_dict = PyscesConfig.ReadConfig(inipath, config=__PyscesConfigDefault)
except Exception as ex:
    print(ex)
    print('Cwd', os.getcwd())
    print('\nWARNING: Cannot read pyscfg.ini using default values\n')
    __config_dict = __PyscesConfigDefault

# Read config
for key in __config_dict:
    if key == 'pitcon':
        pitcon_switch = str2bool(__config_dict[key])
    elif key == 'nleq2':
        nleq2_switch = str2bool(__config_dict[key])
    elif key == 'matplotlib':
        __USE_MATPLOTLIB__ = str2bool(__config_dict[key])
    elif key == 'matplotlib_backend':
        __MATPLOTLIB_BACKEND__ = __config_dict[key]
    elif key == 'gnuplot':
        __USE_GNUPLOT__ = str2bool(__config_dict[key])
    elif key == 'gnuplot_dir':
        GNUPLOT_DIR = __config_dict[key]
        if GNUPLOT_DIR == 'None':
            GNUPLOT_DIR = None
    elif key == 'silentstart':
        __SILENT_START__ = str2bool(__config_dict[key])
    elif key == 'change_dir_on_start':
        __CHGDIR_ON_START__ = str2bool(__config_dict[key])
assert inipath != None, '\nNo configuration file found'

if DEBUG:
    print(time.strftime('3-%H:%M:%S'))

__userdict = None
if os.sys.platform != 'win32':
    if os.path.exists(
        os.path.join(os.path.expanduser('~'), 'Pysces', '.pys_usercfg.ini')
    ):
        __userdict = PyscesConfig.ReadConfig(
            os.path.join(os.path.expanduser('~'), 'Pysces', '.pys_usercfg.ini'),
            PyscesConfig.__DefaultPosixUsr,
        )
    else:
        if not os.path.exists(os.path.join(os.path.expanduser('~'), 'Pysces')):
            os.makedirs(os.path.join(os.path.expanduser('~'), 'Pysces'))
        PyscesConfig.WriteConfig(
            os.path.join(os.path.expanduser('~'), 'Pysces', '.pys_usercfg.ini'),
            config=PyscesConfig.__DefaultPosixUsr,
            section='Pysces',
        )
        __userdict = PyscesConfig.ReadConfig(
            os.path.join(os.path.expanduser('~'), 'Pysces', '.pys_usercfg.ini'),
            PyscesConfig.__DefaultPosixUsr,
        )
else:
    if os.path.exists(
        os.path.join(os.getenv('HOMEDRIVE') + os.path.sep, 'Pysces', '.pys_usercfg.ini')
    ):
        __userdict = PyscesConfig.ReadConfig(
            os.path.join(
                os.getenv('HOMEDRIVE') + os.path.sep, 'Pysces', '.pys_usercfg.ini'
            ),
            PyscesConfig.__DefaultWinUsr,
        )
    elif os.path.exists(
        os.path.join(os.getenv('USERPROFILE'), 'Pysces', '.pys_usercfg.ini')
    ):
        __userdict = PyscesConfig.ReadConfig(
            os.path.join(os.getenv('USERPROFILE'), 'Pysces', '.pys_usercfg.ini'),
            PyscesConfig.__DefaultWinUsr,
        )
    else:
        if not os.path.exists(os.path.join(os.getenv('USERPROFILE'), 'Pysces')):
            os.makedirs(os.path.join(os.getenv('USERPROFILE'), 'Pysces'))
        PyscesConfig.WriteConfig(
            os.path.join(os.getenv('USERPROFILE'), 'Pysces', '.pys_usercfg.ini'),
            config=PyscesConfig.__DefaultWinUsr,
            section='Pysces',
        )
        __userdict = PyscesConfig.ReadConfig(
            os.path.join(os.getenv('USERPROFILE'), 'Pysces', '.pys_usercfg.ini'),
            PyscesConfig.__DefaultWinUsr,
        )
for key in __userdict:
    if key == 'output_dir':
        output_dir = __userdict[key]
        if not os.path.exists(__userdict[key]):
            os.makedirs(__userdict[key])
    elif key == 'model_dir':
        model_dir = __userdict[key]
        if not os.path.exists(__userdict[key]):
            os.makedirs(__userdict[key])
    elif key == 'matplotlib':
        __USE_MATPLOTLIB__ = str2bool(__userdict[key])
    elif key == 'matplotlib_backend':
        __MATPLOTLIB_BACKEND__ = __userdict[key]
    elif key == 'gnuplot':
        __USE_GNUPLOT__ = str2bool(__userdict[key])
    elif key == 'gnuplot_dir':
        GNUPLOT_DIR = __userdict[key]
        if GNUPLOT_DIR == 'None':
            GNUPLOT_DIR = None
    elif key == 'silentstart':
        __SILENT_START__ = str2bool(__userdict[key])
    elif key == 'change_dir_on_start':
        __CHGDIR_ON_START__ = str2bool(__userdict[key])
assert output_dir != None, '\nNo output directory defined'
assert model_dir != None, '\nNo model directory defined'

# following is to get the full path when .pys_usercfg.ini specifies CWD as follows:
# output_dir = ./
backup_dir = os.getcwd()
os.chdir(output_dir)
output_dir = os.getcwd()
os.chdir(backup_dir)
os.chdir(model_dir)
model_dir = os.getcwd()
os.chdir(backup_dir)

del PyscesConfig

if DEBUG:
    print(time.strftime('4-%H:%M:%S'))

# initialise pysces.interface.*
try:
    from . import PyscesInterfaces

    interface = PyscesInterfaces.Core2interfaces()
except Exception as ex:
    print('INFO: pysces.interface.* not available')
    print(ex)
    interface = None

# initialise pysces.plt.*
from . import PyscesPlot2

gplt = None
mplt = None
if __USE_MATPLOTLIB__:
    try:
        mplt = PyscesPlot2.MatplotlibUPI(
            work_dir=output_dir, backend=__MATPLOTLIB_BACKEND__
        )
        if not __SILENT_START__:
            print('Matplotlib interface loaded (pysces.plt.m)')
    except Exception as ex:
        print('Matplotlib interface not available')
        if DEBUG:
            print(ex)
        __USE_MATPLOTLIB__ = False
if __USE_GNUPLOT__:
    if GNUPLOT_DIR == None or not os.path.exists(GNUPLOT_DIR):
        print(
            '''GnuPlot has been enabled but the path to the executable has
not been defined (or does not exist). Please set the "gnuplot_dir" key
in your pyscfg.ini file.
        '''
        )
    else:
        try:
            if DEBUG:
                print(GNUPLOT_DIR)
            gplt = PyscesPlot2.GnuPlotUPI(work_dir=output_dir, gnuplot_dir=GNUPLOT_DIR)
            if not __SILENT_START__:
                print('GnuPlot interface loaded (pysces.plt.g)')
        except Exception as ex:
            print('GnuPlot interface not available')
            if DEBUG:
                print(ex)
            __USE_GNUPLOT__ = False

plt = None
if __USE_MATPLOTLIB__ or __USE_GNUPLOT__:
    plt = PyscesPlot2.PyscesUPI()
if __USE_MATPLOTLIB__ and not __USE_GNUPLOT__:
    plt.p_setInterface('matplotlib', mplt)
elif __USE_GNUPLOT__ and not __USE_MATPLOTLIB__:
    plt.p_setInterface('gnuplot', gplt)
elif __USE_GNUPLOT__ and __USE_MATPLOTLIB__:
    plt.p_setInterface('matplotlib', mplt)
    plt.p_setInterface('gnuplot', gplt)
    plt.p_deactivateInterface('gnuplot')

if DEBUG:
    print(time.strftime('5-%H:%M:%S'))

alt_import = False
alt_import_pitcon = False
alt_import_nleq2 = False
if os.sys.platform == 'win32':
    if pitcon_switch:
        os.sys.path.append(os.path.join(install_dir, 'pitcon'))
        try:
            from .pitcon import pitcon as pitcon

            if not __SILENT_START__:
                print('Continuation routines available')
        except Exception as ex:
            try:
                os.environ['path'] = '%s;%s' % (
                    os.path.join(install_dir, 'win32'),
                    os.environ['path'],
                )
                from .pitcon import pitcon as pitcon

                if not __SILENT_START__:
                    print('Continuation routines available')
            except Exception as ex:
                # print 'Attempting alternate pitcon import ...'
                # alt_import = True
                # alt_import_pitcon = True
                print(ex)
                print('INFO: Pitcon import failed: continuation not available')

    if nleq2_switch:
        os.sys.path.append(os.path.join(install_dir, 'nleq2'))
        try:
            from .nleq2 import nleq2 as nleq2

            if not __SILENT_START__:
                print('NLEQ2 routines available')
        except Exception as ex:
            try:
                os.environ['path'] = '%s;%s' % (
                    os.path.join(install_dir, 'win32'),
                    os.environ['path'],
                )
                from .nleq2 import nleq2 as nleq2

                if not __SILENT_START__:
                    print('NLEQ2 routines available')
            except Exception as ex:
                # print 'Attempting alternate nleq2 import ...'
                # alt_import = True
                # alt_import_nleq2 = True
                print(ex)
                print('INFO: NLEQ2 import failed: option not available')
else:
    if pitcon_switch:
        os.sys.path.append(os.path.join(install_dir, 'pitcon'))
        try:
            from .pitcon import pitcon as pitcon

            if not __SILENT_START__:
                print('Pitcon routines available')
        except Exception as ex:
            # print ex
            alt_import = True
            alt_import_pitcon = True

            print('Attempting alternate pitcon import ...')
    if nleq2_switch:
        os.sys.path.append(os.path.join(install_dir, 'nleq2'))
        try:
            from .nleq2 import nleq2 as nleq2

            if not __SILENT_START__:
                print('NLEQ2 routines available')
        except Exception as ex:
            # print ex
            alt_import = True
            alt_import_nleq2 = True
            print('Attempting alternate nleq2 import ...')

if DEBUG:
    print(time.strftime('6-%H:%M:%S'))

if alt_import:
    savedir = os.getcwd()
    for tpath in os.sys.path:
        if alt_import_pitcon:
            try:
                if (
                    os.path.exists(os.path.join(tpath, 'pysces', 'pitcon'))
                    and tpath != ''
                ):
                    os.chdir(os.path.join(tpath, 'pysces', 'pitcon'))
                    from . import pitcon

                    if not __SILENT_START__:
                        print('Continuation routines available (A)')
            except Exception as ex:
                print(ex)
                print(
                    'INFO: Alternate pitcon import failed: continuation not available'
                )
        if alt_import_nleq2:
            try:
                if (
                    os.path.exists(os.path.join(tpath, 'pysces', 'nleq2'))
                    and tpath != ''
                ):
                    os.chdir(os.path.join(tpath, 'pysces', 'nleq2'))
                    from . import nleq2

                    if not __SILENT_START__:
                        print('NLEQ2 routines available (A)')
            except Exception as ex:
                print(ex)
                nleq2_switch = False
                print('INFO: Alternate NLEQ2 import failed: option not available')
    os.chdir(savedir)

if DEBUG:
    print(time.strftime('7-%H:%M:%S'))

# check for libsbml
try:
    import libsbml as SBML

    if not __SILENT_START__:
        print("SBML support available")
except ImportError as ex:
    SBML = None
    print(ex)
    print("INFO: No SBML library found, SBML support not available")

# This has to come at the end
from .PyscesModel import PysMod as model
from .PyscesModel import ScanDataObj as ScanDataObj

PyscesModel.interface = interface
from .PyscesTest import PyscesTest as test

write = None
try:
    from .PyscesUtils import WriteOutput

    write = WriteOutput()
    del WriteOutput
except ImportError as ex:
    pass

from .PyscesScan import PITCONScanUtils, Scanner

try:
    from .RateChar import RateChar

    if not __SILENT_START__:
        print("RateChar is available")
except Exception as ex:
    RateChar = None
    # print "RateChar not available"

# ParScanner import
try:
    from .PyscesParScan import ParScanner

    if not __SILENT_START__:
        print("Parallel scanner is available")
except ImportError as ex:
    ParScanner = None
    print(ex)
    print("INFO: Parallel scanner not available")

if DEBUG:
    print(time.strftime('9-%H:%M:%S'))

if not __SILENT_START__:
    print('\nPySCeS environment\n******************')
    print('pysces.model_dir = ' + model_dir)
    print('pysces.output_dir = ' + output_dir)

    print('\n\n***********************************************************************')
    print(
        '* Welcome to PySCeS ('
        + __version__
        + ') - Python Simulator for Cellular Systems   *'
    )
    print('*                http://pysces.sourceforge.net                        *')
    ##  print '*                       Somewhere In Time                             *'
    print('* Copyright(C) B.G. Olivier, J.M. Rohwer, J.-H.S. Hofmeyr, 2004-2020  *')
    print('* Triple-J Group for Molecular Cell Physiology                        *')
    print('* Stellenbosch University, ZA and VU University Amsterdam, NL         *')
    print('* PySCeS is distributed under the PySCeS (BSD style) licence, see     *')
    print('* LICENCE.txt (supplied with this release) for details                *')
    ##  print '*                 ** Read about PySCeS **                             *'
    print('* Please cite PySCeS with: doi:10.1093/bioinformatics/bti046          *')
    print('***********************************************************************')

try:
    del os, key, gplt, mplt
except Exception as ex:
    print(ex)
    print('\n\nOops I did it again error ...\n\n')
if DEBUG:
    print(time.strftime('10-%H:%M:%S'))
