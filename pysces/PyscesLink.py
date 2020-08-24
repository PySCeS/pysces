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
          PyscesLink
          ----------

          Interfaces to external software and API's, has replaced the PySCeS contrib classes.
          '''

try:
    input = raw_input  # Py2 compatibility
except NameError:
    pass

# for METATOOLlink
import os, re, io

# for SBWWebLink
try:
    from urllib.request import (
        ProxyHandler,
        build_opener,
        HTTPHandler,
        install_opener,
        urlopen,
    )  # Py 3
except ImportError:
    from urllib2 import (
        ProxyHandler,
        build_opener,
        HTTPHandler,
        install_opener,
        urlopen,
    )  # Py 2

try:
    from urllib.parse import urlencode  # Py 3
except ImportError:
    from urllib import urlencode  # Py 2


class SBWlink(object):
    """Generic access for local SBW services using SBWPython """

    sbw = None
    psbw = None
    sbwModuleProxy = None
    moduleDict = None
    modules = None

    def __init__(self):
        try:
            import SBW as SBW
            import SBW.psbw as psbw

            ##  reload(SBW)
            ##  reload(psbw)
            self.sbw = SBW
            self.psbw = SBW.psbw
            self.sbwModuleProxy = SBW.sbwModuleProxy
            self.moduleDict = SBW.sbwModuleProxy.moduleDict
            self.modules = []
            for m in self.moduleDict:
                if self.moduleDict[m].pythonName not in ['python']:
                    self.SBW_exposeAll(self.moduleDict[m])
                    self.modules.append(self.moduleDict[m].pythonName)
                    setattr(self, self.moduleDict[m].pythonName, self.moduleDict[m])
            print('\nSBWlink established.')
        except Exception as ex:
            print(ex)
            print('\nSBWlink not established.')

    def SBW_exposeAll(self, module):
        for s in module.services:
            s = getattr(module, s)
            for m in s.methods:
                getattr(s, m)

    def SBW_getActiveModules(self):
        idlist = []
        namelst = []
        for id in self.psbw.getModuleIdList():
            idlist.append(id)
            namelst.append(self.psbw.getModuleName(id))
        for id in list(self.moduleDict.keys()):
            if id not in idlist:
                self.moduleDict.pop(id)
        for name in range(len(self.modules) - 1, -1, -1):
            if self.modules[name] not in namelst:
                delattr(self, self.modules[name])
                self.modules.pop(name)
        for name in namelst:
            if name not in self.modules:
                self.SBW_loadModule(name)
        return namelst

    def SBW_loadModule(self, module_name):
        ans = 'Y'
        if module_name[-3:] == 'GUI':
            ans = input(
                'Warning! This may hang the console\n\yPress \'Y\' to continue: '
            )
        if ans == 'Y':
            module_id = self.psbw.SBWGetModuleInstance(module_name)
            assert module_id != None, '\nUnknow module, %s' % module_name
            module = self.sbwModuleProxy.ModuleProxy(module_id)
            self.SBW_exposeAll(module)
            if module_id not in self.moduleDict:
                print(
                    '<PySCeS_SBW> Adding '
                    + module.pythonName
                    + ' to ModuleProxy (id='
                    + str(module_id)
                    + ')'
                )
                self.moduleDict.update({module_id: module})
            if module.pythonName not in self.modules:
                print('<PySCeS_SBW> Adding ' + module.pythonName + ' to SBWlink')
                self.modules.append(module.pythonName)
                setattr(self, module.pythonName, module)
        else:
            print('\nModule %s not loaded' % module_name)


class SBWLayoutWebLink(object):
    """Enables access to DrawNetwork and SBMLLayout web services at www.sys-bio.org"""

    sbwhost = '128.208.17.26'
    sbml = None
    sbmllayout = None
    svg = None
    DEBUGMODE = False
    DEBUGLEVEL = 1
    DRAWNETWORKLOADED = False
    LAYOUTMODULELOADED = False

    def setProxy(self, **kwargs):
        """Set as many proxy settings as you need. You may supply a user name without
        a password in which case you will be prompted to enter one (once) when required
        (NO guarantees, implied or otherwise, on password security AT ALL). Arguments can be:

        user = 'daUser',
        pwd = 'daPassword',
        host = 'proxy.paranoid.net',
        port = 3128
        """
        proxy_info = {}
        for k in list(kwargs.keys()):
            proxy_info.update({k: kwargs[k]})

        if 'user' in proxy_info and 'pwd' not in proxy_info:
            proxy_info.update({'pwd': getpass.getpass()})
        proxy_support = ProxyHandler(
            {"http": "http://%(user)s:%(pwd)s@%(host)s:%(port)d" % proxy_info}
        )
        opener = build_opener(proxy_support, HTTPHandler)
        install_opener(opener)
        del proxy_info, proxy_support

    def loadSBMLFileFromDisk(self, File, Dir=None):
        if Dir != None:
            path = os.path.join(Dir, File)
        else:
            path = File
        if os.path.exists(path):
            self.sbmllayout = None
            self.svg = None
            self.DRAWNETWORKLOADED = False
            self.LAYOUTMODULELOADED = False
            sbmlF = open(path, 'r')
            self.sbml = sbmlF.read()
            sbmlF.close()
            return True
        else:
            print("%s is an invalid path" % path)
            return False

    def loadSBMLFromString(self, str):
        self.sbmllayout = None
        self.svg = None
        self.DRAWNETWORKLOADED = False
        self.LAYOUTMODULELOADED = False
        self.sbml = str
        return True

    def urlGET(self, host, urlpath):
        url = 'http://%s%s' % (host, urlpath)
        con = urlopen(url)
        resp = con.read()
        if self.DEBUGMODE:
            print(con.headers)
        if self.DEBUGMODE and self.DEBUGLEVEL == 2:
            print(resp)
        con.close()
        return resp

    def urlPOST(self, host, urlpath, data):
        assert type(data) == dict, '\nData must be a dictionary'
        url = 'http://%s%s' % (host, urlpath)
        con = urlopen(url, urlencode(data))
        resp = con.read()
        if self.DEBUGMODE:
            print(con.headers)
        if self.DEBUGMODE and self.DEBUGLEVEL == 2:
            print(resp)
        con.close()
        return resp

    def getVersion(self):
        print('Inspector.getVersion()')
        ver = self.urlGET(self.sbwhost, '/generate/Inspector.asmx/getVersion')
        ver = ver.replace('<?xml version="1.0" encoding="utf-8"?>', '')
        ver = ver.replace('<string xmlns="http://www.sys-bio.org/">', '')
        ver = ver.replace('</string>', '')
        return ver

    def drawNetworkLoadSBML(self):
        print('DrawNetwork.loadSBML()')
        assert self.sbml != None, '\nNo SBML file loaded'
        data = {'var0': self.sbml}
        self.DRAWNETWORKLOADED = True
        return self.urlPOST(self.sbwhost, '/generate/DrawNetwork.asmx/loadSBML', data)

    def drawNetworkGetSBMLwithLayout(self):
        print('DrawNetwork.getSBML()')
        assert self.DRAWNETWORKLOADED, '\nSBML not loaded into DrawNetwork module'
        sbml = self.urlGET(self.sbwhost, '/generate/DrawNetwork.asmx/getSBML')
        sbml = sbml.replace('&gt;', '>')
        sbml = sbml.replace('&lt;', '<')
        sbml = sbml.replace(
            '''<string xmlns="http://www.sys-bio.org/"><?xml version="1.0" encoding="utf-8"?>''',
            '',
        )
        sbml = sbml.replace('</string>', '')
        self.sbmllayout = sbml
        return True

    def layoutModuleLoadSBML(self):
        print('SBMLLayoutModule.loadSBML()')
        assert self.sbmllayout != None, '\nNo SBML Layout loaded'
        data = {'var0': self.sbmllayout}
        self.LAYOUTMODULELOADED = True
        return self.urlPOST(
            self.sbwhost, '/generate/SBMLLayoutModule.asmx/loadSBML', data
        )

    def layoutModuleGetSVG(self):
        assert self.LAYOUTMODULELOADED, '\nSBML not loaded into SBMLLayout module'
        svg = self.urlGET(self.sbwhost, '/generate/SBMLLayoutModule.asmx/getSVG')
        svg = svg.replace('&gt;', '>')
        svg = svg.replace('&lt;', '<')
        svg = svg.replace('''<string xmlns="http://www.sys-bio.org/">''', '')
        svg = svg.replace('''<?xml version="1.0" encoding="utf-8"?>''', '')
        svg = svg.replace('</string>', '')
        self.svg = svg
        return True

    def getSBML(self):
        return self.sbml

    def getSBMLlayout(self):
        return self.sbmllayout

    def getSVG(self):
        return self.svg


class METATOOLlink(object):
    """New interface to METATOOL binaries"""

    __metatool_path__ = None
    __mod__ = None
    __emode_exe_int__ = None
    __emode_exe_dbl__ = None
    __emode_intmode__ = 0
    __emode_userout__ = 0
    __emode_file__ = None
    __metatool_file__ = None
    # EModes = ''

    def __init__(self, mod, __metatool_path__=None):
        # Initialise elementary modes
        self.__mod__ = mod
        if __metatool_path__ == None:
            self.__metatool_path__ = os.path.join(mod.__pysces_directory__, 'metatool')
        else:
            self.__metatool_path__ = os.path.join(__metatool_path__, 'metatool')
        assert self.__metatool_path__ != None, '\nPySCeS not found'

        self.__emode_file__ = self.__mod__.ModelFile[:-4] + '_emodes'
        self.__metatool_file__ = self.__mod__.ModelFile[:-4] + '_metatool'
        if os.sys.platform == 'win32':
            self.__emode_exe_int__ = os.path.join(
                self.__metatool_path__, 'meta43_int.exe'
            )
            self.__emode_exe_dbl__ = os.path.join(
                self.__metatool_path__, 'meta43_double.exe'
            )
        else:
            self.__emode_exe_int__ = os.path.join(self.__metatool_path__, 'meta43_int')
            self.__emode_exe_dbl__ = os.path.join(
                self.__metatool_path__, 'meta43_double'
            )

        if os.path.exists(self.__emode_exe_int__):
            print('Using METATOOL int', end=' ')
            self.__emode_intmode__ = True
        else:
            self.__emode_exe_int__ = None
        if os.path.exists(self.__emode_exe_dbl__):
            print('\b\b\b\bdbl')
            self.__emode_intmode__ = False
        else:
            self.__emode_exe_dbl__ = None
        assert (
            self.__emode_exe_dbl__ != None or self.__emode_exe_int__ != None
        ), "\nMETATOOL binaries not available"

    def doEModes(self):
        """
        doEModes()

        Calculate the elementary modes by way of an interface to MetaTool.

        METATOOL is a C program developed from 1998 to 2000 by Thomas Pfeiffer (Berlin)
        in cooperation with Stefan Schuster and Ferdinand Moldenhauer (Berlin) and Juan Carlos Nuno (Madrid).
        http://www.biologie.hu-berlin.de/biophysics/Theory/tpfeiffer/metatool.html

        Arguments:
        None

        """
        print(
            'METATOOL is a C program developed from 1998 to 2000 by Thomas Pfeiffer (Berlin)'
        )
        print(
            'in cooperation with Stefan Schuster and Ferdinand Moldenhauer (Berlin) and Juan Carlos Nuno (Madrid).'
        )
        print(
            'http://www.biologie.hu-berlin.de/biophysics/Theory/tpfeiffer/metatool.html'
        )

        goMode = 0
        fileIn = 'pysces_metatool.dat'
        fileOut = 'pysces_metatool.out'

        goMode = 1
        if goMode == 1:
            # Build MetaTool input file
            File = open(os.path.join(self.__mod__.ModelOutput, fileIn), 'w')
            # Determine type of reaction
            out1 = []
            for key in self.__mod__.__nDict__:
                # print key
                # print self.__mod__.__nDict__[key]['Type']
                out1.append((key, self.__mod__.__nDict__[key]['Type']))
            # print '\nExtracting metatool information from network dictionary ...\n'
            File.write('-ENZREV\n')
            for x in out1:
                if x[1] == 'Rever':
                    File.write(x[0] + ' ')
            File.write('\n\n')
            File.write('-ENZIRREV\n')
            for x in out1:
                if x[1] == 'Irrev':
                    File.write(x[0] + ' ')
            File.write('\n\n')
            File.write('-METINT\n')
            for x in self.__mod__.__species__:
                File.write(x + ' ')
            File.write('\n\n')
            File.write('-METEXT\n')
            for x in self.__mod__.__fixed_species__:
                File.write(x + ' ')
            File.write('\n\n')

            output = []
            allInt = 1
            for x in self.__mod__.__nDict__:
                reList = self.__mod__.__nDict__[x]['Reagents']
                subs = ''
                prods = ''
                # print 'Reaction: ' + x
                for y in reList:
                    if self.__emode_intmode__ == 1:  # use int elementary modes
                        if abs(int(reList[y])) / abs(float(reList[y])) != 1.0:
                            print('INFO: Coefficient not integer = ' + repr(reList[y]))
                            allInt = 0
                        if reList[y] < 0:
                            # print y.replace('self.','') + ' : substrate'
                            if abs(int(reList[y])) != 1:
                                subs += repr(abs(int(reList[y]))) + ' '
                            subs += y.replace('self.', '')
                            subs += ' + '
                        else:
                            # print y.replace('self.','') + ' : product '
                            if abs(int(reList[y])) != 1:
                                prods += repr(abs(int(reList[y]))) + ' '
                            prods += y.replace('self.', '')
                            prods += ' + '
                        # output.append(x + ' : ' + subs[:-3] + ' = ' + prods[:-3] + ' .')
                    else:  # use float/double elementary mode
                        if reList[y] < 0.0:
                            # print y.replace('self.','') + ' : substrate'
                            if abs(float(reList[y])) != 1.0:
                                subs += repr(abs(float(reList[y]))) + ' '
                            subs += y.replace('self.', '')
                            subs += ' + '
                        else:
                            # print y.replace('self.','') + ' : product '
                            if abs(float(reList[y])) != 1.0:
                                prods += repr(abs(float(reList[y]))) + ' '
                            prods += y.replace('self.', '')
                            prods += ' + '
                output.append(x + ' : ' + subs[:-3] + ' = ' + prods[:-3] + ' .')

            File.write('-CAT\n')
            for x in output:
                File.write(x + '\n')
            File.write('\n')

            File.flush()
            File.close()

            if allInt == 1:
                if self.__emode_intmode__ == 1:
                    eModeExe = self.__emode_exe_int__
                else:
                    eModeExe = self.__emode_exe_dbl__

                print('\nMetatool running ...\n')
                ######### UPDATE:
                # Actually works fine on windows and posix - johann 20081128
                print('Generic run')
                os.spawnl(
                    os.P_WAIT,
                    eModeExe,
                    eModeExe,
                    os.path.join(self.__mod__.ModelOutput, fileIn),
                    os.path.join(self.__mod__.ModelOutput, fileOut),
                )
                print('\nMetatool analysis complete\n')

                # Parse MetaTool output file and store the result in a string
                go = 0
                go2 = 0
                result = ''
                end = ''

                try:
                    file2 = open(os.path.join(self.__mod__.ModelOutput, fileOut), 'r')
                    for line in file2:
                        c = re.match('ELEMENTARY MODES', line)
                        d = re.match(' enzymes', line)
                        e = re.match('The elementary mode', line)
                        f = re.match('\n', line)
                        g = re.match('The elementary', line)
                        if c != None:
                            go = 1
                            go2 = 0
                        if d != None:
                            go2 = 1
                        if e != None:
                            go2 = 0
                        if go == 1 and go2 == 1 and f == None:
                            line = line.replace('reversible', '\n  reversible\n')
                            line = line.replace('ir\n  ', '\n  ir')
                            if self.__emode_intmode__ == 1:
                                line = line.replace('] ', ']\n ')
                            else:
                                line = line.replace(') ', ')\n ', 1)
                            result += line
                        if go == 1 and g != None:
                            end += line
                    result += end
                    result += '\n'

                    file2.close()

                    if self.__emode_userout__ == 1:
                        fileo = open(
                            os.path.join(
                                self.__mod__.ModelOutput, self.__metatool_file__
                            )
                            + '.in',
                            'w',
                        )
                        filer = open(
                            os.path.join(self.__mod__.ModelOutput, fileIn), 'r'
                        )
                        for line in filer:
                            fileo.write(line)
                        fileo.write('\n\n')
                        filer.close()
                        fileo.close()
                        filer = open(
                            os.path.join(self.__mod__.ModelOutput, fileOut), 'r'
                        )
                        fileo = open(
                            os.path.join(
                                self.__mod__.ModelOutput, self.__metatool_file__
                            )
                            + '.out',
                            'w',
                        )
                        for line in filer:
                            fileo.write(line)
                        filer.close()
                        fileo.close()
                    os.remove(os.path.join(self.__mod__.ModelOutput, fileIn))
                    os.remove(os.path.join(self.__mod__.ModelOutput, fileOut))
                except Exception as EX:
                    print('doEmode:', EX)
                    print(
                        'WARNING: Unable to open MetaTool output file\nPlease check the MetaTool executables: '
                    )
                    if os.name == 'posix':
                        print(
                            '/MetaTool/meta43_double /MetaTool/meta43_int\nand their permissions'
                        )
                    else:
                        print('/MetaTool/meta43_double.exe /MetaTool/meta43_int.exe')
            else:
                print(
                    '\nINFO: non-integer coefficients\
                \nTry using the double eMode function: self.__emode_intmode__=0'
                )
                result = 'Elementary modes not calculated\n'
        else:
            print('\nNo elementary mode calculation possible - no meta43_xxx.exe')
            result = 'Elementary modes not calculated\n'
        self.EModes = result

    def getEModes(self):
        """
        getEModes()

        Returns the elementary modes as a linked list of fluxes

        """
        try:
            a = self.EModes
            FF = io.StringIO()
            FF.write(self.EModes)
            FF.reset()
            output = []
            for line in FF:
                if (
                    re.match('  ', line)
                    and not re.match('  reversible', line)
                    and not re.match('  irreversible', line)
                ):
                    tmp = [el for el in line.replace('\n', '').split(' ') if el != '']
                    tmpOut = []
                    skip = False
                    for el in range(len(tmp)):
                        if skip:
                            skip = False
                        elif tmp[el][0] != '(':
                            tmpOut.append(tmp[el])
                        elif tmp[el][0] == '(':
                            tmpOut.append(tmp[el] + ')' + tmp[el + 1][:-1])
                            skip = True
                    output.append(tmpOut)
            return output
        except AttributeError as atx:
            print(atx)
            print('\nINFO: Please run doEModes() first\n')

    def showEModes(self, File=None):
        """
        showEModes(File=None)

        Print the results of an elementary mode analysis, generated with doEModes(),
        to screen or file.

        Arguments:
        File [default=None]: Boolean, if True write parsed elementary modes to file

        """
        try:
            if File != None:
                # assert type(File) == file, 'showEmodes() needs an open file object'
                print('\nElementary modes written to file\n')
                f = open(
                    os.path.join(
                        self.__mod__.ModelOutput, self.__emode_file__ + '.out'
                    ),
                    'w',
                )
                f.write('\n## Elementary modes\n')
                f.write(self.EModes)
                f.close()
            else:
                print('\nElementary modes\n')
                print(self.EModes)
        except AttributeError as atx:
            print(atx)
            print('\nINFO: Please run doEModes() first\n')


'''
_HAVE_STOMPY = False
_STOMPY_LOAD_ERROR = ''
try:
    ##  import stompy
    import stochpy as stompy
    _HAVE_STOMPY = True
except Exception, ex:
    _STOMPY_LOAD_ERROR = '%s' % ex
    _HAVE_STOMPY = False


class StomPyInterface(object):
    """
    StomPy interface to PySCeS this may move to pysces.link in the future
    """
    SSA = None
    SSA_REACTIONS = None
    SSA_SPECIES = None
    stompy = None
    _MOD2PSC = None

    TMP_FNAME = None
    TMP_PATH = None
    MODEL_PATH = None
    OUTPUT_PATH = None

    STP_IS_TIME_SIM = False
    STP_METHOD = 'Direct'
    STP_TIMEEND = 1
    STP_TRAJ = 1
    STP_INTERACTIVE = True
    STP_TRACK_PROPENSITIES = True
    STP_WAITING_TIMES = True
    STP_STEPS = 10
    STP_INITIAL_SPECIES = True
    STP_KEEP_PSC_FILES = False

    def __init__(self, model_path, output_path):
        """
        An interface class to the StomPy stochastic simulator

         - *model_path* the default PySCeS model directory
         - *output_path* the default PySCeS output directory
        """
        self.stompy = stompy
        self.OUTPUT_PATH = output_path
        self.MODEL_PATH = model_path
        self.TMP_PATH = os.path.join(model_path, 'orca')
        self._MOD2PSC = interface.writeMod2PSC


    def setProperty(self, **kwargs):
        """
        Sets a StomPy simulation parameter

         - *method* [default='Direct'] select simulation algorithm
         - *trajectories* [default=1]
         - *interactive* [default=True]
         - *track_propensities* [default=True]
         - *steps* [default=10]

        """
        ##  print kwargs

        if kwargs.has_key('method'):
            self.STP_METHOD = kwargs['method']
            ##  print '%s = %s' % ('method', kwargs['method'])
        if kwargs.has_key('trajectories'):
            self.STP_TRAJ = kwargs['trajectories']
            self.STP_TRAJ = 1 # TODO I need to look into this
            ##  print 'Currently only single trajectories are supported via the PySCeS interface'
            ##  print '%s = %s' % ('trajectories', self.STP_TRAJ)
        if kwargs.has_key('interactive'):
            self.STP_INTERACTIVE = kwargs['interactive']
            ##  print '%s = %s' % ('interactive', kwargs['interactive'])
        if kwargs.has_key('track_propensities'):
            self.STP_TRACK_PROPENSITIES = kwargs['track_propensities']
            ##  print '%s = %s' % ('track_propensities', kwargs['track_propensities'])
        if kwargs.has_key('steps'):
            self.STP_STEPS = kwargs['steps']
            ##  print '%s = %s' % ('steps', kwargs['steps'])
        if kwargs.has_key('species_initial'):
            self.STP_INITIAL_SPECIES = kwargs['initial_species']
            ##  print '%s = %s' % ('initial_species', kwargs['initial_species'])
        if kwargs.has_key('keep_psc_files'):
            self.STP_KEEP_PSC_FILES = kwargs['keep_psc_files']
            ##  print '%s = %s' % ('keep_psc_files', kwargs['keep_psc_files'])

    def initModelFromMod(self, pscmod, iValues=False):
        """
        Initialise a StomPy SSA instance from a PySCeS model.

         - *pscmod* an initialised PySCeS model
         - *iValues* [default=False] use initial values (not current)

        """
        self.TMP_FNAME = str(time.time()).replace('.','')+'.psc'

        if self.STP_INITIAL_SPECIES:
            for s in pscmod.species:
                setattr(pscmod, s, pscmod.__sDict__[s]['initial'])

        self._MOD2PSC(pscmod, self.TMP_FNAME, self.TMP_PATH, iValues=iValues)

        self.SSA = self.stompy.SSA(Method=self.STP_METHOD, File=self.TMP_FNAME, dir=self.TMP_PATH, IsInteractive=self.STP_INTERACTIVE)
        self.SSA.Trajectories(self.STP_TRAJ)
        self.SSA_REACTIONS = self.SSA.SSA.rate_names
        self.SSA_SPECIES = self.SSA.SSA.species

        if self.STP_TRACK_PROPENSITIES:
            self.SSA.TrackPropensities()
        try:
            print os.path.join(self.TMP_PATH, self.TMP_FNAME)
            if not self.STP_KEEP_PSC_FILES and self.TMP_PATH != None and self.TMP_FNAME != None:
                os.remove(os.path.join(self.TMP_PATH, self.TMP_FNAME))
        except:
            print 'Could not delete intermediatery StomPy PSC file: %s' % os.path.join(self.TMP_PATH, self.TMP_FNAME)
        self.TMP_FNAME = None
        print 'StomPy model ... initialised.'

    def runTimeSimulation(self, pscmod, endtime=None, method='Direct', iValues=False):
        """
        Run a stochastic simulation

         - *pscmod* and instanitiated PySCeS model
         - *endtime* [default=1] the end time **Note: this could take a long time i.e. generate ahuge amount of steps**
         - *method* [default='Direct'] select the simulation method, one of:

          - *Direct*
          - *FirstReactionMethod*
          - *NextReactionMethod*
          - *TauLeaping*

         - *iValues* [default=False] use initial values (not current)

        """

        if method not in ['Direct','FirstReactionMethod','NextReactionMethod','TauLeaping']:
            print 'Method: %s does not exist using - Direct' % method
            self.STP_METHOD = 'Direct'
        else:
            self.STP_METHOD = method
        if endtime != None:
            self.STP_TIMEEND = endtime
        self.initModelFromMod(pscmod, iValues=iValues)

        ##  self.SSA.Timesteps(self.STP_STEPS)
        self.SSA.Endtime(self.STP_TIMEEND)
        self.SSA.Run()
        self.STP_IS_TIME_SIM = True
        ##  self.SSA.PlotTimeSim()

        print 'StomPy time simulation ... done.'

        # TODO STOCHPY

        ##  if self.SSA.SSA.output[0][-1] == '':
            ##  self.SSA.SSA.output[0][-1] = 0.0
        ##  sim_dat = numpy.array(self.SSA.SSA.output, 'd')

        ##  pscmod.data_stochsim = IntegrationStochasticDataObj()
        ##  pscmod.data_stochsim.setTime(sim_dat[:,0])
        ##  pscmod.data_stochsim.setSpecies(sim_dat[:,1:-1], self.SSA_SPECIES)

        pscmod.data_stochsim = self.SSA.data_stochsim

        if self.STP_WAITING_TIMES:
            wtimes, wt_lbls = self.getWaitingtimesData(reactions=None,lbls=True)
            pscmod.data_stochsim.setWaitingtimes(wtimes, wt_lbls)
        if self.STP_TRACK_PROPENSITIES:
            pscmod.data_stochsim.setPropensities(self.SSA.SSA.propensities_output)
        pscmod.data_stochsim.TYPE_INFO = 'Stochastic'

    def runStepSimulation(self, pscmod, steps=None, method='Direct', iValues=False):
        """
        Run a stochastic simulation

         - *pscmod* and instanitiated PySCeS model
         - *steps* [default=10] the number of steps to simulate
         - *method* [default='Direct'] select the simulation method, one of:

          - *Direct*
          - *FirstReactionMethod*
          - *NextReactionMethod*
          - *TauLeaping*

        - *iValues* [default=False] use initial values (not current)

        """

        if method not in ['Direct','FirstReactionMethod','NextReactionMethod','TauLeaping']:
            print 'Method: %s does not exist using - Direct' % method
            self.STP_METHOD = 'Direct'
        else:
            self.STP_METHOD = method
        if steps != None:
            self.STP_STEPS = steps
        self.initModelFromMod(pscmod, iValues=iValues)

        self.SSA.Timesteps(self.STP_STEPS)
        ##  self.SSA.Endtime(self.STP_TIMEEND)
        self.SSA.Run()
        self.STP_IS_TIME_SIM = False

        print 'StomPy step simulation ... done.'
        ##  print self.SSA.SSA.output[0]
        ##  print self.SSA.SSA.output[1]
        ##  print self.SSA.SSA.output[-1]
        ##  header_line = self.SSA.SSA.output.pop(0)
        ##  if self.SSA.SSA.output[0][-1] == '':
            ##  self.SSA.SSA.output[0][-1] = 0.0
        ##  sim_dat = numpy.array(self.SSA.SSA.output, 'd')

        ##  pscmod.data_stochsim = IntegrationStochasticDataObj()
        ##  pscmod.data_stochsim.setTime(sim_dat[:,0])
        ##  pscmod.data_stochsim.setSpecies(sim_dat[:,1:-1], self.SSA_SPECIES)

        pscmod.data_stochsim = self.SSA.data_stochsim

        if self.STP_WAITING_TIMES:
            wtimes, wt_lbls = self.getWaitingtimesData(reactions=None,lbls=True)
            pscmod.data_stochsim.setWaitingtimes(wtimes, wt_lbls)
        if self.STP_TRACK_PROPENSITIES:
            pscmod.data_stochsim.setPropensities(self.SSA.SSA.propensities_output)
        pscmod.data_stochsim.TYPE_INFO = 'Stochastic'


    def getWaitingtimesData(self,reactions=None,lbls=False):
        """
        Plots the waiting times for each reaction in the model. Makes use of ObtainWaitingtimes to derive the waiting times out of the SSA output.

        Input:

         - *reactions* [default=0] a list of reactions to plot defualts to all reactions
         - *traj* [default=0] trajectory to plot (defaults to first one)
         - *lbls* [default=False] if True return (data_array, column_labels) otherwise just data_array

        This method is derived from StomPy 0.9 (http://stompy.sf.net) Analysis.py

        """
        if self.SSA.IsTauLeaping:
            print 'INFO: waiting times not available when method is Tau Leaping'
            if not lbls:
                return None
            else:
                return None, None
        self.SSA.GetWaitingtimes()
        if reactions == None:
            reactions = self.SSA_REACTIONS
        vect = []
        vect_lbls = []
        for r in reactions:
            if r in self.SSA_REACTIONS:
                vect.append(self.SSA_REACTIONS.index(r)+1)
                vect_lbls.append('wt'+str(r))
            else:
                print "INFO: '%s' is not a valid reaction name" % r
        OUTPUT = []
        ##  for t in range(len(self.SSA.data_stochsim.waiting_times)):
        T_OUTPUT = []
        for i in vect:
            if self.SSA.data_stochsim.waiting_times.has_key(i):
                waiting_time = self.SSA.data_stochsim.waiting_times[i]
                if len(waiting_time) > 1:			# At least 2 waiting times are necessary per reaction
                    T_OUTPUT.append(self.stompy.modules.Analysis.LogBin(waiting_time, 1.5)) 	# Create logarithmic bins
                else:
                    T_OUTPUT.append(None)
            else:
                T_OUTPUT.append(None)
        OUTPUT.append(T_OUTPUT)
        if not lbls:
            return OUTPUT
        else:
            return OUTPUT, vect_lbls


class IntegrationStochasticDataObj(object):
    """
    This class is specifically designed to store the
    results of a stochastic time simulation
    It has methods for setting the Time, Labels, Species and Propensity data and
    getting Time, Species and Rate (including time) arrays. However, of more use:

    - getOutput(\*args) feed this method species/rate labels and it will return
      an array of [time, sp1, r1, ....]
    - getDataAtTime(time) the data generated at time point "time".
    - getDataInTimeInterval(time, bounds=None) more intelligent version of the above
      returns an array of all data points where: time-bounds <= time <= time+bounds

    """
    time = None
    waiting_times = None
    species = None
    propensities = None
    xdata = None
    time_label = 'Time'
    waiting_time_labels = None
    species_labels = None
    propensities_labels = None
    xdata_labels = None
    HAS_SPECIES = False
    HAS_WAITING_TIMES = False
    HAS_PROPENSITIES = False
    HAS_TIME = False
    HAS_XDATA = False
    IS_VALID = True
    TYPE_INFO = 'Stochastic'

    def setLabels(self, species):
        """
        Set the species

         - *species* a list of species labels
        """
        self.species_labels = species

    def setTime(self, time, lbl=None):
        """
        Set the time vector

         - *time* a 1d array of time points
         - *lbl* [default=None] is "Time" set as required

        """
        self.time = time.reshape(len(time), 1)
        self.HAS_TIME = True
        if lbl != None:
            self.time_label = lbl

    def setSpecies(self, species, lbls=None):
        """
        Set the species array

         - *species* an array of species vs time data
         - *lbls* [default=None] a list of species labels
        """
        self.species = species
        self.HAS_SPECIES = True
        if lbls != None:
            self.species_labels = lbls

    def setWaitingtimes(self, waiting_times, lbls=None):
        """
        Set the `waiting_times` this data structure is not an array but a nested list of: waiting time log bins per reaction per trajectory::

        waiting_times = [traj_1, ..., traj_n]
        traj_1 = [wt_J1, ..., wt_Jn] # in order of SSA_REACTIONS
        wt_J1 = (xval, yval, nbin)
        xval =[x_1, ..., x_n]
        yval =[y_1, ..., y_n]
        nbin = n

         - *waiting_times* a list of waiting times
         - *lbls* [default=None] a list of matching reaction names

        """
        self.waiting_times = waiting_times
        self.HAS_WAITING_TIMES = True
        if lbls != None:
            self.waiting_time_labels = lbls


    def setPropensities(self, propensities, lbls=None):
        """
        Sets an array of propensities.

         - *propensities* a list of propensities
         - *lbls* [default=None] a list of matching reaction names

        """

        if lbls == None:
            LB = copy.copy(propensities[0])
            lbls = LB[1:]
            lbls = ['p'+str(r) for r in lbls]
        P_ARR = numpy.zeros((len(propensities), len(propensities[0])-1), 'd')
        P_ARR[-1,:] = numpy.NaN
        for r in range(1, P_ARR.shape[0]):
            P_ARR[r, :] = propensities[r][1:]
        self.propensities = P_ARR
        self.HAS_PROPENSITIES = True
        if lbls != None:
            self.propensities_labels = lbls
        ##  print self.propensities_labels
        ##  print self.propensities


    def setXData(self, xdata, lbls=None):
        """
        Sets an array of extra simulation data

        - *xdata* an array of xdata vs time
        - *lbls* [default=None] a list of xdata labels

        """
        self.xdata = xdata
        self.HAS_XDATA = True
        if lbls != None:
            self.xdata_labels = lbls

    def getTime(self, lbls=False):
        """
        Return the time vector

         - *lbls* [default=False] return only the time array or optionally both the time array and time label

        """
        output = None
        if self.HAS_TIME:
            output = self.time.reshape(len(self.time),)
        if not lbls:
            return output
        else:
            return output, [self.time_label]

    def getSpecies(self, lbls=False):
        """
        Return an array fo time+species

        - *lbls* [default=False] return only the time+species array or optionally both the data array and a list of column label

        """
        output = None
        if self.HAS_SPECIES:
            output = numpy.hstack((self.time, self.species))
            labels = [self.time_label]+self.species_labels
        else:
            output = self.time
            labels = [self.time_label]
        if not lbls:
            return output
        else:
            return output, labels

    def getWaitingtimes(self, lbls=False, traj=[]):
        """
        Return waiting times, time+waiting_time array

         - *lbls* [default=False] return only the time+waiting_time array or optionally both the data array and a list of column label
         - *traj* [default=[0]] return the firs or trajectories defined in this list

        """
        output = None
        labels = None
        if self.HAS_WAITING_TIMES:
            output = []
            if len(traj) == 0:
                traj = range(len(self.waiting_times))
            ##  if len(traj) == 1:
                ##  output = self.waiting_times[0]
            ##  else:
            for t in traj:
                output.append(self.waiting_times[t])
            labels = self.waiting_time_labels
        else:
            output = []
            labels = []
        if not lbls:
            return output
        else:
            return output, labels

    def getPropensities(self, lbls=False):
        """
        Return time+propensity array

         - *lbls* [default=False] return only the time+propensity array or optionally both the data array and a list of column label

        """
        #assert self.propensities != None, "\nNo propensities"
        output = None
        if self.HAS_PROPENSITIES:
            print self.time.shape
            print self.propensities.shape
            output = numpy.hstack((self.time, self.propensities))
            labels = [self.time_label]+self.propensities_labels
        else:
            output = self.time
            labels = [self.time_label]
        if not lbls:
            return output
        else:
            return output, labels

    def getXData(self, lbls=False):
        """
        Return time+xdata array

        - *lbls* [default=False] return only the time+xdata array or optionally both the data array and a list of column label

        """
        output = None
        if self.HAS_XDATA:
            output = numpy.hstack((self.time, self.xdata))
            labels = [self.time_label]+self.xdata_labels
        else:
            output = self.time
            labels = [self.time_label]
        if not lbls:
            return output
        else:
            return output, labels

    def getDataAtTime(self, time):
        """
        Return all data generated at "time"

         - *time* the required exact time point

        """
        #TODO add rate rule data
        t = None
        sp = None
        ra = None
        ru = None
        xd = None
        temp_t = self.time.reshape(len(self.time),)
        for tt in range(len(temp_t)):
            if temp_t[tt] == time:
                t = tt
                if self.HAS_SPECIES:
                    sp = self.species.take([tt], axis=0)
                if self.HAS_PROPENSITIES:
                    ru = self.propensities.take([tt], axis=0)
                if self.HAS_XDATA:
                    xd = self.xdata.take([tt], axis=0)
                break

        output = None
        if t is not None:
            output = numpy.array([[temp_t[t]]])
            if sp is not None:
                output = numpy.hstack((output,sp))
            if ra is not None:
                output = numpy.hstack((output,ra))
            if ru is not None:
                output = numpy.hstack((output,ru))
            if xd is not None:
                output = numpy.hstack((output,xd))

        return output

    def getDataInTimeInterval(self, time, bounds=None):
        """
        Returns an array of all data in interval: time-bounds <= time <= time+bounds
        where bound defaults to stepsize

         - *time* the interval midpoint
         - *bounds* [default=None] interval halfspan defaults to stepsize

        """

        temp_t = self.time.reshape(len(self.time),)
        if bounds == None:
            bounds = temp_t[1] - temp_t[0]
        c1 = (temp_t >= time-bounds)
        c2 = (temp_t <= time+bounds)
        print 'Searching (%s:%s:%s)' % (time-bounds, time, time+bounds)

        t = []
        sp = None
        ra = None

        for tt in range(len(c1)):
            if c1[tt] and c2[tt]:
                t.append(tt)
        output = None
        if len(t) > 0:
            output = self.time.take(t)
            output = output.reshape(len(output),1)
            if self.HAS_SPECIES and self.HAS_TIME:
                output = numpy.hstack((output, self.species.take(t, axis=0)))
            if self.HAS_PROPENSITIES:
                output = numpy.hstack((output, self.propensities.take(t, axis=0)))
            if self.HAS_XDATA:
                output = numpy.hstack((output, self.xdata.take(t, axis=0)))
        return output

    def getAllSimData(self,lbls=False):
        """
        Return an array of time + all available simulation data

         - *lbls* [default=False] return only the data array or (data array, list of labels)

        """
        labels = [self.time_label]
        if self.HAS_SPECIES and self.HAS_TIME:
            output = numpy.hstack((self.time, self.species))
            labels += self.species_labels
        if self.HAS_PROPENSITIES:
            output = numpy.hstack((output, self.propensities))
            labels += self.propensities_labels
        if self.HAS_XDATA:
            output = numpy.hstack((output, self.xdata))
            labels += self.xdata_labels
        if not lbls:
            return output
        else:
            return output, labels

    def getSimData(self, *args, **kwargs):
        """
        Feed this method species/xdata labels and it will return an array of [time, sp1, ....]

         - 'speces_l', 'xdatal' ...
         - *lbls* [default=False] return only the data array or (data array, list of labels)

        """
        output = self.time

        if kwargs.has_key('lbls'):
            lbls = kwargs['lbls']
        else:
            lbls = False
        lout = [self.time_label]
        for roc in args:
            if self.HAS_SPECIES and roc in self.species_labels:
                lout.append(roc)
                output = numpy.hstack((output, self.species.take([self.species_labels.index(roc)], axis=-1)))
            if self.HAS_PROPENSITIES and roc in self.propensities_labels:
                lout.append(roc)
                output = numpy.hstack((output, self.propensities.take([self.propensities_labels.index(roc)], axis=-1)))
            if self.HAS_XDATA and roc in self.xdata_labels:
                lout.append(roc)
                output = numpy.hstack((output, self.xdata.take([self.xdata_labels.index(roc)], axis=-1)))
        if not lbls:
            return output
        else:
            return output, lout


class PysMod:
#STOMPY INSERT START
    def StochSimPlot(self, plot='species', filename=None, title=None, log=None, format='points'):
        """
        Plot the Stochastic simulation results, uses the new UPI pysces.plt interface:

        - *plot* [default='species'] output to plot, can be one of:

         - 'all' species and propensities
         - 'species' species
         - 'waiting_times' waiting_times
         - 'propensities' propensities
         - `['S1', 'R1', ]` a list of model attributes ('species')

        - *filename* [default=None] if defined file is exported to filename
        - *title* [default=None] the plot title
        - *log* [default=None] use log axis for axis 'x', 'y', 'xy'
        - *format* [default='lines'] use UPI or backend specific keys
        """
        data = None
        labels = None
        allowedplots = ['all', 'species', 'propensities','waiting_times']
        ##  allowedplots = ['all', 'species', 'waiting_times']
        if type(plot) != list and plot not in allowedplots:
            raise RuntimeError, '\nPlot must be one of %s not \"%s\"' % (str(allowedplots), plot)
        if plot == 'all':
            ##  data, labels = self.data_stochsim.getSpecies(lbls=True)
            data, labels = self.data_stochsim.getAllSimData(lbls=True)
        elif plot == 'species':
            data, labels = self.data_stochsim.getSpecies(lbls=True)
        elif plot == 'propensities':
            data, labels = self.data_stochsim.getPropensities(lbls=True)
            ##  data, labels = self.data_stochsim.getRates(lbls=True)
        elif plot == 'waiting_times':
            dataLst, labels = self.data_stochsim.getWaitingtimes(lbls=True)
            format='points'
            ##  data, labels = self.data_stochsim.getRates(lbls=True)
        else:
            plot = [at for at in plot if at in self.__species__+[self.data_stochsim.time_label]+self.data_stochsim.propensities_labels]
            kwargs = {'lbls' : True}
            print plot
            if len(plot) > 0:
                data, labels = self.data_stochsim.getSimData(*plot, **kwargs)
        del allowedplots

        xu = 'Time (%(multiplier)s x %(kind)s x 10**%(scale)s)**%(exponent)s' % self.__uDict__['time']

        if plot == 'waiting_times':
            xu = 'Inter-arrival time (%s)' % xu

            xrng_start = 0.1
            xrng_end = 0.1
            yrng_start = 0.1
            yrng_end = 0.1
            for wt in range(len(dataLst)):
                for d in range(len(dataLst[wt])):
                    D = dataLst[wt][d]
                    if plt.__USE_MATPLOTLIB__ and d > 0:
                        plt.m.hold(True)
                    if D != None and len(D[0]) > 0 and len(D[1]) > 0:
                        data = numpy.vstack([D[0], D[1]]).transpose()
                        if min(D[0]) < xrng_start and min(D[0]) > 0.0:
                            xrng_start = min(D[0])
                        if max(D[0]) > xrng_end:
                            xrng_end = max(D[0])
                        if min(D[1]) < yrng_start and min(D[1]) > 0.0:
                            yrng_start = min(D[1])
                        if max(D[1]) > yrng_end:
                            yrng_end = max(D[1])
                        plt.plotLines(data, 0, [1], titles=['Time']+[labels[d]], formats=[format])
            plt.setRange('x', xrng_start*0.8, xrng_end*1.2)
            plt.setRange('y', yrng_start*0.8, yrng_end*1.2)
            if plt.__USE_MATPLOTLIB__:
                plt.m.hold(False)
        else:
            plt.plotLines(data, 0, range(1, data.shape[1]), titles=labels, formats=[format])
            # set the x-axis range so that it is original range + 0.2*sim_end
            # this is a sceintifcally dtermned amount of space that is needed for the title at the
            # end of the line :-) - brett 20040209
            RngTime = self.data_stochsim.getTime()
            end = RngTime[-1] + 0.2*RngTime[-1]
            plt.setRange('x', RngTime[0], end)
            del RngTime

        # For now StochPy results are plotted as Amounts directly from StochPy
        M = 'Amount'
        ##  if self.__KeyWords__['Output_In_Conc']:
            ##  M = 'Concentration'
        ##  else:
            ##  M = 'Amount (%(multiplier)s x %(kind)s x 10**%(scale)s)**%(exponent)s' % self.__uDict__['substance']

        if plot == 'all':
            yl = 'Amount, propensities'
        elif plot == 'propensities':
            yl = 'Propensities'
        elif plot == 'waiting_times':
            yl = 'Frequency'
            if log == None:
                log = 'xy'
        elif plot == 'species':
            yl = '%s' % M
        else:
            yl = 'User defined'
        plt.setAxisLabel('x', xu)
        plt.setAxisLabel('y', yl)
        if log != None:
            plt.setLogScale(log)
        if title == None:
            plt.setGraphTitle('PySCeS/StochPy simulation (' + self.ModelFile + ') ' + time.strftime("%a, %d %b %Y %H:%M:%S"))
        else:
            plt.setGraphTitle(title)
        plt.replot()
        if filename != None:
            plt.export(filename, directory=self.ModelOutput, type='png')


    def doStochSim(self,end=10,mode='steps',method='Direct',trajectories=1):
        """
        doStochSim(end=10, mode='steps', method='Direct')

        Run a stochastic simulation for until `end` is reached. This can be either steps or end time (which could be a *HUGE* number of steps).

        Arguments:

         - *end* [default=10] simulation end (steps or time)
         - *mode* [default='steps'] simulation mode, can be one of:

          - *steps* total number of steps to simulate
          - *time* simulate until time is reached

         - *method* [default='Direct'] stochastic algorithm, can be one of:

          - Direct
          - FirstReactionMethod
          - NextReactionMethod
          - TauLeaping

        """

        if method not in ['Direct', 'FirstReactionMethod','NextReactionMethod','TauLeaping']:
            print 'Method "%s" not recognised using: "Direct"' % method
            method = 'Direct'
        if mode not in ['steps','time']:
            print 'Mode "%s" not recognised using: "steps"' % mode
            mode = 'steps'

        stompy_track_propensities = True
        stompy_keep_psc_files = False
        self.__STOMPY__.setProperty(method=method, trajectories=trajectories, interactive=True, track_propensities=stompy_track_propensities, keep_psc_files=stompy_keep_psc_files)

        if mode == 'time':
            self.__STOMPY__.runTimeSimulation(self, endtime=end, method=method)
        else:
            self.__STOMPY__.runStepSimulation(self, steps=end, method=method)


    def doStochSimPlot(self, end=10.0, mode='steps', method='Direct', plot='species', fmt='points', log=None, filename=None):
        """
        doStochSimPlot(end=10.0, mode='steps', method='Direct', plot='species', fmt='points', log=None, filename=None)

        Run a stochastic simulation for until `end` is reached and plot the results. This can be either steps or end time (which could be a *HUGE* number of steps).

        Arguments:

         - *end* [default=10] simulation end (steps or time)
         - *mode* [default='steps'] simulation mode, can be one of:

          - *steps* total number of 'steps' to simulate
          - *time* simulate until 'time' is reached

         - *method* [default='Direct'] stochastic algorithm, can be one of:

          - Direct
          - FirstReactionMethod
          - NextReactionMethod
          - TauLeaping

        - *plot* [default='species'] output to plot, can be one of:

         - 'all' species and propensities
         - 'species' species
         - 'waiting_times' waiting_times
         - 'propensities' propensities
         - `['S1', 'R1', ]` a list of model attributes ('species')

        - *filename* [default=None] if defined file is exported to filename
        - *title* [default=None] the plot title
        - *log* [default=None] use log axis for axis 'x', 'y', 'xy'
        - *fmt* [default='lines'] use UPI or backend specific keys

        """

        self.doStochSim(end=end, mode=mode, method=method,trajectories=1)
        self.StochSimPlot(plot='species', filename=filename, log=log, format=fmt)

#STOMPY INSERT START

if not _HAVE_STOMPY:
    def nofunc(self, *args, **kwargs):
        print '\nStochastic simulation not available, please download/install *StomPy* from: http://stompy.sf.net\n'
    PysMod.doStochSim = nofunc
    PysMod.doStochSimPlot = nofunc
    PysMod.StochSimPlot = nofunc

'''
