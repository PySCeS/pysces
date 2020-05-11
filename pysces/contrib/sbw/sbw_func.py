"""
This is a prototype PySCeS interface to the System Biology Workbench
(http://www.sys-bio.org). Bits of code in this module have been shamelessly scavenged 
from the SBW/Python inteface modules found in sbwStartup.py and sbwModuleProxy.py.

Essentially what I am trying to do is to get the SBW/Python interface to operate
as a module and not in the global/interactive namespace. Modules added as I 
discover more about SBW ;-) Brett G. Olivier - 20060929

Notes: notification does not work yet, new SBW modules are not autoloaded - the
interface itself checks and, if necessary, loads services as required.

"""
from __future__ import division, print_function
from __future__ import absolute_import
from __future__ import unicode_literals

try:
    import SBW
except Exception as exc:
    print('\nError loading SBW')
    print('Is the Systems Biology Workbench (www.sys-bio.org) installed?')
    print('If SBW is installed and still doesn\'t work try creating an empty file')
    print('called __init__.py and copy it into')
    print('\tpython2.x\\lib\\site-packages\\SBW\\')
    print('if this directory doesn\'t exist copy python2.x\\SBW\\ to')
    print('python2.x\\lib\\site-packages\\SBW\\ add __init__.py and try again')
    print(exc,'\n')

import os
from time import sleep
import SBW.sbwModuleProxy as sbwModuleProxy
import SBW.psbw as psbw       # give user access to psbw in interactive mode

sbwModuleProxy.importAllModules(SBW.__dict__, verbose=1)
sbwModuleProxy.startModuleNotification(SBW.__dict__, verbose=1)

class SBW_base:
    sbml_file = 'The SBML file name'
    sbml_string = 'The SBML file object as a string'
    SBW = 'The inheriting class SBW interface'
    output = 'Holds the output of the analysis'

##    def __init__(self):
##        #Load all modules the generic form of above probably should have its own function someday
##        modules = [ sbwModuleProxy.ModuleProxy(id) for id in psbw.getModuleIdList() ]
##        for module in modules:
##            if not sbwModuleProxy.moduleDict.has_key(module.id):
##                print 'Adding module to ModuleProxy', module.id
##                sbwModuleProxy.moduleDict[module.id] = module
##            if not SBW.__dict__.has_key(module.pythonName):
##                print 'Adding module to SBW namespace', module.pythonName
##                SBW.__dict__[module.pythonName] = module

    def loadSBMLfromFile(self,filename,dir=None):
        """Loads an SBML file where name=filename"""
        if dir != None:
            filename = os.path.join(dir,filename)
        self.sbml_file = filename
        sbml_file = open(filename,'r')
        self.sbml_string = ''
        for line in sbml_file: self.sbml_string += line
        sbml_file.close()

    def Info(self):
        """Display module info"""
        print(' ')
        print('Module: ' + self.SBW.getName() + ' (ver. ' + self.SBW.getVersion() + ')')
        print(' \"' + self.SBW.getDescription() + '\"')
        print('Author: ' + self.SBW.getAuthor() + ' (' + self.SBW.getURL() + ')')

    def SBW_Methods(self):
        """Display available SBW service methods"""
        print('\nMethods available as sbw.'+ self.SBW.getName()+'.SBW.*')
        for x in self.SBW.methods:
            print(x)        
        print('\n')

    def WriteOutputToFile(self,filename,dir=None):
        if dir != None:
            filename = os.path.join(dir,filename)
        Fout = open(filename,'w')
        Fout.write(self.output)
        Fout.flush()
        Fout.close()

class StructAnalysis(SBW_base):
    def __init__(self):
        #Get the module instance and load it if necessary ... i think!?
        #Note to self - find out for sure why this only works in __init__ it is a
        #namespace issue but it would be nice to work around it ;-)  - brett
        struct_module_id = psbw.SBWGetModuleInstance('edu.kgi.StructAnalysis')
        if struct_module_id not in sbwModuleProxy.moduleDict:
            print('<PySCeS_SBW> Adding ' + sbwModuleProxy.ModuleProxy(struct_module_id).pythonName + ' to ModuleProxy (id=' + str(struct_module_id) + ')')
        if sbwModuleProxy.ModuleProxy(struct_module_id).pythonName not in SBW.__dict__:
            print('<PySCeS_SBW> Adding ' + sbwModuleProxy.ModuleProxy(struct_module_id).pythonName + ' to SBW namespace')
            SBW.__dict__[sbwModuleProxy.ModuleProxy(struct_module_id).pythonName] = sbwModuleProxy.ModuleProxy(struct_module_id)

        #wang the StructAnalysis methods into the StructAnalysis instance as SBW.*
        self.SBW = SBW.edu_kgi_StructAnalysis.StructAnalysis

        # MORE COMING SOON TO AN INTERFaCE NEAR YOU ...........

    def _InitStuff(self):
        #reassign some useful stuff for easy access
        self.rank = self.SBW.getRank()

    def loadSBML(self, sbml_str=None):
        if sbml_str == None:
            sbml_str = self.sbml_string
        else:
            self.sbml_string = sbml_str
        try:
            self.output = SBW.edu_kgi_StructAnalysis.StructAnalysis.loadSBML(sbml_str)
        except Exception as ex:
            print(ex)
            self.output = str(ex)
        self._InitStuff()
        
    def Run(self, sbml_str=None):
        if sbml_str == None:
            sbml_str = self.sbml_string
        else:
            self.sbml_string = sbml_str
        try:
            self.output = SBW.edu_kgi_StructAnalysis.StructAnalysis.loadSBML(sbml_str)
        except Exception as ex:
            print(ex)
            self.output = str(ex)
        self._InitStuff()

    def RunWithTests(self, sbml_str=None):
        if sbml_str == None:
            sbml_str = self.sbml_string
        else:
            self.sbml_string = sbml_str
        try:
            self.output = SBW.edu_kgi_StructAnalysis.StructAnalysis.loadSBMLwithTests(sbml_str)
        except Exception as ex:
            self.output = str(ex)
        self._InitStuff()

class DrawNetworkGUI(SBW_base):
    def __init__(self):
        # externally spawn a GUI to avoid locking the console (uses default windows location)
        os.spawnv(os.P_NOWAIT,'c:\\Program Files\\KGI\\SBW\\Layout\\DrawNetwork.exe',[])
        print('Loading ...')
        sleep(3)
        # as we have loaded the thing already, we can just get the Id
        draw_module_id = psbw.getModuleId('DrawNetwork.GUI')
        if draw_module_id not in sbwModuleProxy.moduleDict:
            print('<PySCeS_SBW> Adding ' + sbwModuleProxy.ModuleProxy(draw_module_id).pythonName + ' to ModuleProxy (id=' + str(draw_module_id) + ')')
        if sbwModuleProxy.ModuleProxy(draw_module_id).pythonName not in SBW.__dict__:
            print('<PySCeS_SBW> Adding ' + sbwModuleProxy.ModuleProxy(draw_module_id).pythonName + ' to SBW namespace')
            SBW.__dict__[sbwModuleProxy.ModuleProxy(draw_module_id).pythonName] = sbwModuleProxy.ModuleProxy(draw_module_id)
        self.SBW = SBW.DrawNetwork_GUI.network


    def loadSBML(self,sbml_str=None):
        if sbml_str == None:
            sbml_str = self.sbml_string
        else:
            self.sbml_string = sbml_str
        try:
            self.SBW.loadSBML(sbml_str)
        except Exception as ex:
            print(ex)
            self.output = str(ex)
        print('\n\n\tClick in the DrawNetwork window to refresh display\n\n')

    def getSBMLstringWithLayout(self):
        self.output = self.SBW.getSBML()
        print("SBML+layout stored as output")

##  class Simulate3D(SBW_base):
    ##  def __init__(self):
        ##  # externally spawn a GUI to avoid locking the console (uses default windows location)
        ##  os.spawnv(os.P_NOWAIT,'c:\\Program Files\\KGI\\Simulate3D\\3Dlauncher.exe',[])
        ##  print 'Loading ...'
        ##  sleep(3)
        ##  # as we have loaded the thing already, we can just get the Id
        ##  threeD_module_id = psbw.getModuleId('Simulate3D')
        ##  if not sbwModuleProxy.moduleDict.has_key(threeD_module_id):
            ##  print '<PySCeS_SBW> Adding ' + sbwModuleProxy.ModuleProxy(threeD_module_id).pythonName + ' to ModuleProxy (id=' + str(threeD_module_id) + ')'
        ##  if not SBW.__dict__.has_key(sbwModuleProxy.ModuleProxy(threeD_module_id).pythonName):
            ##  print '<PySCeS_SBW> Adding ' + sbwModuleProxy.ModuleProxy(threeD_module_id).pythonName + ' to SBW namespace'
            ##  SBW.__dict__[sbwModuleProxy.ModuleProxy(threeD_module_id).pythonName] = sbwModuleProxy.ModuleProxy(threeD_module_id)
        ##  self.SBW = SBW.DrawNetwork_GUI.network

    ##  def loadSBML(self,sbml_str=None):
        ##  if sbml_str != None:
            ##  self.sbml_string = sbml_str
            ##  self.SBW.loadSBML(sbml_str)
        ##  else:
            ##  self.SBW.loadSBML(self.sbml_string)
        ##  print '\n\n\tClick in the DrawNetwork window to refresh display\n\n'

    ##  def getSBMLstringWithLayout(self):
        ##  self.output = self.SBW.getSBML()
        ##  print "SBML+layout stored as output"
