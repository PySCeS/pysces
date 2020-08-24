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

__doc__ = '''SBML reading/writing module - now replaced by PySCeS Core2'''

import os, sys
from time import sleep, strftime
from getpass import getuser

try:
    input = raw_input  # Py2 compatibility
except NameError:
    pass

if sys.platform == 'win32':
    try:
        import libsbml as SBML
    except Exception as e:
        print('Windows sbml load error', e)
else:
    try:
        import libsbml as SBML
    except Exception as e:
        print('Posix sbml load error', e)


class PyscesSBML:
    """The PySCeS interface to libSBML and SBML utilities"""

    mode_number_format = '%2.4e'
    sbml_level = 2
    _debug = 0
    SBML = SBML

    def SBML_buildBasicModel(
        self,
        mod,
        filename,
        slvl=2,
        dir=None,
        substance=(1, 0),
        volume=(1, 0),
        time=(1, 0),
        arules=None,
        notes=None,
    ):
        """
        SBML_buildBasicModel(mod,filename,slvl=2,dir=None)

        Create a basic SBML model.
        Arguments:
        =========
        mod: an active PySCeS model object
        filename: the output SBML file name
        slvl [default=2]: the SBML level that should be used
        dir [default=None]: the output directory
        substance [default=(1,0)]: the model substance unit - SBML default is "mole"
        volume [default=(1,0)]: the model volume unit - SBML default is "litre"
        time [default=(1,0)]: the model time unit - SBML default is "second"

        """
        self.SBML_createModel(mod, filename, slvl, dir)
        self.SBML_setCompartment()
        self.SBML_setNotes(txt=notes)
        self.SBML_setUnits(substance=substance, volume=volume, time=time)
        if arules != None:
            self.SBML_setAssignmentRules(arules)
        self.SBML_setSpecies()
        self.SBML_setReactions()
        self.SBML_setModel()

    def write(
        self,
        mod,
        filename,
        slvl=2,
        dir=None,
        substance=(1, 0),
        volume=(1, 0),
        time=(1, 0),
        arules=None,
        notes=None,
    ):
        """
        write(mod,filename,slvl=2,dir=None)

        Write a PySCeS model as an SBML file.

        Arguments:
        =========
        mod: an active PySCeS model object
        filename: the output SBML file name
        slvl [default=2]: the SBML level that should be used
        dir [default=None]: the output directory
        substance [default=(1,0)]: the model substance unit - SBML default is "mole"
        volume [default=(1,0)]: the model volume unit - SBML default is "litre"
        time [default=(1,0)]: the model time unit - SBML default is "second"

        """
        self.SBML_buildBasicModel(
            mod, filename, slvl, dir, substance, volume, time, arules, notes
        )
        self.SBML_writeFile()

    def getSBML_document(
        self, mod, substance=(1, 0), volume=(1, 0), time=(1, 0), arules=None, notes=None
    ):
        """
        Returns an SBML document object

        Arguments:
        =========
        mod: an active PySCeS model object
        substance [default=(1,0)]: the model substance unit - SBML default is "mole"
        volume [default=(1,0)]: the model volume unit - SBML default is "litre"
        time [default=(1,0)]: the model time unit - SBML default is "second"

        """

        filename = 'tempXML'
        slvl = 2
        dir = None
        self.SBML_buildBasicModel(
            mod, filename, slvl, dir, substance, volume, time, arules, notes
        )
        return self.sbml_document

    def getSBML_string(
        self, mod, substance=(1, 0), volume=(1, 0), time=(1, 0), arules=None, notes=None
    ):
        """
        Returns an SBML file as a string

        Arguments:
        =========
        mod: an active PySCeS model object
        substance [default=(1,0)]: the model substance unit - SBML default is "mole"
        volume [default=(1,0)]: the model volume unit - SBML default is "litre"
        time [default=(1,0)]: the model time unit - SBML default is "second"
        """

        filename = 'tempXML'
        slvl = 2
        dir = None
        self.SBML_buildBasicModel(
            mod, filename, slvl, dir, substance, volume, time, arules, notes
        )
        return self.sbml_document.toSBML()

    def __cleanString__(self, s):
        s = s.lstrip()
        s = s.rstrip()
        return s

    def parseForcingFunctions(self):
        self.__forcing_function_dic__ = {}
        ff = self.model_obj._Function_forced.split('\n')
        for f in ff:
            if f != '':
                f = f.split('=')
                f[0] = f[0].replace('self.', '')
                f[1] = f[1].replace('self.', '')
                self.__forcing_function_dic__.setdefault(
                    self.__cleanString__(f[0]), self.__cleanString__(f[1])
                )

    def getAssignmentRules(self):
        self.parseForcingFunctions()
        out = []
        for key in list(self.__forcing_function_dic__.keys()):
            out.append((key, self.__forcing_function_dic__[key]))
        return out

    def SBML_createModel(self, mod, filename, slvl=2, dir=None):
        """
        SBML_createModel(mod,filename,slvl=2,dir=None)

        Set up an SBML document and extract the model NetworkDict

        Arguments:
        =========
        mod: a PySCeS model object
        filename: the output filename
        slvl [default=2]: SBML level required
        dir [default=None]: output directory

        """
        if self._debug:
            print('SBML_createModel')
        self.model_obj = mod
        self.__nDict__ = self.model_obj.__nDict__
        self.model_filename = filename
        if dir == None:
            self.model_dir = self.model_obj.ModelOutput
        else:
            self.model_dir = dir

        self.sbml_level = slvl
        self.sbml_model = self.SBML.Model()
        self.sbml_model.setName(self.model_obj.ModelFile[:-4])
        self.sbml_document = self.SBML.SBMLDocument()
        # new stuff
        self.global_parameters = []
        self.model_compartment_name = None

        # create initdict
        self.__InitStrings__ = [
            s.replace('self.', '') for s in self.model_obj._PysMod__InitStrings
        ]
        self.__InitDict__ = {}
        for ii in self.__InitStrings__:
            l, r = ii.split('=')
            self.__InitDict__.setdefault(
                self.__cleanString__(l), float(self.__cleanString__(r))
            )

        # create forcing function dic
        try:
            self.parseForcingFunctions()
        except:
            print("No pre-defined forcing functions")

        if self.sbml_level == 1:
            if sys.platform == 'win32':
                print(
                    'Due to a bug in self.SBML for Windows writing a lvl 1 file will crash your session writing lvl 2 instead ... sorry'
                )
                self.sbml_document.setLevel(2)
            else:
                self.sbml_document.setLevel(self.sbml_level)
        else:
            self.sbml_document.setLevel(2)

    def SBML_setCompartment(self, name=None, vol=1):
        """
        SBML_setCompartment(name=None,vol=1)

        Initialise SBML compartments (note PySCeS currently utilises a single compartment)

        Arguments:
        =========
        name [default=None]: the compartment name, default is compartment1
        vol [default=1]: the compartment volume

        """
        if self._debug:
            print('SBML_setCompartment')
        comp_def = self.sbml_model.createCompartment()
        if not name:
            self.model_compartment_name = 'Cell'
        else:
            self.model_compartment_name = name
            for char in [' ', '.', '-', '*', '?', '!', '\t', '\n']:
                self.model_compartment_name = self.model_compartment_name.replace(
                    char, '_'
                )
            self.model_compartment_name = name
        comp_def.setId(self.model_compartment_name)
        comp_def.setVolume(vol)

    def SBML_setNotes(self, txt=None):
        notes = '<body xmlns="http://www.w3.org/1999/xhtml">'
        if txt != None:
            notes += '<span style="font-family: Courier New,Courier,monospace;">'
            notes += txt
            notes += '</span>'
        notes += '</body>'
        self.sbml_model.setNotes(notes)

    def SBML_setUnits(self, **kwargs):
        """
        SBML_setUnits(substance=(1,0), volume=(1,0), time=(1,0))

        Set the SBML default units note that the input here is the factor and index multiplying
        the SBML default, so for example the default substance (1,0) is (1*10**0)*mole So e.g.
        if you were specifing default units of millimoles and seconds you would
        set substance=(1,-3) and time=(60,0) i.e. (1*10**-3)*mole and (60*10**0)*seconds

        Arguments:
        =========
        substance [default=(1,0)]: the model substance unit - SBML default is "mole"
        volume [default=(1,0)]: the model volume unit - SBML default is "litre"
        time [default=(1,0)]: the model time unit - SBML default is "second"

        """
        for un in list(kwargs.keys()):
            vdef = self.sbml_model.createUnitDefinition()
            vdef.setId(un)
            vu = self.sbml_model.createUnit()
            if un == 'substance':
                vu.setKind(self.SBML.UnitKind_forName('mole'))
            elif un == 'volume':
                vu.setKind(self.SBML.UnitKind_forName('litre'))
            elif un == 'time':
                vu.setKind(self.SBML.UnitKind_forName('second'))
            vu.setMultiplier(kwargs[un][0])
            vu.setScale(kwargs[un][1])
            vu.setOffset(0)

    def SBML_setSpecies(self):
        """
        SBML_setSpecies()

        Initialise and add species information to the SBML model

        Arguments:
        None

        """
        if self._debug:
            print('SBML_setSpecies')
        reagList = self.model_obj.__species__ + self.model_obj.__fixed_species__
        for reagent in range(len(reagList)):
            s = self.sbml_model.createSpecies()
            s.setId(reagList[reagent])
            s.setName(reagList[reagent])
            s.setCompartment(self.model_compartment_name)
            if reagList[reagent] in self.model_obj.__fixed_species__:
                s.setBoundaryCondition(True)
                s.setConstant(True)
            else:
                s.setBoundaryCondition(False)

            if reagent < len(self.model_obj.__species__):
                reagName = reagList[reagent] + '_init'
            else:
                reagName = reagList[reagent]

            if self.sbml_level == 1:
                s.setInitialAmount(getattr(self.model_obj, reagName))
            else:
                s.setInitialConcentration(getattr(self.model_obj, reagName))

    def SBML_setAssignmentRules(self, rules=[]):

        for rule in rules:
            print(rule)
            self.global_parameters.append(rule[0])
            p = self.sbml_model.createParameter()
            p.setId(rule[0])
            p.setValue(getattr(self.model_obj, rule[0]))
            p.setConstant(False)
            r = self.sbml_model.createAssignmentRule()
            r.setVariable(rule[0])
            r.setFormula(rule[1])
            r.setMathFromFormula()

    def SBML_setReactions(self):
        """
        SBML_setReactions()

        Add kinetic rate laws to the SBMl model

        Arguments:
        None

        """
        if self._debug:
            print('SBML_setReactions')
        # TotSpecies = list(self.model_obj._PysMod__FixedReagents)+list(self.model_obj._PysMod__VarReagents)
        reaction_params = []
        for rxn in self.model_obj._PysMod__ReactionIDs:
            print('Adding reaction:', rxn)
            i = self.sbml_model.createReaction()
            i.setId(rxn)
            ndr = self.model_network_dict[rxn]
            for reagent in ndr['Reagents']:
                stoich = ndr['Reagents'][reagent]
                species = self.SBML.SpeciesReference(
                    reagent.replace('self.', ''), abs(stoich)
                )
                if stoich < 0:
                    i.addReactant(species)
                elif stoich > 0:
                    i.addProduct(species)
                elif stoich == 0:
                    i.addModifier(species)
            # add a volume to convert rate equation to kinetic law
            kineticLaw = ndr['RateEq'].replace('self.', '')
            kineticLaw = kineticLaw.replace('scipy.', '')
            if self.model_compartment_name not in self.model_obj.parameters:
                kineticLaw = self.model_compartment_name + ' * (' + kineticLaw + ')'
            else:
                kineticLaw = kineticLaw
            kineticLaw = self.SBML.KineticLaw(kineticLaw)

            # local parameters retired in favour of globals
            ##  for parameter in ndr['Params']:
            ##  p = parameter.replace('self.','')
            ##  if p not in self.model_obj.__fixed_species__ and p not in self.global_parameters:
            ##  try:
            ##  kineticLaw.addParameter(self.SBML.Parameter(p, getattr(self.model_obj,p)))
            ##  reaction_params.append(p)
            ##  except AttributeError,err :
            ##  print '\n', err
            ##  print "Parameter set error ... are there forcing functions??"
            ##  sleep(0.5)
            i.setKineticLaw(kineticLaw)
            if ndr['Type'] == 'Rever':
                rev = True
            else:
                rev = False
            i.setReversible(rev)

            # Add modifiers to reaction - brett 20050607
            for reac in self.model_obj.__modifiers__:
                if reac[0] == rxn:
                    for x in reac[1]:
                        print(' ' + reac[0] + ' has modifier: ' + x)
                        self.sbml_model.createModifier().setSpecies(x)

        # add extra parameter initialised but not in reactions
        # we have to do this in case the assignment rules are added after we build the model
        hack = list(self.__forcing_function_dic__.keys())

        not_xparams = (
            self.global_parameters
            + reaction_params
            + list(self.model_obj.species)
            + list(self.model_obj.fixed_species)
            + [self.model_compartment_name]
            + hack
        )

        for k in list(self.__InitDict__.keys()):
            if k not in not_xparams:
                print('Adding parameter:', k)
                self.global_parameters.append(k)
                p = self.sbml_model.createParameter()
                p.setId(k)
                p.setValue(getattr(self.model_obj, k))

    def SBML_setModel(self):
        """
        SBML_setModel()

        Add the SBML model to the predefined SBML document

        Arguments:
        None

        """
        if self._debug:
            print('SBML_setModel')
        self.sbml_document.setModel(self.sbml_model)

    def SBML_writeFile(self):
        """
        SBML_writeFile()

        Write the SBML document to predefined output file

        Arguments:
        None

        """

        self.SBML.writeSBML(self.sbml_document, 'pysces_sbml_tmp.xml')
        Fin = open('pysces_sbml_tmp.xml', 'r')
        Fout = open(os.path.join(self.model_dir, self.model_filename + '.xml'), 'w')
        cntr = 0
        try:
            UseR = getuser()
        except:
            UseR = ''
        for line in Fin:
            if cntr == 1:
                Fout.write(
                    '<!-- Created with PySCeS ('
                    + __version__
                    + ') on '
                    + strftime("%a, %d %b %Y %H:%M:%S")
                    + ' by '
                    + UseR
                    + ' -->\n'
                    + line
                )
            else:
                Fout.write(line)
            cntr += 1
        Fout.close()
        Fin.close()

        os.remove('pysces_sbml_tmp.xml')

    def convert2psc(self, filename, dir=None, dirOut=None):
        """
        convert2psc(filename,dir=None,dirOut=None)

        Convert an SBML file into a PySCeS input file

        Arguments:
        =========
        filename: the SBML source file
        dir [default=None]: specify the SBMl file directory
        dirOut [default=None]: the PSC file output directory

        """
        if dir == None:
            dir = os.getcwd()
        File = os.path.join(dir, filename)
        assert os.path.exists(File), "Invalid path"

        self.model_filename = filename

        r = self.SBML.SBMLReader()
        d = r.readSBML(File)
        m = d.getModel()

        def getName(i):
            if d.getLevel() == 1:
                return i.getName()
            else:
                return i.getId()

        reactions = m.getListOfReactions()
        ReactionIDs = []
        for i in reactions:
            ReactionIDs.append(getName(i))

        init_fixed = []
        init_var = []
        init_par = []
        parameters = []

        for i in m.getListOfSpecies():
            parName = getName(i)
            # if a species is a BoundaryCondition or constant it becomes fixed - brett 20050111
            if i.getBoundaryCondition() or i.getConstant():
                if i.getConstant() and not i.getBoundaryCondition():
                    print(
                        parName,
                        ' is set as constant, assuming: BoundaryCondition = True',
                    )
                init_fixed.append((parName, i.getInitialConcentration()))
            else:
                init_var.append((parName, i.getInitialConcentration()))

        NetworkDict = dict(
            [
                (i, dict.fromkeys(['Params', 'RateEq', 'Reagents', 'Type']))
                for i in ReactionIDs
            ]
        )

        for i in reactions:
            rDict = NetworkDict[getName(i)]
            j = i.getKineticLaw()
            par = []
            try:
                for k in j.getListOfParameters():
                    par.append(getName(k))
                    init_par.append((getName(k), k.getValue()))
                    parameters.append(getName(k))
                rDict['Params'] = par
                rDict['RateEq'] = j.getFormula()
                if d.getLevel() == 1:
                    rDict['RateEq'] = rDict['RateEq'].replace(' ', '')
                    rDict['RateEq'] = rDict['RateEq'].replace('^', '**')
            except Exception as err:
                rDict['Params'] = []
                rDict['RateEq'] = ''
                print(err)

            Substrates = []
            Products = []

            for k in i.getListOfReactants():
                species = k.getSpecies()
                stoich = -k.getStoichiometry()
                Substrates.append((species, stoich))

            for k in i.getListOfProducts():
                species = k.getSpecies()
                stoich = k.getStoichiometry()
                Products.append((species, stoich))

            # this is to eliminate zero stoichiometries {0}xyz
            badList = []
            for sub in Substrates:
                if sub[1] == 0:
                    badList.append(sub)
            for bad in badList:
                Substrates.pop(Substrates.index(bad))

            badList = []
            for prod in Products:
                if prod[1] == 0:
                    badList.append(prod)
            for bad in badList:
                Products.pop(Products.index(bad))

            # add source/sink pools to nasty substrate/productless reactions - brett 20050908
            if len(Substrates) == 0:
                Substrates.append(('$pool', -1.0))
            if len(Products) == 0:
                Products.append(('$pool', 1.0))

            # print Substrates
            # print Products

            rDict['Reagents'] = dict(Substrates + Products)
            if i.getReversible() == True:
                t = 'Rever'
            else:
                t = 'Irrev'
            rDict['Type'] = t
            NetworkDict[getName(i)].update(rDict)

        # Add extra model parameters not defined in reactions (apparently)
        if len(m.getListOfParameters()) > 0:
            for x in m.getListOfParameters():
                if getName(x) not in parameters:
                    # print getName(x)
                    init_par.append((getName(x), x.getValue()))

        if dirOut == None:
            self.model_filename = os.path.join(os.getcwd(), self.model_filename)
        else:
            self.model_filename = os.path.join(dirOut, self.model_filename)

        # print 'init_par'
        # print init_par
        # print 'init_var'
        # print init_var
        # print 'init_fixed'
        # print init_fixed

        # sometimes things just work lekker (replaced all the old showS^&t) - brett 20050913
        outFile = open(self.model_filename + '.psc', 'w')
        self.PSC_writeHeader(outFile)
        self.PSC_writeFixedSpeciesList(outFile, init_fixed)
        self.PSC_writeRateEquations(outFile, NetworkDict, number_format='%2.3f')
        self.PSC_writeSpecies(outFile, init_var)
        self.PSC_writeFixedSpecies(outFile, init_fixed)
        self.PSC_writeParameters(outFile, init_par)
        outFile.close()

        # Initialise compartment volumes as a parameter - brett 20050908
        compartmentList = []
        for comp in m.getListOfCompartments():
            # print comp
            compartmentList.append((getName(comp), comp.getVolume()))

        if len(compartmentList) > 1:
            print('\nINFO: PySCeS models are assumed to have a single compartment')

        if len(compartmentList) > 0:
            F = open(self.model_filename + '.psc', 'a')
            F.write('\n## Initialise compartment volumes')
            for comp in compartmentList:
                F.write('\n' + comp[0] + ' = ' + str(comp[1]))
                ##  parameters.append(x[0])
                ##  init_par.append(x)
            F.write('\n')
            F.close()

        # Add assignment rules as forcing functions - brett 20050908
        pscRules = []
        for rule in m.getListOfRules():
            pscRules.append((rule.getVariable(), rule.getFormula()))

        if len(pscRules) > 0:
            F = open(self.model_filename + '.psc', 'a')
            F.write('\n## Assignment rules translated to forcing functions\n')
            for rule in pscRules:
                rule0 = 'self.' + rule[0]
                rule1l = rule[1].split()
                for word in range(len(rule1l)):
                    if rule1l[word].isalnum():
                        if rule1l[word] not in [
                            '1',
                            '2',
                            '3',
                            '4',
                            '5',
                            '6',
                            '7',
                            '8',
                            '9',
                        ]:
                            rule1l[word] = 'self.' + rule1l[word]
                F.write('!F ' + rule0 + ' = ')
                for word in rule1l:
                    F.write(word + ' ')
                F.write('\n')
            F.write('\n')
            F.close()

        if len(m.getNotes()) > 0:
            F = open(self.model_filename + '.psc', 'a')
            F.write('\n## Model notes' + m.getNotes().replace('\n', '\n# ') + '\n\n')
            F.close()

    def PSC_writeHeader(self, File):
        """
        PSC_writeHeader(File)

        Write a PSC file header to an open file object

        Arguments:
        =========
        File: a writable open text file object

        """
        try:
            UseR = getuser()
        except:
            UseR = ''

        header = ''
        # header +=  '############################################################\n'
        header += '# PySCeS (' + __version__ + ') model input file\n'
        header += '# PySCeS can be found at http://pysces.sourceforge.net/\n'
        # header += '###########################################################\n\n'
        header += '# Original input file: ' + File.name.split('\\')[-1][:-4] + '\n'
        header += (
            '# This file generated: '
            + strftime("%a, %d %b %Y %H:%M:%S")
            + ' by '
            + UseR
            + '\n\n'
        )
        File.write(header)
        File.write('\n')

    def PSC_writeSpecies(self, File, species):
        """
        PSC_writeSpecies(File,species)

        Write out model species initiaiisations to file

        Arguments:
        =========
        File: a writable open file object
        species: a list of (species.value) pairs

        """
        out_list = []

        out_list.append('\n## Variable species initial values\n')
        for x in range(len(species)):
            out_list.append(
                species[x][0] + ' = ' + self.mode_number_format % species[x][1] + '\n'
            )
        for x in out_list:
            File.write(x)
        File.write('\n')

    def PSC_writeFixedSpeciesList(self, File, fixed_species):
        """
        PSC_writeFixedSpeciesList(File,fixed_species)

        Write fixed species declaration to a PSC file

        Arguments:
        =========
        File: open, writable file object
        fixed_species: a list of (species,value) pairs

        """
        File.write('## Fixed species\n')
        if len(fixed_species) == 0:
            File.write('#  <none>')
        else:
            File.write('FIX: ')
            for x in fixed_species:
                File.write(x[0] + ' ')
        File.write('\n\n')

    def PSC_writeFixedSpecies(self, File, fixed_species):
        """
        PSC_writeFixedSpecies(File,fixed_species)

        Write fixed species initialisations to a PSC file

        Arguments:
        =========
        File: open, writable file object
        fixed_species: a list of (species,value) pairs

        """
        out_list = []

        out_list.append('\n## Fixed species\n')
        for x in range(len(fixed_species)):
            out_list.append(
                fixed_species[x][0]
                + ' = '
                + self.mode_number_format % fixed_species[x][1]
                + '\n'
            )
        for x in out_list:
            File.write(x)
        File.write('\n')

    def PSC_writeParameters(self, File, parameters):
        """
        PSC_writeParameters(File,parameters)

        Write mode parameter initialisations to a PSC file

        Arguments:
        =========
        File: open, writable file object
        parameters: a list of (parameter,value) pairs

        """
        out_list = []

        out_list.append('\n## Parameters\n')
        for x in range(len(parameters)):
            out_list.append(
                parameters[x][0]
                + ' = '
                + self.mode_number_format % parameters[x][1]
                + '\n'
            )
        for x in out_list:
            File.write(x)
        File.write('\n')

    def PSC_writeRateEquations(self, File, NetworkDict, number_format='%2.3f'):
        """
        PSC_writeRateEquations(File,NetworkDict,number_format='%2.3f')

        Write model rate equations to a PSC file

        Arguments:
        =========
        File: open, writable file object
        NetworkDict: a PySCeS network dictionary
        number_format [default='%2.3f']: number formatting to use in rate laws

        """
        out_list = []

        out_list.append('\n## Reaction stoichiometry and rate equations\n')
        for key in NetworkDict:
            out_list.append(key + ':\n')
            reagL = []
            reagR = []
            for reagent in NetworkDict[key]['Reagents']:
                if NetworkDict[key]['Reagents'][reagent] > 0:
                    if NetworkDict[key]['Reagents'][reagent] == 1.0:
                        reagR.append(reagent.replace('self.', ''))
                    else:
                        reagR.append(
                            '{'
                            + number_format % abs(NetworkDict[key]['Reagents'][reagent])
                            + '}'
                            + reagent.replace('self.', '')
                        )
                elif NetworkDict[key]['Reagents'][reagent] < 0:
                    if NetworkDict[key]['Reagents'][reagent] == -1.0:
                        reagL.append(reagent.replace('self.', ''))
                    else:
                        reagL.append(
                            '{'
                            + number_format % abs(NetworkDict[key]['Reagents'][reagent])
                            + '}'
                            + reagent.replace('self.', '')
                        )
                elif NetworkDict[key]['Reagents'][reagent] == 0:
                    # reagL.append(reagent.replace('self.',''))
                    print(NetworkDict[key]['Reagents'])
                    input('WTF: please contact developers')

            if len(reagL) == 0:
                print('Zero pool substrate', File.name)
                reagL.append('$pool')
            if len(reagR) == 0:
                print('Zero pool product', File.name)
                reagR.append('$pool')

            substring = ''
            count = 0
            for x in reagL:
                if count != 0:
                    substring += ' + '
                substring += x.replace(' ', '')
                count += 1
            prodstring = ''
            count = 0
            for x in reagR:
                if count != 0:
                    prodstring += ' + '
                prodstring += x.replace(' ', '')
                count += 1

            if NetworkDict[key]['Type'] == 'Rever':
                symbol = ' = '
            else:
                symbol = ' > '
            out_list.append('\t' + substring + symbol + prodstring + '\n')
            out_list.append(
                '\t' + NetworkDict[key]['RateEq'].replace('self.', '') + '\n\n'
            )
        for x in out_list:
            File.write(x)
