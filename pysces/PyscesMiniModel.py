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

"""
This module contains simplified methods derived from the Pysces model class
Brett G. Olivier June 2010
"""

try:
    input = raw_input  # Py2 compatibility
except NameError:
    pass
import os, copy, time
import numpy
from pysces import PyscesStoich
from pysces import PyscesParse

mach_spec = numpy.MachAr()
pscParser = PyscesParse.PySCeSParser(debug=0)


class PyscesInputFileParser(object):
    """
    This class contains the PySCeS model loading and Stoichiometric Analysis methods
    """

    ModelDir = None
    ModelFile = None
    ModelOutput = None
    __settings__ = None
    N = None

    def __init__(self, model_file, directory, output_dir=None):
        self.ModelDir = directory
        self.ModelFile = model_file
        if output_dir == None:
            self.ModelOutput = os.getcwd()
        else:
            assert os.path.exists(output_dir), "\n%s is not a valid path" % output_dir
        self.__settings__ = {}
        # Initialize stoichiometric precision
        self.__settings__['stoichiometric_analysis_fp_zero'] = mach_spec.eps * 2.0e4
        self.__settings__['stoichiometric_analysis_lu_precision'] = self.__settings__[
            'stoichiometric_analysis_fp_zero'
        ]
        self.__settings__['stoichiometric_analysis_gj_precision'] = (
            self.__settings__['stoichiometric_analysis_lu_precision'] * 10.0
        )
        self.__settings__['enable_deprecated_attr'] = False
        self.InitialiseInputFile()
        self.N = self.buildN()

    def InitialiseInputFile(self):
        """
        InitialiseInputFile()

        Parse the input file associated with the PySCeS model instance and assign the basic model attributes

        Arguments:
        None

        """
        self.__parseOK = 1  # check that model has parsed ok?
        try:
            if os.path.exists(os.path.join(self.ModelDir, self.ModelFile)):
                pass
            else:
                print(
                    '\nInvalid self.ModelFile: '
                    + os.path.join(self.ModelDir, self.ModelFile)
                )
        except:
            print(
                'WARNING: Problem verifying: '
                + os.path.join(self.ModelDir, self.ModelFile)
            )

        if self.ModelFile[-4:] == '.psc':
            pass
        else:
            print('Assuming extension is .psc')
            self.ModelFile += '.psc'

        print('\nParsing file: %s' % os.path.join(self.ModelDir, self.ModelFile))

        pscParser.ParsePSC(self.ModelFile, self.ModelDir, self.ModelOutput)
        print(' ')

        badlist = pscParser.KeywordCheck(pscParser.ReactionIDs)
        badlist = pscParser.KeywordCheck(pscParser.Inits, badlist)

        if len(badlist) != 0:
            print(
                '\n******************************\nPSC input file contains PySCeS keywords please rename them and reload:'
            )
            for item in badlist:
                print('   --> ' + item)
            print('******************************\n')
            self.__parseOK = 0
            # assert len(badlist) != 0, 'Keyword error, please check input file'

        if self.__parseOK:
            # brett 2008
            self.__nDict__ = pscParser.nDict.copy()
            self.__sDict__ = pscParser.sDict.copy()
            self.__pDict__ = pscParser.pDict.copy()
            self.__uDict__ = pscParser.uDict.copy()

            # model attributes are now initialised here brett2008
            self.__InitDict__ = {}
            # set parameters and add to __InitDict__
            for p in list(self.__pDict__.keys()):
                setattr(self, self.__pDict__[p]['name'], self.__pDict__[p]['initial'])
                self.__InitDict__.update(
                    {self.__pDict__[p]['name']: self.__pDict__[p]['initial']}
                )
            # set species and add to __InitDict__ and set mod.Xi_init
            for s in list(self.__sDict__.keys()):
                setattr(self, self.__sDict__[s]['name'], self.__sDict__[s]['initial'])
                if not self.__sDict__[s]['fixed']:
                    setattr(
                        self,
                        self.__sDict__[s]['name'] + '_init',
                        self.__sDict__[s]['initial'],
                    )
                self.__InitDict__.update(
                    {self.__sDict__[s]['name']: self.__sDict__[s]['initial']}
                )

            # setup keywords
            self.__KeyWords__ = pscParser.KeyWords.copy()
            if self.__KeyWords__['Modelname'] == None:
                self.__KeyWords__['Modelname'] = self.ModelFile.replace('.psc', '')
            if self.__KeyWords__['Description'] == None:
                self.__KeyWords__['Description'] = self.ModelFile.replace('.psc', '')
            # if SpeciesTypes undefined assume []
            if self.__KeyWords__['Species_In_Conc'] == None:
                self.__KeyWords__['Species_In_Conc'] = True
            # if OutputType is undefined assume it is the same as SpeciesType
            if self.__KeyWords__['Output_In_Conc'] == None:
                if self.__KeyWords__['Species_In_Conc']:
                    self.__KeyWords__['Output_In_Conc'] = True
                else:
                    self.__KeyWords__['Output_In_Conc'] = False
            # set the species type in sDict according to 'Species_In_Conc'
            for s in list(self.__sDict__.keys()):
                if not self.__KeyWords__['Species_In_Conc']:
                    self.__sDict__[s]['isamount'] = True
                else:
                    self.__sDict__[s]['isamount'] = False

            # setup compartments
            self.__compartments__ = pscParser.compartments.copy()
            if len(list(self.__compartments__.keys())) > 0:
                self.__HAS_COMPARTMENTS__ = True
            else:
                self.__HAS_COMPARTMENTS__ = False

            # no (self.)
            self.__fixed_species__ = copy.copy(pscParser.fixed_species)
            self.__species__ = copy.copy(pscParser.species)
            self.__parameters__ = copy.copy(pscParser.parameters)
            self.__reactions__ = copy.copy(pscParser.reactions)
            self.__modifiers__ = copy.copy(pscParser.modifiers)
            # Initialize exposed stuff
            self.fixed_species = tuple(pscParser.fixed_species)
            self.species = tuple(pscParser.species)
            self.parameters = tuple(pscParser.parameters)
            self.reactions = tuple(pscParser.reactions)
            self.modifiers = tuple(pscParser.modifiers)

            # Add input file defined fuctions - brett 200500621
            # TODO deprecated
            # self._Function_time = copy.copy(pscParser.TimeFunc)
            self._Function_user = copy.copy(pscParser.UserFunc)
            self._Function_init = pscParser.InitFunc

            self.__functions__ = pscParser.Functions.copy()
            self.__rules__ = pscParser.AssignmentRules.copy()
            self.__InitFuncs__ = pscParser.ModelInit.copy()
            self.__userfuncs__ = pscParser.UserFuncs.copy()
            self.__eDict__ = pscParser.Events.copy()

            ##  if pscParser.ModelUsesNumpyFuncs:
            ##  print 'Numpy functions detected in kinetic laws.\n'
        else:
            print('\nERROR: model parsing error, please check input file.\n')
        # added in a check for model correctness and human error reporting (1=ok, 0=error)
        if len(pscParser.SymbolErrors) != 0:
            print('\nUndefined symbols:\n%s' % self.SymbolErrors)
        if not pscParser.ParseOK:
            print(
                '\n\n*****\nModel parsing errors detected in input file '
                + self.ModelFile
                + '\n*****'
            )
            print('\nInput file errors')
            for error in pscParser.LexErrors:
                print(
                    error[0]
                    + 'in line:\t'
                    + str(error[1])
                    + ' ('
                    + error[2][:20]
                    + ' ...)'
                )
            print('\nParser errors')
            for error in pscParser.ParseErrors:
                try:
                    print(error[0] + '- ' + error[2][:20])
                except:
                    print(error)
            assert pscParser.ParseOK == 1, 'Input File Error'

    def buildN(self):
        """
        buildN()

        Generate the stoichiometric matrix N from the parsed model description.
        Returns a stoichiometric matrix (N)

        Arguments:
        None

        """
        VarReagents = ['self.' + s for s in self.__species__]
        StoicMatrix = numpy.zeros((len(VarReagents), len(self.__reactions__)), 'd')
        for reag in VarReagents:
            for id in self.__reactions__:
                if reag in list(self.__nDict__[id]['Reagents'].keys()):
                    StoicMatrix[VarReagents.index(reag)][
                        self.__reactions__.index(id)
                    ] = self.__nDict__[id]['Reagents'][reag]
        return StoicMatrix

    def Stoichiometry_Init(self, nmatrix):
        """
        Stoichiometry_Init(nmatrix,load=0)

        Initialize the model stoichiometry. Given a stoichiometric matrix N, this method will return an instantiated PyscesStoich instance and status flag.
        and test it's correctness. The status flag indicates 0 = reanalyse stoichiometry or
        1 = complete structural analysis preloaded.

        Arguments:

        nmatrix: The input stoichiometric matrix, N
        load [default=0]: try to load a saved stoichiometry (1)

        """
        # print 'Instantiating new stoichiometry ...'
        stc = PyscesStoich.Stoich(nmatrix)
        status = 0
        return stc, status

    def Stoichiometry_Analyse(self):
        override = 0
        load = 0
        """
        Stoichiometry_Analyse(override=0,load=0)

        Perform a structural analyses. The default behaviour is to construct and analyse the model
        from the parsed model information. Overriding this behaviour analyses the stoichiometry
        based on the current stoichiometric matrix. If load is specified PySCeS tries to load a
        saved stoichiometry, otherwise the stoichiometric analysis is run. The results of
        the analysis are checked for floating point error and nullspace rank consistancy.

        Arguments:

        override [default=0]: override stoichiometric analysis intialisation from parsed data
        load [default=0]: load a presaved stoichiometry

        """
        if not override:
            self.nmatrix = self.buildN()  # Creates the model N
            # print '\nintializing N\n'
        else:
            print('\nStoichiometric override active\n')

        assert len(self.nmatrix) > 0, (
            '\nUnable to generate Stoichiometric Matrix! model has:\n%s reactions\n%s species\nwhat did you have in mind?\n'
            % (len(self.__reactions__), len(self.__species__))
        )

        ##  self.__nmatrix__ = copy.copy(self.nmatrix)
        self.__nmatrix__ = self.nmatrix  # done with caution brett2008
        self.__Nshape__ = self.nmatrix.shape  # Get the shape of N
        ##  self.__Vtemp__ = numpy.zeros((self.__Nshape__[1])) # going going ....

        # get stoich instance and whether it was analysed or loaded - brett 20050830
        self.__structural__, stc_load = self.Stoichiometry_Init(self.nmatrix)

        # if not loaded analyze - brett 20050830
        if not stc_load:
            # technically this means we can define this on the fly - brett #20051013
            self.__structural__.stoichiometric_analysis_fp_zero = self.__settings__[
                'stoichiometric_analysis_fp_zero'
            ]
            self.__structural__.stoichiometric_analysis_lu_precision = self.__settings__[
                'stoichiometric_analysis_lu_precision'
            ]
            self.__structural__.stoichiometric_analysis_gj_precision = self.__settings__[
                'stoichiometric_analysis_gj_precision'
            ]
            self.__structural__.AnalyseL()  # Get all L related stuff
            self.__structural__.AnalyseK()  # Get all K related stuff

        # test matrix values against __settings__['stoichiometric_analysis_lu_precision']
        lsmall, lbig = self.__structural__.MatrixValueCompare(
            self.__structural__.lzeromatrix
        )
        ksmall, kbig = self.__structural__.MatrixValueCompare(
            self.__structural__.kzeromatrix
        )
        SmallValueError = 0
        if (
            abs(lsmall)
            < self.__structural__.stoichiometric_analysis_lu_precision * 10.0
        ):
            print(
                '\nWARNING: values in L0matrix are close to stoichiometric precision!'
            )
            print(
                'Stoichiometric LU precision:',
                self.__structural__.stoichiometric_analysis_lu_precision,
            )
            print('L0 smallest abs(value)', abs(lsmall))
            print('Machine precision:', mach_spec.eps)
            SmallValueError = 1
        if (
            abs(ksmall)
            < self.__structural__.stoichiometric_analysis_lu_precision * 10.0
        ):
            print(
                '\nWARNING: values in K0matrix are close to stoichiometric precision!'
            )
            print(
                'Stoichiometric precision:',
                self.__structural__.stoichiometric_analysis_lu_precision,
            )
            print('K0 smallest abs(value)', abs(ksmall))
            print('Machine precision:', mach_spec.eps)
            SmallValueError = 1
        if SmallValueError:
            input(
                '\nStructural Analysis results may not be reliable!!!.\n\nTry change <mod>.__settings__["stoichiometric_analysis_lu_precision"] (see reference manual for details)\n\n\t press any key to continue: '
            )

        # cross check that rank is consistant between K0 and L0
        if (
            self.__structural__.kzeromatrix.shape[0]
            != self.__structural__.lzeromatrix.shape[1]
        ):
            print(
                '\nWARNING: the rank calculated by the Kand L analysis methods are not the same!'
            )
            print(
                '\tK analysis calculates the rank as: '
                + repr(self.__structural__.kzeromatrix.shape[0])
            )
            print(
                '\tL analysis calculates the rank as: '
                + repr(self.__structural__.lzeromatrix.shape[1])
            )
            print('This is not good! Structural Analysis results are not reliable!!!\n')
            assert (
                self.__structural__.kzeromatrix.shape[0]
                == self.__structural__.lzeromatrix.shape[1]
            ), '\nStructuralAnalysis Error: rank mismatch'

        self.__HAS_FLUX_CONSERVATION__ = self.__structural__.info_flux_conserve
        self.__HAS_MOIETY_CONSERVATION__ = self.__structural__.info_moiety_conserve

        if self.__settings__['enable_deprecated_attr']:
            self.nmatrix_row = self.__structural__.nmatrix_row
            self.nmatrix_col = self.__structural__.nmatrix_col

            self.kmatrix = self.__structural__.kmatrix
            self.kmatrix_row = self.__structural__.kmatrix_row
            self.kmatrix_col = self.__structural__.kmatrix_col
            self.kzeromatrix = self.__structural__.kzeromatrix
            self.kzeromatrix_row = self.__structural__.kzeromatrix_row
            self.kzeromatrix_col = self.__structural__.kzeromatrix_col

            self.lmatrix = self.__structural__.lmatrix
            self.lmatrix_row = self.__structural__.lmatrix_row
            self.lmatrix_col = self.__structural__.lmatrix_col
            self.lzeromatrix = self.__structural__.lzeromatrix
            self.lzeromatrix_row = self.__structural__.lzeromatrix_row
            self.lzeromatrix_col = self.__structural__.lzeromatrix_col
            self.conservation_matrix = self.__structural__.conservation_matrix
            self.conservation_matrix_row = self.__structural__.conservation_matrix_row
            self.conservation_matrix_col = self.__structural__.conservation_matrix_col
            self.nrmatrix = self.__structural__.nrmatrix
            self.nrmatrix_row = self.__structural__.nrmatrix_row
            self.nrmatrix_col = self.__structural__.nrmatrix_col

            self.__kmatrix__ = copy.copy(self.kmatrix)
            self.__kzeromatrix__ = copy.copy(self.kzeromatrix)
            self.__lmatrix__ = copy.copy(self.lmatrix)
            self.__lzeromatrix__ = copy.copy(self.lzeromatrix)
            self.__nrmatrix__ = copy.copy(self.nrmatrix)
        # switch that is set if the stoichiometric analysis is up to date

        self.__structural__.species = self.species
        self.__structural__.reactions = self.reactions
        self.Nmatrix = PyscesStoich.StructMatrix(
            self.__structural__.nmatrix,
            self.__structural__.nmatrix_row,
            self.__structural__.nmatrix_col,
        )
        self.Nmatrix.setRow(self.species)
        self.Nmatrix.setCol(self.reactions)

        self.Nrmatrix = PyscesStoich.StructMatrix(
            self.__structural__.nrmatrix,
            self.__structural__.nrmatrix_row,
            self.__structural__.nrmatrix_col,
        )
        self.Nrmatrix.setRow(self.species)
        self.Nrmatrix.setCol(self.reactions)

        self.Kmatrix = PyscesStoich.StructMatrix(
            self.__structural__.kmatrix,
            self.__structural__.kmatrix_row,
            self.__structural__.kmatrix_col,
        )
        self.Kmatrix.setRow(self.reactions)
        self.Kmatrix.setCol(self.reactions)

        self.K0matrix = PyscesStoich.StructMatrix(
            self.__structural__.kzeromatrix,
            self.__structural__.kzeromatrix_row,
            self.__structural__.kzeromatrix_col,
        )
        self.K0matrix.setRow(self.reactions)
        self.K0matrix.setCol(self.reactions)

        self.Lmatrix = PyscesStoich.StructMatrix(
            self.__structural__.lmatrix,
            self.__structural__.lmatrix_row,
            self.__structural__.lmatrix_col,
        )
        self.Lmatrix.setRow(self.species)
        self.Lmatrix.setCol(self.species)

        self.L0matrix = PyscesStoich.StructMatrix(
            self.__structural__.lzeromatrix,
            self.__structural__.lzeromatrix_row,
            self.__structural__.lzeromatrix_col,
        )
        self.L0matrix.setRow(self.species)
        self.L0matrix.setCol(self.species)

        if self.__structural__.info_moiety_conserve:
            self.Consmatrix = PyscesStoich.StructMatrix(
                self.__structural__.conservation_matrix,
                self.__structural__.conservation_matrix_row,
                self.__structural__.conservation_matrix_col,
            )
            self.Consmatrix.setRow(self.species)
            self.Consmatrix.setCol(self.species)
        else:
            self.Consmatrix = None
        self.__StoichOK = 1
        print(' ')


if __name__ == '__main__':
    ModelFile = 'pysces_test_linear1.psc'
    ModelDir = '/home/bgoli/Pysces/psc'
    mod = PyscesInputFileParser(ModelFile, ModelDir)
    # ~ mod.Stoichiometry_Analyse()
