"""
PySCeS - Python Simulator for Cellular Systems (http://pysces.sourceforge.net)

Copyright (C) 2004-2024 B.G. Olivier, J.M. Rohwer, J.-H.S Hofmeyr all rights reserved,

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
    input = raw_input
except NameError:
    pass
import os, copy, time
import numpy
from pysces import PyscesStoich
from pysces import PyscesParse
mach_eps = numpy.finfo(float).eps
pscParser = PyscesParse.PySCeSParser(debug=0)


class PyscesInputFileParser(object):
    """
    This class contains the PySCeS model loading and Stoichiometric Analysis methods.

    Attributes
    ----------
    ModelDir : str
        The directory where the model file is located.
    ModelFile : str
        The name of the model file to be parsed.
    ModelOutput : str
        The directory where the output will be saved. Defaults to the current working directory if not set.
    __settings__ : dict
        Internal settings used for stoichiometric analysis, including floating point zero and precision settings.
    N : numpy.ndarray
        The stoichiometric matrix generated from the parsed model.
    
    Methods
    -------
    InitialiseInputFile()
        Parses the input file associated with the PySCeS model instance and assigns basic model attributes.
    buildN()
        Generates the stoichiometric matrix N from the parsed model description.
    Stoichiometry_Init(nmatrix, load=0)
        Initializes the model stoichiometry with a given stoichiometric matrix and tests its correctness.
    Stoichiometry_Analyse(override=0, load=0)
        Performs structural analyses based on parsed model information or current stoichiometric matrix, and checks results for precision issues.
    """
    ModelDir = None
    ModelFile = None
    ModelOutput = None
    __settings__ = None
    N = None

    def __init__(self, model_file, directory, output_dir=None):
        """
    Initializes a PyscesInputFileParser instance.

    This constructor sets up the model file and directory paths, initializes 
    the output directory, prepares internal settings for stoichiometric analysis,
    and initializes the input file and stoichiometric matrix. If the output 
    directory is not provided, it defaults to the current working directory.

    Parameters
    ----------
    model_file : str
        The name of the model file to be processed.
    directory : str
        The directory where the model file is located.
    output_dir : str, optional
        The directory where output files will be stored. Defaults to the current 
        working directory if not specified.

    Raises
    ------
    AssertionError
        If the specified output directory does not exist.

    Attributes
    ----------
    ModelDir : str
        Directory where the model file is located.
    ModelFile : str
        Name of the model file.
    ModelOutput : str
        Directory where output files are stored.
    __settings__ : dict
        A dictionary containing settings for stoichiometric analysis, including:
        - 'stoichiometric_analysis_fp_zero' : float
        - 'stoichiometric_analysis_lu_precision' : float
        - 'stoichiometric_analysis_gj_precision' : float
        - 'enable_deprecated_attr' : bool
    N : matrix
        The stoichiometric matrix associated with the model file.

    Examples
    --------
    >>> parser = PyscesInputFileParser("model.psc", "/path/to/models")
    >>> print(parser.ModelDir)
    /path/to/models
    >>> print(parser.ModelOutput)
    /current/working/directory
    """
        self.ModelDir = directory
        self.ModelFile = model_file
        if output_dir is None:
            self.ModelOutput = os.getcwd()
        else:
            assert os.path.exists(output_dir
                ), '\n%s is not a valid path' % output_dir
        self.__settings__ = {}
        self.__settings__['stoichiometric_analysis_fp_zero'
            ] = mach_eps * 20000.0
        self.__settings__['stoichiometric_analysis_lu_precision'
            ] = self.__settings__['stoichiometric_analysis_fp_zero']
        self.__settings__['stoichiometric_analysis_gj_precision'
            ] = self.__settings__['stoichiometric_analysis_lu_precision'
            ] * 10.0
        self.__settings__['enable_deprecated_attr'] = False
        self.InitialiseInputFile()
        self.N = self.buildN()

    def InitialiseInputFile(self):
        """
    Parse the input file associated with the PySCeS model instance.

    This method reads the specified model input file and assigns the necessary
    basic model attributes to the instance. It validates the file's existence and
    extension, and performs various parsing checks to ensure the input is
    defined correctly. The parsed data is then used to initialize internal
    dictionaries and attributes representing model entities like species, parameters,
    reactions, and more.

    Notes
    -----
    - If the input file does not have the '.psc' extension, it will be automatically appended.
    - The existence of the `ModelFile` in the specified `ModelDir` is verified.
    - If parsing finds PySCeS keywords as variable names, you will be prompted
      to rename them.
    - Ensures the file can be parsed correctly and initializes all necessary model 
      attributes.

    Examples
    --------
    Initialize a `PyscesInputFileParser` instance and parse a model file:

    >>> parser = PyscesInputFileParser(model_file='example_model', directory='/path/to/models')
    >>> parser.InitialiseInputFile()

    This would parse the 'example_model.psc' file located in '/path/to/models'
    and initialize the model attributes within the parser instance.

    Raises
    ------
    AssertionError
        If the input file cannot be parsed correctly, an assertion error is raised
        indicating an input file error.

    """
        self.__parseOK = 1
        try:
            if os.path.exists(os.path.join(self.ModelDir, self.ModelFile)):
                pass
            else:
                print('\nInvalid self.ModelFile: ' + os.path.join(self.
                    ModelDir, self.ModelFile))
        except:
            print('WARNING: Problem verifying: ' + os.path.join(self.
                ModelDir, self.ModelFile))
        if self.ModelFile[-4:] == '.psc':
            pass
        else:
            print('Assuming extension is .psc')
            self.ModelFile += '.psc'
        print('\nParsing file: %s' % os.path.join(self.ModelDir, self.
            ModelFile))
        pscParser.ParsePSC(self.ModelFile, self.ModelDir, self.ModelOutput)
        print(' ')
        badlist = pscParser.KeywordCheck(pscParser.ReactionIDs)
        badlist = pscParser.KeywordCheck(pscParser.Inits, badlist)
        if len(badlist) != 0:
            print(
                """
******************************
PSC input file contains PySCeS keywords please rename them and reload:"""
                )
            for item in badlist:
                print('   --> ' + item)
            print('******************************\n')
            self.__parseOK = 0
        if self.__parseOK:
            self.__nDict__ = pscParser.nDict.copy()
            self.__sDict__ = pscParser.sDict.copy()
            self.__pDict__ = pscParser.pDict.copy()
            self.__uDict__ = pscParser.uDict.copy()
            self.__InitDict__ = {}
            for p in list(self.__pDict__.keys()):
                setattr(self, self.__pDict__[p]['name'], self.__pDict__[p][
                    'initial'])
                self.__InitDict__.update({self.__pDict__[p]['name']: self.
                    __pDict__[p]['initial']})
            for s in list(self.__sDict__.keys()):
                setattr(self, self.__sDict__[s]['name'], self.__sDict__[s][
                    'initial'])
                if not self.__sDict__[s]['fixed']:
                    setattr(self, self.__sDict__[s]['name'] + '_init', self
                        .__sDict__[s]['initial'])
                self.__InitDict__.update({self.__sDict__[s]['name']: self.
                    __sDict__[s]['initial']})
            self.__KeyWords__ = pscParser.KeyWords.copy()
            if self.__KeyWords__['Modelname'] is None:
                self.__KeyWords__['Modelname'] = self.ModelFile.replace('.psc',
                    '')
            if self.__KeyWords__['Description'] is None:
                self.__KeyWords__['Description'] = self.ModelFile.replace(
                    '.psc', '')
            if self.__KeyWords__['Species_In_Conc'] is None:
                self.__KeyWords__['Species_In_Conc'] = True
            if self.__KeyWords__['Output_In_Conc'] is None:
                if self.__KeyWords__['Species_In_Conc']:
                    self.__KeyWords__['Output_In_Conc'] = True
                else:
                    self.__KeyWords__['Output_In_Conc'] = False
            for s in list(self.__sDict__.keys()):
                if not self.__KeyWords__['Species_In_Conc']:
                    self.__sDict__[s]['isamount'] = True
                else:
                    self.__sDict__[s]['isamount'] = False
            self.__compartments__ = pscParser.compartments.copy()
            if len(list(self.__compartments__.keys())) > 0:
                self.__HAS_COMPARTMENTS__ = True
            else:
                self.__HAS_COMPARTMENTS__ = False
            self.__fixed_species__ = copy.copy(pscParser.fixed_species)
            self.__species__ = copy.copy(pscParser.species)
            self.__parameters__ = copy.copy(pscParser.parameters)
            self.__reactions__ = copy.copy(pscParser.reactions)
            self.__modifiers__ = copy.copy(pscParser.modifiers)
            self.fixed_species = tuple(pscParser.fixed_species)
            self.species = tuple(pscParser.species)
            self.parameters = tuple(pscParser.parameters)
            self.reactions = tuple(pscParser.reactions)
            self.modifiers = tuple(pscParser.modifiers)
            self._Function_user = copy.copy(pscParser.UserFunc)
            self._Function_init = pscParser.InitFunc
            self.__functions__ = pscParser.Functions.copy()
            self.__rules__ = pscParser.AssignmentRules.copy()
            self.__InitFuncs__ = pscParser.ModelInit.copy()
            self.__userfuncs__ = pscParser.UserFuncs.copy()
            self.__eDict__ = pscParser.Events.copy()
        else:
            print('\nERROR: model parsing error, please check input file.\n')
        if len(pscParser.SymbolErrors) != 0:
            print('\nUndefined symbols:\n%s' % self.SymbolErrors)
        if not pscParser.ParseOK:
            print('\n\n*****\nModel parsing errors detected in input file ' +
                self.ModelFile + '\n*****')
            print('\nInput file errors')
            for error in pscParser.LexErrors:
                print(error[0] + 'in line:\t' + str(error[1]) + ' (' +
                    error[2][:20] + ' ...)')
            print('\nParser errors')
            for error in pscParser.ParseErrors:
                try:
                    print(error[0] + '- ' + error[2][:20])
                except:
                    print(error)
            assert pscParser.ParseOK == 1, 'Input File Error'

    def buildN(self):
        """
    Generate the stoichiometric matrix N from the parsed model description.

    This method constructs the stoichiometric matrix `N`, representing the
    relationships between species and reactions in a biochemical model. Each
    element of the matrix reflects the stoichiometric coefficient of a species
    in a reaction.

    Returns
    -------
    numpy.ndarray
        A 2D numpy array representing the stoichiometric matrix `N`, where
        rows correspond to species and columns correspond to reactions.

    Examples
    --------
    >>> parser = PyscesInputFileParser('model.psc', '/path/to/models')
    >>> parser.InitialiseInputFile()
    >>> stoichiometric_matrix = parser.buildN()
    >>> print(stoichiometric_matrix)
    array([[ 1.0, -1.0, 0.0],
           [ 0.0,  1.0, -1.0],
           [ 1.0,  0.0,  1.0]])

    Here, the `stoichiometric_matrix` is generated based on the species and
    reactions defined in the 'model.psc' file loaded by the `PyscesInputFileParser`.
    """
        VarReagents = [('self.' + s) for s in self.__species__]
        StoicMatrix = numpy.zeros((len(VarReagents), len(self.__reactions__
            )), 'd')
        for reag in VarReagents:
            for id in self.__reactions__:
                if reag in list(self.__nDict__[id]['Reagents'].keys()):
                    StoicMatrix[VarReagents.index(reag)][self.__reactions__
                        .index(id)] = self.__nDict__[id]['Reagents'][reag]
        return StoicMatrix

    def Stoichiometry_Init(self, nmatrix):
        """
    Initialize the model stoichiometry.

    This method initializes the model stoichiometry by creating an instance of 
    the `PyscesStoich` class with the provided stoichiometric matrix `N`. It 
    also tests the correctness of the stoichiometric matrix. Based on the `load` 
    parameter, it attempts to load a pre-saved stoichiometry configuration. This 
    function returns a tuple containing the `PyscesStoich` instance and a status 
    flag that indicates the outcome of the initialization:
    - 0: reanalyse stoichiometry
    - 1: the complete structural analysis is preloaded

    Parameters
    ----------
    nmatrix : numpy.ndarray
        The input stoichiometric matrix, `N`.
    load : int, optional
        A flag to attempt loading a saved stoichiometry configuration. Default 
        is 0, meaning no attempt to load.

    Returns
    -------
    tuple
        A tuple containing the instantiated `PyscesStoich` object and the status 
        flag.

    Examples
    --------
    >>> parser = PyscesInputFileParser()
    >>> nmatrix = parser.buildN()
    >>> stoich_instance, status = parser.Stoichiometry_Init(nmatrix)
    >>> print(stoich_instance)
    <PyscesStoich instance>
    >>> print(status)
    0
    """
        stc = PyscesStoich.Stoich(nmatrix)
        status = 0
        return stc, status

    def Stoichiometry_Analyse(self):
        """
    Stoichiometry_Analyse(override=0, load=0)

    Perform a structural analysis of the model's stoichiometry. The method's default
    behavior constructs and analyzes the model from parsed model information.
    By setting `override`, the stoichiometry analysis is performed directly based on
    the current stoichiometric matrix. If `load` is specified, the method attempts
    to load a saved stoichiometry; otherwise, it performs a new stoichiometric
    analysis. The results of the analysis are checked for floating-point errors and
    null-space rank consistency.

    The analysis involves calculating and validating the L and K matrices, checking
    values for stoichiometric precision, and ensuring reliable structural analysis
    results.

    Parameters
    ----------
    override : int, optional
        Override stoichiometric analysis initialization from parsed data (default is 0).
    load : int, optional
        Load a pre-saved stoichiometry instead of performing analysis (default is 0).

    Examples
    --------
    Analyze the model's stoichiometry using the default settings:

    >>> parser = PyscesInputFileParser()
    >>> parser.Stoichiometry_Analyse()
    
    Override stoichiometry parsing with existing matrix:

    >>> parser.Stoichiometry_Analyse(override=1)
    
    Load a pre-saved stoichiometry:

    >>> parser.Stoichiometry_Analyse(load=1)

    Notes
    -----
    If floating point values in the L0 or K0 matrices are close to the specified stoichiometric precision,
    the results may not be reliable. Ensure that stoichiometric precision settings are correct when encountering warnings.
    """
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
            self.nmatrix = self.buildN()
        else:
            print('\nStoichiometric override active\n')
        assert len(self.nmatrix
            ) > 0, """
Unable to generate Stoichiometric Matrix! model has:
%s reactions
%s species
what did you have in mind?
""" % (
            len(self.__reactions__), len(self.__species__))
        self.__nmatrix__ = self.nmatrix
        self.__Nshape__ = self.nmatrix.shape
        self.__structural__, stc_load = self.Stoichiometry_Init(self.nmatrix)
        if not stc_load:
            self.__structural__.stoichiometric_analysis_fp_zero = (self.
                __settings__['stoichiometric_analysis_fp_zero'])
            self.__structural__.stoichiometric_analysis_lu_precision = (self
                .__settings__['stoichiometric_analysis_lu_precision'])
            self.__structural__.stoichiometric_analysis_gj_precision = (self
                .__settings__['stoichiometric_analysis_gj_precision'])
            self.__structural__.AnalyseL()
            self.__structural__.AnalyseK()
        lsmall, lbig = self.__structural__.MatrixValueCompare(self.
            __structural__.lzeromatrix)
        ksmall, kbig = self.__structural__.MatrixValueCompare(self.
            __structural__.kzeromatrix)
        SmallValueError = 0
        if abs(lsmall
            ) < self.__structural__.stoichiometric_analysis_lu_precision * 10.0:
            print(
                '\nWARNING: values in L0matrix are close to stoichiometric precision!'
                )
            print('Stoichiometric LU precision:', self.__structural__.
                stoichiometric_analysis_lu_precision)
            print('L0 smallest abs(value)', abs(lsmall))
            print('Machine precision:', mach_eps)
            SmallValueError = 1
        if abs(ksmall
            ) < self.__structural__.stoichiometric_analysis_lu_precision * 10.0:
            print(
                '\nWARNING: values in K0matrix are close to stoichiometric precision!'
                )
            print('Stoichiometric precision:', self.__structural__.
                stoichiometric_analysis_lu_precision)
            print('K0 smallest abs(value)', abs(ksmall))
            print('Machine precision:', mach_eps)
            SmallValueError = 1
        if SmallValueError:
            input(
                """
Structural Analysis results may not be reliable!!!.

Try change <mod>.__settings__["stoichiometric_analysis_lu_precision"] (see reference manual for details)

	 press any key to continue: """
                )
        if self.__structural__.kzeromatrix.shape[0
            ] != self.__structural__.lzeromatrix.shape[1]:
            print(
                '\nWARNING: the rank calculated by the Kand L analysis methods are not the same!'
                )
            print('\tK analysis calculates the rank as: ' + repr(self.
                __structural__.kzeromatrix.shape[0]))
            print('\tL analysis calculates the rank as: ' + repr(self.
                __structural__.lzeromatrix.shape[1]))
            print(
                'This is not good! Structural Analysis results are not reliable!!!\n'
                )
            assert self.__structural__.kzeromatrix.shape[0] == self.__structural__.lzeromatrix.shape[1], """
StructuralAnalysis Error: rank mismatch"""
        self.__HAS_FLUX_CONSERVATION__ = self.__structural__.info_flux_conserve
        self.__HAS_MOIETY_CONSERVATION__ = (self.__structural__.
            info_moiety_conserve)
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
            self.conservation_matrix_row = (self.__structural__.
                conservation_matrix_row)
            self.conservation_matrix_col = (self.__structural__.
                conservation_matrix_col)
            self.nrmatrix = self.__structural__.nrmatrix
            self.nrmatrix_row = self.__structural__.nrmatrix_row
            self.nrmatrix_col = self.__structural__.nrmatrix_col
            self.__kmatrix__ = copy.copy(self.kmatrix)
            self.__kzeromatrix__ = copy.copy(self.kzeromatrix)
            self.__lmatrix__ = copy.copy(self.lmatrix)
            self.__lzeromatrix__ = copy.copy(self.lzeromatrix)
            self.__nrmatrix__ = copy.copy(self.nrmatrix)
        self.__structural__.species = self.species
        self.__structural__.reactions = self.reactions
        self.Nmatrix = PyscesStoich.StructMatrix(self.__structural__.
            nmatrix, self.__structural__.nmatrix_row, self.__structural__.
            nmatrix_col)
        self.Nmatrix.setRow(self.species)
        self.Nmatrix.setCol(self.reactions)
        self.Nrmatrix = PyscesStoich.StructMatrix(self.__structural__.
            nrmatrix, self.__structural__.nrmatrix_row, self.__structural__
            .nrmatrix_col)
        self.Nrmatrix.setRow(self.species)
        self.Nrmatrix.setCol(self.reactions)
        self.Kmatrix = PyscesStoich.StructMatrix(self.__structural__.
            kmatrix, self.__structural__.kmatrix_row, self.__structural__.
            kmatrix_col)
        self.Kmatrix.setRow(self.reactions)
        self.Kmatrix.setCol(self.reactions)
        self.K0matrix = PyscesStoich.StructMatrix(self.__structural__.
            kzeromatrix, self.__structural__.kzeromatrix_row, self.
            __structural__.kzeromatrix_col)
        self.K0matrix.setRow(self.reactions)
        self.K0matrix.setCol(self.reactions)
        self.Lmatrix = PyscesStoich.StructMatrix(self.__structural__.
            lmatrix, self.__structural__.lmatrix_row, self.__structural__.
            lmatrix_col)
        self.Lmatrix.setRow(self.species)
        self.Lmatrix.setCol(self.species)
        self.L0matrix = PyscesStoich.StructMatrix(self.__structural__.
            lzeromatrix, self.__structural__.lzeromatrix_row, self.
            __structural__.lzeromatrix_col)
        self.L0matrix.setRow(self.species)
        self.L0matrix.setCol(self.species)
        if self.__structural__.info_moiety_conserve:
            self.Consmatrix = PyscesStoich.StructMatrix(self.__structural__
                .conservation_matrix, self.__structural__.
                conservation_matrix_row, self.__structural__.
                conservation_matrix_col)
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
