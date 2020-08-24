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

__doc__ = "PySCeS JWS parser module -- uses  PLY 1.5 or newer"

try:
    input = raw_input  # Py2 compatibility
except NameError:
    pass

import os, copy
from .lib import lex
from .lib import yacc
from getpass import getuser
from time import sleep, strftime
from scipy import MachAr

MyMachArr = MachAr()


class JWSParser:
    """JWSParser written by Johann, based on Jannie's lexparse and integrated into PySCeS by brett  -- converts PySCeS (.psc) files to JWS Online (jws) files"""

    ReactionIDs = []  # List of reaction names
    Names = []  # List of all reagent, parameter and function names
    LexErrors = []  # List of lexing errors

    NetworkDict = {}  # Dictionary containing all reaction information
    InitStrings = []  # Initialisation strings
    InitParStrings = []  # Initialisation strings for parameters -- johann new
    InitVarStrings = []  # Initialisation strings for variables -- johann new
    Inits = []  # Initialised entities
    Reagents = []  # All reagents found during parsing of reactions
    VarReagents = []  # Variable reagents that occur in reactions
    FixedReagents = []  # Fixed reagents
    ReacParams = []  # Temporary list of reaction parameters
    InitParams = []  # Initialised parameters
    ParseErrors = []

    mach_spec = MyMachArr
    AllRateEqsGiven = 1  # Flag to check that all rate equations have been given
    Debug = 0

    ##############
    # Build the lexer
    ##############

    # elementary regular expressions used as building blocks
    Int = r'\d+'  # Integer
    Dec = Int + '\.' + Int  # Decimal

    # List of token names
    tokens = (
        'FIXDEC',
        'IRREV',
        #'REAL', # johann -- now build up real in a p function since we want to make exponent explicit
        'INT',
        'DEC',  # johann -- added explicitly since we no longer have REAL token
        'PLUS',
        'MINUS',
        'TIMES',
        'DIVIDE',
        'POWER',
        'LPAREN',
        'RPAREN',
        'EQUALS',
        'COMMA',
        'REACTION_ID',
        'STOICH_COEF',
        'NAME',
        'EXP',
    )  # johann -- new EXP token

    # Simple tokens
    t_IRREV = r'>'
    # t_REAL = Real # johann -- no longer have a real token, now a p function
    t_INT = Int
    t_DEC = Dec  # new DEC token
    t_PLUS = r'\+'
    t_MINUS = r'-'
    t_TIMES = r'\*'
    t_DIVIDE = r'/'
    t_POWER = '\*\*'
    t_LPAREN = r'\('
    t_RPAREN = r'\)'
    t_EQUALS = r'='
    t_COMMA = r','
    t_ignore = ' \t\r'  # Ignore spaces and tabs --- and windows return - brett 20040229

    def t_comment(self, t):
        r'\#.+\n'  # Match from # to newline
        t.lineno += 1  # Increment line number

    def t_newline(self, t):
        r'\n+'  # Match newline
        t.lineno += len(t.value)  # Increment with number of consecutive newlines

    def t_EXP(self, t):  # johann -- need separate EXP token to replace for Mathematica
        r'\d+\.?\d*[E|e][\+|\-]?'  # define EXP token merely as digits[.]digits[E|e][+|-]
        t.type = 'EXP'  # parse final integer separately in 'Real' p-function to remove leading zeros
        t.value = t.value.replace('e', ' 10^')
        t.value = t.value.replace('E', ' 10^')
        return t

    def t_FIXDEC(self, t):
        r'FIX:'
        t.type = 'FIXDEC'
        t.value = 'FIX:'
        return t

    def t_REACTION_ID(self, t):
        r'[a-zA-Z]\w*:'  # Match any letter followed by zero or more characters
        # in [a-zA-Z0-9_] up to a colon
        t.type = 'REACTION_ID'
        if t.value[0] == 'v' and len(t.value) > 1:
            t.value = t.value[
                1:
            ]  # remove initial 'v' if present to avoid constructions like 'v[vR1]'
        t.value = (
            'v[' + t.value[:-1] + ']'
        )  # remove the colon and add v[] for JWS -- johann
        if t.value in self.ReactionIDs:
            self.LexErrors.append(('Duplicate ReactionID ', t.lineno, t.value, t.type))
        else:
            self.ReactionIDs.append(t.value)
        return t

    def t_STOICH_COEF(self, t):
        r'\{\d+\}|\{\d+\.\d+\}'
        t.type = 'STOICH_COEF'
        t.value = t.value[1:-1]
        return t

    def t_NAME(self, t):
        r'[a-zA-Z][\w]*'  # Match any letter followed by zero or characters in the set [a-zA-Z0-9_]
        if (t.value + '[t]' not in self.Names) and (
            t.value not in self.FuncNames
        ):  # Only add to list if absent in list
            # self.Names.append(t.value)
            self.Names.append(t.value + '[t]')  # -- johann
            # print self.Names[-1]
        # hack! - brett
        if (
            t.value not in self.FuncNames
        ):  # make class attributes, ignore function names
            # print 't value before', t.value
            gt = t.value + '[t]'
            t.value = gt
        # print 't value after', t.value
        t.type = 'NAME'
        return t

    def t_error(self, t):
        self.LexErrors.append(('Lexer error ', t.lineno, t.value, t.type))
        print('Illegal character, Line ' + str(t.lineno) + ' :' + str(t.value[0]))
        t.skip(1)

    ##############
    # The parser #
    ##############

    FuncNames = (
        'acos',
        'asin',
        'atan',
        'atan2',
        'ceil',
        'cos',
        'cosh',
        'exp',
        'fabs',
        'floor',
        'fmod',
        'frexp',
        'hypot',
        'ldexp',
        'log',
        'log10',
        'modf',
        'pow',
        'sin',
        'sinh',
        'sqrt',
        'tan',
        'tanh',
    )

    precedence = (
        ('left', 'PLUS', 'MINUS'),
        ('left', 'TIMES', 'DIVIDE'),
        ('left', 'POWER'),
        ('right', 'UMINUS'),
    )

    def Show(self, name, tok):
        if self.Debug:
            print(name, tok)

    def p_error(self, t):
        self.ParseErrors.append(('Syntax error ', t.lineno, t.value, t.type))
        print('Syntax error, Line ' + str(t.lineno) + ' : ' + str(t.value))
        tok = yacc.token()
        while tok and tok.type != 'REACTION_ID':
            tok = yacc.token()
        return tok

    def p_model(self, t):
        '''Model : Statement
                 | Model Statement '''

        self.Show('Model', t[0])

    def p_statement(self, t):
        '''Statement : Fixed
                     | ReactionLine
                     | Initialise'''
        self.Show('Statement', t[0])

    def p_fixed(self, t):
        '''Fixed : FIXDEC FixList'''

        self.Show('Fixed:', t[0])

    def p_fixedreagents(self, t):
        '''FixList : NAME
                   | NAME FixList'''
        if t[1] != None:
            self.FixedReagents.append(t[1][:-3])  # johann -- remove [t] off end
        t[0] = [t[1]]
        try:
            t[0] += t[2]
        except:
            pass
        self.Show('FixList', t[0])

    def p_initialise(self, t):
        '''Initialise : NAME EQUALS Expression'''
        t[1] = t[1][:-3] + '[0]'  # johann 20050302 -- Mathematica initialisation
        t[0] = t[1] + t[2] + t[3]
        ##    del temp
        self.InitStrings.append(t[0].replace('=', ' = '))
        self.Inits.append(t[1])
        self.Show('Initialisation', t[0])

    def p_reaction_line(self, t):
        '''ReactionLine : REACTION_ID ReactionEq
                        | REACTION_ID ReactionEq Expression'''

        # global self.AllRateEqsGiven, ReacParams
        ReacID = t[1]
        if ReacID in self.NetworkDict:
            self.ParseErrors.append(('Duplicate Reaction ', t.lineno, ReacID, None))
        self.NetworkDict[ReacID] = {}  # Reaction dictionary for ReacID
        self.NetworkDict[ReacID]['Reagents'] = {}  # Reagent dictionary within ReacID

        # brett: if an index exists sum the coefficients instead of adding a new one
        # this seems to deal with multiple definitions like X + X > Y and 2{X} + Y > Z + X
        for i in t[2][
            0
        ]:  # First tuple member of ReactionEq contains list of (name,stoichcoef)
            if i[0] in self.NetworkDict[ReacID]['Reagents']:
                self.NetworkDict[ReacID]['Reagents'][i[0]] = (
                    self.NetworkDict[ReacID]['Reagents'][i[0]] + i[1]
                )
            else:
                self.NetworkDict[ReacID]['Reagents'][i[0]] = i[
                    1
                ]  # Key for reagent with stoichcoef value
        killList = []
        # brett: however for the case of X + Y > Y + Z where the sum of the coefficients
        # is zero we can delete the key (Y) out of the reaction list altgether (hopefully!)
        for i in self.NetworkDict[ReacID]['Reagents']:
            if (
                abs(self.NetworkDict[ReacID]['Reagents'][i])
                < self.mach_spec.eps * 100.0
            ):
                killList.append(i)
                # print self.mach_spec.eps*100.0, self.NetworkDict[ReacID]['Reagents']
        # print killList, self.NetworkDict[ReacID]['Reagents']
        # brett: and the easiest way of doing this is putting the zero keys in a list
        # and deleting them out of the dictionary
        if len(killList) != 0:
            for i in killList:
                del self.NetworkDict[ReacID]['Reagents'][i]
        # print killList, self.NetworkDict[ReacID]['Reagents']

        self.NetworkDict[ReacID]['Type'] = t[2][
            1
        ]  # Second tuple member of ReactionEq contains type
        try:  # Save rate equation and create parameter list
            self.NetworkDict[ReacID]['RateEq'] = t[3]
            self.NetworkDict[ReacID]['Params'] = self.ReacParams
            self.ReacParams = []  # Reset global self.ReacParams list
        except:
            self.NetworkDict[ReacID]['RateEq'] = ''
            self.NetworkDict[ReacID]['Params'] = []
            self.AllRateEqsGiven = 0  # Set global flag to false
        self.Show('ReactionLine', t[0])
        self.Show('t1', t[1])
        self.Show('t2', t[2])
        self.Show('t3', t[3])

    def p_reaction_eq(self, t):
        '''ReactionEq : LeftHalfReaction EQUALS RightHalfReaction
                      | LeftHalfReaction IRREV  RightHalfReaction'''

        ReacType = ''
        if t[2] == '=':
            ReacType = 'Rever'
        elif t[2] == '>':
            ReacType = 'Irrev'
        t[0] = (t[1] + t[3], ReacType)
        self.Show('ReactionEq', t[0])

    def p_left_half_reaction(self, t):
        ''' LeftHalfReaction : SubstrateTerm
                             | SubstrateTerm PLUS LeftHalfReaction'''

        # Make a list of substrate terms
        t[0] = [t[1]]

        try:
            t[0] += t[3]
        except:
            pass
        # brett
        # print "lhr ", t[0]
        self.Show('LeftHalfReaction', t[0])

    def p_right_half_reaction(self, t):
        ''' RightHalfReaction : ProductTerm
                              | ProductTerm PLUS RightHalfReaction'''

        # Make a list of product terms
        t[0] = [t[1]]
        try:
            t[0] += t[3]
        except:
            pass
        # brett
        # print "rhr ", t[0]
        self.Show('RightHalfReaction', t[0])

    def p_substrate_term(self, t):
        '''SubstrateTerm : STOICH_COEF NAME
                         | NAME'''

        # Make tuple of NAME and stoichiometric coefficient
        # (< 0 because substrate)
        try:
            t[0] = (t[2], -float(t[1]))
            if t[2] not in self.Reagents:
                self.Reagents.append(t[2])
        except:
            t[0] = (t[1], -1.0)
            if t[1] not in self.Reagents:
                self.Reagents.append(t[1])
        self.Show('SubstrateTerm', t[0])

    def p_product_term(self, t):
        '''ProductTerm : STOICH_COEF NAME
                       | NAME'''

        # Make tuple of NAME and stoichiometric coefficient
        # (> 0 because product)
        try:
            t[0] = (t[2], float(t[1]))
            if t[2] not in self.Reagents:
                self.Reagents.append(t[2])
        except:
            t[0] = (t[1], 1.0)
            if t[1] not in self.Reagents:
                self.Reagents.append(t[1])
        self.Show('ProductTerm', t[0])

    def p_rate_eq(self, t):
        '''Expression : Expression PLUS Expression
                      | Expression MINUS Expression
                      | Expression TIMES Expression
                      | Expression DIVIDE Expression
                      | Power
                      | Number
                      | Func'''
        # |UMINUS : add if the
        #  alternative for p_uminus is used

        if len(t.slice) == 4:
            t[0] = t[1] + t[2] + t[3]
        else:
            t[0] = t[1]

    def p_power(self, t):
        '''Power : Expression POWER Expression'''

        t[0] = (
            'Power[' + t[1] + ',' + t[3] + ']'
        )  # changed to Mathematica notation -- johann

    def p_uminus(self, t):
        '''Expression : MINUS Expression %prec UMINUS'''
        # Alternative '''UMINUS : MINUS Expression'''

        t[0] = t[1] + t[2]

    def p_number(self, t):
        '''Number : Real
                  | INT
                  | DEC
                  | NAME'''

        # Build list of entities
        try:
            float(t[1])  # check for a number
        except:
            if (
                (t[1] not in self.FuncNames)
                and (t[1] not in self.ReacParams)
                and (' 10^' not in t[1])
            ):
                # ignore function names, duplications and exponentials
                self.ReacParams.append(t[1])
                # self.ReacParams.append('self.' + t[1])
        t[0] = t[1]

    def p_real(self, t):
        '''Real : EXP INT'''
        loop = 1
        while loop == 1:  # remove leading zeros from exponent
            if t[2][0] == '0' and len(t[2]) > 1:
                t[2] = t[2][1:]
            else:
                loop = 0
        t[0] = t[1] + t[2]

    def p_function(self, t):
        '''Func : LPAREN ArgList RPAREN
                | NAME LPAREN ArgList RPAREN'''

        try:
            t[0] = t[1] + t[2] + t[3] + t[4]
        except:
            t[0] = t[1] + t[2] + t[3]

    def p_arglist(self, t):
        '''ArgList : Expression
                   | Expression COMMA Expression'''

        t[0] = t[1]
        try:
            t[0] += t[2] + t[3]
        except:
            pass

    ############################################
    # end of lexer and parser definitions
    ############################################

    def psc2jws(self, File, indir=None, outdir=None, quiet=1, debug=0):
        """
        psc2jws(File,indir=None,outdir=None,quiet=1,debug=0)
        Convert a PySCeS (.psc) file to a JWS Online (.jws) file. Call with the input file name, note the input (indir) and output (outdir) can optionally be specified.

        Arguments:
        =========
        File: PSC input file
        indir [default=None]: directory of PSC file
        outdir [default=None]: output directory for JWS file
        quiet [default=1]: turn lex/parse noise on/off
        debug [default=0]:  optionally output debug information

        """
        if indir == None:
            indir = os.getcwd()
        if outdir == None:
            outdir = os.getcwd()
        if os.path.exists(os.path.join(indir, File)) and File[-4:] == '.psc':
            go = 1
        else:
            print('\nIgnoring non-PySCeS model file: ' + os.path.join(indir, File))
            go = 0

        if go == 1:
            # clean up the modules
            reload(lex)  # brett's bugbear code these have to be here ALWAYS!!
            reload(yacc)
            # clean up the instance
            self.ReactionIDs = []  # List of reaction names
            self.Names = []  # List of all reagent, parameter and function names
            self.LexErrors = []  # List of lexing errors
            self.NetworkDict = {}  # Dictionary containing all reaction information
            self.InitStrings = []  # Initialisation strings

            self.Inits = []  # Initialised entities
            self.Reagents = []  # All reagents found during parsing of reactions
            self.FixedReagents = []  # Fixed reagents
            self.ReacParams = []  # Temporary list of reaction parameters
            self.ParseErrors = []

            self.InitParStrings = (
                []
            )  # Initialisation strings for parameters -- johann new
            self.InitVarStrings = (
                []
            )  # Initialisation strings for variables -- johann new
            self.VarReagents = []  # Variable reagents that occur in reactions
            self.InitParams = []  # Initialised parameters

            print('\nParsing file: ' + os.path.join(indir, File))

            Data = open(os.path.join(indir, File), 'r')
            Model = Data.read()
            Data.close()

            self.Debug = debug
            self.AllRateEqsGiven = (
                1  # Flag to check that all rate equations have been given
            )

            # try and find a temporary workspace or use cwd
            if 'TMP' in os.environ:
                tempDir = os.environ['TMP']
            elif 'TEMP' in os.environ:
                tempDir = os.environ['TEMP']
            else:
                tempDir = os.getcwd()

            os.chdir(tempDir)

            # fix filenames for intermediary files - brett
            if not File[:-4].isalnum():
                FileL = list(File)
                FileT = ''
                for let in FileL:
                    if let.isalnum():
                        FileT += let

                # instantiate the lexer and parser
                self.debugfile = '_jws' + FileT[:-3] + ".dbg"
                self.tabmodule = '_jws' + FileT[:-3] + "_" + "parsetab"
            else:
                self.debugfile = '_jws' + File[:-4] + ".dbg"
                self.tabmodule = '_jws' + File[:-4] + "_" + "parsetab"

            if self.Debug:
                print(self.tabmodule)
                print(self.debugfile)

            lex.lex(module=self, debug=self.Debug)
            lex.input(Model)
            yacc.yacc(
                module=self,
                debug=self.Debug,
                debugfile=self.debugfile,
                tabmodule=self.tabmodule,
            )

            os.chdir(outdir)

            while 1:
                tok = lex.token()
                if not tok:
                    break
            if self.LexErrors != []:
                print('self.LexErrors = ', self.LexErrors, '\n')

            while 1:
                p = yacc.parse(Model)
                if not p:
                    break

            # we have the dictionary get rid of this stuff
            del Model, p

            # Create list of variable reagents and remove '[t]' from fixed reagents
            for i in range(
                len(self.Reagents)
            ):  # johann -- new construction otherwise list elements not replaced
                if self.Reagents[i][:-3] not in self.FixedReagents:
                    self.VarReagents.append(self.Reagents[i])
                if self.Reagents[i][:-3] in self.FixedReagents:
                    self.Reagents[i] = self.Reagents[i][:-3]

            # Create list of initialised parameters
            for i in range(len(self.Inits)):  # johann -- reworked extensively
                if self.Inits[i][:-3] + '[t]' not in self.VarReagents:
                    self.InitStrings[i] = self.InitStrings[i].replace('[0]', '')
                    self.InitStrings[i] = self.InitStrings[i].replace(
                        '[t]', ''
                    )  # capture params initialised i.t.o. other params
                    self.Inits[i] = self.Inits[i][:-3]
                    self.InitParams.append(self.Inits[i])
                    self.InitParStrings.append(self.InitStrings[i])
                elif self.Inits[i][:-3] + '[t]' in self.VarReagents:
                    self.InitVarStrings.append(self.InitStrings[i])

            # In self.NetworkDict, clean rate equation parameter list of variables that occur in that reaction
            # Add FixedReagent to Params even if not a parameter in rate eqn (requirement to add '$' below)
            for id in list(self.NetworkDict.keys()):
                for reag in self.VarReagents:
                    if reag in self.NetworkDict[id]['Params']:
                        self.NetworkDict[id]['Params'].remove(reag)
                for reag in self.FixedReagents:
                    if (
                        reag + '[t]' in list(self.NetworkDict[id]['Reagents'].keys())
                    ) and (reag not in self.NetworkDict[id]['Params']):
                        self.NetworkDict[id]['Params'].append(reag + '[t]')

            # Warn if no reagents have been fixed
            if self.FixedReagents == []:
                print('Warning: No reagents have been fixed')
            else:  # Warn if a fixed reagent does not occur in a reaction equation
                for reag in self.FixedReagents:
                    if reag not in self.Reagents:
                        print(
                            'Warning: '
                            + reag
                            + ' (fixed) does not occur in any reaction'
                        )

            # Check whether all parameters have been initialised
            # johann -- remove [t] from params
            for id in list(self.NetworkDict.keys()):
                for i in range(len(self.NetworkDict[id]['Params'])):
                    self.NetworkDict[id]['Params'][i] = self.NetworkDict[id]['Params'][
                        i
                    ][:-3]
                    if self.NetworkDict[id]['Params'][i] not in self.InitParams:
                        print(
                            'Warning: Parameter '
                            + self.NetworkDict[id]['Params'][i]
                            + ' has not been initialised'
                        )

            # Check whether all variable reagents have been initialised
            for reag in self.VarReagents:
                if reag[:-3] + '[0]' not in self.Inits:
                    print('Warning: Variable ' + reag + ' has not been initialised')

            # Check that all initialised parameters actually occur in self.Inits
            known = 0
            for param in self.InitParams:
                for id in list(self.NetworkDict.keys()):
                    if param in self.NetworkDict[id]['Params']:
                        known = 1
                        break
                    else:
                        known = 0
                if not known:
                    print(
                        'Warning: '
                        + param
                        + ' has been initialised but does not occur in any rate equation'
                    )

            # clean up rate equations in self.NetworkDict to remove [t] for Params
            # clean up Reagents to remove [t] and add $ for fixed
            for id in list(self.NetworkDict.keys()):
                for param in self.NetworkDict[id]['Params']:
                    self.NetworkDict[id]['RateEq'] = self.NetworkDict[id][
                        'RateEq'
                    ].replace(param + '[t]', param)
                for reag in list(self.NetworkDict[id]['Reagents'].keys()):
                    if reag[:-3] in self.NetworkDict[id]['Params']:
                        saveval = self.NetworkDict[id]['Reagents'].pop(reag)
                        self.NetworkDict[id]['Reagents']['$' + reag[:-3]] = saveval
                    else:
                        saveval = self.NetworkDict[id]['Reagents'].pop(reag)
                        self.NetworkDict[id]['Reagents'][reag[:-3]] = saveval

            # output errors
            if self.ParseErrors != []:
                print('Parse errors occurred: ', self.ParseErrors)

            # debugging
            if debug:
                print('\n\n\n')
                print('\nself.ReactionIDs: ', self.ReactionIDs)
                print('\nself.NetworkDict: ', self.NetworkDict)
                print('\nself.Names: ', self.Names)
                print('\nself.Inits: ', self.Inits)
                print('\nself.InitStrings: ', self.InitStrings)
                print('\nself.InitParStrings: ', self.InitParStrings)
                print('\nself.InitVarStrings: ', self.InitVarStrings)
                print('\nself.InitParams: ', self.InitParams)
                print('\nself.Reagents: ', self.Reagents)
                print('\nself.FixedReagents: ', self.FixedReagents)
                print('\nself.VarReagents: ', self.VarReagents)
                print('\nParseErrors: ', self.ParseErrors)

            # now write the jws output file
            filename = File[:-4]
            filename = self.chkjws(filename)
            go = 0
            loop = 0
            filex = ''
            while loop == 0:
                try:
                    filex = os.path.join(outdir, filename)
                    f = open(filex, 'r')
                    f.close()
                    inp = input('\nFile "' + filex + '" exists.\nOverwrite? ([y]/n) ')
                    if inp == 'y' or inp == '':
                        go = 1
                        loop = 1
                    elif inp == 'n':
                        filename = input(
                            '\nFile "' + filename + '" exists. Enter a new filename: '
                        )
                        go = 1
                        filex = os.path.join(outdir, filename)
                        filename = self.chkjws(filename)
                    else:
                        print('\nInvalid input')
                except:
                    print('\nFile "' + filex + '" does not exist, proceeding...')
                    loop = 1
                    go = 1
            if go == 1:
                try:
                    UseR = getuser()
                except:
                    UseR = ''

                outFile = open(filex, 'w')
                header = ''
                # header +=  '############################################################\n'
                header += '# JWS model input file \n'
                header += (
                    '# Generated by PySCeS ('
                    + __version__
                    + ') (http://pysces.sourceforge.net) \n'
                )
                header += '# Pysces input file: ' + File + '\n'
                header += (
                    '# This file generated: '
                    + strftime("%a, %d %b %Y %H:%M:%S")
                    + ' by '
                    + UseR
                    + ' \n'
                )
                header += (
                    '###########################################################\n\n'
                )
                outFile.write(header)

                # modelname
                modelname = File[:-4]
                outFile.write('begin name\n' + modelname + '\nend name\n\n')

                # reactions and rate equations
                reaction_list = []
                rateeq_list = []

                nd = self.NetworkDict
                reaclist = copy.copy(
                    list(nd.keys())
                )  # johann -- to sort self.ReactionIDs neatly ;-)
                reaclist.sort()
                for key in reaclist:  # key = reaction name
                    reagL = []
                    reagR = []
                    Req = copy.copy(nd[key]['RateEq'])
                    for reagent in nd[key]['Reagents']:
                        if nd[key]['Reagents'][reagent] > 0:
                            reagR.append(
                                '{'
                                + str(abs(nd[key]['Reagents'][reagent]))
                                + '}'
                                + reagent
                            )
                        elif nd[key]['Reagents'][reagent] < 0:
                            reagL.append(
                                '{'
                                + str(abs(nd[key]['Reagents'][reagent]))
                                + '}'
                                + reagent
                            )
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
                    symbol = ' = '
                    reaction_list.append(key + '\t' + substring + symbol + prodstring)
                    rateeq_list.append(key + ' = ' + Req)
                outFile.write('begin reactions\n')
                for x in reaction_list:
                    outFile.write(x + '\n')
                outFile.write('end reactions\n\n')
                outFile.write('begin rate equations\n')
                for x in rateeq_list:
                    outFile.write(x + '\n')
                outFile.write('end rate equations\n\n')

                # parameters
                outFile.write('begin parameters\n')
                for x in self.InitParStrings:
                    outFile.write(x + '\n')
                outFile.write('end parameters\n\n')

                # species initial values
                outFile.write('begin initial conditions\n')
                for x in self.InitVarStrings:
                    outFile.write(x + '\n')
                outFile.write('end initial conditions\n\n')

                # close output file
                outFile.close()

                # print to stdout if quiet is set to zero
                if quiet == 0:
                    print('\nModel name: ' + modelname)

                    print("\nReactions:")
                    for x in reaction_list:
                        print(x)

                    print("\nRate Equations:")
                    for x in rateeq_list:
                        print(x)

                    print('\nParameters:')
                    for x in self.InitParStrings:
                        print(x)

                    print('\nSpecies Initial Values:')
                    for x in self.InitVarStrings:
                        print(x)

    def chkjws(self, File):
        """
        chkjws(File)

        Checks if a filename has a .jws extension and adds one to the returned filename if needed

        Arguments:
        =========
        File: the filename to check

        """
        try:
            if File[-4:] == '.jws':
                pass
            else:
                print('Assuming extension is .jws')
                File += '.jws'
        except:
            print('Chkjws error')
        return File


if __name__ == '__main__':
    import os, sys
    from time import sleep

    inDiR = 'c://mypysces//pscmodels'
    outDiR = 'c://mypysces//jws'

    jwp = JWSParser()
    for mod in os.listdir(inDiR):
        jwp.psc2jws(mod, indir=inDiR, outdir=outDiR, quiet=1, debug=0)

    psp = PySCeSParser(debug=0)
