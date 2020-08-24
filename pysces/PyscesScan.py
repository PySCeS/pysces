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

__doc__ = """
          PyscesScan
          ----------

          PySCeS classes for continuations and multi-dimensional parameter scans
          """

import scipy
from pysces.PyscesUtils import TimerBox

__psyco_active__ = 0


class Scanner(object):
    """
    Arbitrary dimension generic scanner. This class is initiated with
    a loaded PySCeS model and then allows the user to define scan parameters
    see self.addScanParameter() and user output see self.addUserOutput().
    Steady-state results are always stored in self.SteadyStateResults while
    user output can be found in self.UserOutputResults - brett 2007.
    """

    genOn = True
    quietRun = False
    _MODE_ = 'state'
    HAS_USER_OUTPUT = False
    HAS_STATE_OUTPUT = True
    nan_on_bad_state = True
    MSG_PRINT_INTERVAL = 500
    __AnalysisModes__ = ('state', 'elasticity', 'mca', 'stability', 'null')
    invalid_state_list = None
    invalid_state_list_idx = None

    def __init__(self, mod):
        try:
            N = mod.nmatrix
        except Exception as ex:
            print(
                '\nSCANNER: Please load a model <i.e. mod. doLoad()> before trying to use it.'
            )

        self.GenDict = {}
        self.GenOrder = []
        self.ScanSpace = []
        self.mod = mod
        self.SteadyStateResults = []
        self.UserOutputList = []
        self.UserOutputResults = []
        self.scanT = TimerBox()

    def testInputParameter(self, name):
        """
        This tests whether a str(name) is an attribute of the model
        """
        if hasattr(self.mod, name):
            return True
        else:
            return False

    def resetInputParameters(self):
        """
        Just remembered what this does, I think it resets the input model
        parameters after a scan run.
        """
        for key in self.GenDict:
            setattr(self.mod, key, self.GenDict[key][6])

    def addUserOutput(self, *kw):
        """
        Add output parameters to the scanner as a collection of one or more
        string arguments ('O1','O2','O3', 'On'). These are evaluated at
        each iteration of the scanner and stored in the self.UserOutputResults
        array. The list of output is stored in self.UserOutputList.
        """
        output = [attr for attr in kw if type(attr) == str]
        # print output
        self.HAS_USER_OUTPUT = True
        ModeList = list(self.__AnalysisModes__)
        MaxMode = 0
        for attr in output:
            if attr[:2] == 'ec' and self._MODE_ != 'null':
                cM = ModeList.index('elasticity')
                if cM > MaxMode:
                    self._MODE_ = 'elasticity'
                    MaxMode = cM
            elif attr[:2] == 'cc' and self._MODE_ != 'null':
                cM = ModeList.index('mca')
                if cM > MaxMode:
                    self._MODE_ = 'mca'
                    MaxMode = cM
            elif attr[:6] == 'lambda' and self._MODE_ != 'null':
                cM = ModeList.index('stability')
                if cM > MaxMode:
                    self._MODE_ = 'stability'
                    MaxMode = cM
        print('MaxMode', MaxMode)
        self.UserOutputList = output

    def addScanParameter(self, name, start, end, points, log=False, follower=False):
        """
        Add a parameter to scan (an axis if you like) input is:

        - str(name) = model parameter name
        - float(start) = lower bound of scan
        - float(end) = upper bound of scan
        - int(points) = number of points in scan range
        - bool(log) = Use a logarithmic (base10) range
        - bool(follower) = Scan parameters can be leaders i.e. an independent axis or a "follower" which moves synchronously with the previously defined parameter range.

        The first ScanParameter cannot be a follower.
        """
        offset = 1
        assert self.testInputParameter(name), (
            '\nSCANNER: Model does not have an attribute \"%s\" \n' % name
        )
        assert (
            name not in self.GenOrder
        ), '\nSCANNER: This operation is currently not allowed\n'
        if not len(self.GenOrder) == 0:
            for el in self.GenOrder:
                if self.GenDict[el][4] == False:  # test if not follower
                    offset = offset * self.GenDict[el][2]  # increment offset
            if follower == True:
                prevpar = self.GenOrder[-1]  # previous parameter, i.e. leader
                offset = self.GenDict[prevpar][5]  # don't increment for follower
                if points != self.GenDict[prevpar][2]:
                    print(
                        'SCANNER: Follower parameter needs to iterate over same number of\npoints as leader...resetting points.'
                    )
                    points = self.GenDict[prevpar][2]
        else:
            if follower == True:
                follower = False
                print(
                    'SCANNER: Inner range cannot be a follower ... resetting to leader'
                )
        self.GenDict.setdefault(
            name, (start, end, points, log, follower, offset, getattr(self.mod, name))
        )
        setattr(self, name + '_test', self.stepGen(offset))
        setattr(self, name, self.rangeGen(name, start, end, points, log))
        if name not in self.GenOrder:
            self.GenOrder.append(name)

    def makeRange(self, start, end, points, log):
        """
        Should be pretty self evident it defines a range:

        - float(start)
        - float(end)
        - int(points)
        - bool(log)
        """
        if log:
            rng = scipy.logspace(scipy.log10(start), scipy.log10(end), points)
        else:
            rng = scipy.linspace(start, end, points)
        return rng

    def rangeGen(self, name, start, end, points, log):
        """
        This is where things get more interesting. This function creates
        a cycling generator which loops over a parameter range.

        - *parameter* name
        - *start* value
        - *end* value
        - *points*
        - *log* scale
        """
        range = self.makeRange(start, end, points, log)
        cnt = 0
        while self.genOn:
            if next(getattr(self, name + '_test')):
                if cnt >= len(range) - 1:
                    cnt = 0
                else:
                    cnt += 1
            yield range[cnt]

    def stepGen(self, offset):
        """
        Another looping generator function. The idea here is to create
        a set of generators for the scan parameters. These generators then
        all fire together and determine whether the range generators should
        advance or not. Believe it or not this dynamically creates the matrix
        of parameter values to be evaluated.
        """
        val = 0
        while self.genOn:
            if val >= offset:
                val = 0
                yield True
            else:
                yield False
            val += 1

    def __nextStep__(self):
        """
        Fire all step generators
        """
        return [next(getattr(self, genN + '_test')) for genN in self.GenOrder]

    def __nextValue__(self):
        """
        Fire all range generators
        """
        return [next(getattr(self, genN)) for genN in self.GenOrder]

    def setModValue(self, name, value):
        """
        An easy one, assign value to name of the instantiated PySCeS model attribute
        """
        assert hasattr(self.mod, name), (
            '\nModel does not have an attribute: %s \n' % name
        )
        setattr(self.mod, name, float(value))

    def RunAgain(self):
        """
        While it is impossible to change the generator/range structure of a scanner
        (just build another one) you can 'in principle' change the User Output and run
        it again.
        """
        self.Run(ReRun=True)

    def Run(self, ReRun=False):
        """
        Run the parameter scan
        """
        self.scanT.normal_timer('RUN')
        self._MODE_ = self._MODE_.lower()
        assert self._MODE_ in self.__AnalysisModes__, (
            '\nSCANNER: \"%s\" is not a valid analysis mode!' % self._MODE_
        )
        if self.quietRun:
            self.mod.SetQuiet()
        else:
            self.mod.SetLoud()
        if ReRun:
            self.ScanSpace = []
            self.UserOutputResults = []
            self.SteadyStateResults = []
        self.invalid_state_list = []
        self.invalid_state_list_idx = []

        if self.nan_on_bad_state:
            self.mod.__settings__["mode_state_nan_on_fail"] = True
            # self.mod.mode_state_nan_on_fail = True
        Tsteps = 1
        for gen in self.GenOrder:
            if self.GenDict[gen][4] == False:  # don't increase Tsteps for follower
                Tsteps *= self.GenDict[gen][2]
        print(next(self.scanT.RUN))
        analysis_counter = 0
        print('SCANNER: Tsteps', Tsteps)
        pcntr = 0
        for step in range(Tsteps):
            pars = self.__nextValue__()
            # if pars not in self.ScanSpace:
            self.ScanSpace.append(pars)
            for par in range(len(self.GenOrder)):
                self.setModValue(self.GenOrder[par], pars[par])
            analysis_counter += 1
            self.Analyze()
            if not self.mod.__StateOK__:
                self.invalid_state_list.append(pars)
                self.invalid_state_list_idx.append(analysis_counter)
            self.StoreData()
            pcntr += 1
            if pcntr >= self.MSG_PRINT_INTERVAL:
                print('\t', analysis_counter, next(self.scanT.RUN))
                pcntr = 0
        self.ScanSpace = scipy.array(self.ScanSpace)
        self.SteadyStateResults = scipy.array(self.SteadyStateResults)
        self.UserOutputResults = scipy.array(self.UserOutputResults)
        self.resetInputParameters()
        if self.nan_on_bad_state:
            self.mod.mode_state_nan_on_fail = False
        print("\nSCANNER: %s states analysed\n" % analysis_counter)
        if len(self.invalid_state_list) > 0:
            print('Bad steady states encountered at:\n')
            print(self.invalid_state_list)
            print(self.invalid_state_list_idx)

    def Analyze(self):
        """
        The analysis method, the mode is automatically set by the self.addUserOutput()
        method but can be reset by the user.
        """
        if self._MODE_ == 'state':
            self.mod.doState()
        elif self._MODE_ == 'elasticity':
            self.mod.doElas()
        elif self._MODE_ == 'mca':
            self.mod.doMca()
        elif self._MODE_ == 'stability':
            self.mod.doEigenMca()
        elif self._MODE_ == 'null':
            self.HAS_STATE_OUTPUT = False
        self.mod.User_Function()

    def StoreData(self):
        """
        Internal function which concatenates and stores the data generated by Analyze.
        """
        if self.HAS_STATE_OUTPUT and self._MODE_ != 'null':
            self.SteadyStateResults.append(
                scipy.hstack((self.mod.state_species, self.mod.state_flux))
            )
        if self.HAS_USER_OUTPUT and self._MODE_ != 'null':
            self.UserOutputResults.append(
                scipy.array([getattr(self.mod, res) for res in self.UserOutputList])
            )

    def getOutput(self):
        """
        Will be the new output function.
        """
        pass

    def getResultMatrix(self, stst=False, lbls=False):
        """
        Returns an array of result data. I'm keepin this for backwards compatibility but
        it will be replaced by a getOutput() method when this scanner is updated to use
        the new data_scan object.

        - *stst* add steady-state data to output array
        - *lbls* return a tuple of (array, column_header_list)

        If *stst* is True output has dimensions [scan_parameters]+[state_species+state_flux]+[Useroutput]
        otherwise [scan_parameters]+[Useroutput].
        """
        output_array = None
        labels = []
        if stst:
            if self.HAS_USER_OUTPUT:
                output_array = scipy.hstack(
                    [self.ScanSpace, self.SteadyStateResults, self.UserOutputResults]
                )
                labels = (
                    self.GenOrder
                    + list(self.mod.species)
                    + list(self.mod.reactions)
                    + self.UserOutputList
                )
            else:
                output_array = scipy.hstack([self.ScanSpace, self.SteadyStateResults])
                labels = (
                    self.GenOrder + list(self.mod.species) + list(self.mod.reactions)
                )
        else:
            output_array = scipy.hstack([self.ScanSpace, self.UserOutputResults])
            labels = self.GenOrder + self.UserOutputList
        if lbls:
            return output_array, labels
        else:
            return output_array


class PITCONScanUtils(object):
    """
    Static Bifurcation Scanning utilities using PITCON, call with loaded model object.
    Hopefully nobody else was trying to use the older class as it was horrible. This new
    one is is leaner, meaner and pretty cool ;-)
    """

    model = None
    timer = None
    pitcon_range_low = None
    pitcon_range_high = None
    pitcon_scan_density = None
    pitcon_scan_parameter = None
    pitcon_scan_parameter_3d = None
    pitcon_res = None

    user_output = None
    res_user = None

    res_idx = None
    res_metab = None
    res_flux = None
    res_eigen = None

    def __init__(self, model):
        """
        Instantiate PITCONScanUtils with a loaded PySCeS model instance.
        """
        self.model = model
        self.timer = TimerBox()

    def setUserOuput(self, *args):
        """
        Set the user output required as n string arguments.
        """
        self.user_output = [str(arg) for arg in args]

    def runContinuation(
        self, parameter, low, high, density, par3d=None, logrange=True, runQuiet=True
    ):
        """
        Run the continuation using the following parameters:

        Args:

        - parameter = str(the parameter to be scanned)
        - low = float(lower bound)
        - high = float(upper bound)
        - density = int(the number of initial points)
        - par3d = float(extra 3d parameter to insert into the output array) this parameter is not set ONLY used in output
        - logrange = boolean [default = True], if True generate the result using logspace(log10(low), log10(high), density) otherwise use a linear range
        - runQuiet = boolean [default = True], if True do not display intermediate results to screen, disable for debugging

        After running the continuation the results are stored in numpy arrays

        - mod.res_idx  = scan parameter values (and optionally par3d)
        - mod.res_metab = steady-state species concentrations
        - mod.res_flux = steady-state flux values

        """

        self.pitcon_scan_density = density
        self.pitcon_scan_parameter = parameter
        self.pitcon_scan_parameter_3d = par3d
        if logrange:
            self.pitcon_range_low = scipy.log10(low)
            self.pitcon_range_high = scipy.log10(high)
            self.model.pitcon_par_space = scipy.logspace(
                self.pitcon_range_low, self.pitcon_range_high, self.pitcon_scan_density
            )
        else:
            self.pitcon_range_low = low
            self.pitcon_range_high = high
            self.model.pitcon_par_space = scipy.linspace(
                self.pitcon_range_low, self.pitcon_range_high, self.pitcon_scan_density
            )

        self.model.pitcon_flux_gen = 1
        if runQuiet:
            self.model.SetQuiet()
        else:
            self.model.SetLoud()

        if self.pitcon_scan_parameter_3d != None:
            self.pitcon_res = self.model.PITCON(
                self.pitcon_scan_parameter, self.pitcon_scan_parameter_3d
            )
            self.res_idx = self.pitcon_res[:, :2]
            self.res_metab = self.pitcon_res[:, 2 : len(self.model.species) + 2 :]
            self.res_flux = self.pitcon_res[:, len(self.model.species) + 2 :]
        else:
            self.pitcon_res = self.model.PITCON(self.pitcon_scan_parameter)
            self.res_idx = self.pitcon_res[:, 0]
            self.res_idx = self.res_idx.reshape(self.res_idx.shape[0], 1)
            self.res_metab = self.pitcon_res[:, 1 : len(self.model.species) + 1]
            self.res_flux = self.pitcon_res[:, len(self.model.species) + 1 :]

        print('\n\tContinuation complete\n')

    def analyseData(self, analysis='elas'):
        """
        Performs "analysis" on the PITCON generated set of steady-state results where analysis is:

        - 'elasv' = variable elasticities
        - 'elasp' = parameter elasticities
        - 'elas'  = all elasticities
        - 'mca'   = control coefficients
        - 'resp'  = response coefficients
        - 'eigen' = eigen values
        - 'all'   = all of the above

        Higher level analysis types automatically enable the lower level analysis needed e.g.
        selecting 'mca' implies 'elasv' etc. User output defined with mod.setUserOutput() is
        stored in the mod.res_user array.
        """
        # reset result lists
        self.res_eigen = []
        self.res_user = []

        step_cntr = 0
        self.timer.step_timer('st1', self.res_idx.shape[0])
        for row in range(self.res_idx.shape[0]):
            self.timer.st1_cstep += 1
            step_cntr += 1
            if self.timer.st1_cstep == 1 or step_cntr == 200:
                print('\n', next(self.timer.st1))
                step_cntr = 0

            PitParVal = self.res_idx[row, -1]
            setattr(self.model, self.pitcon_scan_parameter, PitParVal)

            self.model.state_species = self.res_metab[row, :]
            self.model.state_flux = self.res_flux[row, :]
            self.model.SetStateSymb(self.model.state_flux, self.model.state_species)
            self.model.Forcing_Function()

            if analysis in ['elasv', 'elas', 'mca', 'resp', 'eigen', 'all']:
                self.model.EvalEvar(self.model.state_species, self.model.state_flux)
            if analysis in ['elasp', 'elas', 'resp', 'all']:
                self.model.EvalEpar(self.model.state_species, self.model.state_flux)
            if analysis in ['mca', 'resp', 'all']:
                self.model.EvalCC()
            if analysis in ['resp', 'all']:
                self.model.EvalRC()
            if analysis in ['eigen', 'all']:
                self.model.EvalEigen()
                self.res_eigen.append(self.model.eigen_values)
            state_out = []
            if self.user_output != None and len(self.user_output) > 0:
                for attr in self.user_output:
                    try:
                        state_out.append(getattr(self.model, attr))
                    except Exception as ex:
                        print(ex)
                self.res_user.append(state_out)
        self.res_user = scipy.array(self.res_user)
        print('\n\t%s evaluation complete\n' % analysis)

    def getArrayListAsArray(self, array_list):
        """
        Stack (concatenate) the list of arrays into a single array.
        """
        output = None
        if len(array_list) > 1:
            output = scipy.vstack(array_list)
        elif len(array_list) == 1:
            output = array_list[0]
        return output
