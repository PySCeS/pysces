"""
PySCeS - Python Simulator for Cellular Systems (http://pysces.sourceforge.net)

Copyright (C) 2004-2024 B.G. Olivier, J.M. Rohwer, J.-H.S Hofmeyr all rights reserved,

Johann Rohwer (jrohwer@users.sourceforge.net)
Triple-J Group for Molecular Cell Physiology
Stellenbosch University, South Africa.

Permission to use, modify, and distribute this software is given under the
terms of the PySceS (BSD style) license. See LICENSE.txt that came with
this distribution for specifics.

NO WARRANTY IS EXPRESSED OR IMPLIED.  USE AT YOUR OWN RISK.
Johann M. Rohwer
"""
from __future__ import division, print_function
from __future__ import absolute_import
from __future__ import unicode_literals
from pysces.version import __version__
__doc__ = """
          PyscesParScan
          -------------

          PySCeS class distributed multi-dimensional parameter scans with IPython
          """
import numpy as np
from pysces.PyscesUtils import TimerBox
from pysces.PyscesScan import Scanner
import sys, os
import dill as pickle
flush = sys.stdout.flush
from time import sleep, time
import subprocess
try:
    import ipyparallel
except ImportError as ex:
    print()
    print(ex)
    raise ImportError(
        'PARSCANNER: Requires ipyparallel version >=4.0 (http://ipython.org).')
MULTISCANFILE = __file__.replace('PyscesParScan', '_multicorescan')
__psyco_active__ = 0


class ParScanner(Scanner):
    """
    Arbitrary dimension generic distributed scanner.
    Subclassed from *pysces.PyscesScan.Scanner*.
    This class is initiated with a loaded PySCeS model and then allows
    the user to define scan parameters, see ``self.addScanParameter()``
    and user output, see ``self.addUserOutput()``.
    Steady-state results are always stored in ``self.SteadyStateResults`` while
    user output can be found in ``self.UserOutputResults``.
    Distributed (parallel) execution is achieved with the clustering capability
    of IPython. See *ipcluster --help*.

    The optional 'engine' argument specifies the parallel engine to use.

    - *'multiproc'* -- multiprocessing (default)
    - *'ipcluster'* -- IPython cluster
    """
    genOn = True
    _MODE_ = 'state'
    HAS_USER_OUTPUT = False
    nan_on_bad_state = True
    MSG_PRINT_INTERVAL = 500
    __AnalysisModes__ = 'state', 'elasticity', 'mca', 'stability'
    invalid_state_list = None
    invalid_state_list_idx = None
    scans_per_run = 100

    def __init__(self, mod, engine='multiproc'):
        """
		
        Initialize the ParScanner with a specified model and parallel engine.

        Parameters
        ----------
        mod : model object
            The model to be used for scanning.
        engine : {'multiproc', 'ipcluster'}, optional
            The parallel engine to use. Defaults to 'multiproc'. 
            'multiproc' uses multiprocessing, while 'ipcluster' requires a running IPython cluster.

        Raises
        ------
        OSError
            If the engine is set to 'ipcluster' and there is no active IPython cluster running.
        UserWarning
            If the provided engine is not valid.

        Examples
        --------
        >>> from pysces.PyscesParScan import Analyze, setModValue
        >>> mod = SomeModelClass()
        >>> scanner = ParScanner(mod, engine='multiproc')
        parallel engine: multiproc

        >>> mod = AnotherModelClass()
        >>> scanner = ParScanner(mod, engine='ipcluster')
        parallel engine: ipcluster
        Requires a running IPython cluster. See "ipcluster --help".
        
        
PARSCANNER: Requires a running IPython cluster. See "ipcluster --help".
		
		"""
        self.engine = engine
        if engine == 'multiproc':
            print('parallel engine: multiproc')
        elif engine == 'ipcluster':
            print('parallel engine: ipcluster')
            from ipyparallel import Client
            try:
                rc = Client()
                self.rc = rc
            except OSError as ex:
                raise OSError(str(ex) +
                    """
PARSCANNER: Requires a running IPython cluster. See "ipcluster --help".
"""
                    )
            dv = rc[:]
            dv.use_dill()
            lv = rc.load_balanced_view()
            self.dv = dv
            self.lv = lv
            dv.execute('from pysces.PyscesParScan import Analyze, setModValue')
        else:
            raise UserWarning(engine + ' is not a valid parallel engine!')
        from ipyparallel.serialize import codeutil
        self.GenDict = {}
        self.GenOrder = []
        self.ScanSpace = []
        self.mod = mod
        self.SteadyStateResults = []
        self.UserOutputList = []
        self.UserOutputResults = []
        self.scanT = TimerBox()

    def genScanSpace(self):
        """
		
        Generates the parameter scan space, partitioned according to self.scans_per_run.

        This method creates the scan space for the parameters defined in the model, and partitions
        the scan space into manageable chunks based on the `scans_per_run` attribute.

        The method follows these steps:
        1. It calculates the total number of steps (`Tsteps`) required for the scan.
        2. It generates the parameter values for each step and appends them to the `ScanSpace`.
        3. It partitions the `ScanSpace` and sequence array (`SeqArray`), based on `scans_per_run`.

        Attributes
        ----------
        ScanSpace : numpy.ndarray
            Array containing the parameter values for each scan step.
        SeqArray : numpy.ndarray
            Array containing the sequence numbers for each scan step.
        ScanPartition : list of numpy.ndarray
            List containing partitions of the `ScanSpace`, each partition of size `scans_per_run`.
        SeqPartition : list of numpy.ndarray
            List containing partitions of the `SeqArray`, each partition of size `scans_per_run`.
        Tsteps : int
            Total number of steps in the scan.

        Examples
        --------
        >>> mod = SomeModelClass()
        >>> scanner = ParScanner(mod, engine='multiproc')
        parallel engine: multiproc
        >>> scanner.GenOrder = ['param1', 'param2']
        >>> scanner.GenDict = {
        ...     'param1': [None, None, 10, None, False],
        ...     'param2': [None, None, 5, None, False]
        ... }
        >>> scanner.scans_per_run = 20
        >>> scanner.genScanSpace()
        >>> print(scanner.ScanSpace)
        [[...], [...], ...]
        >>> print(scanner.ScanPartition)
        [array([...]), array([...]), ...]

        		
		"""
        spr = self.scans_per_run
        Tsteps = 1
        for gen in self.GenOrder:
            if self.GenDict[gen][4] == False:
                Tsteps *= self.GenDict[gen][2]
        for step in range(Tsteps):
            pars = self.__nextValue__()
            self.ScanSpace.append(pars)
        self.ScanSpace = np.array(self.ScanSpace)
        self.SeqArray = np.arange(1, Tsteps + 1)
        if Tsteps % spr == 0:
            numparts = Tsteps // spr
        else:
            numparts = Tsteps // spr + 1
        self.ScanPartition = [self.ScanSpace[n * spr:(n + 1) * spr] for n in
            range(numparts)]
        self.SeqPartition = [self.SeqArray[n * spr:(n + 1) * spr] for n in
            range(numparts)]
        self.Tsteps = Tsteps

    def Prepare(self, ReRun=False):
        """
		
        Internal method to prepare the parameters and generate ScanSpace.

        This method performs the necessary preparations for parameter scanning, including:
        - Validating the analysis mode.
        - Resetting the scan space and results if a re-run is requested.
        - Generating the parameter scan space.
        - Printing preparation details and timing information.

        Parameters
        ----------
        ReRun : bool, optional
            If True, the scan space and results will be reset before generating the scan space. Defaults to False.

        Raises
        ------
        AssertionError
            If the analysis mode specified in `self._MODE_` is not in the list of valid analysis modes (`self.__AnalysisModes__`).

        Examples
        --------
        >>> mod = SomeModelClass()
        >>> scanner = ParScanner(mod, engine='multiproc')
        parallel engine: multiproc
        >>> scanner._MODE_ = 'default'
        >>> scanner.__AnalysisModes__ = ['default', 'advanced']
        >>> scanner.GenOrder = ['param1', 'param2']
        >>> scanner.GenDict = {
        ...     'param1': [None, None, 10, None, False],
        ...     'param2': [None, None, 5, None, False]
        ... }
        >>> scanner.scans_per_run = 20
        >>> scanner.Prepare(ReRun=True)
        
        PREPARATION
        -----------
        Generated ScanSpace: 0.1234
        PARSCANNER: Tsteps 50
        
        		
		"""
        self._MODE_ = self._MODE_.lower()
        assert self._MODE_ in self.__AnalysisModes__, """
SCANNER: "%s" is not a valid analysis mode!""" % self._MODE_
        if ReRun:
            self.ScanSpace = []
            self.UserOutputResults = []
            self.SteadyStateResults = []
        self.invalid_state_list = []
        self.invalid_state_list_idx = []

    def Run(self, ReRun=False):
        """
		
        Run the parameter scan with a load balancing task client.

        This method executes the parameter scan based on the chosen parallel engine.
        It supports both `multiproc` and `ipcluster` engines for parallel execution.
        The scan results are gathered and stored once the computation is complete.

        Parameters
        ----------
        ReRun : bool, optional
            If True, the scan space and results will be reset before generating the scan space. Defaults to False.

        Examples
        --------
        >>> mod = SomeModelClass()
        >>> scanner = ParScanner(mod, engine='multiproc')
        parallel engine: multiproc
        >>> scanner._MODE_ = 'default'
        >>> scanner.__AnalysisModes__ = ['default', 'advanced']
        >>> scanner.GenOrder = ['param1', 'param2']
        >>> scanner.GenDict = {
        ...     'param1': [None, None, 10, None, False],
        ...     'param2': [None, None, 5, None, False]
        ... }
        >>> scanner.scans_per_run = 20
        >>> scanner.Run(ReRun=True)
        
        PREPARATION
        -----------
        Generated ScanSpace: 0.1234
        PARSCANNER: Tsteps 50
        Preparation completed: 0.5678
        
        PARALLEL COMPUTATION
        --------------------
        Tasks to go...  10
        Tasks to go...  5
        Parallel calculation completed: 1.2345
        
        GATHER RESULT
        -------------
        
        		
		"""
        arl = []
        if self.engine == 'multiproc':
            fN = str(time()).split('.')[0]
            with open(fN, 'wb') as F:
                pickle.dump((self.mod, self.ScanPartition, self.
                    SeqPartition, self.GenOrder, self.UserOutputList), F,
                    protocol=-1)
            fN = os.path.abspath(fN)
            print('Preparation completed:', next(self.scanT.PREP))
            self.scanT.normal_timer('RUN')
            subprocess.call([sys.executable, MULTISCANFILE, self._MODE_, fN])
            with open(fN, 'rb') as F:
                res_list = pickle.load(F)
            os.remove(fN)
            for result in res_list:
                self.StoreData(result)
        elif self.engine == 'ipcluster':
            for i in range(len(self.ScanPartition)):
                arl.append(self.lv.apply(Analyze, self.ScanPartition[i],
                    self.SeqPartition[i], self.GenOrder, self._MODE_, self.
                    UserOutputList, self.mod))
            print('Submitted tasks:', len(arl))
            print('Preparation completed:', next(self.scanT.PREP))
            print('\nPARALLEL COMPUTATION\n--------------------')
            flush()
            self.scanT.normal_timer('RUN')
            while self.lv.queue_status()['unassigned'] > 0:
                sleep(5)
                print('Tasks to go... ', self.lv.queue_status()['unassigned'])
            self.lv.wait()
            flush()
            print('\nGATHER RESULT\n-------------')
            flush()
            for ar in arl:
                result = ar.get()
                self.StoreData(result)

    def RunScatter(self, ReRun=False):
        """
		
        Run the parameter scan by using scatter and gather for the ScanSpace.
        Not load balanced, equal number of scan runs per node.

        This method executes the parameter scan in parallel by scattering the scan space
        across multiple nodes and gathering the results once the computation is finished.
        It only supports the `ipcluster` parallel engine.

        Parameters
        ----------
        ReRun : bool, optional
            If True, the scan space and results will be reset before generating the scan space. Defaults to False.

        Raises
        ------
        RuntimeError
            If the engine is not 'ipcluster'.

        Examples
        --------
        >>> mod = SomeModelClass()
        >>> scanner = ParScanner(mod, engine='ipcluster')
        parallel engine: ipcluster
        >>> scanner._MODE_ = 'default'
        >>> scanner.__AnalysisModes__ = ['default', 'advanced']
        >>> scanner.GenOrder = ['param1', 'param2']
        >>> scanner.GenDict = {
        ...     'param1': [None, None, 10, None, False],
        ...     'param2': [None, None, 5, None, False]
        ... }
        >>> scanner.scans_per_run = 20
        >>> scanner.RunScatter(ReRun=True)
        
        RunScatter() only supported with ipcluster!

        PREPARATION
        -----------
        Generated ScanSpace: 0.1234
        PARSCANNER: Tsteps 50
        Scattered ScanSpace on number of engines: 4
        Preparation completed: 0.5678

        PARALLEL COMPUTATION
        --------------------
        Parallel calculation completed: 1.2345

        GATHER RESULT
        -------------
        
        		
		"""
        if self.engine != 'ipcluster':
            print('RunScatter() only supported with ipcluster!')
            return
        ar = self.dv.gather('y')
        results = ar.get()
        results = [results[i:i + 5] for i in range(0, len(results), 5)]
        for result in results:
            self.StoreData(result)

    def StoreData(self, result):
        """
		
        Internal function which concatenates and stores single result generated by Analyze.

        This method takes the result produced by the `Analyze` function, and stores it appropriately
        within the `ParScanner` instance. The result typically includes steady-state data,
        user-defined output, and invalid states.

        Parameters
        ----------
        result : object
            The result object from an IPython client, generated by the `Analyze` function. It is expected
            to be a list or tuple where:
            - result[0] is the steady-state data.
            - result[1] is additional steady-state data to be concatenated.
            - result[2] is user-defined output data (if present).
            - result[3] is a list of invalid states.
            - result[4] is a list of indices corresponding to invalid states.

        Examples
        --------
        >>> mod = SomeModelClass()
        >>> scanner = ParScanner(mod, engine='multiproc')
        parallel engine: multiproc
        >>> result = [
        ...     [1.0, 2.0, 3.0],  # Steady-state data
        ...     [4.0, 5.0, 6.0],  # Additional steady-state data
        ...     [7.0, 8.0, 9.0],  # User-defined output data
        ...     ['invalid_state_1'],  # Invalid states
        ...     [4]  # Indices of invalid states
        ... ]
        >>> scanner.HAS_USER_OUTPUT = True
        >>> scanner.StoreData(result)
        >>> print(scanner.SteadyStateResults)
        [array([1., 2., 3., 4., 5., 6.])]
        >>> print(scanner.UserOutputResults)
        [array([7., 8., 9.])]
        >>> print(scanner.invalid_state_list)
        ['invalid_state_1']
        >>> print(scanner.invalid_state_list_idx)
        [4]
        
        		
		"""
        if self.HAS_USER_OUTPUT:
            self.UserOutputResults.append(np.array(result[2]))
        self.invalid_state_list += result[3]
        self.invalid_state_list_idx += result[4]

    def GatherScanResult(self):
        """
		
        Concatenates and combines output result fragments from the parallel scan.

        This method aggregates the results of the parameter scan by concatenating the steady-state
        results and user output results (if applicable) into single arrays. It also resets the input 
        parameters and provides summary information about the scan, including any invalid states encountered.

        Examples
        --------
        >>> mod = SomeModelClass()
        >>> scanner = ParScanner(mod, engine='multiproc')
        parallel engine: multiproc
        >>> scanner.SteadyStateResults = [
        ...     np.array([1.0, 2.0, 3.0]),
        ...     np.array([4.0, 5.0, 6.0])
        ... ]
        >>> scanner.UserOutputResults = [
        ...     np.array([7.0, 8.0, 9.0]),
        ...     np.array([10.0, 11.0, 12.0])
        ... ]
        >>> scanner.HAS_USER_OUTPUT = True
        >>> scanner.invalid_state_list = ['invalid_state_1']
        >>> scanner.invalid_state_list_idx = [4]
        >>> scanner.scanT = TimerBox()
        >>> scanner.GatherScanResult()
        
        PARSCANNER: 2 states analysed
        Total time taken:  0.1234
        
        Bad steady states encountered at:
        Sequence:  [4]
        Parameters:  ['invalid_state_1']
        
        >>> print(scanner.SteadyStateResults)
        [[1. 2. 3.]
         [4. 5. 6.]]
        >>> print(scanner.UserOutputResults)
        [[ 7.  8.  9.]
         [10. 11. 12.]]

        		
		"""
        self.SteadyStateResults = np.vstack([i for i in self.
            SteadyStateResults])
        if self.HAS_USER_OUTPUT:
            self.UserOutputResults = np.vstack([i for i in self.
                UserOutputResults])
        if len(self.invalid_state_list) > 0:
            print('\nBad steady states encountered at:\n')
            print('Sequence: ', self.invalid_state_list_idx)
            print('Parameters: ', self.invalid_state_list)


def setModValue(mod, name, value):
    """
		
    Sets a specified attribute of a model to a given value.

    This function checks if the given model object `mod` has an attribute with the 
    name `name`. If the attribute exists, it sets the attribute to the specified 
    `value`. The value is converted to a float before setting the attribute.

    Parameters
    ----------
    mod : object
        The model object whose attribute is to be set.
    name : str
        The name of the attribute to be set.
    value : float
        The value to set the attribute to.

    Raises
    ------
    AssertionError
        If the model object does not have an attribute with the given name.

    Examples
    --------
    >>> class SomeModelClass:
    ...     def __init__(self):
    ...         self.param1 = 0.0
    >>> mod = SomeModelClass()
    >>> setModValue(mod, 'param1', 42)
    >>> print(mod.param1)
    42.0

    	
	"""
    assert hasattr(mod, name), 'Model does not have an attribute: %s ' % name


def Analyze(partition, seqarray, GenOrder, mode, UserOutputList, mod):
    """
		
    Analyzes the model's states and outputs based on a partitioned parameter space.

    This function runs various analysis modes on a model object `mod` using a set of 
    parameters defined by `partition`. It sets the model parameters, executes the 
    specified analysis mode, and collects the state and flux results. It also collects 
    user-defined outputs and keeps track of any invalid states encountered.

    Parameters
    ----------
    partition : list of array-like
        The sets of parameter values to be applied to the model.
    seqarray : array-like
        The sequence numbers corresponding to the parameter sets.
    GenOrder : list of str
        The order of parameters to be set in the model.
    mode : {'state', 'elasticity', 'mca', 'stability'}
        The type of analysis to perform:
        - 'state' for steady-state analysis.
        - 'elasticity' for elasticity analysis.
        - 'mca' for metabolic control analysis.
        - 'stability' for stability analysis.
    UserOutputList : list of str
        The list of user-defined outputs to collect from the model.
    mod : object
        The model object on which analyses are to be performed.

    Returns
    -------
    tuple
        A tuple containing:
        - state_species : list
            The steady-state species concentrations for each parameter set.
        - state_flux : list
            The steady-state fluxes for each parameter set.
        - user_output_results : list
            The user-defined outputs for each parameter set.
        - invalid_state_list : list
            The list of parameter sets that resulted in invalid states.
        - invalid_state_list_idx : list
            The sequence numbers corresponding to the invalid states.

    Raises
    ------
    AssertionError
        If the analysis mode is not recognized.

    Examples
    --------
    >>> class SomeModelClass:
    ...     def __init__(self):
    ...         self.param1 = 0.0
    ...         self.__StateOK__ = True
    ...         self.state_species = []
    ...         self.state_flux = []
    ...     def doState(self):
    ...         pass
    ...     def doElas(self):
    ...         pass
    ...     def doMca(self):
    ...         pass
    ...     def doEigenMca(self):
    ...         pass
    ...     def User_Function(self):
    ...         pass
    >>> mod = SomeModelClass()
    >>> partition = [[1.0, 2.0], [3.0, 4.0]]
    >>> seqarray = [1, 2]
    >>> GenOrder = ['param1', 'param2']
    >>> mode = 'state'
    >>> UserOutputList = ['some_output']
    >>> result = Analyze(partition, seqarray, GenOrder, mode, UserOutputList, mod)
    >>> state_species, state_flux, user_output_results, invalid_state_list, invalid_state_list_idx = result

    	
	"""
    state_species = []
    state_flux = []
    user_output_results = []
    invalid_state_list = []
    invalid_state_list_idx = []
    for i in range(len(partition)):
        pars = partition[i]
        for par in range(len(GenOrder)):
            setModValue(mod, GenOrder[par], pars[par])
        if mode == 'state':
            mod.doState()
        elif mode == 'elasticity':
            mod.doElas()
        elif mode == 'mca':
            mod.doMca()
        elif mode == 'stability':
            mod.doEigenMca()
        mod.User_Function()
        if not mod.__StateOK__:
            invalid_state_list.append(pars)
            invalid_state_list_idx.append(seqarray[i])
        state_species.append(mod.state_species)
        state_flux.append(mod.state_flux)
        user_output_results.append([getattr(mod, res) for res in
            UserOutputList])
    return (state_species, state_flux, user_output_results,
        invalid_state_list, invalid_state_list_idx)
