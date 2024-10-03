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

    Subclass of `pysces.PyscesScan.Scanner`. This class is designed for scanning
    through parameters of a loaded PySCeS model in a distributed manner. It allows 
    users to define scan parameters and customize user output. The results of 
    steady-state computations are stored in `self.SteadyStateResults`, and user-defined 
    output results are stored in `self.UserOutputResults`.

    Parallel execution of scans is facilitated by using IPython's clustering 
    capabilities, providing options for different parallel engines.

    Parameters
    ----------
    mod : object
        A loaded PySCeS model to be used for scanning.
    engine : {'multiproc', 'ipcluster'}, optional
        Specifies the parallel engine to use for distributed execution. 
        Options are:
        - *'multiproc'*: Use Python's multiprocessing (default).
        - *'ipcluster'*: Use IPython cluster.
    
    Attributes
    ----------
    genOn : bool
        Indicates if the generic scanner is active. Default is `True`.
    _MODE_ : str
        Current mode of analysis, with initial value 'state'.
    HAS_USER_OUTPUT : bool
        Indicates whether user output is enabled. Default is `False`.
    nan_on_bad_state : bool
        Determines if NaN should be assigned on encountering a bad state. 
        Default is `True`.
    MSG_PRINT_INTERVAL : int
        Interval for printing progress messages during scan execution.
    __AnalysisModes__ : tuple
        Supported analysis modes: ('state', 'elasticity', 'mca', 'stability').
    invalid_state_list : list
        Stores a list of invalid states encountered during the scan.
    invalid_state_list_idx : list
        Stores the indices of invalid states.
    scans_per_run : int
        Number of scan partitions to execute per run.
    engine : str
        The engine used for parallel execution. Set during initialization.
    rc : ipyparallel.Client, optional
        IPython Client object, initialized if using 'ipcluster'.
    dv : ipyparallel.DirectView or None
        Direct view of IPython engines for executing commands in parallel. 
        Initialized if using 'ipcluster'.
    lv : ipyparallel.LoadBalancedView or None
        A load balanced view for dispatching tasks to engines. Initialized 
        if using 'ipcluster'.
    GenDict : dict
        Dictionary storing scan generator details.
    GenOrder : list
        List maintaining the order of scan generators.
    ScanSpace : list
        List storing the generated parameter scan space.
    mod : object
        Reference to the loaded PySCeS model.
    SteadyStateResults : list
        Accumulates the results of steady-state computations.
    UserOutputList : list
        Contains user-defined output parameters.
    UserOutputResults : list
        Stores results corresponding to user-defined outputs.
    scanT : TimerBox
        Timing utility for measuring preparation and running phases of the scan.

    Notes
    -----
    - Ensure an IPython cluster is running before using 'ipcluster' as the engine. 
      Run `ipcluster --help` for more information on starting an IPython cluster.
    - User must define scan parameters and user output via `self.addScanParameter()` 
      and `self.addUserOutput()` respectively before executing scans.

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
    Initializes the ParScanner with the specified model and parallel engine.

    This method sets up the parallel computation environment for performing
    parameter scans on the specified model. It can use either 'multiproc' for
    multiprocessing or 'ipcluster' for IPython parallel computing. 

    Parameters
    ----------
    mod : object
        The model object to be used for parameter scanning.
    engine : str, optional
        The type of parallel engine to use, by default 'multiproc'.
        Accepted values are 'multiproc' for multiprocessing and 'ipcluster' 
        for using an IPython cluster.

    Attributes
    ----------
    engine : str
        Stores the selected parallel engine.
    rc : ipyparallel.Client, optional
        An IPython parallel Client object if the 'ipcluster' engine is used.
    dv : ipyparallel.DirectView, optional
        A direct view on the engines that can execute tasks in parallel.
    lv : ipyparallel.LoadBalancedView, optional
        A load-balanced view that schedules tasks for parallel execution.
    GenDict : dict
        A dictionary to store generation results.
    GenOrder : list
        A list to store the order of the generations.
    ScanSpace : list
        A list to hold parameters to be scanned.
    mod : object
        The model object used for parameter scanning.
    SteadyStateResults : list
        A list to store the results of steady state computations.
    UserOutputList : list
        A list of user-defined outputs to be captured.
    UserOutputResults : list
        A list to store results of user-defined outputs.
    scanT : TimerBox
        A timing utility for managing scan durations.

    Raises
    ------
    OSError
        If the 'ipcluster' engine is selected but no IPython cluster is running.
    UserWarning
        If an invalid engine type is specified.

    Example
    -------
    # Initialize with default 'multiproc' engine
    >>> scanner = ParScanner(mod)

    # Initialize with 'ipcluster' engine
    >>> scanner = ParScanner(mod, engine='ipcluster')
    
    Notes
    -----
    If using 'ipcluster', ensure that an IPython cluster is running prior to 
    initialization. 'multiproc' is suitable for simpler parallel tasks not 
    requiring distributed computing.

    See also
    --------
    ipyparallel.Client : Documentation for the IPython parallel Client.
    ipyparallel.DirectView : API for directly controlling engines.
    ipyparallel.LoadBalancedView : API for load-balanced task distribution.
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
    Generates the parameter scan space, partitioned according to `self.scans_per_run`.
    
    The method computes a parameter scan space based on the configuration specified 
    in `self.GenDict` and `self.GenOrder`. It partitions the scan space into segments 
    determined by `self.scans_per_run` and stores the partitioned scan space and 
    sequence arrays in `self.ScanPartition` and `self.SeqPartition`, respectively.

    Attributes
    ----------
    ScanSpace : np.ndarray
        An array containing the complete parameter scan space.
    SeqArray : np.ndarray
        An array containing sequence numbers corresponding to each element in the scan space.
    ScanPartition : list of np.ndarray
        A list of arrays, each representing a partition of the scan space.
    SeqPartition : list of np.ndarray
        A list of arrays, each representing a partition of the sequence array.
    Tsteps : int
        Total number of steps in the generated scan space.

    Examples
    --------
    Suppose you have a `ParScanner` instance `scanner` properly initialized with 
    parameter generation settings:

    >>> scanner = ParScanner(mod=your_model_instance)
    >>> scanner.genScanSpace()
    >>> print(scanner.ScanSpace)
    [[1.0, 2.0, 3.0],
     [1.5, 2.5, 3.5],
     ...]
    >>> print(scanner.ScanPartition)
    [array([[1.0, 2.0, 3.0],
            [1.5, 2.5, 3.5]]),
     array([[2.0, 3.0, 4.0],
            [2.5, 3.5, 4.5]])]
    
    This will generate the scan space and partition it based on `scans_per_run` value.
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

    This method configures the initial settings for a parameter scan process
    by confirming the analysis mode and generating the parameter scan space.
    It also resets specific internal state variables if the `ReRun` flag is set.
    Additionally, it validates the current analysis mode against known modes,
    ensuring it is correctly configured before proceeding with the scan generation.

    Parameters
    ----------
    ReRun : bool, optional
        If True, the method will reset the `ScanSpace`, `UserOutputResults`,
        and `SteadyStateResults` before generating a new scan space. Default
        is False.

    Raises
    ------
    AssertionError
        If the current analysis mode (`self._MODE_`) is not one of the valid
        options contained in `self.__AnalysisModes__`.

    Notes
    -----
    This method is typically invoked as part of the setup process prior to
    executing a parameter scan. It uses the `genScanSpace` method to populate
    `self.ScanSpace` with the various parameter configurations.

    See Also
    --------
    ParScanner.genScanSpace : Generates the parameter scan space.

    Examples
    --------
    >>> par_scanner = ParScanner(mod=some_model)
    >>> par_scanner._MODE_ = 'steady_state'
    >>> par_scanner.Prepare(ReRun=True)
    PREPARATION
    -----------
    Generated ScanSpace: ...
    PARSCANNER: Tsteps N

    In the example above, `Prepare` is called with `ReRun=True` to regenerate
    the scan space from scratch, resetting prior results if any exist.
    """
        print('\nPREPARATION\n-----------')
        self.scanT.normal_timer('PREP')
        self._MODE_ = self._MODE_.lower()
        assert self._MODE_ in self.__AnalysisModes__, """
SCANNER: "%s" is not a valid analysis mode!""" % self._MODE_
        if ReRun:
            self.ScanSpace = []
            self.UserOutputResults = []
            self.SteadyStateResults = []
        self.invalid_state_list = []
        self.invalid_state_list_idx = []
        self.genScanSpace()
        print('Generated ScanSpace:', next(self.scanT.PREP))
        print('PARSCANNER: Tsteps', self.Tsteps)
        flush()

    def Run(self, ReRun=False):
        """
    Run the parameter scan with a load balancing task client.

    This method orchestrates the parameter scanning process using a specified
    computational engine. It prepares the scan environment, divides the tasks,
    and manages the execution using either multiprocessing or an IPython
    cluster. After execution, results are stored and processed.

    Parameters
    ----------
    ReRun : bool, optional
        If True, resets the ScanSpace and prepares for a re-run of the scan.
        Default is False.

    Raises
    ------
    AssertionError
        If the specified analysis mode (`self._MODE_`) is invalid.

    Notes
    -----
    The method relies on either the 'multiproc' or 'ipcluster' engine as specified
    in the `self.engine` attribute. It performs the following steps:
    
    1. Prepares the scan environment by invoking `self.Prepare(ReRun)`.
    2. For 'multiproc' engine, handles task management using subprocesses.
    3. For 'ipcluster' engine, uses a load balancing task client to distribute tasks.
    4. Gathers and stores computed results with `self.StoreData(result)`.

    Example
    -------
    >>> scanner = ParScanner()
    >>> scanner.Run(ReRun=True)

    Here, `ParScanner` is assumed to be a pre-defined class that encapsulates
    the parameter scanning configuration and logic.
    """
        self.Prepare(ReRun)
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
        print('Parallel calculation completed:', next(self.scanT.RUN))
        self.GatherScanResult()

    def RunScatter(self, ReRun=False):
        """
    Run the parameter scan by using scatter and gather for the ScanSpace.
    Not load balanced, equal number of scan runs per node.

    This method performs a parallel parameter scan using a scatter and gather
    approach. It distributes (scatters) parts of the `ScanSpace` and 
    `SeqArray` to multiple compute nodes (engines) and collects (gathers) the
    results after computation. The distribution is not load balanced, meaning
    each node processes the same number of scan runs.

    Parameters
    ----------
    ReRun : bool, optional
        If `True`, the scan is prepared for rerunning, which clears the 
        existing results. Default is `False`.

    Notes
    -----
    - This method is only supported when using the 'ipcluster' engine.
    - The parallel computation and result gathering is facilitated via 
      IPython's parallel computing tools.

    See Also
    --------
    ParScanner.Prepare : Prepares parameters and generates ScanSpace.
    ParScanner.StoreData : Stores individual scan results.
    ParScanner.GatherScanResult : Combines and finalizes gathered results.

    Examples
    --------
    Create an instance of `ParScanner` and run the parameter scan using 
    `RunScatter`.

    >>> scanner = ParScanner(engine='ipcluster', GenOrder=...)
    >>> scanner.RunScatter(ReRun=False)

    This will distribute the scan tasks across the engines and gather the results after completion.
    """
        if self.engine != 'ipcluster':
            print('RunScatter() only supported with ipcluster!')
            return
        self.Prepare(ReRun)
        self.dv.push({'GenOrder': self.GenOrder, 'mode': self._MODE_, 'mod':
            self.mod, 'UserOutputList': self.UserOutputList})
        self.dv.scatter('partition', self.ScanSpace)
        self.dv.scatter('seqarray', self.SeqArray)
        print('Scattered ScanSpace on number of engines:', len(self.dv))
        print('Preparation completed:', next(self.scanT.PREP))
        print('\nPARALLEL COMPUTATION\n--------------------')
        flush()
        self.scanT.normal_timer('RUN')
        self.dv.execute(
            'y=Analyze(partition,seqarray,GenOrder,mode,UserOutputList,mod)',
            block=True)
        print('Parallel calculation completed:', next(self.scanT.RUN))
        flush()
        print('\nGATHER RESULT\n-------------')
        flush()
        ar = self.dv.gather('y')
        results = ar.get()
        results = [results[i:i + 5] for i in range(0, len(results), 5)]
        for result in results:
            self.StoreData(result)
        print('Parallel calculation completed:', next(self.scanT.RUN))
        self.GatherScanResult()

    def StoreData(self, result):
        """
    Internal function which concatenates and stores a single result generated by Analyze.

    This method processes and stores the result from a parameter scan. It concatenates
    steady-state results and maintains lists of user-defined outputs and any invalid states.

    Parameters
    ----------
    result : list
        An IPython client result object. It is expected to be a list where:
        - `result[0]` and `result[1]` contain data to be concatenated into steady-state results.
        - `result[2]` contains user-defined outputs, if applicable.
        - `result[3]` and `result[4]` contain invalid state information.

    Attributes Modified
    -------------------
    SteadyStateResults : list
        Appends the concatenated steady-state results.
    UserOutputResults : list
        Appends user-defined outputs, if `HAS_USER_OUTPUT` is True.
    invalid_state_list : list
        Extends with invalid state data from `result[3]`.
    invalid_state_list_idx : list
        Extends with indices of invalid states from `result[4]`.

    Examples
    --------
    >>> par_scanner = ParScanner()
    >>> some_result = [
    ...     [1.0, 2.0, 3.0],  # Result part 0
    ...     [4.0, 5.0, 6.0],  # Result part 1
    ...     [7.0, 8.0, 9.0],  # Result part 2 (optional user output)
    ...     ['invalid_state_1'],  # Result part 3 (invalid state list)
    ...     [0]               # Result part 4 (invalid state indices)
    ... ]
    >>> par_scanner.StoreData(some_result)
    >>> par_scanner.SteadyStateResults
    [array([1., 2., 3., 4., 5., 6.])]
    >>> par_scanner.UserOutputResults
    [array([7., 8., 9.])]
    >>> par_scanner.invalid_state_list
    ['invalid_state_1']
    >>> par_scanner.invalid_state_list_idx
    [0]
    """
        self.SteadyStateResults.append(np.hstack((np.array(result[0]), np.
            array(result[1]))))
        if self.HAS_USER_OUTPUT:
            self.UserOutputResults.append(np.array(result[2]))
        self.invalid_state_list += result[3]
        self.invalid_state_list_idx += result[4]

    def GatherScanResult(self):
        """
    Concatenates and combines output result fragments from the parallel scan.

    This method aggregates the results obtained from parallel computations
    into unified structures, handling both the steady-state results and any
    user-defined output results if applicable. It also prints a summary of 
    the number of states analyzed and any issues encountered with invalid 
    steady states.

    After collecting and processing the results, it resets input parameters 
    to prepare for any subsequent scans.

    Examples
    --------
    >>> scanner = ParScanner()
    >>> scanner.RunScatter()
    >>> scanner.GatherScanResult()
    PARSCANNER: 150 states analysed
    Total time taken: ...
    Bad steady states encountered at:
    Sequence:  [2, 48]
    Parameters:  [[1.0, 0.5, 0.2], [0.9, 0.4, 0.1]]

    Notes
    -----
    - This function assumes that `SteadyStateResults`, `UserOutputResults`, 
      and the lists for invalid states (`invalid_state_list`, 
      `invalid_state_list_idx`) have been populated during the parallel scan.
    - The timing details are managed using internal variable `scanT` which 
      tracks preparation and execution durations.

    See Also
    --------
    ParScanner.RunScatter : Initiates the parallel scan by scattering data
    ParScanner.StoreData : Handles storage of individual result fragments
    """
        self.SteadyStateResults = np.vstack([i for i in self.
            SteadyStateResults])
        if self.HAS_USER_OUTPUT:
            self.UserOutputResults = np.vstack([i for i in self.
                UserOutputResults])
        self.resetInputParameters()
        print('\nPARSCANNER: %s states analysed' % len(self.SteadyStateResults)
            )
        print('Total time taken: ', next(self.scanT.PREP))
        self.scanT.PREP.close()
        self.scanT.RUN.close()
        if len(self.invalid_state_list) > 0:
            print('\nBad steady states encountered at:\n')
            print('Sequence: ', self.invalid_state_list_idx)
            print('Parameters: ', self.invalid_state_list)


def setModValue(mod, name, value):
    """
    Sets the value of a specified attribute in a given model object to a float.

    This function checks if the given model object `mod` has an attribute named `name`.
    If the attribute exists, it sets its value to the provided `value`, converted to a float.
    If the attribute does not exist, it raises an `AssertionError`.

    Parameters
    ----------
    mod : object
        The model object whose attribute needs to be set.
    name : str
        The name of the attribute to be set within the model object.
    value : int or float or str
        The value to set for the specified attribute. It will be converted to a float.

    Raises
    ------
    AssertionError
        If `mod` does not have an attribute named `name`.

    Examples
    --------
    >>> class Model:
    ...     def __init__(self):
    ...         self.parameter = 10.0
    ...
    >>> model = Model()
    >>> setModValue(model, 'parameter', 20)
    >>> print(model.parameter)
    20.0

    >>> setModValue(model, 'parameter', '30')
    >>> print(model.parameter)
    30.0

    >>> setModValue(model, 'nonexistent', 40)
    AssertionError: Model does not have an attribute: nonexistent
    """
    assert hasattr(mod, name), 'Model does not have an attribute: %s ' % name
    setattr(mod, name, float(value))


def Analyze(partition, seqarray, GenOrder, mode, UserOutputList, mod):
    """
    Analyzes system states or parameters based on a given partition and mode, 
    updating the model accordingly and collecting results.

    This function iteratively configures a model using a partition, performs
    specific analyses (state, elasticity, MCA, stability), and gathers 
    results for states, fluxes, and user-defined outputs. It also tracks any 
    invalid states encountered during the analysis.

    Parameters
    ----------
    partition : list of list of float
        A list of parameter sets, each defining a configuration for the model.
    seqarray : list of int
        An array indicating the sequence index for each partition set.
    GenOrder : list of str
        A list of model attribute names that defines the order in which 
        the parameters are set.
    mode : {'state', 'elasticity', 'mca', 'stability'}
        Specifies the type of analysis to perform on the model.
        - 'state': Performs state analysis.
        - 'elasticity': Calculates elasticity.
        - 'mca': Conducts metabolic control analysis.
        - 'stability': Executes eigenvalue stability analysis.
    UserOutputList : list of str
        A list of attribute names for which user-defined outputs are extracted 
        from the model after analysis.
    mod : object
        The model object to be analyzed. It must have certain attributes and 
        methods including those specified in `GenOrder` and 
        `UserOutputList`, as well as `doState()`, `doElas()`, `doMca()`, 
        `doEigenMca()`, and a `User_Function()` method.

    Returns
    -------
    tuple
        A tuple containing:
        - state_species (list of array-like): Collected states for each analysis.
        - state_flux (list of array-like): Collected fluxes for each analysis.
        - user_output_results (list of list): Results from user-defined outputs.
        - invalid_state_list (list of list): Parameters corresponding to 
          unsuccessful analyses (invalid states).
        - invalid_state_list_idx (list of int): Sequence indices corresponding 
          to invalid states.

    Raises
    ------
    AssertionError
        If the model does not have an attribute specified in `GenOrder`.
    
    Examples
    --------
    >>> class MockModel:
    ...     def __init__(self):
    ...         self.param1 = 0.0
    ...         self.param2 = 0.0
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
    ...
    >>> partition = [[1.0, 2.0], [3.0, 4.0]]
    >>> seqarray = [0, 1]
    >>> GenOrder = ['param1', 'param2']
    >>> mode = 'state'
    >>> UserOutputList = []
    >>> mod = MockModel()
    >>> results = Analyze(partition, seqarray, GenOrder, mode, UserOutputList, mod)
    >>> print(results[0])  # state_species
    [[], []]
    >>> print(results[1])  # state_flux
    [[], []]
    >>> print(results[3])  # invalid_state_list
    []

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
