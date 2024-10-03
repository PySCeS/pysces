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
from pysces.version import __version__
import os
import copy
import time
import contextlib
import re
import pickle
import numpy
import scipy
import scipy.linalg
import scipy.integrate
from getpass import getuser
from . import model_dir as MODEL_DIR
from . import output_dir as OUTPUT_DIR
from . import install_dir as INSTALL_DIR
from . import PyscesRandom as random
from .PyscesScan import Scanner
from .core2.InfixParser import MyInfixParser
from .core2.PyscesCore2 import NewCoreBase, NumberBase
from . import nleq2_switch, pitcon_switch, plt, gplt, interface, PyscesStoich, PyscesParse, __SILENT_START__, __CHGDIR_ON_START__, __CUSTOM_DATATYPE__, SED, _checkPandas
if nleq2_switch:
    from . import nleq2
if pitcon_switch:
    from . import pitcon
"""
TODO: Parameter elasticities wrt the compartments
"""
__doc__ = """
            PyscesModel
            -----------

            This module contains the core PySCeS classes which
            create the model and associated data objects

            """
try:
    input = raw_input
except NameError:
    pass
ScipyDerivative = None
HAVE_SCIPY_DERIV = False
try:
    ScipyDerivative = scipy.derivative
    HAVE_SCIPY_DERIV = True
except AttributeError:
    try:
        import scipy.misc
        ScipyDerivative = scipy.misc.derivative
        HAVE_SCIPY_DERIV = True
    except ImportError as AttributeError:
        pass
if not HAVE_SCIPY_DERIV:
    raise RuntimeError('\nSciPy derivative function not available')
if __CHGDIR_ON_START__:
    CWD = OUTPUT_DIR
else:
    CWD = os.getcwd()
del (random.setstate, random.getstate)
del random.randrange, random.Random, random.choice
del random.shuffle, random.jumpahead
del random.SystemRandom, random.WichmannHill, random.triangular
if int(scipy.__version__.split('.')[0]) < 1:
    print('\nINFO: Your version of SciPy (' + scipy.version.version +
        """) might be too old
	Version 1.0.x or newer is strongly recommended
"""
        )
elif not __SILENT_START__:
    print('You are using NumPy ({}) with SciPy ({})'.format(numpy.
        __version__, scipy.__version__))
_HAVE_ASSIMULO = False
_ASSIMULO_LOAD_ERROR = ''
try:
    with contextlib.redirect_stderr(open(os.devnull, 'w')):
        from assimulo.solvers import CVode
        from assimulo.problem import Explicit_Problem
    _HAVE_ASSIMULO = True
    if not __SILENT_START__:
        print('Assimulo CVode available')
except Exception as ex:
    _ASSIMULO_LOAD_ERROR = '{}'.format(ex)
    _HAVE_ASSIMULO = False
if _HAVE_ASSIMULO:


    class EventsProblem(Explicit_Problem):
        """
    A class that models problems with state events using an explicit problem solver.

    This class extends the `Explicit_Problem` class and manages events in
    a simulation model compliant with the behavior defined in the `mod`
    object. Events determine specific conditions during the simulation
    that trigger changes in the state of the model.

    Parameters
    ----------
    mod : object
        An object that contains the model description, including events,
        parameters, and rules. It is expected to have specific attributes
        such as `__KeyWords__` for metadata and `__events__` for event
        definitions.

    **kwargs : dict
        Additional keyword arguments used for the initialization of the
        superclass `Explicit_Problem`.

    Attributes
    ----------
    mod : object
        The model object passed during initialization that contains the
        simulation details.

    name : str
        The name of the model, extracted from `mod.__KeyWords__['Modelname']`.

    parvals : list
        A list that stores the parameter values during the simulation at
        various states.

    event_times : list
        A list that records the times at which events occur during the
        simulation.

    Methods
    -------
    state_events(t, y, sw)
        Defines the state events and checks for conditions or changes in
        state at each timestep during the simulation.

    handle_event(solver, event_info)
        Manages the actions to be taken when a state event occurs, such as
        updating parameter values or executing assignments.

    setSequence(event_list)
        Arranges a given list of events based on their priority, returned as
        a processed sequence of events to be handled.
    """

        def __init__(self, mod, **kwargs):
            """
    Initialize an EventsProblem instance.

    This constructor sets up an EventsProblem by initializing its parent class
    `Explicit_Problem` and storing model-specific data such as parameters and event times.

    Parameters
    ----------
    mod : object
        A model object containing the necessary attributes for setting up the problem.
        The model must have a `__KeyWords__` dictionary with a 'Modelname' key and
        a list of `parameters`.
        
    **kwargs
        Additional keyword arguments passed to the `Explicit_Problem` initialization.

    Attributes
    ----------
    mod : object
        Stores the model object provided with the problem setup.
    name : str
        The name of the model, extracted from the `mod.__KeyWords__` dictionary.
    parvals : list of list
        A list containing a list of parameter values obtained from the `mod.parameters`.
    event_times : list
        An initially empty list meant to store event times pertinent to the problem.

    Examples
    --------
    Here's a basic example demonstrating how to initialize an `EventsProblem` instance:

    >>> class MockModel:
    ...     __KeyWords__ = {'Modelname': 'SampleModel'}
    ...     parameters = ['param1', 'param2']
    ...     param1 = 10
    ...     param2 = 20
    ...
    >>> mock_model = MockModel()
    >>> problem = EventsProblem(mock_model)
    >>> problem.name
    'SampleModel'
    >>> problem.parvals
    [[10, 20]]
    """
            Explicit_Problem.__init__(self, **kwargs)
            self.mod = mod
            self.name = self.mod.__KeyWords__['Modelname']
            self.parvals = []
            self.parvals.append([getattr(self.mod, p) for p in self.mod.
                parameters])
            self.event_times = []

        def state_events(self, t, y, sw):
            """
    Evaluate and process state events for the given time and state.

    This method checks each of the state events defined in the model at the 
    specified time `t`. It returns a numpy array where each element corresponds 
    to an event, indicating whether the event condition has been met (0 if 
    true, 1 if false). If the model includes rate rules, those are executed 
    as part of the processing.

    Parameters
    ----------
    t : float
        The current time point at which the state events are being evaluated.
    y : array-like
        The current state of the system, typically representing the state variables.
    sw : array-like
        A status or switch array, usually provided by the system or solver 
        to manage the state event handling.

    Returns
    -------
    numpy.ndarray
        An array of zeros and ones. Each entry corresponds to a specific 
        event; 0 if the event's condition is met (i.e., the event has occurred), 
        or 1 if the condition is not met.

    Notes
    -----
    - The method assumes that the model has the attribute `__events__` 
      which is a list of callable event functions.
    - The attribute `__HAS_RATE_RULES__` is used to determine if there 
      are rate rules that need to be executed during the event processing.

    Example
    -------
    Assume you have an instance of `EventsProblem` initialized with a model `mod`
    which defines state events.

    >>> problem = EventsProblem(mod)
    >>> current_time = 5.0
    >>> current_state = [1.0, 0.5, 0.2]
    >>> switch_array = [0, 1, 0]
    >>> event_statuses = problem.state_events(current_time, current_state, switch_array)
    >>> print(event_statuses)
    array([0., 1., 1.])

    Here, the method checks the model's events at `time = 5.0` and evaluates the 
    current system state, returning an array indicating the status of each event.
    """
            self.events = self.mod.__events__
            eout = numpy.zeros(len(self.events))
            self.mod._TIME_ = t
            for ev in range(len(self.events)):
                if self.events[ev](t):
                    eout[ev] = 0
                else:
                    eout[ev] = 1
            if self.mod.__HAS_RATE_RULES__:
                self.mod.ExecRateRules()
            return eout

        def handle_event(self, solver, event_info):
            """
    Handles events triggered during the simulation by executing the defined
    event assignments and updating the state of the model accordingly.

    Parameters
    ----------
    solver : object
        The solver instance holding the current state of the model and time point.
    event_info : tuple
        Information about the event, including the state information array
        which indicates any state changes (-1 for events triggered).

    Notes
    -----
    - The method updates the model's state when events occur, according to
      the assignments defined in each event.
    - Handles special cases where the model uses rate rules or assignments 
      to dependent species.
    - Appends the parameter values to `self.parvals` after executing events.

    Examples
    --------
    Suppose `solver` is an instance of a solver class with current simulation
    time `t` and `event_info` is a tuple representing the event state information.

    >>> events_problem = EventsProblem(mod)
    >>> solver = MySolver()
    >>> solver.t = 10.0
    >>> event_info = ([-1, 0, -1],)  # Assume two events are triggered

    Handle the events at the current time step:

    >>> events_problem.handle_event(solver, event_info)
    executing Event1
    executing Event3

    The method will print the names of the executed events and update the
    state of the species involved in the events.
    """
            self.event_times.append(solver.t)
            state_info = event_info[0]
            state_info = numpy.array(state_info)
            idx = numpy.where(state_info == -1)
            event_list = [self.events[j] for j in idx[0]]
            sequence = self.setSequence(event_list)
            for ev in sequence:
                if ev._assign_now:
                    print('executing', ev.name)
                    for ass in ev.assignments:
                        ass.evaluateAssignment()
                        if (ass.variable in self.mod.L0matrix.getLabels()[1
                            ] or self.mod.mode_integrate_all_odes and ass.
                            variable in self.mod.__species__):
                            assVal = ass.getValue()
                            assIdx = self.mod.__species__.index(ass.variable)
                            if self.mod.__KeyWords__['Species_In_Conc']:
                                solver.y[assIdx] = assVal * getattr(self.
                                    mod, self.mod.__CsizeAllIdx__[assIdx])
                                setattr(self.mod, ass.variable, assVal *
                                    getattr(self.mod, self.mod.
                                    __CsizeAllIdx__[assIdx]))
                            else:
                                solver.y[assIdx] = assVal
                                setattr(self.mod, ass.variable, assVal)
                        elif not self.mod.mode_integrate_all_odes and ass.variable in self.mod.L0matrix.getLabels(
                            )[0]:
                            print(
                                """Event assignment to dependent species!
Consider setting "mod.mode_integrate_all_odes = True\""""
                                )
                            setattr(self.mod, ass.variable, ass.value)
                        elif self.mod.__HAS_RATE_RULES__ and ass.variable in self.mod.__rate_rules__:
                            assert self.mod.mode_integrate_all_odes, """Assigning events to RateRules requires integrating all ODEs.
 Set "mod.mode_integrate_all_odes = True\""""
                            assVal = ass.getValue()
                            rrIdx = self.mod.__rate_rules__.index(ass.variable)
                            self.mod.__rrule__[rrIdx] = assVal
                            solver.y[self.mod.Nmatrix.shape[0] + rrIdx
                                ] = assVal
                            setattr(self.mod, ass.variable, assVal)
                        else:
                            ass()
                            setattr(self.mod, ass.variable, ass.value)
                self.parvals.append([getattr(self.mod, p) for p in self.mod
                    .parameters])

        def setSequence(self, event_list):
            """
    Determine the sequence of events based on their priority.

    This method takes a list of events and orders them first by priority, if
    specified. Events without specified priorities are inserted randomly
    into the sequence.

    Parameters
    ----------
    event_list : list of Event
        A list of `Event` instances, each representing an event that needs
        to be sequenced based on its priority attribute.

    Returns
    -------
    sequence : list of Event
        The list of events ordered first by priority and then randomly 
        for those without priorities.

    Examples
    --------
    >>> ev1 = Event(name="event1", priority=2)
    >>> ev2 = Event(name="event2", priority=1)
    >>> ev3 = Event(name="event3", priority=None)
    >>> ev4 = Event(name="event4", priority=3)

    >>> events_problem = EventsProblem()
    >>> sequence = events_problem.setSequence([ev1, ev2, ev3, ev4])

    >>> for ev in sequence:
    >>>     print(ev.name)

    This might print, for instance:
    event4
    event1
    event2
    event3
    where `event4`, `event1`, and `event2` are sorted by their priorities, 
    and `event3` is inserted at a random position.
    """
            no_prio = [ev for ev in event_list if ev.priority is None]
            prio = [ev for ev in event_list if ev.priority is not None]
            random_prio_sample = random.sample(prio, len(prio))
            prio_sorted = sorted(random_prio_sample, reverse=True, key=lambda
                x: x.priority)
            if not no_prio:
                sequence = prio_sorted
            elif not prio_sorted:
                sequence = random.sample(no_prio, len(no_prio))
            else:
                sequence = prio_sorted
                for i in no_prio:
                    sequence.insert(random.randint(0, len(prio_sorted)), i)
            return sequence
_HAVE_VPYTHON = False
_VPYTHON_LOAD_ERROR = ''
__psyco_active__ = 0
"""
try:
    import visual
    _HAVE_VPYTHON = True
    vpyscene = visual.display(x=150, y=50, title='PySCeS '+__version__+' - C Brett G. Olivier, Stellenbosch 2022', width=640, height=480,    center=(0,0,0), background=(0,0,0))
    vpyscene.select()
    # vpyscene.autocenter = True
    l1 = visual.label(pos=(0,0,0), text="Welcome to the PySCeS Visualisation Terminal", color=visual.color.red, box=0)
    l2 = visual.label(pos=(0,0.5,0.5), text="It just works!", color=visual.color.green, box=0)
except Exception, ex:
    _VPYTHON_LOAD_ERROR = '%s' % ex
    _HAVE_VPYTHON = False
try:
    import psyco
    psyco.profile()
    __psyco_active__ = 1
    print('PySCeS Model module is now PsycoActive!')
except:
    __psyco_active__ = 0
"""
mach_eps = numpy.finfo(float).eps
pscParser = PyscesParse.PySCeSParser(debug=0)


def chkpsc(File):
    """
    Checks whether the filename has a '.psc' extension and appends one if it does not.

    Ensures that the provided filename string ends with the '.psc' extension.
    If the filename does not already have this extension, '.psc' is appended to it.

    Parameters
    ----------
    File : str
        The filename string to be checked and potentially modified.

    Returns
    -------
    str
        The corrected filename with a '.psc' extension.

    Examples
    --------
    >>> chkpsc("model")
    'model.psc'

    >>> chkpsc("data.psc")
    'data.psc'

    Notes
    -----
    If an unexpected error occurs during execution, 'Chkpsc error' is printed.
    
    """
    try:
        if File[-4:] == '.psc':
            pass
        else:
            print('Assuming extension is .psc')
            File += '.psc'
    except:
        print('Chkpsc error')
    return File


def chkmdir():
    """
    Import and retrieve `pysces.model_dir`.

    This function is used to import the `pysces` module and access the 
    `model_dir` attribute, which indicates the directory path where `pysces` 
    models are stored. Note that this function does not take any arguments and 
    does not return any value. It is meant to demonstrate the process of 
    importing a module and using its attributes.

    Examples
    --------
    To use `chkmdir`, simply call the function as shown below. Assuming that 
    the `pysces` module is installed and the `model_dir` attribute is accessible,
    this example demonstrates how you might use the function to retrieve the 
    model directory.

    >>> import pysces
    >>> model_directory = chkmdir()
    >>> print(model_directory)

    Note: Ensure that the `pysces` module is installed in your Python environment 
    before attempting to import it.
    """
    pass


class BagOfStuff(object):
    """
    A collection of attributes defined by row and column lists, used by response coefficients.

    This class manages a matrix of values (`matrix`) alongside lists of row and
    column names (`row` and `col`). These lists are used to access entries in the matrix
    and create attributes dynamically, which can be queried or listed.

    Attributes
    ----------
    matrix : numpy.ndarray
        A two-dimensional array of values.
    row : list of str
        A list of row name strings.
    col : list of str
        A list of column name strings.

    Methods
    -------
    __init__(matrix, row, col):
        Initializes the BagOfStuff instance with a matrix, row names, and column names.
    load():
        Attaches matrix values as attributes named by `row` and `col` combinations.
    get(attr1, attr2):
        Returns the attribute value for the specified row and column combination.
    select(attr, search='a'):
        Selects and returns a dictionary of attributes associated with a given attribute.
    list():
        Lists all attributes with their corresponding values as a dictionary.
    listAllOrdered(order='descending', absolute=True):
        Returns an ordered list of (attribute, value) tuples.
    """
    matrix = None
    row = None
    col = None

    def __init__(self, matrix, row, col):
        """
        Initialize the BagOfStuff class by loading attributes from arrays.

        This method sets the `matrix`, `row`, and `col` attributes using the provided input
        parameters. It also asserts that the dimensions of `row` and `col` match the dimensions
        of `matrix`.

        Parameters
        ----------
        matrix : np.ndarray
            A 2D NumPy array representing the matrix.
        
        row : array-like
            An iterable containing elements corresponding to the rows of the matrix.
            The length of `row` should match the number of rows in `matrix`.
        
        col : array-like
            An iterable containing elements corresponding to the columns of the matrix.
            The length of `col` should match the number of columns in `matrix`.

        Raises
        ------
        AssertionError
            If the length of `row` does not match the number of rows in `matrix`.
            If the length of `col` does not match the number of columns in `matrix`.

        Examples
        --------
        >>> import numpy as np
        >>> matrix = np.array([[1, 2, 3], [4, 5, 6]])
        >>> row = ['row1', 'row2']
        >>> col = ['col1', 'col2', 'col3']
        >>> bag = BagOfStuff(matrix, row, col)
        >>> bag.matrix
        array([[1, 2, 3],
               [4, 5, 6]])
        >>> bag.row
        ['row1', 'row2']
        >>> bag.col
        ['col1', 'col2', 'col3']
        """
        assert len(row) == matrix.shape[0], '\nRow dimension mismatch'
        assert len(col) == matrix.shape[1], '\nCol dimension mismatch'
        self.matrix = matrix
        self.row = list(row)
        self.col = list(col)

    def load(self):
        """
    Populate object attributes with matrix values.

    This method dynamically creates attributes on the object for each 
    combination of rows and columns of the matrix. The attribute names 
    are formed by concatenating the string representations of each 
    row and column label, joined by an underscore. The values of these 
    attributes are set to the corresponding values from the matrix.

    Raises
    ------
    AttributeError
        If there is any issue setting an attribute on the object.
    
    Examples
    --------
    >>> import numpy as np
    >>> matrix = np.array([[1, 2], [3, 4]])
    >>> rows = ['row1', 'row2']
    >>> cols = ['col1', 'col2']
    >>> bag = BagOfStuff(matrix, rows, cols)
    >>> bag.load()
    >>> bag.row1_col1
    1
    >>> bag.row1_col2
    2
    >>> bag.row2_col1
    3
    >>> bag.row2_col2
    4
    
    """
        for r in range(self.matrix.shape[0]):
            for c in range(self.matrix.shape[1]):
                setattr(self, str(self.row[r]) + '_' + str(self.col[c]),
                    self.matrix[r, c])

    def get(self, attr1, attr2):
        """
    Retrieve the value of an attribute specified by `attr1_attr2` pattern.

    Given two attribute identifiers, `attr1` and `attr2`, this method 
    returns the value of the combined attribute `attr1_attr2` if it exists 
    in the instance. If the attribute does not exist, None is returned, 
    and an exception message is printed.

    Parameters
    ----------
    attr1 : str
        The first part of the attribute name.
    attr2 : str
        The second part of the attribute name.

    Returns
    -------
    any or None
        The value of the attribute `attr1_attr2` if it exists, otherwise 
        None.

    Examples
    --------
    >>> instance = BagOfStuff(matrix=[[1, 2], [3, 4]], row=['row1', 'row2'], col=['col1', 'col2'])
    >>> instance.load()
    >>> instance.get('row1', 'col1')
    1
    >>> instance.get('row2', 'col3')  # Non-existing attribute
    'row2_col3'
    None
    """
        try:
            return getattr(self, attr1 + '_' + attr2)
        except Exception as ex:
            print(ex)
            return None

    def select(self, attr, search='a'):
        """
    Selects and returns a dictionary of attributes in the format "<attr>_<name>" 
    or "<name>_<attr>" with their corresponding values based on the specified 
    search criteria.

    If `attr` is found as an index in both `self.row` and `self.col`, the 
    search parameter determines which attributes are selected:
        - search='a': Selects from both left (row-based) and right (column-based) attributes.
        - search='l': Selects from left (row-based) attributes only.
        - search='r': Selects from right (column-based) attributes only.

    Parameters
    ----------
    attr : str
        The attribute name to search for in `self.row` and/or `self.col`.
    search : {'a', 'l', 'r'}, optional
        Specifies the attribute search strategy. Default is 'a'.

    Returns
    -------
    dict
        A dictionary with keys in the format "<attr>_<name>" or "<name>_<attr>" 
        and their corresponding matrix values. Returns an empty dictionary if 
        no matching attributes are found.

    Examples
    --------
    >>> bag = BagOfStuff()
    >>> bag.row = ['row1', 'row2']
    >>> bag.col = ['col1', 'col2']
    >>> bag.matrix = np.array([[1, 2], [3, 4]])
    >>> bag.select(attr='row1', search='a')
    {'row1_col1': 1, 'row1_col2': 2}

    >>> bag.select(attr='col1', search='r')
    {'row1_col1': 1, 'row2_col1': 3}

    >>> bag.select(attr='row2', search='l')
    {'row2_col1': 3, 'row2_col2': 4}
    """
        output = {}
        if attr in self.row and attr in self.col:
            if search in ['a', 'l']:
                row = self.row.index(attr)
                for col in range(self.matrix.shape[1]):
                    output.setdefault(attr + '_' + self.col[col], self.
                        matrix[row, col])
            if search in ['a', 'r']:
                col = self.col.index(attr)
                for row in range(self.matrix.shape[0]):
                    output.setdefault(self.row[row] + '_' + attr, self.
                        matrix[row, col])
        elif attr in self.row:
            row = self.row.index(attr)
            for col in range(self.matrix.shape[1]):
                output.setdefault(attr + '_' + self.col[col], self.matrix[
                    row, col])
        elif attr in self.col:
            col = self.col.index(attr)
            for row in range(self.matrix.shape[0]):
                output.setdefault(self.row[row] + '_' + attr, self.matrix[
                    row, col])
        return output

    def list(self):
        """
    Return all attributes as an attr:val dictionary.

    This method iterates over the rows and columns of a matrix attribute
    of the instance and constructs a dictionary that maps attribute names,
    formed by combining row and column identifiers, to their corresponding
    values in the matrix.

    Returns
    -------
    dict
        A dictionary where keys are strings in the format `<row>_<col>`
        and values are the corresponding elements from the matrix.

    Examples
    --------
    >>> bag = BagOfStuff()
    >>> bag.row = ['row1', 'row2']
    >>> bag.col = ['col1', 'col2']
    >>> bag.matrix = np.array([[1, 2], [3, 4]])
    >>> bag.list()
    {'row1_col1': 1, 'row1_col2': 2, 'row2_col1': 3, 'row2_col2': 4}
    """
        output = {}
        for row in range(self.matrix.shape[0]):
            for col in range(self.matrix.shape[1]):
                output.setdefault(self.row[row] + '_' + self.col[col], self
                    .matrix[row, col])
        return output

    def listAllOrdered(self, order='descending', absolute=True):
        """
    Return an ordered list of (attr, value) tuples.

    Parameters
    ----------
    order : str, optional
        The order to return the tuples in. Can be `'descending'` or `'ascending'`.
        The default is `'descending'`.
    absolute : bool, optional
        Determines whether to order by the absolute value of the elements. The default is `True`.

    Returns
    -------
    list of tuple
        A list of tuples where each tuple contains an attribute name and its corresponding value,
        ordered as specified by the `order` and `absolute` parameters.

    Examples
    --------
    Assume `bag` is an instance of `BagOfStuff` with `matrix`, `row`, and `col` appropriately initialized.

    >>> bag = BagOfStuff()
    >>> bag.matrix = numpy.array([[3, -4], [2, -1]])
    >>> bag.row = ['attr1', 'attr2']
    >>> bag.col = ['val1', 'val2']
    >>> bag.listAllOrdered()
    [('attr1_val2', -4), ('attr1_val1', 3), ('attr2_val1', 2), ('attr2_val2', -1)]

    >>> bag.listAllOrdered(order='ascending')
    [('attr2_val2', -1), ('attr2_val1', 2), ('attr1_val1', 3), ('attr1_val2', -4)]

    >>> bag.listAllOrdered(absolute=False)
    [('attr1_val2', -4), ('attr1_val1', 3), ('attr2_val1', 2), ('attr2_val2', -1)]

    >>> bag.listAllOrdered(order='ascending', absolute=False)
    [('attr2_val2', -1), ('attr2_val1', 2), ('attr1_val1', 3), ('attr1_val2', -4)]
    """
        if order != 'ascending':
            order = True
        else:
            order = False
        output_v = []
        output_n = []
        for row in range(self.matrix.shape[0]):
            for col in range(self.matrix.shape[1]):
                output_v.append(self.matrix[row, col])
                output_n.append(self.row[row] + '_' + self.col[col])
        if absolute:
            new_idx = numpy.argsort([abs(v) for v in output_v])
        else:
            new_idx = numpy.argsort(output_v)
        output = []
        if order:
            new_idx = new_idx.tolist()
            new_idx.reverse()
        for i_ in new_idx:
            output.append((output_n[i_], output_v[i_]))
        return output


class ScanDataObj(object):
    """
    A class designed to store parameter scan data, utilizing the StateDataObj class.

    This class is intended for the collection and management of data resulting from
    parameter scans, including various forms of biological, chemical, or physical
    state data such as species concentrations, fluxes, rules, and additional
    experimental data (xdata). It facilitates the structured storage of these data
    sets and provides methods to access them.

    Attributes
    ----------
    species : list or None
        A list to store species data; initialized as an empty list, or None if not used.
    fluxes : list or None
        A list to store fluxes data; initialized as an empty list, or None if not used.
    rules : list or None
        A list to store rules data; initialized as an empty list, or None if not used.
    xdata : list or None
        A list to store additional xdata; initialized as an empty list, or None if not used.
    mod_data : list or None
        A list to store modified data; initialized as an empty list, or None if not used.
    flux_labels : list or None
        Labels for fluxes, if applicable.
    species_labels : list or None
        Labels for species, if applicable.
    rules_labels : list or None
        Labels for rules, if applicable.
    xdata_labels : list or None
        Labels for xdata, if applicable.
    mod_data_labels : list or None
        Labels for modified data, if applicable.
    invalid_states : list or None
        A list to store invalid state inputs.
    parameters : list or None
        A list to track parameters used in the scan.
    parameter_labels : list
        Labels for the parameters.
    HAS_FLUXES : bool
        Indicates whether flux data is present.
    HAS_SPECIES : bool
        Indicates whether species data is present.
    HAS_RULES : bool
        Indicates whether rules data is present.
    HAS_XDATA : bool
        Indicates whether xdata is present.
    HAS_MOD_DATA : bool
        Indicates whether modified data is present.
    HAS_SET_LABELS : bool
        Indicates whether the labels have been set from an initial state data object.
    ALL_VALID : bool
        Indicates whether all states in the scan are valid.
    OPEN : bool
        Indicates whether the scan session is open for data addition.

    Methods
    -------
    __init__(par_label)
        Initializes a ScanDataObj instance with a given parameter label.

    setLabels(ssdata)
        Sets the labels for species, fluxes, rules, and xdata based on an input state data object.

    addPoint(ipar, ssdata)
        Adds a point to the scan, validating the state data object.

    addModData(mod, *args)
        Adds modified data attributes to the current scan session.

    closeScan()
        Finalizes the scan, preventing further additions, and converts lists to numpy arrays.

    getSpecies(lbls=False)
        Returns species data, optionally with labels.

    getFluxes(lbls=False)
        Returns fluxes data, optionally with labels.

    getRules(lbls=False)
        Returns rules data, optionally with labels.

    getXData(lbls=False)
        Returns xdata, optionally with labels.

    getModData(lbls=False)
        Returns modified data, optionally with labels.

    getAllScanData(lbls=False)
        Returns all scan data, optionally with labels.

    getScanData(*args, **kwargs)
        Returns specific scan data based on input labels, optionally with labels.
    """
    species = None
    fluxes = None
    rules = None
    xdata = None
    mod_data = None
    flux_labels = None
    species_labels = None
    rules_labels = None
    xdata_labels = None
    mod_data_labels = None
    invalid_states = None
    parameters = None
    parameter_labels = None
    HAS_FLUXES = False
    HAS_SPECIES = False
    HAS_RULES = False
    HAS_XDATA = False
    HAS_MOD_DATA = False
    HAS_SET_LABELS = False
    ALL_VALID = True
    OPEN = True

    def __init__(self, par_label):
        """
        Initialize a ScanDataObj instance.

        This constructor sets up a ScanDataObj with several lists initialized 
        to empty. The `parameter_labels` attribute is initialized based on 
        the provided `par_label`. If `par_label` is a list, it is directly 
        used to set `parameter_labels`; otherwise, it's converted into a 
        single-element list.

        Parameters
        ----------
        par_label : str or list of str
            The label(s) for the parameters. If a single string is provided,
            it is converted into a list containing one string. If a list is
            provided, it is used directly.

        Attributes
        ----------
        species : list
            A list to store species data.
        fluxes : list
            A list to store flux data.
        rules : list
            A list to store rules data.
        xdata : list
            A list to store x-axis data points for plotting or analysis.
        mod_data : list
            A list to store modified data.
        invalid_states : list
            A list to store states that are considered invalid under certain
            conditions.
        parameters : list
            A list to store parameters associated with the scan.
        parameter_labels : list of str
            A list to store labels for the parameters.

        Examples
        --------
        Create a ScanDataObj with a single parameter label:

        >>> scan_obj = ScanDataObj('Parameter1')
        >>> print(scan_obj.parameter_labels)
        ['Parameter1']

        Create a ScanDataObj with multiple parameter labels:

        >>> scan_obj = ScanDataObj(['Parameter1', 'Parameter2'])
        >>> print(scan_obj.parameter_labels)
        ['Parameter1', 'Parameter2']
        """
        self.species = []
        self.fluxes = []
        self.rules = []
        self.xdata = []
        self.mod_data = []
        self.invalid_states = []
        self.parameters = []
        if isinstance(par_label, list):
            self.parameter_labels = par_label
        else:
            self.parameter_labels = [par_label]

    def setLabels(self, ssdata):
        """
    Set labels for species, fluxes, rules, and other attributes in the ScanDataObj.

    This method extracts label information from the `ssdata` object and sets them
    on the ScanDataObj. It also sets flags indicating the presence of each type of
    label.

    Parameters
    ----------
    ssdata : object
        An object containing various label attributes that need to be set on the 
        ScanDataObj instance. The object should have attributes such as 
        `species_labels`, `flux_labels`, `rules_labels`, `xdata_labels`, and
        boolean flags like `HAS_SPECIES`, `HAS_FLUXES`, `HAS_RULES`, and `HAS_XDATA`
        indicating the existence of each label type.

    Attributes
    ----------
    species_labels : list
        Labels corresponding to species in the ScanDataObj. Present only if 
        `ssdata.HAS_SPECIES` is True.
    flux_labels : list
        Labels corresponding to fluxes in the ScanDataObj. Present only if 
        `ssdata.HAS_FLUXES` is True.
    rules_labels : list
        Labels corresponding to rules in the ScanDataObj. Present only if 
        `ssdata.HAS_RULES` is True.
    xdata_labels : list
        Labels corresponding to xdata in the ScanDataObj. Present only if 
        `ssdata.HAS_XDATA` is True.
    HAS_SET_LABELS : bool
        A flag set to True after `setLabels` is successfully executed.

    Notes
    -----
    Before calling this method, ensure that `ssdata` contains the relevant
    label attributes and that the respective boolean flags are correctly set.

    Example
    -------
    >>> ssdata = SomeDataClass()  # An instance with appropriate label attributes
    >>> ssdata.HAS_SPECIES = True
    >>> ssdata.species_labels = ['species1', 'species2']
    >>> sd_obj = ScanDataObj(parameter_labels=['param1'])
    >>> sd_obj.setLabels(ssdata)
    >>> sd_obj.species_labels
    ['species1', 'species2']

    """
        if ssdata.HAS_SPECIES:
            self.species_labels = ssdata.species_labels
            self.HAS_SPECIES = True
        if ssdata.HAS_FLUXES:
            self.flux_labels = ssdata.flux_labels
            self.HAS_FLUXES = True
        if ssdata.HAS_RULES:
            self.rules_labels = ssdata.rules_labels
            self.HAS_RULES = True
        if ssdata.HAS_XDATA:
            self.xdata_labels = ssdata.xdata_labels
            self.HAS_XDATA = True
        self.HAS_SET_LABELS = True

    def addPoint(self, ipar, ssdata):
        """
    Adds a data point consisting of input parameter values and associated `ssdata` object.

    This method updates the internal state of the `ScanDataObj` instance by appending parameter values,
    species data, fluxes, rules, and xdata to the corresponding attributes. If the labels have not been
    set, they will be initialized using the provided `ssdata` object. If any of the states in `ssdata`
    are invalid, it marks the scan as not fully valid.

    Parameters
    ----------
    ipar : list or array
        A list or array of input parameter values to be added to the scan.
    
    ssdata : object
        An object representing the current state of the scan data, which must have methods:
        - `getSpecies()`: Returns the species data.
        - `getFluxes()`: Returns the fluxes data.
        - `getRules()`: Returns the rules data.
        - `getXData()`: Returns the xdata.
        It should also have the following attributes:
        - `HAS_SPECIES`, `HAS_FLUXES`, `HAS_RULES`, `HAS_XDATA` (boolean flags).
        - `IS_VALID` (boolean indicating the validity of the data).

    Raises
    ------
    AssertionError
        If the scan has been finalized, indicated by the `OPEN` attribute being `False`.

    Notes
    -----
    - This method assumes that the ssdata methods and attributes (`getSpecies()`, `getFluxes()`, etc.)
      are correctly implemented and available.
    - The validity of the entire scan may be affected if any of the added points are invalid.

    Examples
    --------
    >>> scan_data = ScanDataObj(['param1'])
    >>> ssdata = SomeSSDataObject()  # Assuming this is a valid ssdata object with necessary methods
    >>> parameters = [1.2, 3.4, 5.6]
    >>> scan_data.OPEN = True  # Open the scan for adding data
    >>> scan_data.HAS_SET_LABELS = False  # Labels need to be set initially
    >>> scan_data.addPoint(parameters, ssdata)
    >>> print(scan_data.parameters)
    [[1.2, 3.4, 5.6]]
    >>> print(scan_data.species)
    # Output from ssdata.getSpecies() for this point

    """
        assert self.OPEN, '\nScan has been finalised no new data may be added'
        if not self.HAS_SET_LABELS:
            self.setLabels(ssdata)
        if hasattr(ipar, '__iter__'):
            self.parameters.append(ipar)
        else:
            self.parameters.append([ipar])
        if self.HAS_SPECIES:
            self.species.append(ssdata.getSpecies())
        if self.HAS_FLUXES:
            self.fluxes.append(ssdata.getFluxes())
        if self.HAS_RULES:
            self.rules.append(ssdata.getRules())
        if self.HAS_XDATA:
            self.xdata.append(ssdata.getXData())
        if not ssdata.IS_VALID:
            self.ALL_VALID = False
            self.invalid_states.append(ipar)

    def addModData(self, mod, *args):
        """
    Add modified data from a model object to the scan data.

    This method extracts attributes from the given model object `mod`, 
    defined by `args`, and appends them as a tuple to the `mod_data` list. 
    If the `mod_data_labels` are not already set, they will be initialized 
    with the provided attributes from `args`.

    Parameters
    ----------
    mod : object
        The model object from which to extract data attributes. The model 
        should have attributes specified in `args` that will be extracted 
        and stored.
    *args : str
        A variable length argument list of attribute names to be extracted 
        from the model object `mod`.

    Attributes
    ----------
    HAS_MOD_DATA : bool
        Indicates whether the modified data has been previously added. 
        It's set to `True` if not already.
    mod_data_labels : list of str
        Lists the attribute names that are being tracked and extracted from 
        the model object `mod`.
    mod_data : list of tuples
        A collection of tuples where each tuple contains the extracted 
        attribute values from `mod`.

    Example
    -------
    Let's assume a model object with attributes `a` and `b`:

    >>> class Model:
    ...     def __init__(self):
    ...         self.a = 10
    ...         self.b = 20
    >>> mod = Model()
    >>> scan_data_obj = ScanDataObj()
    >>> scan_data_obj.addModData(mod, 'a', 'b')
    >>> print(scan_data_obj.mod_data)
    [(10, 20)]
    """
        if not self.HAS_MOD_DATA:
            self.HAS_MOD_DATA = True
            self.mod_data_labels = [attr for attr in args if hasattr(mod, attr)
                ]
        self.mod_data.append(tuple([getattr(mod, attr) for attr in self.
            mod_data_labels]))

    def closeScan(self):
        """
    Finalize the scan, converting all collected data to numpy arrays.

    This method finalizes the current scan by converting all the stored
    data lists into numpy arrays. Once the scan is closed, no new data
    can be added. A diagnostic message is printed to confirm closure.

    Attributes
    ----------
    self.species : list
        If `HAS_SPECIES` is True, converts the `species` list to a numpy array.
    self.fluxes : list
        If `HAS_FLUXES` is True, converts the `fluxes` list to a numpy array.
    self.rules : list
        If `HAS_RULES` is True, converts the `rules` list to a numpy array.
    self.xdata : list
        If `HAS_XDATA` is True, converts the `xdata` list to a numpy array.
    self.mod_data : list
        If `HAS_MOD_DATA` is True, converts the `mod_data` list to a numpy array.
    self.parameters : list
        Converts the `parameters` list to a numpy array.
    self.OPEN : bool
        Sets `OPEN` to False to indicate that the scan has been closed.

    Examples
    --------
    >>> scan_data_obj = ScanDataObj()
    >>> scan_data_obj.addPoint([1.0, 2.0, 3.0], ssdata)
    >>> scan_data_obj.closeScan()
    INFO: closing scan no new data may be added.
    >>> print(scan_data_obj.OPEN)
    False
    >>> print(type(scan_data_obj.parameters))
    <class 'numpy.ndarray'>
    """
        print('\nINFO: closing scan no new data may be added.')
        if self.HAS_SPECIES:
            self.species = numpy.array(self.species)
        if self.HAS_FLUXES:
            self.fluxes = numpy.array(self.fluxes)
        if self.HAS_RULES:
            self.rules = numpy.array(self.rules)
        if self.HAS_XDATA:
            self.xdata = numpy.array(self.xdata)
        if self.HAS_MOD_DATA:
            self.mod_data = numpy.array(self.mod_data)
        self.parameters = numpy.array(self.parameters)
        self.OPEN = False

    def getSpecies(self, lbls=False):
        """
    Retrieve the species data, optionally with labels.

    This method returns the species data, combined with the parameters.
    If `lbls` is set to `True`, it also returns the labels for the parameters
    and the species. If the scan is still open, it will first be closed by
    calling `closeScan()` to ensure all data is prepared for retrieval.

    Parameters
    ----------
    lbls : bool, optional
        If True, return a tuple containing both the species data and the
        corresponding labels. Default is False, which means only the species
        data will be returned without labels.

    Returns
    -------
    output : numpy.ndarray or None
        The combined array of parameters and species data. It returns `None`
        if `HAS_SPECIES` is False.
    labels : list of str, optional
        A list of strings representing the labels for parameters and species.
        This is only returned if `lbls` is True.

    Examples
    --------
    >>> scan_data_obj = ScanDataObj()
    >>> scan_data_obj.parameters = numpy.array([[1, 2], [3, 4]])
    >>> scan_data_obj.species = numpy.array([[5, 6], [7, 8]])
    >>> scan_data_obj.parameter_labels = ['param1', 'param2']
    >>> scan_data_obj.species_labels = ['species1', 'species2']
    >>> scan_data_obj.HAS_SPECIES = True
    >>> scan_data_obj.getSpecies()
    array([[1, 2, 5, 6],
           [3, 4, 7, 8]])
    >>> data, lbls = scan_data_obj.getSpecies(lbls=True)
    >>> data
    array([[1, 2, 5, 6],
           [3, 4, 7, 8]])
    >>> lbls
    ['param1', 'param2', 'species1', 'species2']
    """
        if self.OPEN:
            self.closeScan()
        output = None
        labels = None
        if self.HAS_SPECIES:
            output = numpy.hstack((self.parameters, self.species))
            labels = self.parameter_labels + self.species_labels
        if lbls:
            return output, labels
        else:
            return output

    def getFluxes(self, lbls=False):
        """
    Retrieve the flux data from the scan, optionally including labels.

    Checks if the scan is open and closes it if necessary. If flux data
    is available, it combines the flux data with parameters and returns
    the resulting array. Optionally, it can also return associated labels.

    Parameters
    ----------
    lbls : bool, optional
        If True, returns a tuple containing both the data array and 
        the labels. Otherwise, only returns the data array. Default is False.

    Returns
    -------
    numpy.ndarray or tuple of (numpy.ndarray, list)
        If `lbls` is False, returns the flux data combined with parameters
        as a NumPy array. If `lbls` is True, also returns the labels as a list.

    Examples
    --------
    >>> scan_data = ScanDataObj()
    >>> flux_data = scan_data.getFluxes()
    >>> print(flux_data)
    array([[1.5, 0.3, 0.2, ...],
           [1.6, 0.4, 0.1, ...],
           ...])

    >>> flux_data, labels = scan_data.getFluxes(lbls=True)
    >>> print(labels)
    ['param1', 'param2', 'flux1', ...]
    >>> print(flux_data)
    array([[1.5, 0.3, 0.2, ...],
           [1.6, 0.4, 0.1, ...],
           ...])
    """
        if self.OPEN:
            self.closeScan()
        output = None
        labels = None
        if self.HAS_FLUXES:
            output = numpy.hstack((self.parameters, self.fluxes))
            labels = self.parameter_labels + self.flux_labels
        if lbls:
            return output, labels
        else:
            return output

    def getRules(self, lbls=False):
        """
    Retrieve the rules data along with associated parameters.

    This method retrieves a concatenated array of parameters and rules
    data if the object has rules data available. Optionally, it can also
    return the labels for these data.

    Parameters
    ----------
    lbls : bool, optional
        If set to True, the method will return both the data array 
        and its corresponding labels. If False, only the data array 
        is returned. Defaults to False.

    Returns
    -------
    numpy.ndarray or tuple
        If `lbls` is False, returns a numpy array containing the 
        parameters and rules data. If `lbls` is True, returns a tuple
        where the first element is the numpy array of data and the 
        second element is a list of corresponding labels.

    Notes
    -----
    If there are no rules data (`HAS_RULES` is False), `None` will 
    be returned.

    Examples
    --------
    >>> scan_data = ScanDataObj()
    >>> data = scan_data.getRules()
    >>> data_with_labels = scan_data.getRules(lbls=True)
    """
        if self.OPEN:
            self.closeScan()
        output = None
        labels = None
        if self.HAS_RULES:
            output = numpy.hstack((self.parameters, self.rules))
            labels = self.parameter_labels + self.rules_labels
        if lbls:
            return output, labels
        else:
            return output

    def getXData(self, lbls=False):
        """
    Retrieve the X data from the scan, optionally with labels.

    This method checks if the scan is currently open and closes it if necessary.
    It then retrieves the X data from the scan, appending it to any associated
    parameters. If the `HAS_XDATA` attribute is True, the X data are combined
    with parameters using a horizontal stack. If `lbls` is set to True, the method
    also returns corresponding labels.

    Parameters
    ----------
    lbls : bool, optional
        A flag indicating whether to include labels with the output. If True,
        the method returns a tuple containing the output data and the labels.
        Defaults to False.

    Returns
    -------
    output : numpy.ndarray or None
        A stacked array of scan parameters and X data. Returns None if `HAS_XDATA`
        is False.
    labels : list of str, optional
        A list of labels corresponding to the output data. Returned only if
        `lbls` is True.

    Examples
    --------
    >>> scan_data_obj = ScanDataObj()
    >>> xdata = scan_data_obj.getXData()
    >>> xdata_with_labels = scan_data_obj.getXData(lbls=True)
    >>> data, labels = xdata_with_labels
    """
        if self.OPEN:
            self.closeScan()
        output = None
        labels = None
        if self.HAS_XDATA:
            output = numpy.hstack((self.parameters, self.xdata))
            labels = self.parameter_labels + self.xdata_labels
        if lbls:
            return output, labels
        else:
            return output

    def getModData(self, lbls=False):
        """
    Retrieve the modified data associated with this scan, along with optional labels.

    If the scan is currently open, it will be closed before retrieving the data.

    Parameters
    ----------
    lbls : bool, optional
        If True, returns both the modified data and their corresponding labels.
        If False, only the modified data is returned. Default is False.

    Returns
    -------
    numpy.ndarray or tuple of (numpy.ndarray, list of str)
        If `lbls` is False, returns a numpy array containing the modified data.
        If `lbls` is True, returns a tuple containing:
        - a numpy array with the modified data
        - a list of strings representing the labels of the modified data

    Notes
    -----
    - The method automatically closes the scan if it is open.
    - The modified data includes the merged parameters and modifiable data fields.

    Example
    -------
    >>> scan_data_obj = ScanDataObj()
    >>> mod_data = scan_data_obj.getModData(lbls=False)
    >>> mod_data_with_labels = scan_data_obj.getModData(lbls=True)
    >>> print(mod_data_with_labels[0])  # Output modified data
    >>> print(mod_data_with_labels[1])  # Output labels associated with the data
    """
        if self.OPEN:
            self.closeScan()
        output = None
        labels = None
        if self.HAS_MOD_DATA:
            output = numpy.hstack((self.parameters, self.mod_data))
            labels = self.parameter_labels + self.mod_data_labels
        if lbls:
            return output, labels
        else:
            return output

    def getAllScanData(self, lbls=False):
        """
    Retrieves all available scan data from the object, including parameters,
    species, fluxes, rules, x-data, and modified data. Optionally returns associated labels.

    Parameters
    ----------
    lbls : bool, optional
        If True, returns a tuple with the data and corresponding labels;
        otherwise, only the data is returned. Default is False.

    Returns
    -------
    numpy.ndarray or tuple of (numpy.ndarray, list of str)
        Returns a NumPy array containing all scan data. If `lbls` is True,
        also returns a list of labels corresponding to the columns in the
        data array.

    Examples
    --------
    >>> scan_obj = ScanDataObj()
    >>> data = scan_obj.getAllScanData()
    >>> print(data)
    [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], ...]  # Example output

    >>> data, labels = scan_obj.getAllScanData(lbls=True)
    >>> print(data)
    [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], ...]  # Example output
    >>> print(labels)
    ['param1', 'param2', 'species1', 'flux1', ...]  # Example output
    """
        if self.OPEN:
            self.closeScan()
        output = self.parameters
        labels = self.parameter_labels
        if self.HAS_SPECIES:
            output = numpy.hstack((output, self.species))
            labels = labels + self.species_labels
        if self.HAS_FLUXES:
            output = numpy.hstack((output, self.fluxes))
            labels = labels + self.flux_labels
        if self.HAS_RULES:
            output = numpy.hstack((output, self.rules))
            labels = labels + self.rules_labels
        if self.HAS_XDATA:
            output = numpy.hstack((output, self.xdata))
            labels = labels + self.xdata_labels
        if self.HAS_MOD_DATA:
            output = numpy.hstack((output, self.mod_data))
            labels = labels + self.mod_data_labels
        if lbls:
            return output, labels
        else:
            return output

    def getScanData(self, *args, **kwargs):
        """
    Fetches the scan data corresponding to the given labels for species, flux, rules, or mod data.
    
    It returns an array containing parameters and the selected types of data, which are specified via
    the labels provided in the arguments. The method supports optional retrieval of labels alongside
    the data if requested.

    Parameters
    ----------
    *args : str
        Variable-length argument list. Each argument should be a label corresponding to species, flux,
        rules, or mod data you want to retrieve.
    **kwargs : dict, optional
        lbls : bool
            If set to True, the method returns a tuple consisting of the output data and the corresponding labels.
            Defaults to False.

    Returns
    -------
    numpy.ndarray or tuple
        If `lbls` is False, returns a numpy array of the parameters and the selected data.
        If `lbls` is True, returns a tuple containing:
        - output : numpy.ndarray
            An array of the parameters and the selected data.
        - lout : list of str
            The list of labels for each element in the output array.

    Examples
    --------
    >>> scan_obj = ScanDataObj()
    >>> scan_obj.getScanData('species_label_1', 'flux_label_2')
    array([...])

    >>> scan_obj.getScanData('species_label_1', 'flux_label_2', lbls=True)
    (array([...]), ['param_label_1', 'species_label_1', 'flux_label_2'])
    """
        if 'lbls' in kwargs:
            lbls = kwargs['lbls']
        else:
            lbls = False
        output = self.parameters
        lout = self.parameter_labels
        for roc in args:
            if self.HAS_SPECIES and roc in self.species_labels:
                lout.append(roc)
                output = numpy.hstack((output, self.species.take([self.
                    species_labels.index(roc)], axis=-1)))
            if self.HAS_FLUXES and roc in self.flux_labels:
                lout.append(roc)
                output = numpy.hstack((output, self.fluxes.take([self.
                    flux_labels.index(roc)], axis=-1)))
            if self.HAS_RULES and roc in self.rules_labels:
                lout.append(roc)
                output = numpy.hstack((output, self.rules.take([self.
                    rules_labels.index(roc)], axis=-1)))
            if self.HAS_XDATA and roc in self.xdata_labels:
                lout.append(roc)
                output = numpy.hstack((output, self.xdata.take([self.
                    xdata_labels.index(roc)], axis=-1)))
            if self.HAS_MOD_DATA and roc in self.mod_data_labels:
                lout.append(roc)
                output = numpy.hstack((output, self.mod_data.take([self.
                    mod_data_labels.index(roc)], axis=-1)))
        if not lbls:
            return output
        else:
            return output, lout


class StateDataObj(object):
    """
    A class to store and manage steady-state data for a variety of biochemical
    simulations or analyses. This class facilitates the storage and retrieval
    of different categories of data, such as species concentrations, fluxes,
    rules, and additional simulation data.

    Attributes
    ----------
    fluxes : numpy.ndarray or None
        An array storing the steady-state fluxes of the system.
    species : numpy.ndarray or None
        An array storing the steady-state concentrations of species.
    rules : numpy.ndarray or None
        An array storing the results of rate rules.
    xdata : numpy.ndarray or None
        An array storing additional simulation data.
    flux_labels : list of str or None
        A list of labels corresponding to the flux data.
    species_labels : list of str or None
        A list of labels corresponding to the species data.
    rules_labels : list of str or None
        A list of labels corresponding to the rules data.
    xdata_labels : list of str or None
        A list of labels corresponding to the xdata.
    HAS_FLUXES : bool
        A flag indicating whether flux data is available.
    HAS_SPECIES : bool
        A flag indicating whether species data is available.
    HAS_RULES : bool
        A flag indicating whether rules data is available.
    HAS_XDATA : bool
        A flag indicating whether additional simulation xdata is available.
    IS_VALID : bool
        A flag indicating the validity of the state data object.
    
    Methods
    -------
    setSpecies(species, lbls=None)
        Sets the species array and optionally labels for these species.
    setFluxes(fluxes, lbls=None)
        Sets the flux array and optionally labels for these fluxes.
    setRules(rules, lbls=None)
        Sets the rules array and optionally labels for these rules.
    setXData(xdata, lbls=None)
        Sets the additional simulation xdata and optionally labels.
    getSpecies(lbls=False)
        Returns the species array, optionally with labels.
    getFluxes(lbls=False)
        Returns the flux array, optionally with labels.
    getRules(lbls=False)
        Returns the rules array, optionally with labels.
    getXData(lbls=False)
        Returns the xdata array, optionally with labels.
    getAllStateData(lbls=False)
        Returns all available data as a combined array, optionally with labels.
    getStateData(*args, **kwargs)
        Returns an array of data for specified species or rate labels, optionally with labels.
    """
    fluxes = None
    species = None
    rules = None
    xdata = None
    flux_labels = None
    species_labels = None
    rules_labels = None
    xdata_labels = None
    HAS_FLUXES = False
    HAS_SPECIES = False
    HAS_RULES = False
    HAS_XDATA = False
    IS_VALID = True

    def setSpecies(self, species, lbls=None):
        """
    Set the species array and optionally set species labels.

    This method assigns the provided `species` array to the instance and
    sets the `HAS_SPECIES` attribute to `True`. If `lbls` is provided,
    it updates the `species_labels` attribute with the provided labels
    and creates attributes on the instance corresponding to each label
    with values from the `species` array.

    Parameters
    ----------
    species : array-like
        The array representing species data to be associated with the object.
    
    lbls : list of str, optional
        A list of labels to be used as attribute names for each species element.
        If not provided, no labels or attributes are set.

    Examples
    --------
    >>> state_obj = StateDataObj()
    >>> species_data = [0.1, 0.2, 0.3]
    >>> species_labels = ['species1', 'species2', 'species3']
    >>> state_obj.setSpecies(species_data, species_labels)
    >>> print(state_obj.species)
    [0.1, 0.2, 0.3]
    >>> print(state_obj.species_labels)
    ['species1', 'species2', 'species3']
    >>> print(state_obj.species1)
    0.1
    """
        self.species = species
        self.HAS_SPECIES = True
        if lbls != None:
            self.species_labels = lbls
            for s in range(len(self.species_labels)):
                setattr(self, self.species_labels[s], self.species[s])

    def setFluxes(self, fluxes, lbls=None):
        """Set the flux array.

    This method initializes the `fluxes` attribute of the `StateDataObj` 
    instance with the provided array of flux values. It also sets the 
    `HAS_FLUXES` flag to `True`.

    Optionally, if a list of labels is provided, these labels are used to 
    create additional attributes on the instance that directly map to the 
    respective flux values.

    Parameters
    ----------
    fluxes : array_like
        An array that contains the flux values to set for the object.

    lbls : list of str, optional
        A list of labels corresponding to each flux value. If provided, 
        attributes will be dynamically added to the instance, with each 
        attribute name being the label and the attribute value being the 
        corresponding flux.

    Examples
    --------
    >>> obj = StateDataObj()
    >>> fluxes = [0.1, 0.2, 0.3]
    >>> labels = ['flux1', 'flux2', 'flux3']
    >>> obj.setFluxes(fluxes, labels)
    >>> obj.fluxes
    [0.1, 0.2, 0.3]
    >>> obj.flux1
    0.1
    >>> obj.flux2
    0.2
    >>> obj.flux3
    0.3
    """
        self.fluxes = fluxes
        self.HAS_FLUXES = True
        if lbls != None:
            self.flux_labels = lbls
            for f in range(len(self.flux_labels)):
                setattr(self, self.flux_labels[f], self.fluxes[f])

    def setRules(self, rules, lbls=None):
        """
    Set the results of rate rules.

    This method sets the `rules` attribute with the provided rule values
    and updates the `HAS_RULES` attribute to indicate that rules have
    been defined. Optionally, if labels are provided, it assigns each rule 
    value to a dynamically created attribute named after each label.

    Parameters
    ----------
    rules : list
        A list of rule values to be stored.
    lbls : list of str, optional
        A list of label strings corresponding to each rule value in 
        `rules`. If provided, attributes will be dynamically assigned to 
        the instance using these labels.

    Attributes
    ----------
    HAS_RULES : bool
        A flag that indicates whether rules are set for this instance.
    rules : list
        Stores the provided rule values.
    rules_labels : list of str, optional
        Stores the provided labels if any.

    Examples
    --------
    >>> data = StateDataObj()
    >>> data.setRules([0.1, 0.2, 0.3], lbls=['rule1', 'rule2', 'rule3'])
    >>> data.rules
    [0.1, 0.2, 0.3]
    >>> data.HAS_RULES
    True
    >>> data.rule1
    0.1
    >>> data.rule2
    0.2
    >>> data.rule3
    0.3
    """
        self.rules = rules
        self.HAS_RULES = True
        if lbls != None:
            self.rules_labels = lbls
            for r in range(len(self.rules_labels)):
                setattr(self, self.rules_labels[r], self.rules[r])

    def setXData(self, xdata, lbls=None):
        """
    Sets extra simulation data.

    This method allows you to set additional simulation data for the `StateDataObj` instance. 
    If provided, labels for the data entries will be set as attributes of the instance, 
    allowing for easier access to individual data elements.

    Parameters
    ----------
    xdata : list or array-like
        A list or array containing the extra simulation data.
    lbls : list of str, optional
        A list of labels for each element in `xdata`. If provided, each label will be set 
        as an attribute of the instance referring to the corresponding element in `xdata`.

    Attributes
    ----------
    xdata : list or array-like
        Stores the extra simulation data.
    HAS_XDATA : bool
        Flag indicating whether extra simulation data has been set. It's set to True after calling this method.
    xdata_labels : list of str, optional
        Stores the labels for the extra simulation data, if provided.

    Examples
    --------
    >>> state_obj = StateDataObj()
    >>> xdata = [0.5, 1.2, 3.3]
    >>> labels = ['rate1', 'rate2', 'rate3']
    >>> state_obj.setXData(xdata, labels)
    >>> state_obj.rate1
    0.5
    >>> state_obj.rate2
    1.2
    >>> state_obj.rate3
    3.3
    """
        self.xdata = xdata
        self.HAS_XDATA = True
        if lbls != None:
            self.xdata_labels = lbls
            for x in range(len(self.xdata_labels)):
                setattr(self, self.xdata_labels[x], self.xdata[x])

    def getSpecies(self, lbls=False):
        """
    Return species array.

    Depending on the value of `lbls`, this method returns either the
    species array or a tuple containing the species array and its labels.

    Parameters
    ----------
    lbls : bool, optional
        If `True`, return a tuple containing both the species array and
        its labels. If `False` (default), return only the species array.

    Returns
    -------
    species : array or None
        The array of species. Returns `None` if `HAS_SPECIES` is `False`.
    tuple of (array, list) or None
        A tuple consisting of the species array and a list of species labels,
        if `lbls` is `True`. Returns `None` if `HAS_SPECIES` is `False`.

    Notes
    -----
    This method assumes that the `species` and `species_labels` attributes
    exist if `HAS_SPECIES` is `True`.

    Example
    -------
    >>> state_data = StateDataObj()
    >>> state_data.HAS_SPECIES = True
    >>> state_data.species = ['species_1', 'species_2', 'species_3']
    >>> state_data.species_labels = ['label_1', 'label_2', 'label_3']
    >>> species = state_data.getSpecies()
    >>> print(species)
    ['species_1', 'species_2', 'species_3']

    >>> species_with_labels = state_data.getSpecies(lbls=True)
    >>> print(species_with_labels)
    (['species_1', 'species_2', 'species_3'], ['label_1', 'label_2', 'label_3'])
    """
        output = None
        if self.HAS_SPECIES:
            output = self.species
        if not lbls:
            return output
        else:
            return output, self.species_labels

    def getFluxes(self, lbls=False):
        """
    Return the flux array.

    Parameters
    ----------
    lbls : bool, optional
        If True, the method also returns the flux labels. The default is False.

    Returns
    -------
    output : list or None
        The flux array if `HAS_FLUXES` is True, otherwise None.
    flux_labels : list, optional
        The corresponding labels for the flux array, only returned if `lbls`
        is set to True.

    Examples
    --------
    >>> state_data = StateDataObj()
    >>> # Assuming state_data.HAS_FLUXES is True and state_data.fluxes is defined
    >>> fluxes = state_data.getFluxes()
    >>> print(fluxes)
    [1.2, 2.3, 3.4]

    >>> # To get both fluxes and their labels
    >>> fluxes, labels = state_data.getFluxes(lbls=True)
    >>> print(fluxes)
    [1.2, 2.3, 3.4]
    >>> print(labels)
    ['flux1', 'flux2', 'flux3']
    """
        output = None
        if self.HAS_FLUXES:
            output = self.fluxes
        if not lbls:
            return output
        else:
            return output, self.flux_labels

    def getRules(self, lbls=False):
        """Return rule array.

    This method retrieves the array of rules associated with the `StateDataObj` instance.
    If the `lbls` parameter is set to `True`, it returns a tuple containing the rules
    and their corresponding labels.

    Parameters
    ----------
    lbls : bool, optional
        If set to `True`, the method returns both the rules and their labels as a tuple.
        The default is `False`, which means only the rules will be returned.

    Returns
    -------
    output : array-like or tuple
        Returns an array of rules if `lbls` is `False`. If `lbls` is `True`, it returns
        a tuple containing the rule array and an array of rule labels.

    Examples
    --------
    >>> state_data = StateDataObj()
    >>> rules = state_data.getRules()
    >>> print(rules)
    [Rule1, Rule2, Rule3]

    >>> rules_with_labels = state_data.getRules(lbls=True)
    >>> print(rules_with_labels)
    ([Rule1, Rule2, Rule3], ['Label1', 'Label2', 'Label3'])
    """
        output = None
        if self.HAS_RULES:
            output = self.rules
        if not lbls:
            return output
        else:
            return output, self.rules_labels

    def getXData(self, lbls=False):
        """Return xdata array.

    Retrieves the xdata array associated with the `StateDataObj` instance. If `lbls` is set to True,
    the method will also return the labels corresponding to the xdata array.

    Parameters
    ----------
    lbls : bool, optional
        A flag indicating whether to return xdata labels along with the xdata array. Default is False.

    Returns
    -------
    ndarray or tuple of (ndarray, array-like)
        If `HAS_XDATA` is True and `lbls` is False, returns the xdata array.
        If `HAS_XDATA` is True and `lbls` is True, returns a tuple containing the xdata array and the xdata labels.
        If `HAS_XDATA` is False, returns None (or optionally (None, None) if `lbls` is also True).

    Raises
    ------
    AttributeError
        If the `HAS_XDATA` attribute is not defined on the instance.

    Examples
    --------
    >>> state_data = StateDataObj()
    >>> state_data.HAS_XDATA = True
    >>> state_data.xdata = [1.0, 2.0, 3.0]
    >>> state_data.xdata_labels = ['x1', 'x2', 'x3']

    >>> # Retrieve only xdata
    >>> xdata = state_data.getXData()
    >>> print(xdata)
    [1.0, 2.0, 3.0]

    >>> # Retrieve both xdata and labels
    >>> xdata, labels = state_data.getXData(lbls=True)
    >>> print(xdata)
    [1.0, 2.0, 3.0]
    >>> print(labels)
    ['x1', 'x2', 'x3']
    """
        output = None
        if self.HAS_XDATA:
            output = self.xdata
        if not lbls:
            return output
        else:
            return output, self.xdata_labels

    def getAllStateData(self, lbls=False):
        """
    Return all available data as well as labels across species, fluxes, rules, and xdata.

    This method combines all available data from species, fluxes, rules, and xdata
    into a single array. If `lbls` is set to `True`, it will return a tuple containing 
    the data array and corresponding labels. Otherwise, it only returns the data array.

    Parameters
    ----------
    lbls : bool, optional
        A flag indicating whether to return labels alongside the data array. 
        Default is `False`.

    Returns
    -------
    output : numpy.ndarray
        The combined array of all available data.
    (output, labels) : tuple of (numpy.ndarray, list)
        If `lbls` is `True`, returns a tuple where `output` is the combined array 
        of all available data and `labels` is a list of the corresponding labels.

    Notes
    -----
    This method assumes that the instance has boolean attributes `HAS_SPECIES`, 
    `HAS_FLUXES`, `HAS_RULES`, and `HAS_XDATA` that indicate the availability 
    of each respective data type. If data is not available for any type, 
    it will be skipped in the aggregation.

    Examples
    --------
    >>> state_data_obj = StateDataObj()
    >>> # Example when lbls is False
    >>> data = state_data_obj.getAllStateData(lbls=False)
    >>> print(data)
    array([...])  # Combined array of species, fluxes, and rules data

    >>> # Example when lbls is True
    >>> data, labels = state_data_obj.getAllStateData(lbls=True)
    >>> print(data)
    array([...])  # Combined array of species, fluxes, rules, and xdata
    >>> print(labels)
    ['species_label_1', 'species_label_2', ..., 'flux_label_1', ..., 'rule_label_1', ..., 'xdata_label_1', ...]
    """
        labels = []
        output = None
        if self.HAS_SPECIES:
            output = self.species
            labels += self.species_labels
        if self.HAS_FLUXES:
            if output is None:
                output = self.fluxes
            else:
                output = numpy.hstack((output, self.fluxes))
            labels += self.flux_labels
        if self.HAS_RULES:
            if output is None:
                output = self.rules
            else:
                output = numpy.hstack((output, self.rules))
            labels += self.rules_labels
        if self.HAS_XDATA:
            if output is None:
                output = self.xdata
            else:
                output = numpy.hstack((output, self.xdata))
            labels += self.xdata_labels
        if not lbls:
            return output
        else:
            return output, labels

    def getStateData(self, *args, **kwargs):
        """
    Retrieve specified species or rate data.

    This method takes a variable number of arguments, each representing a label
    for species or rate data. It returns an array containing time and the requested
    data points. If the `lbls` keyword argument is set to `True`, it returns a tuple
    containing the data array and a list of the corresponding labels.

    Parameters
    ----------
    *args : str
        Labels for species or rate data to be retrieved.
    **kwargs : keyword arguments
        Optional arguments:
        - lbls : bool, optional
            If `True`, return the labels associated with the data. Default is `False`.

    Returns
    -------
    numpy.ndarray or tuple
        An array containing the time and requested data. If `lbls` is `True`, a tuple
        of the data array and a list of labels is returned.

    Notes
    -----
    If a label is not available in the dataset, it will be ignored, and a message
    will be printed to standard output.

    Examples
    --------
    >>> state_data_obj = StateDataObj()
    >>> output = state_data_obj.getStateData('sp1', 'r1', lbls=True)
    >>> print(output)
    (array([...]), ['sp1', 'r1'])

    >>> output = state_data_obj.getStateData('sp1', 'non_existent_label')
    I don't have an attribute non_existent_label ... ignoring.
    >>> print(output)
    [...]
    """
        if 'lbls' in kwargs:
            lbls = kwargs['lbls']
        else:
            lbls = False
        lout = []
        output = []
        for roc in args:
            if self.HAS_SPECIES and roc in self.species_labels:
                lout.append(roc)
                output.append(self.species[self.species_labels.index(roc)])
            elif self.HAS_FLUXES and roc in self.flux_labels:
                lout.append(roc)
                output.append(self.fluxes[self.flux_labels.index(roc)])
            elif self.HAS_RULES and roc in self.rules_labels:
                lout.append(roc)
                output.append(self.rules[self.rules_labels.index(roc)])
            elif self.HAS_XDATA and roc in self.xdata_labels:
                lout.append(roc)
                output.append(self.xdata[self.xdata_labels.index(roc)])
            else:
                print("I don't have an attribute {} ... ignoring.".format(roc))
        if not lbls:
            return output
        else:
            return numpy.array(output), lout


class IntegrationDataObj(object):
    """
    IntegrationDataObj is a class designed to store and manage the results of a time-based simulation. 
    It offers methods for setting and retrieving time, species, rate, and additional simulation data. 
    This class provides utility functions for obtaining specific slices of data over time, making 
    it useful for analyzing the progress and results of simulations.

    Attributes
    ----------
    time : numpy.ndarray
        An array representing the time points of the simulation.
    rates : numpy.ndarray
        An array storing the rate data at each time point.
    species : numpy.ndarray
        An array storing the species data at each time point.
    rules : numpy.ndarray
        An array storing the results of rate rules at each time point.
    xdata : numpy.ndarray
        An array storing additional simulation data.
    time_label : str
        A label for the time data, default is 'Time'.
    rate_labels : list of str
        A list of labels for the rate data.
    species_labels : list of str
        A list of labels for the species data.
    rules_labels : list of str
        A list of labels for the rules data.
    xdata_labels : list of str
        A list of labels for the additional simulation data.
    HAS_SPECIES : bool
        Indicates if species data is available.
    HAS_RATES : bool
        Indicates if rate data is available.
    HAS_RULES : bool
        Indicates if rules data is available.
    HAS_TIME : bool
        Indicates if time data is available.
    HAS_XDATA : bool
        Indicates if additional simulation data is available.
    IS_VALID : bool
        Indicates if the data object is valid. Default is True.
    TYPE_INFO : str
        Describes the type of simulation, default is 'Deterministic'.

    Methods
    -------
    setLabels(species=None, rates=None, rules=None)
        Set the label lists for species, rates, and rules.
    setTime(time, lbl=None)
        Set the time vector for the simulation.
    setSpecies(species, lbls=None)
        Set the data array for species.
    setRates(rates, lbls=None)
        Set the data array for rates.
    setRules(rules, lbls=None)
        Set the data array for rules.
    setXData(xdata, lbls=None)
        Set the array for additional simulation data.
    getTime(lbls=False)
        Return the time array; optionally includes label.
    getSpecies(lbls=False)
        Return a combined time and species data array; optionally includes labels.
    getRates(lbls=False)
        Return a combined time and rate data array; optionally includes labels.
    getRules(lbls=False)
        Return a combined time and rules data array; optionally includes labels.
    getXData(lbls=False)
        Return a combined time and xdata array; optionally includes labels.
    getDataAtTime(time)
        Returns the data at a specified time point.
    getDataInTimeInterval(time, bounds=None)
        Returns data points within a specified time interval.
    getOutput(*args, **kwargs)
        An alias for getSimData(); returns selected simulation data arrays.
    getAllSimData(lbls=False)
        Returns all available data arrays with time; optionally includes labels.
    getSimData(*args, **kwargs)
        Returns an array of time plus selected simulation data based on species/rate labels.
    """
    time = None
    rates = None
    species = None
    rules = None
    xdata = None
    time_label = 'Time'
    rate_labels = None
    species_labels = None
    rules_labels = None
    xdata_labels = None
    HAS_SPECIES = False
    HAS_RATES = False
    HAS_RULES = False
    HAS_TIME = False
    HAS_XDATA = False
    IS_VALID = True
    TYPE_INFO = 'Deterministic'

    def setLabels(self, species=None, rates=None, rules=None):
        """
    Set the species, rate, and rule label lists.

    This method allows updating the internal label lists for species, rates, and rules.
    Each argument is optional; if provided, it replaces the current label list for the corresponding entity.

    Parameters
    ----------
    species : list of str, optional
        A list of labels for species. If `None`, the species label list remains unchanged.
    rates : list of str, optional
        A list of labels for rates. If `None`, the rate label list remains unchanged.
    rules : list of str, optional
        A list of labels for rules. If `None`, the rule label list remains unchanged.

    Examples
    --------
    >>> obj = IntegrationDataObj()
    >>> obj.setLabels(
    ...     species=['H2O', 'CO2', 'O2'],
    ...     rates=['rate1', 'rate2'],
    ...     rules=['ruleA', 'ruleB']
    ... )

    After calling `setLabels`, the object `obj` will have its species, rate, and rule label lists
    updated to the specified values.

    >>> obj.species_labels
    ['H2O', 'CO2', 'O2']
    >>> obj.rate_labels
    ['rate1', 'rate2']
    >>> obj.rules_labels
    ['ruleA', 'ruleB']
    """
        if species != None:
            self.species_labels = species
        if rates != None:
            self.rate_labels = rates
        if rules != None:
            self.rules_labels = rules

    def setTime(self, time, lbl=None):
        """
    Set the time vector for the integration data object.

    This method assigns a given time vector to the object and sets an
    indicator that the time data is available. Optionally, a label for
    the time vector can be assigned.

    Parameters
    ----------
    time : numpy.ndarray
        A 1-dimensional array representing the time points for the integration
        data. This will be reshaped into a column vector internally.
    lbl : str, optional
        The label for the time vector. Defaults to None if no label is provided.

    Attributes
    ----------
    time : numpy.ndarray
        The reshaped time vector stored as a column vector in the object.
    HAS_TIME : bool
        A flag set to True indicating that the time data has been initialized.
    time_label : str, optional
        The label for the time vector, if provided.

    Examples
    --------
    >>> import numpy as np
    >>> obj = IntegrationDataObj()
    >>> time_points = np.array([0, 1, 2, 3, 4])
    >>> obj.setTime(time_points, lbl='Simulation Time')
    >>> obj.time
    array([[0],
           [1],
           [2],
           [3],
           [4]])
    >>> obj.HAS_TIME
    True
    >>> obj.time_label
    'Simulation Time'
    """
        self.time = time.reshape(len(time), 1)
        self.HAS_TIME = True
        if lbl != None:
            self.time_label = lbl

    def setSpecies(self, species, lbls=None):
        """Set the species array.

    This method assigns the provided species matrix to the object's attribute
    and marks that the species data is set. Optionally, it can also set labels
    for each species if provided.

    Parameters
    ----------
    species : numpy.ndarray
        The array representing the species data. It should be structured
        such that its rows correspond to different time points and its columns
        to different species.

    lbls : list of str, optional
        A list of labels corresponding to the species. The length of this list
        should match the number of columns in the species array.

    Attributes
    ----------
    species : numpy.ndarray
        The stored species data for the object.
        
    HAS_SPECIES : bool
        A flag indicating whether the species data has been set.

    species_labels : list of str
        The labels for the species, if provided.

    Examples
    --------
    Set the species data without labels:

    >>> data_obj = IntegrationDataObj()
    >>> species_data = np.array([[0.1, 0.2], [0.3, 0.4]])
    >>> data_obj.setSpecies(species_data)

    Set the species data with labels:

    >>> species_labels = ['Species A', 'Species B']
    >>> data_obj.setSpecies(species_data, lbls=species_labels)
    """
        self.species = species
        self.HAS_SPECIES = True
        if lbls != None:
            self.species_labels = lbls

    def setRates(self, rates, lbls=None):
        """Set the rate array.

    This method assigns the passed rate data to the instance's `rates`
    attribute and marks the instance as having rates by setting
    `HAS_RATES` to True. Optionally, it also assigns provided labels
    for the rates.

    Parameters
    ----------
    rates : array_like
        The rates data to be set, typically a NumPy array or similar
        structure.
    lbls : array_like, optional
        An array of labels corresponding to each element in the `rates`
        data. If provided, it should have the same length as the number
        of rates, and it will be used to set the `rate_labels` attribute.

    Attributes
    ----------
    rates : array_like
        Holds the rates data provided by the user.
    HAS_RATES : bool
        A boolean flag indicating whether the rate data has been set.
    rate_labels : array_like, optional
        If the `lbls` parameter is provided, this attribute will store
        the label information for each rate.

    Examples
    --------
    >>> integration_data = IntegrationDataObj()
    >>> rates = np.array([0.1, 0.2, 0.3])
    >>> rate_labels = ['reaction1_rate', 'reaction2_rate', 'reaction3_rate']
    >>> integration_data.setRates(rates, rate_labels)
    >>> integration_data.rates
    array([0.1, 0.2, 0.3])
    >>> integration_data.HAS_RATES
    True
    >>> integration_data.rate_labels
    ['reaction1_rate', 'reaction2_rate', 'reaction3_rate']
    """
        self.rates = rates
        self.HAS_RATES = True
        if lbls != None:
            self.rate_labels = lbls

    def setRules(self, rules, lbls=None):
        """Set the results of rate rules.

    This method updates the `rules` attribute of the `IntegrationDataObj` instance, indicating that rate rules have been defined for the object.
    If labels are provided, they are stored in the `rules_labels` attribute.

    Parameters
    ----------
    rules : array-like
        An array or list containing the results of rate rules to be assigned to the instance.
    lbls : array-like, optional
        An array or list of labels corresponding to the rate rules. If provided, these labels will
        be stored in the `rules_labels` attribute.

    Attributes
    ----------
    rules : array-like
        Stores the rate rules results.
    HAS_RULES : bool
        Flag to indicate that the rate rules have been set.
    rules_labels : array-like, optional
        Stores the labels corresponding to the rate rules, if provided.

    Examples
    --------
    >>> integration_data_obj = IntegrationDataObj()
    >>> rate_rules = [0.1, 0.2, 0.3]
    >>> rule_labels = ['rule1', 'rule2', 'rule3']
    >>> integration_data_obj.setRules(rate_rules, rule_labels)
    >>> integration_data_obj.rules
    [0.1, 0.2, 0.3]
    >>> integration_data_obj.rules_labels
    ['rule1', 'rule2', 'rule3']
    >>> integration_data_obj.HAS_RULES
    True
    """
        self.rules = rules
        self.HAS_RULES = True
        if lbls != None:
            self.rules_labels = lbls

    def setXData(self, xdata, lbls=None):
        """Sets extra simulation data.

    This method assigns additional simulation data provided in `xdata` to the
    `IntegrationDataObj` instance. The `HAS_XDATA` attribute is set to `True` 
    to indicate that extra data is present. Optionally, labels for the data 
    can be set using `lbls`.

    Parameters
    ----------
    xdata : array-like
        An array or list containing the additional simulation data.
    lbls : array-like, optional
        An optional array or list of labels corresponding to each element in
        `xdata`. If provided, these labels will be associated with the 
        respective data in `xdata`.

    Attributes
    ----------
    xdata : array-like
        Stores the provided extra simulation data.
    HAS_XDATA : bool
        A boolean flag that indicates whether extra data has been set.
    xdata_labels : array-like, optional
        Stores the provided labels for the extra simulation data if `lbls` is not `None`.

    Examples
    --------
    >>> obj = IntegrationDataObj()
    >>> extra_data = [0.5, 1.2, 3.4]
    >>> labels = ['param1', 'param2', 'param3']
    >>> obj.setXData(extra_data, labels)
    >>> print(obj.xdata)
    [0.5, 1.2, 3.4]
    >>> print(obj.xdata_labels)
    ['param1', 'param2', 'param3']
    >>> print(obj.HAS_XDATA)
    True
    """
        self.xdata = xdata
        self.HAS_XDATA = True
        if lbls != None:
            self.xdata_labels = lbls

    def getTime(self, lbls=False):
        """Return the time vector.

    Parameters
    ----------
    lbls : bool, optional
        If True, returns a tuple containing the time vector and its label. Default is False.

    Returns
    -------
    numpy.ndarray or tuple
        If `lbls` is False, returns a 1-D numpy array representing the time vector.
        If `lbls` is True, returns a tuple with the first element being the time vector and the second element
        a list containing the time label.

    Examples
    --------
    >>> integration_data = IntegrationDataObj()
    >>> integration_data.time = np.array([[0], [1], [2], [3]])
    >>> integration_data.time_label = 'Time (s)'
    >>> integration_data.HAS_TIME = True
    
    >>> # Retrieve time vector only
    >>> time_vector = integration_data.getTime()
    >>> print(time_vector)
    [0 1 2 3]

    >>> # Retrieve time vector with label
    >>> time_with_label = integration_data.getTime(lbls=True)
    >>> print(time_with_label)
    (array([0, 1, 2, 3]), ['Time (s)'])
    """
        output = None
        if self.HAS_TIME:
            output = self.time.reshape(len(self.time))
        if not lbls:
            return output
        else:
            return output, [self.time_label]

    def getSpecies(self, lbls=False):
        """Return a concatenated array of time and species data.

    This method returns an array that combines time and species data if species 
    data is available. If species data is not available, it returns only the time 
    data. Additionally, it can return the labels associated with the time and 
    species data.

    Parameters
    ----------
    lbls : bool, optional
        If True, also return the labels for the time and species data. 
        The default is False.

    Returns
    -------
    output : numpy.ndarray
        The concatenated time and species data array. If only time data is 
        available, returns only the time data array.
    labels : list of str, optional
        The labels for the time and species data. Returned only if `lbls` is True.

    Examples
    --------
    >>> # Assume `data_obj` is an instance of IntegrationDataObj with relevant data
    >>> time_species_data = data_obj.getSpecies()
    >>> print(time_species_data)
    array([...])

    >>> # To retrieve both data and labels
    >>> time_species_data, labels = data_obj.getSpecies(lbls=True)
    >>> print(time_species_data)
    array([...])
    >>> print(labels)
    ['time_label', 'species_label1', 'species_label2', ...]
    """
        output = None
        if self.HAS_SPECIES:
            output = numpy.hstack((self.time, self.species))
            labels = [self.time_label] + self.species_labels
        else:
            output = self.time
            labels = [self.time_label]
        if not lbls:
            return output
        else:
            return output, labels

    def getRates(self, lbls=False):
        """Return the time and rate array.

    This method retrieves an array that combines the time vector with the rate data. 
    If the `lbls` parameter is set to True, it also returns the labels for the time and rates.

    Parameters
    ----------
    lbls : bool, optional
        If True, returns a tuple containing the array and the list of labels. 
        If False, only returns the array. Default is False.

    Returns
    -------
    numpy.ndarray
        An array consisting of the time vector and rate data. If rate data is 
        unavailable, returns only the time vector.
    
    tuple of (numpy.ndarray, list of str)
        When `lbls` is True, returns a tuple with the array and a list of labels 
        for the time and rates. If rate data is unavailable, the label list will 
        contain only the time label.

    Examples
    --------
    >>> integration_data = IntegrationDataObj()
    >>> rates = integration_data.getRates()
    >>> rates_with_labels = integration_data.getRates(lbls=True)
    >>> print(rates_with_labels)
    (array([0. , 0.1, 0.2, ..., 1.0]), ['time', 'rate1', 'rate2', ...])
    """
        output = None
        if self.HAS_RATES:
            output = numpy.hstack((self.time, self.rates))
            labels = [self.time_label] + self.rate_labels
        else:
            output = self.time
            labels = [self.time_label]
        if not lbls:
            return output
        else:
            return output, labels

    def getRules(self, lbls=False):
        """
    Return a concatenated array of time and rule data.

    This method retrieves the time data and, if available, appends it with rule data to form a single array.
    Optionally, it can return both the data array and their corresponding labels.

    Parameters
    ----------
    lbls : bool, optional
        A flag indicating whether to return the labels alongside the data array. 
        If set to `True`, the method returns a tuple of (array, labels). 
        Default is `False`.

    Returns
    -------
    ndarray
        If `lbls` is `False`, returns an array combining time and rule data.
    tuple
        If `lbls` is `True`, returns a tuple containing:
        - ndarray: The combined array of time and rule data.
        - list: A list of labels for the data array.

    Notes
    -----
    The method will only append the rule data if the `HAS_RULES` attribute is set to `True`.

    Examples
    --------
    >>> data_obj = IntegrationDataObj()
    >>> time_rule_array = data_obj.getRules()
    >>> print(time_rule_array)
    array([...])  # Example output format, actual data may vary

    >>> time_rule_array, labels = data_obj.getRules(lbls=True)
    >>> print(labels)
    ['time', 'rule1', 'rule2', ...]  # Example label format, actual labels may vary
    >>> print(time_rule_array)
    array([...])  # Example output format, actual data may vary
    """
        output = None
        if self.HAS_RULES:
            output = numpy.hstack((self.time, self.rules))
            labels = [self.time_label] + self.rules_labels
        else:
            output = self.time
            labels = [self.time_label]
        if not lbls:
            return output
        else:
            return output, labels

    def getXData(self, lbls=False):
        """
    Return the time and xdata array.

    This method returns an array combining time data with xdata if available.
    If `lbls` is set to True, it also returns labels associated with the data.

    Parameters
    ----------
    lbls : bool, optional
        If True, the method will return a tuple containing both the data
        array and a list of labels, by default False.

    Returns
    -------
    numpy.ndarray
        A numpy array containing time and, if available, xdata.
    tuple of (numpy.ndarray, list of str)
        If `lbls` is True, returns a tuple where the first element is a numpy
        array containing time and xdata, and the second element is a list of
        labels corresponding to the data columns.

    Examples
    --------
    >>> data_obj = IntegrationDataObj()
    >>> data_obj.time = numpy.array([0, 1, 2])
    >>> data_obj.xdata = numpy.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]])
    >>> data_obj.time_label = "Time"
    >>> data_obj.xdata_labels = ["X1", "X2"]
    >>> data_obj.HAS_XDATA = True
    >>> data_obj.getXData()
    array([[0, 1, 2],
           [0.1, 0.2, 0.3],
           [0.4, 0.5, 0.6]])
    
    >>> data_obj.getXData(lbls=True)
    (array([[0, 1, 2],
            [0.1, 0.2, 0.3],
            [0.4, 0.5, 0.6]]),
     ['Time', 'X1', 'X2'])
    """
        output = None
        if self.HAS_XDATA:
            output = numpy.hstack((self.time, self.xdata))
            labels = [self.time_label] + self.xdata_labels
        else:
            output = self.time
            labels = [self.time_label]
        if not lbls:
            return output
        else:
            return output, labels

    def getDataAtTime(self, time):
        """
    Return all data generated at a specified time.

    This method retrieves all relevant datasets associated with a given
    time point from the stored data of the `IntegrationDataObj` object. 
    The returned array may include species, rates, rules, and x-data, 
    depending on what is available.

    Parameters
    ----------
    time : float
        The specific time point at which the data is to be retrieved.

    Returns
    -------
    numpy.ndarray or None
        An array where the first column contains the time point and subsequent
        columns contain available data such as species, rates, rules, or
        x-data. If no data exists for the specified time, returns None.

    Notes
    -----
    The method assumes that the `time` attribute of the object is a NumPy array 
    of times and that additional relevant datasets (species, rates, rules, 
    x-data) are stored within the object.

    Example
    -------
    >>> integration_data = IntegrationDataObj()
    >>> time_point = 5.0
    >>> data_at_time = integration_data.getDataAtTime(time_point)
    >>> print(data_at_time)
    array([[5.0, species_data, rate_data, rule_data, xdata]])

    In this example, `species_data`, `rate_data`, `rule_data`, and
    `xdata` are placeholders for the actual data associated with the 
    given time point.

    """
        t = None
        sp = None
        ra = None
        ru = None
        xd = None
        temp_t = self.time.reshape(len(self.time))
        for tt in range(len(temp_t)):
            if temp_t[tt] == time:
                t = tt
                if self.HAS_SPECIES:
                    sp = self.species.take([tt], axis=0)
                if self.HAS_RATES:
                    ra = self.rates.take([tt], axis=0)
                if self.HAS_RULES:
                    ru = self.rules.take([tt], axis=0)
                if self.HAS_XDATA:
                    xd = self.xdata.take([tt], axis=0)
                break
        output = None
        if t is not None:
            output = numpy.array([[temp_t[t]]])
            if sp is not None:
                output = numpy.hstack((output, sp))
            if ra is not None:
                output = numpy.hstack((output, ra))
            if ru is not None:
                output = numpy.hstack((output, ru))
            if xd is not None:
                output = numpy.hstack((output, xd))
        return output

    def getDataInTimeInterval(self, time, bounds=None):
        """
    Returns an array of all data points where `time-bounds <= time <= time+bounds`.

    Parameters
    ----------
    time : float
        The target time around which to retrieve data points.
    bounds : float, optional
        The time interval bounds, defaulting to the step size between the first
        two time points in the data.

    Returns
    -------
    numpy.ndarray
        An array containing all data points within the specified time interval. 
        The array includes the time points and, if available, species, rates,
        rules, and any additional data.

    Notes
    -----
    The method finds and returns all data points that fall within the
    time interval `[time - bounds, time + bounds]`. If `bounds` is not provided,
    it defaults to the difference between the first two recorded time points.

    Examples
    --------
    Assume an instance of `IntegrationDataObj` is available as `data_obj`.

    >>> data_obj = IntegrationDataObj()
    >>> time_point = 5.0
    >>> bounds = 1.0
    >>> data_interval = data_obj.getDataInTimeInterval(time_point, bounds)
    >>> print(data_interval)
    array([[4.0, species_1, species_2, ..., rate_1, rule_1, ...],
           [5.0, species_1, species_2, ..., rate_1, rule_1, ...],
           [6.0, species_1, species_2, ..., rate_1, rule_1, ...]])

    In this example, the returned data includes all data points recorded at 
    times 4, 5, and 6, with respective species, rates, and rules data included.
    """
        temp_t = self.time.reshape(len(self.time))
        if bounds is None:
            bounds = temp_t[1] - temp_t[0]
        c1 = temp_t >= time - bounds
        c2 = temp_t <= time + bounds
        print('Searching ({}:{}:{})'.format(time - bounds, time, time + bounds)
            )
        t = []
        sp = None
        ra = None
        for tt in range(len(c1)):
            if c1[tt] and c2[tt]:
                t.append(tt)
        output = None
        if len(t) > 0:
            output = self.time.take(t)
            output = output.reshape(len(output), 1)
            if self.HAS_SPECIES and self.HAS_TIME:
                output = numpy.hstack((output, self.species.take(t, axis=0)))
            if self.HAS_RATES:
                output = numpy.hstack((output, self.rates.take(t, axis=0)))
            if self.HAS_RULES:
                output = numpy.hstack((output, self.rules.take(t, axis=0)))
            if self.HAS_XDATA:
                output = numpy.hstack((output, self.xdata.take(t, axis=0)))
        return output

    def getOutput(self, *args, **kwargs):
        """
    Returns simulation data for specified species or rate labels.

    This method is an alias for `getSimData()`. By providing species or
    rate labels as arguments, it returns an array composed of time and
    the values of the specified species and rates, in the form
    `[time, sp1, r1, ...]`.

    Parameters
    ----------
    *args : list
        Species or rate labels to include in the output array.
    **kwargs : dict
        Additional keyword arguments to pass to the `getSimData()` method.

    Returns
    -------
    numpy.ndarray
        An array containing the time, followed by the requested species 
        and/or rates for each time point.

    See Also
    --------
    getSimData : Core function that retrieves simulation data.

    Notes
    -----
    This method serves as a backward-compatible alias to `getSimData()`
    and is maintained for legacy support.

    Examples
    --------
    >>> # Assuming an instance of IntegrationDataObj named `data_obj`:
    >>> species_labels = ['species1', 'species2']
    >>> rate_labels = ['rate1']
    >>> output_data = data_obj.getOutput(*species_labels, *rate_labels)
    >>> print(output_data)
    array([[time1, sp1_value1, sp2_value1, r1_value1],
           [time2, sp1_value2, sp2_value2, r1_value2],
           ...])
    """
        return self.getSimData(*args, **kwargs)

    def getAllSimData(self, lbls=False):
        """
    Returns all available data as a combination of time, species, rates, and rules.

    This method aggregates all simulation data into a single array. By default, it 
    returns the data array. If the `lbls` parameter is set to `True`, it returns a 
    tuple containing both the array and the associated labels for each column in the array.

    Parameters
    ----------
    lbls : bool, optional
        If `True`, the function returns a tuple containing the data array and the labels. 
        If `False` (default), it returns only the data array.

    Returns
    -------
    output : numpy.ndarray
        An array containing all available simulation data, including time, species, rates, 
        rules, and optionally other data.
    
    labels : list of str, optional
        A list of labels corresponding to each column in the output array. Only returned 
        if `lbls` is set to `True`.

    Examples
    --------
    Assuming `obj` is an instance of `IntegrationDataObj`:

    >>> data = obj.getAllSimData()
    >>> data.shape  # Returns the shape of the data array

    To include labels:

    >>> data, labels = obj.getAllSimData(lbls=True)
    >>> labels  # Returns the list of labels such as ['time', 'species1', 'rate1', ...]
    """
        labels = [self.time_label]
        if self.HAS_SPECIES and self.HAS_TIME:
            output = numpy.hstack((self.time, self.species))
            labels += self.species_labels
        if self.HAS_RATES:
            output = numpy.hstack((output, self.rates))
            labels += self.rate_labels
        if self.HAS_RULES:
            output = numpy.hstack((output, self.rules))
            labels += self.rules_labels
        if self.HAS_XDATA:
            output = numpy.hstack((output, self.xdata))
            labels += self.xdata_labels
        if not lbls:
            return output
        else:
            return output, labels

    def getSimData(self, *args, **kwargs):
        """
    Retrieve simulation data for specified species and rate labels.

    This method allows you to specify the species and rate labels for 
    which you want to retrieve the simulation data. The returned output 
    is an array where the first column is the time and the subsequent 
    columns correspond to the specified species and rates.

    Parameters
    ----------
    *args : str
        Species and rate labels for which data is to be retrieved.
    **kwargs : dict, optional
        Additional arguments, where 'lbls' key can be set to `True` 
        to include the labels in the returned output.

    Returns
    -------
    numpy.ndarray
        Concatenated array of simulation data. The array will have 
        columns representing time and the specified species and rates.
    tuple of (numpy.ndarray, list of str)
        If `lbls=True` is provided as a keyword argument, returns a 
        tuple where the first element is the concatenated array and 
        the second element is the list of labels corresponding to 
        the columns in the array.

    Example
    -------
    >>> integration_data = IntegrationDataObj()
    >>> time_species_rate_data = integration_data.getSimData('species1', 'rate1')
    >>> print(time_species_rate_data)
    array([[ 0. ,  1.2,  0.5],
           [ 0.1,  1.5,  0.7],
           ...
           [10. ,  3.2,  2.1]])
    
    >>> time_species_rate_data, labels = integration_data.getSimData('species1', 
                                                                     'rate1', 
                                                                     lbls=True)
    >>> print(labels)
    ['time', 'species1', 'rate1']
    >>> print(time_species_rate_data)
    array([[ 0. ,  1.2,  0.5],
           [ 0.1,  1.5,  0.7],
           ...
           [10. ,  3.2,  2.1]])
    """
        output = self.time
        if 'lbls' in kwargs:
            lbls = kwargs['lbls']
        else:
            lbls = False
        lout = [self.time_label]
        for roc in args:
            if self.HAS_SPECIES and roc in self.species_labels:
                lout.append(roc)
                output = numpy.hstack((output, self.species.take([self.
                    species_labels.index(roc)], axis=-1)))
            if self.HAS_RATES and roc in self.rate_labels:
                lout.append(roc)
                output = numpy.hstack((output, self.rates.take([self.
                    rate_labels.index(roc)], axis=-1)))
            if self.HAS_RULES and roc in self.rules_labels:
                lout.append(roc)
                output = numpy.hstack((output, self.rules.take([self.
                    rules_labels.index(roc)], axis=-1)))
            if self.HAS_XDATA and roc in self.xdata_labels:
                lout.append(roc)
                output = numpy.hstack((output, self.xdata.take([self.
                    xdata_labels.index(roc)], axis=-1)))
        if not lbls:
            return output
        else:
            return output, lout


class ReactionObj(NewCoreBase):
    """
    Represents a reaction within a model, defined by a kinetic law, formula, and name. Each 
    reaction is bound to a specific model instance.

    Attributes
    ----------
    formula : str
        The mathematical representation of the reaction formula.
    mod : object
        The model instance to which this reaction is bound.
    rate : float
        The computed rate of the reaction after evaluating the kinetic law.
    xcode : code object
        The compiled code object that, when executed, computes the reaction's rate.
    code_string : str
        The string representation of the Python code that calculates the reaction rate.
    symbols : list of str
        List of symbols involved in the reaction's kinetic law.
    compartment : optional
        The compartment (if any) where the reaction takes place.
    piecewises : list
        List of piecewise function definitions parsed from the kinetic law.

    Methods
    -------
    __call__(*args)
        Executes the reaction's compiled code (xcode) to compute and return the rate.
        
    setKineticLaw(kl, klrepl='self.')
        Parses the kinetic law string, compiles it, and prepares it for execution.
    """
    formula = None
    mod = None
    rate = None
    xcode = None
    code_string = None
    symbols = None
    compartment = None
    piecewises = None

    def __init__(self, mod, name, kl, klrepl='self.'):
        """
    Initialize a ReactionObj instance.

    Parameters
    ----------
    mod : object
        The model instance this reaction belongs to.
    name : str
        The name of the reaction.
    kl : str
        The kinetic law associated with the reaction.
    klrepl : str, optional
        The replacement string for the KineticLaw, by default 'self.'

    Attributes
    ----------
    mod : object
        The model instance associated with this reaction.
    name : str
        The name of the reaction.
    
    Methods
    -------
    setKineticLaw(kl, klrepl='self.')
        Sets the kinetic law for the reaction.
    
    Example
    -------
    >>> model_instance = SomeModelClass()
    >>> reaction_name = "Reaction1"
    >>> kinetic_law = "k1 * A * B"
    >>> reaction = ReactionObj(model_instance, reaction_name, kinetic_law)
    >>> print(reaction.name)
    Reaction1
    >>> print(reaction.mod)
    <SomeModelClass instance>
    """
        self.mod = mod
        self.name = name
        self.setKineticLaw(kl, klrepl='self.')

    def __call__(self, *args):
        """
    Execute the pre-defined reaction code and return the rate of the reaction.

    This method allows an instance of `ReactionObj` to be called as if it were a function. 
    When called, it executes a block of code defined by `self.xcode` and returns the 
    computed rate of the reaction.

    Parameters
    ----------
    *args : tuple
        Arguments passed to the method. These arguments are not used in the 
        current implementation, but can be utilized if required modifications are 
        made to the `self.xcode` execution process.

    Returns
    -------
    float
        The calculated rate of the reaction as determined by `self.xcode`.

    Example
    -------
    >>> reaction = ReactionObj(mod, "reaction1", kl)
    >>> reaction.xcode = "self.rate = 0.1 * 42"  # Example xcode to define the rate
    >>> rate = reaction()
    >>> print(rate)
    4.2
    """
        exec(self.xcode)
        return self.rate

    def setKineticLaw(self, kl, klrepl='self.'):
        """
    Parses and sets the kinetic law for the reaction, compiling it into executable code.

    Parameters
    ----------
    kl : str
        A string representing the kinetic law of the reaction.
    klrepl : str, optional
        A string to replace within the kinetic law, default is 'self.'.

    Notes
    -----
    This method utilizes an `InfixParser` to parse the kinetic law string, replacing
    specified parts with an executable format, and compiles it into Python bytecode.
    The method also updates internal attributes that represent the symbols and 
    piecewise functions involved in the formula.

    The kinetic law string, after processing, is stored in the `self.formula` attribute,
    and the compiled code is stored in `self.xcode`.

    Example
    -------
    >>> reaction = ReactionObj(mod=model_instance, name="reaction1", kl="rate*self.concentration")
    >>> reaction.setKineticLaw(kl="rate*selfA", klrepl="self.")
    >>> print(reaction.formula)
    'rate*selfA'
    """
        InfixParser.setNameStr('self.mod.', '')
        InfixParser.parse(kl.replace(klrepl, ''))
        self.symbols = InfixParser.names
        self.piecewises = InfixParser.piecewises
        formula = InfixParser.output
        for pw in InfixParser.piecewises:
            formula = formula.replace('self.mod.{}'.format(pw),
                'self.mod.{}()'.format(pw))
        self.code_string = 'self.rate={}'.format(formula)
        self.xcode = compile(self.code_string, 'Req: {}'.format(self.name),
            'exec')
        self.formula = kl.replace(klrepl, '').replace('numpy.', '').replace(
            'operator.', '')


InfixParser = MyInfixParser()
InfixParser.buildlexer()
InfixParser.buildparser(debug=0, debugfile='infix.dbg', tabmodule=
    'infix_tabmodule', outputdir=OUTPUT_DIR)
InfixParser.setNameStr('self.', '')
os.chdir(CWD)


class Function(NewCoreBase):
    """
    Function class ported from Core2 to enable the use of functions in PySCeS.

    This class allows the representation and execution of mathematical functions
    within the PySCeS framework. It translates a provided formula into executable
    Python code using an infix parser, enabling dynamic evaluation of functions 
    with symbolic arguments.

    Attributes
    ----------
    formula : str
        The mathematical formula representing the function.
    code_string : str
        The generated Python code that evaluates the function.
    xcode : code
        The compiled code object for the function evaluation.
    value : float
        The result of the function after execution.
    symbols : list of str
        The symbols (variables) used in the formula.
    argsl : list
        The list of arguments used within the function.
    mod : object
        The module within which the function is contextualized.
    piecewises : list of str
        Piecewise functions used within the formula, if any.
    functions : list of str
        The functions detected within the formula.
    
    Notes
    -----
    The class simplifies the process of defining functions that can be executed 
    as part of larger modeling activities using PySCeS. It leverages `InfixParser` 
    to parse formulae and compile them into executable form.

    See Also
    --------
    InfixParser : The parser used for converting the formula into executable code.
    NewCoreBase : The base class from which `Function` inherits.
    """
    formula = None
    code_string = None
    xcode = None
    value = None
    symbols = None
    argsl = None
    mod = None
    piecewises = None
    functions = None

    def __init__(self, name, mod):
        """
Initialize a new Function instance.

This constructor sets up a Function with a given name and module, preparing
it for further configuration and use.

Parameters
----------
name : str
    The name of the function.
mod : Any
    The module or context in which the function operates, allowing the function to access various attributes or methods specific to this context.

Attributes
----------
name : str
    Stores the name of the function for identification and operational purposes.
argsl : list
    An empty list initialized to hold argument specifications or details needed for the function's operation.
functions : list
    An empty list initialized to potentially store other functions or operations associated with this function.
mod : Any
    The module or context associated with the function, offering access to context-specific methods or attributes.

Examples
--------
>>> from some_module import Function
>>> function_instance = Function(name='Transformation', mod=my_module)
>>> print(function_instance.name)
'Transformation'
>>> print(function_instance.mod)
<module 'my_module' from '/path/to/my_module.py'>
"""
        self.name = name
        self.argsl = []
        self.functions = []
        self.mod = mod

    def __call__(self, *args):
        """
    Executes the compiled code with the given arguments and returns the result.

    The `__call__` method allows instances of the `Function` class to be used as callable objects.
    It sets the given arguments to the attributes in `argsl` and executes the pre-compiled expression
    stored in `xcode`. The result of the execution is returned.

    Parameters
    ----------
    *args : tuple
        Variable length argument list. Each element corresponds to an argument expected by the function.
        The number of arguments provided should match the number of expected arguments.

    Returns
    -------
    result : variable type
        The result of executing the pre-compiled expression for the given arguments.

    Example
    -------
    Here's an example of how you might use the `__call__` method:

    >>> my_func = Function('my_function', mod)
    >>> my_func.argsl = ['arg1', 'arg2']
    >>> my_func.xcode = compile('self.value = arg1 + arg2', '<string>', 'exec')
    >>> result = my_func(5, 3)
    >>> print(result)
    8

    In this example, `my_func` is a callable object where the function 'arg1 + arg2' is executed
    with the provided arguments `5` and `3`, resulting in an output of `8`.
    """
        for ar in range(len(args)):
            self.__setattr__(self.argsl[ar], args[ar])
        exec(self.xcode)
        return self.value

    def setArg(self, var, value=None):
        """
    Sets an argument for the function by assigning a specified value and tracks the argument name.

    Parameters
    ----------
    var : str
        The name of the argument to be set for the function.
    value : optional
        The value to assign to the argument. If not provided, the default is `None`.

    Examples
    --------
    >>> f = Function(name="example_function", mod="example_module")
    >>> f.setArg("x", 42)
    >>> f.setArg("y")
    >>> f.argsl
    ['x', 'y']
    >>> f.x
    42
    >>> f.y
    None
    """
        self.__setattr__(var, value)
        self.argsl.append(var)

    def addFormula(self, formula):
        """
    Adds a mathematical formula to the function.

    This method processes the input formula using the `InfixParser`, 
    converting it into executable code. It updates the function's 
    internal representation of piecewise expressions, symbols, and functions 
    based on the parsed formula, and sets up a code string for evaluation.

    Parameters
    ----------
    formula : str
        A string representing the mathematical formula to be added. The 
        formula can include standard mathematical symbols and notations.

    Notes
    -----
    This method uses `InfixParser` to interpret the formula. It replaces 
    global placeholders such as `_TIME_` with specific module references 
    before parsing the formula. The resulting computations are stored as 
    executable code string, ready for execution when the object is called.

    Example
    -------
    >>> func = Function()
    >>> func.name = "example_func"
    >>> formula = "x + y + _TIME_"
    >>> func.addFormula(formula)
    >>> print(func.code_string)
    'self.value=<parsed_output_from_InfixParser>'

    The `code_string` generated will assign the calculation outcome to 
    `self.value`, using the parsed result of the formula.

    """
        self.formula = formula
        InfixParser.setNameStr('self.', '')
        InfixParser.SymbolReplacements = {'_TIME_': 'mod._TIME_'}
        InfixParser.parse(formula)
        self.piecewises = InfixParser.piecewises
        self.symbols = InfixParser.names
        self.functions = InfixParser.functions
        self.code_string = 'self.value={}'.format(InfixParser.output)
        self.xcode = compile(self.code_string, 'Func: {}'.format(self.name),
            'exec')


class EventAssignment(NumberBase):
    """
    EventAssignment class serves as a mechanism for handling event-triggered assignments
    within the PySCeS event handling framework. It facilitates the evaluation and execution
    of assignments when corresponding events occur.

    Attributes
    ----------
    variable : str or None
        The name of the variable to which the formula result is assigned.
    symbols : list or None
        List of symbols used within the formula.
    formula : str or None
        The string representation of the formula to be executed when the event is triggered.
    code_string : str or None
        The string representation of the Python code that computes the formula result.
    xcode : code object or None
        A compiled code object for the formula that can be executed within Python.
    mod : object
        A reference to the module or component associated with the event assignment.
    piecewises : list or None
        A list of piecewise expressions derived from the formula, if any.
    __DEBUG__ : bool
        A flag for enabling or disabling debug output during event assignment.

    Methods
    -------
    __call__():
        Updates the module variable with the computed value and optionally prints debug information.
    __init__(name, mod):
        Initializes an EventAssignment instance with a given name and module reference.
    setVariable(var):
        Sets the target variable for the assignment.
    setFormula(formula):
        Parses and compiles the formula into executable code.
    evaluateAssignment():
        Executes the compiled code to perform the assignment.

    Notes
    -----
    This class uses an `InfixParser` to parse mathematical formulae. The parser replaces certain
    symbols and compiles the string representation into runnable Python code for dynamic execution.

    """
    variable = None
    symbols = None
    formula = None
    code_string = None
    xcode = None
    mod = None
    piecewises = None
    __DEBUG__ = False

    def __call__(self):
        """
    Represents an assignment operation that can be executed on a model object.
    
    The `__call__` method is used to update the attribute `variable` of the object `mod` with the value `value`.
    If debugging is enabled, the operation is printed to the console.
    
        Executes the event assignment by setting an attribute on the model object.

        This method assigns the specified `value` to the `variable` attribute of the
        `mod` object. If the instance variable `__DEBUG__` is set to True, it prints
        the assignment operation for debugging purposes.

        Returns
        -------
        bool
            Returns True after performing the assignment operation.

        Example
        -------
        >>> class Model:
        ...     pass
        ...
        >>> mod = Model()
        >>> assignment = EventAssignment()
        >>> assignment.mod = mod
        >>> assignment.variable = 'attribute_name'
        >>> assignment.value = 42
        >>> assignment.__DEBUG__ = True
        >>> assignment()
        Assigning attribute_name = 42
        True
        >>> mod.attribute_name
        42

        Notes
        -----
        - Ensure that the `mod` object has an attribute named `variable` before executing.
        - The debugging output shows the variable name and assigned value when debugging is enabled.
        """
        setattr(self.mod, self.variable, self.value)
        if self.__DEBUG__:
            print('\tAssigning {} = {}'.format(self.variable, self.value))
        return True

    def __init__(self, name, mod):
        """
    Initialize an EventAssignment instance.

    This constructor method initializes an `EventAssignment` object with a
    specific name and module. The name is typically used to identify the
    assignment operation, while the module (`mod`) represents the target
    module where the assignment will take place.

    Parameters
    ----------
    name : str
        A string representing the name of the event assignment. This is often
        used as an identifier within the context where this instance is being
        used.
    mod : object
        An object that represents the module or context in which the event
        assignment operation is to occur. This should be the object that
        holds the variable to be modified by this event assignment.

    Examples
    --------
    >>> class Module:
    ...     pass
    ...
    >>> mod_instance = Module()
    >>> event_assignment = EventAssignment(name='SampleEvent', mod=mod_instance)
    """
        self.setName(name)
        self.mod = mod

    def setVariable(self, var):
        """
    Sets the variable for the event assignment.

    This method assigns a new variable to the event assignment instance, which 
    will be used when the event is triggered.

    Parameters
    ----------
    var : str
        The name of the variable to be set for the event assignment. This is 
        typically a string representing a variable name within the module 
        `mod`.

    Examples
    --------
    >>> event_assignment = EventAssignment("event1", some_module)
    >>> event_assignment.setVariable("temperature")
    >>> print(event_assignment.variable)
    temperature
    """
        self.variable = var

    def setFormula(self, formula):
        """
    Set the formula for the event assignment and update internal state.

    The method takes a formula string, processes it using the `InfixParser`,
    and updates the internal state, including piecewise structures, symbol names,
    and a code string that can be executed to evaluate the formula.

    Parameters
    ----------
    formula : str
        A string representing the mathematical formula to be used for the event assignment.

    Examples
    --------
    >>> from event_assignment import EventAssignment
    >>> event = EventAssignment("SampleEvent", mod_instance)
    >>> event.setFormula("a * b + c / d")
    
    After setting the formula, the internal state is updated to include parsed
    piecewise expressions, symbol replacements, and a compiled code object 
    ready for execution.
    """
        self.formula = formula
        InfixParser.setNameStr('self.mod.', '')
        InfixParser.SymbolReplacements = {'_TIME_': 'self.mod._TIME_'}
        InfixParser.parse(formula)
        self.piecewises = InfixParser.piecewises
        self.symbols = InfixParser.names
        self.code_string = 'self.value={}'.format(InfixParser.output)
        self.xcode = compile(self.code_string, 'EvAs: {}'.format(self.name),
            'exec')

    def evaluateAssignment(self):
        """
    Executes the compiled code for the event assignment evaluation.

    This method performs the evaluation of the event assignment using 
    previously compiled Python code (stored in `self.xcode`). The execution 
    of this code determines the value of the event assignment. It ensures 
    the event's conditions are processed as expected, allowing for the 
    dynamic adjustment of model variables during simulation.

    Examples
    --------
    >>> # Assuming 'event_assignment' is an instance of EventAssignment,
    >>> # and 'setFormula' has been called to define 'self.xcode'.
    >>> event_assignment.evaluateAssignment()
    >>> # This execution sets the evaluated output to the instance's
    >>> # 'self.value' attribute, applying the transformations or calculations
    >>> # defined by the provided formula.
    """
        exec(self.xcode)


class Event(NewCoreBase):
    """
    Represents an event with a trigger condition that, when met, executes 
    a set of event assignments. This class handles the evaluation of the 
    trigger condition, the application of delays, and the persistence of 
    events. Events are ported from the Core2 system and include debug 
    capabilities.

    Attributes
    ----------
    trigger : str, optional
        The condition that triggers the event. Initialized during trigger setup.
    delay : float
        The delay time after the trigger condition is met before the event 
        assignments are executed. Default is 0.0.
    formula : str, optional
        The string representation of the trigger formula.
    code_string : str, optional
        The compiled string of the trigger's logic, ready for execution.
    xcode : code, optional
        The compiled bytecode representation of the trigger logic.
    state0 : bool
        Previous state of the event; checks if the event was triggered.
    state : bool
        Current state of the event; checks if the event is active.
    assignments : list of EventAssignment
        A list of event assignments that are executed when the event is triggered.
    _TIME_ : float
        The current simulation time.
    _ASS_TIME_ : float
        The scheduled execution time for the assignments.
    _need_action : bool
        Flag indicating whether the event needs action (i.e., assignments to be
        executed) at the current time.
    _assign_now : bool
        Flag indicating whether assignments should be performed immediately.
    symbols : list, optional
        A list of symbol names used in the trigger formula.
    _time_symbol : str, optional
        Custom symbol used for replacing time-related variables in the formula.
    piecewises : list, optional
        Parsed piecewise functions from the formula.
    mod : Module
        The module associated with this event providing context and dependencies.
    priority : any
        Priority value for determining event execution order.
    persistent : bool
        Flag indicating whether the event will remain active once triggered.
    __DEBUG__ : bool
        Enables or disables debug printing inside the class methods.

    Methods
    -------
    __init__(name, mod)
        Initializes the event with a given name and associated module.
    __call__(time)
        Evaluates the trigger at the given time and checks for actions.
    setTrigger(formula, delay=0.0)
        Sets the trigger condition and delay for the event.
    setAssignment(var, formula)
        Associates an assignment formula with a given variable.
    setPriority(priority)
        Sets the priority of the event.
    setPersistent(persistent)
        Sets the persistence of the event.
    reset()
        Resets the state of the event to its initial conditions, clearing the 
        effect of previous triggers.
    """
    trigger = None
    delay = 0.0
    formula = None
    code_string = None
    xcode = None
    state0 = False
    state = False
    assignments = None
    _TIME_ = 0.0
    _ASS_TIME_ = 0.0
    _need_action = False
    _assign_now = False
    symbols = None
    _time_symbol = None
    piecewises = None
    mod = None
    priority = None
    persistent = True
    __DEBUG__ = True

    def __init__(self, name, mod):
        """
    Initialize an Event instance.

    This method sets the name of the event and initializes its
    assignments list. It also sets the `mod` attribute for the event.

    Parameters
    ----------
    name : str
        The name of the event to be set.
    
    mod : object
        The module or model associated with the event.

    Examples
    --------
    >>> model = SomeModel()
    >>> event = Event(name='start', mod=model)
    >>> print(event.name)
    start
    >>> print(event.mod)
    <SomeModel object at 0x...>

    Notes
    -----
    The method `setName` must be implemented to handle the setting
    of the event's name.
    """
        self.setName(name)
        self.assignments = []
        self.mod = mod

    def __call__(self, time):
        """
    Evaluates the event at a given time and determines whether actions need
    to be executed based on the current state and whether certain conditions 
    have been met according to the event's logic.

    Parameters
    ----------
    time : float
        The current time at which the event is being evaluated.

    Returns
    -------
    bool
        Returns `True` if the event has been triggered and either needs to 
        take action or has completed the required action. Otherwise, 
        returns `False`.

    Notes
    -----
    The method works by executing embedded code and checking the state 
    transitions to determine if any actions need to be performed. It uses 
    attributes like self.state, self.state0, self.delay, and others to control 
    the logic flow.

    Attributes
    ----------
    _TIME_ : float
        The time at which the event is considered for execution.
    _ASS_TIME_ : float
        The computed time after which the event's assigned action should be 
        performed.
    _need_action : bool
        Boolean flag indicating whether the event requires action to be taken.
    _assign_now : bool
        Boolean flag indicating whether the event's actions need to be 
        executed immediately.
    persistent : bool
        Flag indicating whether the event remains active beyond the initial 
        trigger.
    __DEBUG__ : bool
        A debug flag to control the verbosity of print statements for debugging 
        output.

    Example
    -------
    >>> event = Event(...)  # Initialize with appropriate parameters
    >>> current_time = 10.0
    >>> if event(current_time):
    ...     print("Event actions are being executed.")

    In this example, if the Event condition is met at the time of 10.0, the 
    event is executed, and a message is printed. The implementation assumes 
    the existence of the attributes like `self.state`, `self.state0`, etc., 
    which should be defined elsewhere within the `Event` class.
    """
        self._TIME_ = time
        exec(self.xcode)
        ret = False
        if self.state0 and not self.state:
            self.state0 = self.state
        if not self.state0 and self.state:
            self.state0 = self.state
            self._need_action = True
            self._ASS_TIME_ = time + self.delay
            if self.__DEBUG__:
                print('\nevent {} is evaluating at {}'.format(self.name, time))
            ret = True
        if self._need_action and self._TIME_ >= self._ASS_TIME_:
            if not self.persistent and not self.state:
                self._assign_now = False
                if self.__DEBUG__:
                    print('event {} NOT assigning after delay...'.format(
                        self.name))
                    print(
                        '    ...non-persistent event, trigger ({}) no longer valid'
                        .format(self.formula))
            else:
                self._assign_now = True
                if self.__DEBUG__:
                    print('event {} is assigning at {} (delay={})'.format(
                        self.name, time, self.delay))
            self._need_action = False
            ret = True
        return ret

    def setTrigger(self, formula, delay=0.0):
        """
    Sets the trigger condition and delay for the event.

    This method configures the event's condition, defined by the `formula`, 
    and an optional delay `delay`. It parses the formula and sets up the 
    internal state for event evaluation.

    Parameters
    ----------
    formula : str
        The logical expression that defines when the event should be triggered.
    delay : float, optional
        The delay in time units after the trigger condition is met 
        before the event executes, by default 0.0.

    Example
    -------
    >>> event = Event('Sample Event', 'mod')
    >>> event.setTrigger('x > 10', delay=5.0)
    >>> print(event.formula)
    x > 10
    >>> print(event.delay)
    5.0

    Notes
    -----
    This method utilizes the `InfixParser` for parsing the provided formula. 
    The parsed result is used to dynamically compile Python code that is 
    executed to determine the event's state.
    """
        self.formula = formula
        self.delay = delay
        InfixParser.setNameStr('self.mod.', '')
        if self._time_symbol != None:
            InfixParser.SymbolReplacements = {self._time_symbol: '_TIME_'}
        else:
            InfixParser.SymbolReplacements = {'_TIME_': '_TIME_'}
        InfixParser.parse(formula)
        self.piecewises = InfixParser.piecewises
        self.symbols = InfixParser.names
        self.code_string = 'self.state={}'.format(InfixParser.output)
        self.xcode = compile(self.code_string, 'Ev: {}'.format(self.name),
            'exec')

    def setAssignment(self, var, formula):
        """
    Sets up an event assignment for a specific variable with a given formula.

    This method creates an `EventAssignment` object for a specified variable 
    and formula, adds the assignment to the event's list of assignments, and 
    binds it to the event instance.

    Parameters
    ----------
    var : str
        The variable to be assigned when the event is triggered.
    formula : str
        The formula represented as a string that computes the value to be 
        assigned to `var` when the event occurs.

    Examples
    --------
    >>> event = Event(mod=model)
    >>> event.setAssignment('concentration', 'initial_concentration * exp(-k * time)')
    >>> # This will create an assignment for 'concentration' that, when the 
    >>> # event is triggered, assigns it a value based on its initial value,
    >>> # a decay rate `k`, and the elapsed `time`.
    """
        ass = EventAssignment(var, mod=self.mod)
        ass.setVariable(var)
        ass.setFormula(formula)
        self.assignments.append(ass)
        self.__setattr__('_' + var, ass)

    def setPriority(self, priority):
        """
Set the priority level of the event.

This method assigns a priority level to the event, where lower numbers
indicate higher priority. The priority level can be used to determine
the order in which events should be processed, with the most prioritized
events being handled first.

Parameters
----------
priority : int
    An integer representing the priority level of the event. Lower numbers
    correspond to higher priority.

Examples
--------
>>> event = Event()
>>> event.setPriority(1)
>>> print(event.priority)
1

In the above example, an Event instance is created, and its priority is
set to 1, indicating a high priority level.
"""
        self.priority = priority

    def setPersistent(self, persistent):
        """
    Sets the persistent attribute of the event.

    This method determines whether the event's triggered state should remain
    true after the triggering condition becomes false, allowing the event to 
    continue acting or to return to an untriggered state.

    Parameters
    ----------
    persistent : bool
        Specifies whether the event is persistent. If `True`, the event remains
        triggered even after the trigger condition no longer holds. If `False`,
        the event will become untriggered once the condition is not met.

    Examples
    --------
    Create an event and set it to be persistent:

    >>> event = Event()
    >>> event.setPersistent(True)
    >>> print(event.persistent)
    True

    Set the event to be non-persistent:

    >>> event.setPersistent(False)
    >>> print(event.persistent)
    False
    """
        self.persistent = persistent

    def reset(self):
        """
    Resets the state of the Event object to its initial conditions.

    This method sets the `state0` and `state` attributes to `False`, and 
    resets the `_TIME_` and `_ASS_TIME_` attributes to `0.0`. It is useful 
    for reinitializing the event's state in scenarios that require restarting 
    the event's lifecycle.

    Examples
    --------
    >>> event = Event()
    >>> event.state0 = True
    >>> event.state = True
    >>> event._TIME_ = 10.0
    >>> event._ASS_TIME_ = 5.0
    >>> event.reset()
    >>> event.state0
    False
    >>> event.state
    False
    >>> event._TIME_
    0.0
    >>> event._ASS_TIME_
    0.0
    """
        self.state0 = False
        self.state = False
        self._TIME_ = 0.0
        self._ASS_TIME_ = 0.0


class PieceWise(NewCoreBase):
    """
    A generic piecewise class adapted from Core2, designed to generate
    a compiled Python code block for evaluating piecewise functions of
    arbitrary length. This class enables the specification of piecewise
    statements within assignment rules using a syntax such as:
    `piecewise(<Piece>, <Conditional>, <OtherValue>)`, where there can
    be an arbitrary number of `<Piece>, <Conditional>` pairs.

    Attributes
    ----------
    name : any
        A placeholder attribute for storing variable name or relevant data.
    value : any
        The evaluated result of the piecewise function.
    formula : str
        The string representation of the piecewise logic used in evaluation.
    code_string : str
        The actual Python code representing the piecewise logic.
    xcode : code
        The compiled Python code ready for execution.
    _names : list
        The list of variable names involved in the piecewise function.
    _TIME_ : any
        An attribute typically used for time-based evaluation.

    Parameters
    ----------
    pwd : dict
        A dictionary of piecewise information parsed by the InfixParser,
        containing pieces and their corresponding conditions.
    mod : object
        A module or object context to resolve variable references within
        the piecewise function.

    Notes
    -----
    Piecewise functions are used extensively in mathematical modeling to handle
    functions that have different expressions based on certain conditions. This
    class automates the translation of these logical expressions into executable
    Python code that can be evaluated efficiently at runtime.

    See Also
    --------
    InfixParser : A parser utilized to interpret and transform piecewise logic
                  into executable Python code.

    """
    name = None
    value = None
    formula = None
    code_string = None
    xcode = None
    _names = None
    _TIME_ = None

    def __init__(self, pwd, mod):
        """"""
        pwd = pwd.copy()
        self.mod = mod
        if pwd['other'] != None:
            other = 'self.value = {}'.format(pwd.pop('other').replace(
                'self.', 'self.mod.'))
        else:
            other = 'pass'
            pwd.pop('other')
        InfixParser.setNameStr('self.mod.', '')
        self._names = []
        formula = pwd[0][0]
        InfixParser.parse(formula)
        for n in InfixParser.names:
            if n not in self._names and n != '_TIME_':
                self._names.append(n)
        formula = InfixParser.output
        thenStat = pwd[0][1].replace('self.', 'self.mod.')
        if len(list(pwd.keys())) == 1:
            self.code_string = ('if {}:\n    self.value = {}\nelse:\n    {}'
                .format(formula, thenStat, other))
            self.formula = self.code_string.replace('self.', '')
        else:
            self.code_string = 'if {}:\n    self.value = {}\n'.format(formula,
                thenStat)
            pwd.pop(0)
            for p in pwd:
                formula = pwd[p][0]
                InfixParser.parse(formula)
                for n in InfixParser.names:
                    if n not in self._names and n != '_TIME_':
                        self._names.append(n)
                formula = InfixParser.output
                thenStat = pwd[p][1].replace('self.', 'self.mod.')
                self.code_string += 'elif {}:\n    self.value = {}\n'.format(
                    formula, thenStat)
            self.code_string += 'else:\n    {}'.format(other)
            self.formula = self.code_string.replace('self.', '')
        self.xcode = compile(self.code_string, 'PieceWise', 'exec')

    def __call__(self):
        """
    Execute the compiled piecewise function and return the resulting value.

    The `__call__` method is responsible for running the piecewise logic 
    that has been compiled and stored in `self.xcode`. The method executes 
    the code in the current context and returns the value that results from 
    the execution of the piecewise function.

    Returns
    -------
    object
        The resulting value after executing the compiled piecewise logic.

    Examples
    --------
    >>> # Assume `PieceWise` is set up to handle specific conditions.
    >>> pw = PieceWise(pwd=[{'0': ('x > 0', '10'), '1': ('x < 0', '-10'), 'other': '0'}], mod=my_model)
    >>> pw.x = 5
    >>> result = pw()  # This will execute PieceWise logic, resulting in value 10
    >>> print(result)
    10

    >>> pw.x = -3
    >>> result = pw()  # This will now result in value -10
    >>> print(result)
    -10

    >>> pw.x = 0
    >>> result = pw()  # Since x is 0 and not covered by above conditions, it results in 'other' value 0
    >>> print(result)
    0
    """
        exec(self.xcode)
        return self.value


class PysMod(object):
    """
    Create a model object and instantiate a PySCeS model for further analyses. PySCeS model descriptions can be loaded from files or strings (see the *loader* argument for details).

    Parameters
    ----------
    File : str or None, optional
        Name of the PySCeS input file. If not explicit, a `.psc` extension is assumed.
    dir : str or None, optional
        Path to the input file. If not specified, the default PySCeS model directory (defined in the `pys_config.ini` file) is assumed.
    loader : str, optional
        Loader type. Defaults to 'file'. If set to 'string', an input file can be supplied as the *fString* argument.
    fString : str or None, optional
        A string containing a PySCeS model file (use with *loader='string'*). The *File* argument now specifies the new input file name.
    autoload : bool, optional
        Autoload the model, pre version 0.7.1, `mod.doLoad()` was required. Default is True.

    Attributes
    ----------
    __settings__ : dict
        Stores the model's configuration settings.
    WorkDir : str
        Working directory for the model.
    mode_integrator : str
        Default integrator for simulations. Can be 'LSODA' or 'CVODE'.
    random : module
        Reference to the random module, used within the class.
    sim_start : float
        Simulation start time.
    sim_end : float
        Simulation end time.
    sim_points : int
        Number of data points in the simulation.

    Methods
    -------
    enableDataPandas(var=True)
        Toggle output data type between numpy recarray and pandas DataFrame.
    ModelLoad(stoich_load=0)
        Load and instantiate a PySCeS model for further analyses.
    doLoad(stoich_load=0)
        Load and instantiate a PySCeS model. This function is retained for backwards compatibility and is being replaced by the ModelLoad() method.
    reLoad(stoich_load=0)
        Reloads the PySCeS model for further analyses.
    """
    __version__ = __version__
    __pysces_directory__ = INSTALL_DIR
    __settings__ = None
    random = random

    def __init__(self, File=None, dir=None, loader='file', fString=None,
        autoload=True):
        """
    Initializes a PysMod object and instantiates a PySCeS model for analysis. 
    The model descriptions can be loaded from files or strings based on the 
    specified loader.

    Parameters
    ----------
    File : str, optional
        The name of the PySCeS input file. If no extension is provided, 
        a default *.psc is assumed.
    dir : str, optional
        The directory path to the input file. If unspecified, defaults to 
        the PyscesModel directory defined in the pys_config.ini file.
    loader : {'file', 'string'}, optional
        Specifies the method of loading the PySCeS model. The default is 
        'file', in which case the model is loaded from a PSC file. If set to 
        'string', the model is loaded from a string provided as *fString*.
    fString : str, optional
        A string containing a PySCeS model file content when using 
        *loader='string'*. The *File* argument should then specify the file 
        name.
    autoload : bool, optional
        Whether to autoload the model upon initialization. Default is True. 
        If False, the user must manually call 'mod.doLoad()' for versions 
        prior to 0.7.1.

    Example
    -------
    >>> # Example of loading a model from a file
    >>> mod = PysMod(File='my_model.psc', dir='/models/', autoload=True)
    >>> 
    >>> # Example of loading a model from a string
    >>> model_string = 
    >>> mod = PysMod(File='new_model.psc', loader='string', fString=model_string)

    Notes
    -----
    Post initialization, the `ModelLoad()` method will be called if `autoload` 
    is set to True. If the model data type is set to 'pandas', it checks for 
    the availability of pandas and enables its use accordingly.
    """
        self.__settings__ = {}
        self.__settings__.update({'enable_deprecated_attr': True})
        self.WorkDir = CWD
        if loader == 'file':
            self.LoadFromFile(File, dir)
        elif loader == 'string':
            self.LoadFromString(File, fString)
        else:
            self.LoadFromFile(File, dir)
        self.__settings__['mode_substitute_assignment_rules'] = False
        self.__settings__['cvode_track_assignment_rules'] = True
        self.__settings__['display_compartment_warnings'] = False
        self.__settings__['custom_datatype'] = __CUSTOM_DATATYPE__
        self._TIME_ = 0.0
        self.piecewise_functions = []
        self.__piecewises__ = {}
        self.__HAS_PIECEWISE__ = False
        if autoload:
            self.ModelLoad()
            self.__PSC_auto_load = True
        else:
            self.__PSC_auto_load = False
        if __CUSTOM_DATATYPE__ == 'pandas':
            if _checkPandas():
                self.enableDataPandas(True)
            else:
                self.enableDataPandas(False)

    def enableDataPandas(self, var=True):
        """
    Toggle custom data type for `mod.sim` and `mod.scan` objects.

    If `var` is set to `True`, this method attempts to initialize the pandas library for 
    data processing, changing the data type for simulation and scan results to `pandas.DataFrame`.
    If pandas is not installed, it prints instructions for installation and defaults to not setting 
    a custom data type. If `var` is set to `False`, the method uses `numpy.recarray` as the data 
    type for results.

    Parameters
    ----------
    var : bool, optional
        If `True`, return results as a pandas DataFrame.
        If `False`, return results as a numpy recarray. 
        Default is `True`.

    Examples
    --------
    Enable pandas DataFrame output for simulation and scan results:

    >>> model = PysMod(File='example.psc')
    >>> model.enableDataPandas(True)
    Pandas output enabled.

    Disable pandas DataFrame output, revert to numpy recarray:

    >>> model.enableDataPandas(False)
    Pandas output disabled.
    
        Pandas is not installed. Install with:
            pip install pandas        or
            conda install pandas
        Unsetting custom datatype!
                """
        if var:
            try:
                import pandas
                self.__pandas = pandas
                self.__settings__['custom_datatype'] = 'pandas'
                print('Pandas output enabled.')
            except ModuleNotFoundError:
                print(
                    """
        Pandas is not installed. Install with:
            pip install pandas        or
            conda install pandas
        Unsetting custom datatype!
                    """
                    )
                self.__settings__['custom_datatype'] = None
        else:
            self.__settings__['custom_datatype'] = None
            print('Pandas output disabled.')

    def ModelLoad(self, stoich_load=0):
        """
    Load and instantiate a PySCeS model so that it can be used for further analyses. This function
    replaces the pre-0.7.1 doLoad() method.

    Parameters
    ----------
    stoich_load : int, optional
        Try to load a structural analysis saved with `Stoichiometry_Save_Serial()`. The default is 0,
        which means no structural analysis will be loaded.

    Raises
    ------
    AssertionError
        If there is an error in the input file, an assertion error is raised indicating parsing could
        not complete successfully.

    Examples
    --------
    Create and load a PySCeS model with structural analysis:

    >>> model = PysMod('example_model.psc')
    >>> model.ModelLoad(stoich_load=1)

    Load a PySCeS model without structural analysis (default behavior):

    >>> model = PysMod('example_model.psc')
    >>> model.ModelLoad()
    """
        self.InitialiseInputFile()
        assert self.__parseOK, '\nError in input file, parsing could not complete'
        self.Stoichiometry_Analyse(override=0, load=stoich_load)
        self.InitialiseFunctions()
        self.InitialiseCompartments()
        self.InitialiseRules()
        self.InitialiseEvents()
        self.InitialiseModel()
        self.InitialiseRuleChecks()
        self.InitialiseOldFunctions()

    def doLoad(self, stoich_load=0):
        """
    Load and instantiate a PySCeS model for further analysis. This function is
    being replaced by the `ModelLoad()` method.

    Parameters
    ----------
    stoich_load : int, optional
        Attempts to load a structural analysis saved with `Stoichiometry_Save_Serial()`.
        Default is 0, which indicates no loading of saved structural analysis.

    Notes
    -----
    The `doLoad()` method has been deprecated in favor of `ModelLoad()`. If auto-loading
    is enabled (`__PSC_auto_load`), the model is loaded automatically upon instantiation
    of the model object. If you prefer not to auto-load the model, pass `autoload=False`
    to the constructor. To reload the model explicitly, use the `reLoad()` method.

    Examples
    --------
    >>> model = PysMod(autoload=False)
    >>> model.doLoad(stoich_load=1)

    This will load the model and attempt to include any existing structural analysis data
    saved previously. If the model is already set to auto-load, a message will be printed
    explaining the auto-load behavior and suggesting alternatives.
    """
        if not self.__PSC_auto_load:
            self.ModelLoad(stoich_load=stoich_load)
        else:
            print(
                'PySCeS now automatically loads the model on model object instantiation. If you do not want this behaviour pass the autoload=False argument to the constructor, if you really want to reload the model, run reLoad().'
                )

    def reLoad(self, stoich_load=0):
        """
    Re-loads and instantiates a PySCeS model so that it can be used for further analyses. 
    This method serves as a convenience call to the `ModelLoad` method.

    Parameters
    ----------
    stoich_load : int, optional
        This parameter attempts to load a structural analysis that was saved using 
        `Stoichiometry_Save_Serial()`. The default value is 0, which indicates that 
        no structural analysis will be loaded unless specified otherwise.

    Notes
    -----
    This function acts as a simple wrapper around the `ModelLoad` method, providing 
    an easy way to re-instantiate the model without needing to directly invoke 
    `ModelLoad`. It allows users to quickly refresh their model analyses.

    Examples
    --------
    Basic usage example of `reLoad`:

    >>> model = PysMod()
    >>> model.reLoad(stoich_load=1)
    This will load the PySCeS model and attempt to load a previously saved 
    structural analysis if available.
    """
        self.ModelLoad(stoich_load=stoich_load)

    def LoadFromString(self, File=None, fString=None):
        """
    Load a PySCeS model from a string and save it as a file for further analyses.

    This method takes a model definition provided as a string and writes it to a file
    in the designated model directory. The model is prepared for subsequent operations,
    including analyses and simulations.

    Parameters
    ----------
    File : str, optional
        The name of the file to which the model will be saved. This file is stored
        within the subdirectory 'orca' of the `MODEL_DIR`. If not specified, a default
        filename or handling mechanism should be used.
    fString : str, optional
        The model's definition as a string. This string is expected to contain the
        entire specification of the model to be saved.

    Raises
    ------
    Exception
        If there are issues while writing the string to the file or if the directory
        structure does not support file creation.

    Notes
    -----
    This function automatically initializes certain settings related to the model's
    output directory and precision configurations for stoichiometric analysis.

    Example
    -------
    >>> model_definition = "A -> B; k1*A - k2*B"
    >>> model = PysMod()
    >>> model.LoadFromString(File="example_model.psc", fString=model_definition)
    Using model directory: /path/to/MODEL_DIR
    Using file: example_model.psc
    /path/to/MODEL_DIR/orca/example_model.psc loading .....
    Done.
    """
        chkmdir()
        mdir = MODEL_DIR
        File = chkpsc(File)
        print('Using model directory: ' + mdir)
        if not os.path.isdir(os.path.join(MODEL_DIR, 'orca')):
            os.mkdir(os.path.join(MODEL_DIR, 'orca'))
        mdir = os.path.join(MODEL_DIR, 'orca')
        try:
            outFile = open(os.path.join(mdir, File), 'w')
            outFile.write(fString)
            outFile.close()
            print('Using file: ' + File)
            if os.path.exists(os.path.join(mdir, File)):
                print(os.path.join(mdir, File) + ' loading .....', end=' ')
                self.ModelDir = mdir
                self.ModelFile = File
            else:
                print(os.path.join(mdir, File) + ' does not exist')
                print('Please set with ModelDir and ModelFile .....', end=' ')
                self.ModelFile = 'None'
                self.ModelDir = mdir
        except Exception as e:
            print(e)
            print(os.path.join(mdir, File) +
                ' does not exist please re-instantiate model .....', end=' ')
        self.__settings__['display_debug'] = 0
        self.ModelOutput = OUTPUT_DIR
        self.__settings__['serial_dir'] = os.path.join(self.ModelOutput,
            'pscdat')
        self.__settings__['stoichiometric_analysis_fp_zero'
            ] = mach_eps * 20000.0
        self.__settings__['stoichiometric_analysis_lu_precision'
            ] = self.__settings__['stoichiometric_analysis_fp_zero']
        self.__settings__['stoichiometric_analysis_gj_precision'
            ] = self.__settings__['stoichiometric_analysis_lu_precision'
            ] * 10.0
        """
        # Initialise elementary modes
        if os.sys.platform == 'win32':
            self.eModeExe_int = os.path.join(self.__metatool,'meta43_int.exe')
                self.eModeExe_dbl = os.path.join(self.__metatool,'meta43_double.exe')
        else:
            self.eModeExe_int = os.path.join(self.__metatool,'meta43_int')
            self.eModeExe_dbl = os.path.join(self.__metatool,'meta43_double')
        print 'Done.'
        """

    def LoadFromFile(self, File=None, dir=None):
        """
    Initialise a PySCeS model object with a PSC file that can be found in an optional directory.
    
    This method initializes the PySCeS model by loading a specified PSC file. If the file is not 
    specified, it displays the contents of the default model directory, allowing the user to 
    input the desired model name interactively. The loading process can be exited at any point 
    by pressing `<ctrl>+C`.

    Parameters
    ----------
    File : str, optional
        The name of the PySCeS input file. If not supplied, the user will be prompted to input 
        a model file name. Default is `None`.
    dir : str, optional
        The directory in which to look for the PSC file. If not specified, the function defaults 
        to the `pysces.model_dir` directory.

    Examples
    --------
    Initialise a PySCeS model by specifying the model file and directory:
    
    >>> model = PysMod()
    >>> model.LoadFromFile(File='example_model.psc', dir='/path/to/models')
    
    Initialise a PySCeS model interactively by specifying only the directory:
    
    >>> model = PysMod()
    >>> model.LoadFromFile(dir='/path/to/models')
    
    If no directory is specified, the default `pysces.model_dir` is used:
    
    >>> model = PysMod()
    >>> model.LoadFromFile(File='example_model.psc')

    Notes
    -----
    If neither `File` nor `dir` is specified, the user will be prompted to choose a model from 
    the default model directory. The user needs to have interactive console access to utilize this 
    feature effectively.
    """
        if dir != None:
            if os.path.isdir(dir):
                pass
            else:
                chkmdir()
                dir = MODEL_DIR
        else:
            chkmdir()
            dir = MODEL_DIR
        mfgo = 0
        if File is None:
            mfgo = 1
        try:
            if File != None:
                File = chkpsc(File)
            if not os.path.exists(os.path.join(dir, File)):
                mfgo = 1
        except:
            mfgo = 1
        while mfgo == 1:
            print('Models available in your model_dir: \n************')
            cntr = 0
            namelen = 0
            while len(os.listdir(dir)) == 0:
                print(
                    "No models available in model directory, please set using pys_usercfg.ini or call                with pysces.model('file.psc',dir='path\\to\\models')"
                    )
                dir = input(
                    '\nPlease enter full path to models <CTRL+C> exits: ')
            dirList = os.listdir(dir)
            for x in range(len(dirList) - 1, -1, -1):
                if dirList[x][-4:] != '.psc':
                    a = dirList.pop(x)
            for x in dirList:
                if len(x) > namelen:
                    namelen = len(x)
            for x in dirList:
                if cntr < 2:
                    print(x + ' ' * (namelen - len(x)), end=' ')
                    cntr += 1
                else:
                    print(x)
                    cntr = 0
            print('\n************\n')
            print('\nYou need to specify a valid model file ...\n')
            File = input('\nPlease enter filename: ')
            try:
                File = chkpsc(File)
                if os.path.exists(os.path.join(dir, File)):
                    mfgo = 0
            except:
                mfgo = 1
        print('Using model directory: ' + dir)
        try:
            if os.path.exists(os.path.join(dir, File)):
                print(os.path.join(dir, File) + ' loading .....', end=' ')
                self.ModelDir = dir
                self.ModelFile = File
            else:
                print(os.path.join(dir, File) + ' does not exist')
                print('Please set with ModelDir and ModelFile .....', end=' ')
                self.ModelFile = 'None'
                self.ModelDir = dir
        except:
            print(os.path.join(dir, File) +
                ' does not exist please re-instantiate model .....', end=' ')
        self.__settings__['display_debug'] = 0
        self.ModelOutput = OUTPUT_DIR
        self.__settings__['serial_dir'] = os.path.join(self.ModelOutput,
            'pscdat')
        self.__settings__['stoichiometric_analysis_fp_zero'
            ] = mach_eps * 20000.0
        self.__settings__['stoichiometric_analysis_lu_precision'
            ] = self.__settings__['stoichiometric_analysis_fp_zero']
        self.__settings__['stoichiometric_analysis_gj_precision'
            ] = self.__settings__['stoichiometric_analysis_lu_precision'
            ] * 10.0

    def __ParsePiecewiseFunctions__(self, piecewises):
        """
    Parses a dictionary of piecewise functions and creates corresponding piecewise objects,
    substituting the object names in the formula string. This method is highly experimental.

    Parameters
    ----------
    piecewises : dict
        A dictionary of piecewise definitions created by the InfixParser, where each key is
        the name of the piecewise function and the value is its definition.

    Notes
    -----
    This function will update the `__piecewises__` attribute with new piecewise objects and
    add them to the `piecewise_functions` list. Additionally, it will set attributes on the
    instance for each piecewise object.

    Examples
    --------
    Consider having a class `PysMod` with an instance `model`:

    >>> piecewise_dict = {
    ...     "pw1": "x < 5 ? 2 * x : 3 * x",
    ...     "pw2": "y > 10 ? y - 10 : y + 1"
    ... }
    >>> model.__ParsePiecewiseFunctions__(piecewise_dict)
    Info: adding piecewise object: pw1
    Info: adding piecewise object: pw2

    After execution, you can access the piecewise objects on the `model` object:

    >>> model.pw1  # Access the piecewise object named 'pw1'
    <PieceWise object at ...>
    >>> model.pw2  # Access the piecewise object named 'pw2'
    <PieceWise object at ...>
    """
        if piecewises != None and len(list(piecewises.keys())) > 0:
            self.__HAS_PIECEWISE__ = True
            for p in piecewises:
                print('Info: adding piecewise object: {}'.format(p))
                self.__piecewises__.update({p: piecewises[p]})
                P = PieceWise(piecewises[p], self)
                P.setName(p)
                self.piecewise_functions.append(P)
                setattr(self, p, P)

    def InitialiseInputFile(self):
        """
    Parses the input file associated with the PySCeS model instance and
    assigns the basic model attributes.

    The function checks for the existence of the model file in the specified
    directory, verifies its extension, and attempts to parse it to initialize
    the model attributes, which include species, parameters, and reactions.
    It raises errors in cases where parsing is unsuccessful due to invalid 
    symbols or syntax errors.

    If the input file does not have a `.psc` extension, this method assumes 
    the extension to be `.psc`. It then uses the `pscParser` to parse the file 
    and set up the model attributes. In the event of errors in parsing, it 
    prints the details of the errors.

    Raises
    ------
    AssertionError
        If the model parsing is unsuccessful.

    Attributes Set
    --------------
    __nDict__ : dict
        Dictionary of numerical identifiers for the entities in the model.
    __sDict__ : dict
        Dictionary of species attributes in the model.
    __pDict__ : dict
        Dictionary of parameter attributes in the model.
    __compartments__ : dict
        Dictionary storing compartment information.
    __functions__ : dict
        Dictionary containing user-defined functions in the model.
    __rules__ : dict
        Dictionary of assignment rules in the model.
    __eDict__ : dict
        Dictionary of events defined in the model.

    Example
    -------
    >>> model = PysMod()
    >>> model.ModelFile = 'example.psc'
    >>> model.ModelDir = '/path/to/model_directory'
    >>> model.InitialiseInputFile()
    Parsing file: /path/to/model_directory/example.psc
    Model file parsed successfully.
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
        print('\nParsing file: {}'.format(os.path.join(self.ModelDir, self.
            ModelFile)))
        pscParser.ParsePSC(self.ModelFile, self.ModelDir, self.WorkDir)
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
            InfixParser.__pwcntr__ = 0
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
            if self.__KeyWords__['ModelType'] is None:
                self.__KeyWords__['ModelType'] = ['Deterministic']
            else:
                self.__KeyWords__['ModelType'] = tuple([t.strip() for t in
                    self.__KeyWords__['ModelType'].split(',')])
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
            self.__cbm_fluxbounds__ = pscParser.cbm_FluxBounds.copy()
            self.__cbm_objfuncs__ = pscParser.cbm_ObjectiveFunctions.copy()
            self.__cbm_userfluxconstraints__ = (pscParser.
                cbm_UserFluxConstraints.copy())
        else:
            print('\nERROR: model parsing error, please check input file.\n')
        if len(pscParser.SymbolErrors) != 0:
            print('\nUndefined symbols:\n{}'.format(pscParser.SymbolErrors))
        if not pscParser.ParseOK:
            print('\n\n*****\nModel parsing errors detected in input file ' +
                self.ModelFile + '\n*****')
            print('\nInput file errors')
            for error in pscParser.LexErrors:
                print(error)
            print('\nParser errors')
            for error in pscParser.ParseErrors:
                print(error)
            assert pscParser.ParseOK == 1, 'Input File Error'

    def InitialiseRules(self):
        """
    Initializes the rules for a PySCeS model, categorizing them into 
    assignment rules and rate rules, and preparing them for further 
    computation. 

    The function determines the presence of forced functions (assignment 
    rules) and rate rules within the model, processes each type, and 
    compiles necessary Python code strings for evaluating the rules 
    during simulation.

    Attributes
    ----------
    __HAS_FORCED_FUNCS__ : bool
        Indicates if the model has assignment rules.
    __HAS_RATE_RULES__ : bool
        Indicates if the model has rate rules.
    _CVODE_extra_output : list
        Stores additional CVODE output requirements if needed.
    _CVODE_XOUT : bool
        A CVODE output flag, initially set to False.
    _NewRuleXCode : dict
        Contains code strings for new assignment rules.
    _Function_forced : str
        Code string for forced assignment functions.
    __CODE_forced : code object
        Compiled code for forced assignment rules.
    __rate_rules__ : list
        List to store identified rate rules.
    _CODE_raterule : code object
        Compiled code for rate rule evaluation.
    _CODE_raterule_map : code object
        Compiled code for mapping rate rule results.

    Example
    -------
    >>> model = PysMod()
    >>> model.InitialiseRules()
    Assignment rule(s) detected.
    Rate rule(s) detected.

    This method should be called after parsing and setting up the 
    model's input structure and before simulating the model to ensure 
    all rules are accurately compiled and ready for execution.
    """
        self.__HAS_FORCED_FUNCS__ = False
        self.__HAS_RATE_RULES__ = False
        self._CVODE_extra_output = []
        self._CVODE_XOUT = False
        rate_rules = {}
        assignment_rules = {}
        for ar in self.__rules__:
            if self.__rules__[ar]['type'] == 'assignment':
                self.__HAS_FORCED_FUNCS__ = True
                assignment_rules.update({ar: self.__rules__[ar]})
            elif self.__rules__[ar]['type'] == 'rate':
                self.__HAS_RATE_RULES__ = True
                rate_rules.update({ar: self.__rules__[ar]})
        InfixParser.setNameStr('self.', '')
        self._NewRuleXCode = {}
        if self.__HAS_FORCED_FUNCS__ and not self.__settings__[
            'mode_substitute_assignment_rules']:
            code_string = ''
            all_names = []
            for ar in assignment_rules:
                InfixParser.setNameStr('self.', '')
                InfixParser.parse(assignment_rules[ar]['formula'])
                assignment_rules[ar]['symbols'] = InfixParser.names
                formula = InfixParser.output
                for pw in InfixParser.piecewises:
                    formula = formula.replace('self.{}'.format(pw),
                        'self.{}()'.format(pw))
                self.__ParsePiecewiseFunctions__(InfixParser.piecewises)
                assignment_rules[ar]['code_string'] = formula
                all_names += InfixParser.names
            keep = []
            rules = list(assignment_rules.keys())
            dep = rules
            while len(dep) > 0:
                dep = []
                indep = []
                for ar in rules:
                    if ar in all_names:
                        indep.append(ar)
                    else:
                        dep.append(ar)
                keep += dep
                rules = indep
            for ar in (indep + keep):
                if self.__settings__['cvode_track_assignment_rules']:
                    self._CVODE_extra_output.append(assignment_rules[ar][
                        'name'])
                evalCode = 'self.{} = {}\n'.format(assignment_rules[ar][
                    'name'], assignment_rules[ar]['code_string'])
                self._NewRuleXCode.update({assignment_rules[ar]['name']:
                    evalCode})
                code_string += evalCode
            print('\nAssignment rule(s) detected.')
            self._Function_forced = code_string
        elif self.__HAS_FORCED_FUNCS__ and self.__settings__[
            'mode_substitute_assignment_rules']:
            InfixParser.setNameStr('self.', '')
            symbR = {}
            for ass in assignment_rules:
                InfixParser.setNameStr('self.', '')
                InfixParser.parse(assignment_rules[ass]['formula'])
                formula = InfixParser.output
                for pw in InfixParser.piecewises:
                    formula = formula.replace('self.{}'.format(pw),
                        'self.{}()'.format(pw))
                self.__ParsePiecewiseFunctions__(InfixParser.piecewises)
                symbR.update({assignment_rules[ass]['name']: formula})
            for ass in assignment_rules:
                InfixParser.setNameStr('self.', '')
                InfixParser.FunctionReplacements = symbR
                InfixParser.parse(assignment_rules[ass]['formula'])
                assignment_rules[ass]['code_string'] = InfixParser.output
            self._Function_forced = 'pass\n'
            print('Assignment rule(s) detected and substituted.')
        else:
            self._Function_forced = 'pass\n'
        self.__CODE_forced = compile(self._Function_forced, 'AssignRules',
            'exec')
        self.__rate_rules__ = []
        self.__rrule__ = None
        rr_code_block = ''
        rr_map_block = ''
        self._CODE_raterule = None
        self._NewRateRuleXCode = {}
        if self.__HAS_RATE_RULES__:
            self.__rrule__ = numpy.ones(len(list(rate_rules.keys())), 'd')
            cntr = 0
            rr_keys = list(rate_rules.keys())
            rr_keys.sort()
            for ar in rr_keys:
                name = rate_rules[ar]['name']
                InfixParser.setNameStr('self.', '')
                InfixParser.parse(rate_rules[ar]['formula'])
                formula = InfixParser.output
                rate_rules[ar]['symbols'] = InfixParser.names
                for pw in InfixParser.piecewises:
                    formula = formula.replace('self.{}'.format(pw),
                        'self.{}()'.format(pw))
                self.__ParsePiecewiseFunctions__(InfixParser.piecewises)
                rate_rules[ar]['code_string'] = formula
                self.__rate_rules__.append(name)
                rr_code_block += 'self.__rrule__[{}] = {}\n'.format(cntr,
                    formula)
                rr_map_block += 'self.{} = self.__rrule__[{}]\n'.format(name,
                    cntr)
                cntr += 1
                setattr(self, '{}_init'.format(name), getattr(self, name))
            print('Rate rule(s) detected.')
        else:
            rr_code_block = 'pass\n'
            rr_map_block = 'pass\n'
        self._CODE_raterule = compile(rr_code_block, 'RateRules', 'exec')
        self._CODE_raterule_map = compile(rr_map_block, 'RateRuleMap', 'exec')
        del rate_rules, assignment_rules

    def ExecRateRules(self):
        """
    Execute the precompiled rate rule code.

    This method runs the rate rules that have been precompiled into bytecode
    and stored in the `_CODE_raterule` attribute. The execution of these rules
    updates the internal variables of the model based on the defined rate
    rules.

    Parameters
    ----------
    None

    Notes
    -----
    This method assumes that the rate rules have been initialized and compiled
    using the `InitialiseRules` method. If `InitialiseRules` has not been run,
    this method may result in an error or undefined behavior.

    Examples
    --------
    >>> model = PysMod()
    >>> model.InitialiseRules()  # Make sure rules are initialized and compiled
    >>> model.ExecRateRules()    # Execute the rate rules
    """
        exec(self._CODE_raterule)

    def InitialiseRuleChecks(self):
        """
    Execute and validate the rate and assignment rules, handling exceptions.

    This method initializes and checks rate and assignment rules by executing
    precompiled code objects. It handles exceptions such as `ZeroDivisionError`
    and other generic exceptions during execution while attempting to resolve
    issues such as division by zero in subsequent executions.

    Raises
    ------
    ZeroDivisionError
        If a division by zero occurs during the execution of rate rules or
        assignment rules.
    Exception
        For any other exception that may occur during rule processing.

    Warns
    -----
    Prints a warning message if a `ZeroDivisionError` or any other exception
    is encountered during execution. A message will indicate if zero division
    elimination has failed, or if assignment rules are disabled due to an error.

    Example
    -------
    >>> p = PysMod()
    >>> p.InitialiseRules() # Assume this method sets up the rules
    >>> p.InitialiseRuleChecks()
    WARNING: Assignment RateRule ZeroDivision on initialisation (continuing)
    WARNING: RateRule initialisation error
    [Exception message]
    WARNING: Assignment rule ZeroDivision on initialisation
    WARNING: Assignment rule error (disabling all rules)
    [Exception message]
    """
        try:
            exec(self._CODE_raterule)
        except ZeroDivisionError:
            print(
                'WARNING: Assignment RateRule ZeroDivision on initialisation (continuing)'
                )
        except Exception as ex:
            print('WARNING: RateRule initialisation error\n', ex)
        zeroDivErr = []
        for k in self._NewRuleXCode:
            try:
                exec(compile(self._NewRuleXCode[k], 'NewRuleXCode', 'exec'))
            except Exception as ex:
                zeroDivErr.append(k)
        exit = len(list(self._NewRuleXCode.keys()))
        while len(zeroDivErr) > 0 and exit > 0:
            zeroDivErr2 = []
            for kk in zeroDivErr:
                try:
                    exec(compile(self._NewRuleXCode[kk], 'NewRuleXCode',
                        'exec'))
                    if self.__HAS_RATE_RULES__:
                        exec(self._CODE_raterule)
                        exec(self._CODE_raterule_map)
                except Exception as ex:
                    zeroDivErr2.append(kk)
            exit -= 1
            zeroDivErr = zeroDivErr2
            if exit < 0:
                print('WARNING: ZeroDivision elimination failed')
        try:
            exec(self.__CODE_forced)
        except ZeroDivisionError:
            print('WARNING: Assignment rule ZeroDivision on intialisation')
        except Exception as ex:
            print('WARNING: Assignment rule error (disabling all rules)\n', ex)
            self.__CODE_forced = compile('pass\n', 'AssignRules', 'exec')

    def Stoichiometry_Init(self, nmatrix, load=0):
        """
    Initialize the model stoichiometry.

    This method initializes the stoichiometry for a given model by processing a stoichiometric matrix `nmatrix`. 
    If `load` is enabled, the method attempts to load a previously saved stoichiometric analysis using 
    `Stoichiometry_Save_Serial` and verifies its correctness. The method returns an instantiated `PyscesStoich` 
    instance and a status flag indicating whether the stoichiometry needs reanalysis or if a complete structural 
    analysis is preloaded.

    Parameters
    ----------
    nmatrix : array-like
        The input stoichiometric matrix, N.
    load : int, optional
        Indicates whether to attempt loading a saved stoichiometry. Default is 0.
        If set to 1, the method will try to load and verify previously saved stoichiometric data.

    Returns
    -------
    PyscesStoich
        An instance of PyscesStoich representing the stoichiometry of the model.
    status : int
        A status flag where 0 means the stoichiometry needs to be reanalyzed, and 
        1 indicates a complete structural analysis was successfully loaded.

    Notes
    -----
    This method checks the consistency and compatibility of previously saved stoichiometric analysis 
    if `load` is specified, ensuring that the species and reaction orders match. Any detected discrepancies 
    will result in a reanalysis of the stoichiometry.

    Examples
    --------
    Initialize a new stoichiometry analysis without loading:

    >>> model = PysMod()
    >>> nmatrix = [[1, -1, 0], [0, 1, -1], [-1, 0, 1]]  # Example stoichiometric matrix
    >>> stc, status = model.Stoichiometry_Init(nmatrix)
    >>> print(status)  # Output: 0, since it reanalyzed the stoichiometry

    Attempt to load a previously saved stoichiometry:

    >>> stc, status = model.Stoichiometry_Init(nmatrix, load=1)
    >>> print(status)  # Output: 1 if loading and verification were successful, otherwise 0
    """
        if load:
            print('Loading saved stoichiometry ...')
            try:
                stc = self.Stoichiometry_Load_Serial()
                go = 1
            except Exception as slx:
                print(slx)
                go = 0
            row, col = nmatrix.shape
            if go:
                badList = []
                for x in range(row):
                    if stc.__species__[x] != self.__species__[x]:
                        badList.append((stc.__species__[x], self.
                            __species__[x]))
                        go = 0
                    for y in range(col):
                        if stc.__reactions__[y] != self.__reactions__[y]:
                            badList.append((stc.__reactions__[y], self.
                                __reactions__[y]))
                            go = 0
                        if abs(nmatrix[x, y] - stc.nmatrix[x, y]
                            ) > stc.stoichiometric_analysis_fp_zero:
                            badList.append(((x, y), nmatrix[x, y], stc.
                                nmatrix[x, y]))
                            go = 0
            if not go:
                print('\nProblem loading stoichiometry, reanalysing ...')
                stc = PyscesStoich.Stoich(nmatrix)
                status = 0
            else:
                print('Stoichiometry verified ... we have liftoff')
                status = 1
        else:
            stc = PyscesStoich.Stoich(nmatrix)
            status = 0
        return stc, status

    def Stoichiometry_Save_Serial(self):
        """
    Serialize and save a Stoichiometric instance to a binary pickle file.

    This method serializes and saves the current model's stoichiometry to a file named
    `<model>_stoichiometry.pscdat`. The file is stored in the directory specified by
    `mod.__settings__['serial_dir']`, which defaults to `mod.model_output/pscdat`.

    Parameters
    ----------
    None

    Examples
    --------
    Assuming `mod` is an instance of the `PysMod` class and the model's stoichiometry
    has been defined, you can save the stoichiometry as follows:

    >>> mod = PysMod()
    >>> mod.Stoichiometry_Save_Serial()

    This will save the stoichiometric data to a file in the specified directory.

    """
        self.__structural__.__species__ = copy.copy(self.__species__)
        self.__structural__.__reactions__ = copy.copy(self.__reactions__)
        self.SerialEncode(self.STOICH, self.ModelFile[:-4] + '_stoichiometry')

    def Stoichiometry_Load_Serial(self):
        """
    Load a saved stoichiometry instance.

    This method loads a previously saved stoichiometry analysis that was stored 
    using the `Stoichiometry_Save_Serial` method. It returns a `PyscesStoich` 
    instance representing the stoichiometry of the model.

    Returns
    -------
    PyscesStoich
        An instance of the `PyscesStoich` class that contains the loaded 
        stoichiometry data.

    Notes
    -----
    The stoichiometry data is expected to be stored in a file named 
    `<model>_stoichiometry.pscdat`, located in the directory specified by 
    `mod.__settings__['serial_dir']`.

    Raises
    ------
    FileNotFoundError
        If the stoichiometry file does not exist.
    Exception
        If there is an error during the deserialization process.

    Examples
    --------
    Load a previously saved stoichiometry:

    >>> model = PysMod()
    >>> stoichiometry_instance = model.Stoichiometry_Load_Serial()
    >>> print(stoichiometry_instance)
    <PyscesStoich instance with loaded data>
    """
        stc = self.SerialDecode(self.ModelFile[:-4] + '_stoichiometry')
        return stc

    def Stoichiometry_Analyse(self, override=0, load=0):
        """
    Perform a structural analysis of the stoichiometry of the model.

    The default behavior is to construct and analyze the model from the parsed 
    model information. When `override` is set to 1, the analysis will be based 
    on the current stoichiometric matrix instead. If `load` is specified, 
    PySCeS attempts to load a saved stoichiometry; otherwise, it performs a 
    stoichiometric analysis. The results of the analysis are checked for 
    floating point error and nullspace rank consistency.

    Parameters
    ----------
    override : int, optional
        Override stoichiometric analysis initialization from parsed data. 
        Default is 0.
    load : int, optional
        Load a pre-saved stoichiometry. Default is 0.

    Examples
    --------
    >>> mod = PysMod()
    >>> mod.Stoichiometry_Analyse()
    Performing structural analysis using parsed model information...

    >>> mod.Stoichiometry_Analyse(override=1)
    Stoichiometric override active
    Performing structural analysis based on current stoichiometric matrix...

    >>> mod.Stoichiometry_Analyse(load=1)
    Attempting to load saved stoichiometry... successful.

    Notes
    -----
    The analysis checks for:
    - Floating point errors in the L0 and K0 matrices.
    - Consistency in the rank between K and L analysis methods.
    If any inconsistencies or potential issues are detected, warnings will be
    displayed, and the user should consider adjusting analysis precision settings.
    """
        if not override:
            self.__nmatrix__ = self.__initmodel__()
        else:
            print('\nStoichiometric override active\n')
        assert len(self.__nmatrix__
            ) > 0, """
Unable to generate Stoichiometric Matrix! model has:
{} reactions
{} species
what did you have in mind?
""".format(
            len(self.__reactions__), len(self.__species__))
        self.__Nshape__ = self.__nmatrix__.shape
        self.__structural__, stc_load = self.Stoichiometry_Init(self.
            __nmatrix__, load=load)
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
            self.nmatrix = self.__nmatrix__
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

    def Stoichiometry_ReAnalyse(self):
        """
    Reanalyse the stoichiometry using the current N matrix.

    This method reanalyzes the stoichiometry by forcing the use of the current 
    stoichiometric (N) matrix. It effectively sets `override=1` to ensure that 
    the analysis is based on the most recent modifications to the N matrix,
    which can be useful after changes made with `mod.Stoich_matrix_SetValue`.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Example
    -------
    >>> mod = PysMod()
    >>> mod.Stoich_matrix_SetValue(...)
    >>> mod.Stoichiometry_ReAnalyse()

    This will update the stoichiometry analysis to reflect the current 
    stoichiometric changes made via `Stoich_matrix_SetValue`.
    """
        self.Stoichiometry_Analyse(override=1)
        self.InitialiseConservationArrays()

    def Stoich_nmatrix_SetValue(self, species, reaction, value):
        """
    Change a stoichiometric coefficient's value in the N matrix.

    This method updates the value of a stoichiometric coefficient in the N matrix for a specified
    species and reaction. Only the magnitude of the coefficient may be set, meaning its sign
    (negative, positive, or zero) must remain unchanged. After modifying a coefficient, it is
    necessary to reanalyse the stoichiometry to ensure consistency.

    Parameters
    ----------
    species : str
        The name of the species for which the stoichiometric coefficient will be set, e.g., 's0'.
    reaction : str
        The name of the reaction in which the stoichiometric coefficient will be set, e.g., 'R4'.
    value : float
        The new value for the stoichiometric coefficient. Must adhere to the original sign of the
        coefficient.

    Examples
    --------
    To set a new stoichiometric coefficient value of -1.5 for species 's1' in reaction 'R3', use:

    >>> model = PysMod()  # Assume model is an instance of PysMod
    >>> model.Stoich_nmatrix_SetValue(species='s1', reaction='R3', value=-1.5)
    >>> model.Stoichiometry_ReAnalyse()  # Reanalyse to reflect changes in the N matrix

    Notes
    -----
    It is crucial to reanalyse the stoichiometry after changing any coefficients to maintain
    the accuracy and reliability of the model's analysis.
    """
        index = self.__Stoich_nmatrix_FindIndex__(species, reaction)
        if self.__Stoich_nmatrix_CheckValue__(index, value):
            self.__Stoich_nmatrix_UpdateValue__(index, value)
            self.__StoichOK = 0

    def __Stoich_nmatrix_FindIndex__(self, species, reaction):
        """
    Return the N matrix coordinates of the coefficient referenced by species and reaction name.

    Parameters
    ----------
    species : str
        The name of the species, e.g., 's0'.
    reaction : str
        The name of the reaction, e.g., 'R4'.

    Returns
    -------
    tuple
        A tuple containing the row index and column index corresponding to the specified species and
        reaction in the N matrix. Returns (None, None) if either species or reaction is not found.

    Example
    -------
    Suppose you have a model with species ['s0', 's1'] and reactions ['R1', 'R2', 'R4'] and you want to
    find the position of the coefficient for species 's0' in reaction 'R4':

    >>> pys_mod = PysMod()
    >>> index = pys_mod.__Stoich_nmatrix_FindIndex__('s0', 'R4')
    >>> print(index)
    (0, 2)
    
    Notes
    -----
    This method assumes that `self.species` and `self.reactions` are lists that contain the names of the
    species and reactions in the order they appear in the N matrix.
    """
        rval = None
        cval = None
        for x in enumerate(self.species):
            if species == x[1]:
                rval = x[0]
        for y in enumerate(self.reactions):
            if reaction == y[1]:
                cval = y[0]
        return rval, cval

    def __Stoich_nmatrix_UpdateValue__(self, xxx_todo_changeme, val):
        """
    Updates the N matrix at the specified coordinates with the given value.

    Parameters
    ----------
    xxx_todo_changeme : tuple of int
        A tuple (x, y) representing the row and column coordinates in the matrix.
    val : float
        The new value to set at the specified coordinates in the matrix.

    Examples
    --------
    >>> pys_mod = PysMod()
    >>> x, y = 2, 3
    >>> val = 4.5
    >>> pys_mod.__Stoich_nmatrix_UpdateValue__((x, y), val)
    >>> print(pys_mod.nmatrix[x, y])
    4.5
    >>> print(pys_mod.__nmatrix__[x, y])
    4.5
    >>> print(pys_mod.Nmatrix.array[x, y])
    4.5

    Notes
    -----
    This is a private method meant to update the N matrix, and its usage is typically internal to the class.
    
    """
        x, y = xxx_todo_changeme
        self.nmatrix[x, y] = val
        self.__nmatrix__[x, y] = val
        self.Nmatrix.array[x, y] = val

    def __Stoich_nmatrix_CheckValue__(self, xxx_todo_changeme1, val):
        """
    __Stoich_nmatrix_CheckValue__((x,y), val)

    Check the validity of a new coefficient value against the existing N matrix coefficient 
    at position (x, y) for existence and/or sign. Returns `1` if valid, otherwise `0`.

    Parameters
    ----------
    (x, y) : tuple of int
        The coordinates (row, column) within the N matrix to check.
    val : float
        The new coefficient value to be validated against the existing value in the N matrix.

    Returns
    -------
    int
        `1` if the new coefficient value is valid, otherwise `0`.

    Raises
    ------
    ValueError
        If the indices (x, y) do not match the expected species or reaction in the N matrix.

    Notes
    -----
    - Prints error messages describing the mismatch if the indices are out of expected bounds, 
      or if there is a coefficient sign mismatch.
    - It assumes that valid N matrix coefficients are expected to have specific signs (positive, 
      negative, or zero) and flags any violation of these rules.

    Examples
    --------
    >>> pysmod_instance = PysMod()
    >>> result = pysmod_instance.__Stoich_nmatrix_CheckValue__((2, 3), -0.5)
    >>> if result:
    ...     print("The coefficient value is valid.")
    ... else:
    ...     print("The coefficient value is not valid.")
    The coefficient value is valid.

    """
        x, y = xxx_todo_changeme1
        go = 1
        if x is None or y is None:
            print('\nSpecies (nmatrix) index  =', x)
            print('Reaction (nmatrix) index =', y)
            print(
                "\nI'm confused, perhaps you entered an incorrect species or reaction name"
                )
            print('or they are in the wrong order?')
            go = 0
        elif abs(self.__nmatrix__[x, y]) == 0.0:
            if val != 0.0:
                go = 0
                print('\nZero coefficient violation')
                print('  nmatrix[' + repr(x) + ',' + repr(y) +
                    '] can only be = 0.0 (input ' + str(val) + ')')
        elif self.__nmatrix__[x, y] > 0.0:
            if val <= 0.0:
                go = 0
                print('\nPositive coefficient violation')
                print('  nmatrix[' + repr(x) + ',' + repr(y) +
                    '] can only be > 0.0 (input ' + str(val) + ')')
        elif self.__nmatrix__[x, y] < 0.0:
            if val >= 0.0:
                go = 0
                print('\nNegative coefficient violation')
                print('  nmatrix[' + repr(x) + ',' + repr(y) +
                    '] can only be < 0.0 (input ' + str(val) + ')')
        if go:
            return 1
        else:
            return 0

    def SerialEncode(self, data, filename):
        """
    SerialEncode(data, filename)

    Serialise and save a Python object using a binary pickle to a file. The serialised object
    is saved as <filename>.pscdat in the directory defined by `mod.model_serial`.

    Parameters
    ----------
    data : any
        A pickleable Python object to be serialised and saved.
    filename : str
        The base name of the output file without extension. The saved file will have a `.pscdat` extension.

    Raises
    ------
    Exception
        If an error occurs during the pickling process, the exception is caught, and an error message 
        is printed to the console.

    Notes
    -----
    The method will check if the directory specified in `self.__settings__['serial_dir']` exists. 
    If it does not exist, the method will create the directory before saving the file.

    Example
    -------
    >>> mod = PysMod()
    >>> my_data = {'key1': 'value1', 'key2': [1, 2, 3]}
    >>> mod.SerialEncode(my_data, 'my_serialised_data')
    
    After execution, a file named `my_serialised_data.pscdat` is created in the directory 
    specified by `mod.model_serial`.
    """
        filename = str(filename) + '.pscdat'
        if os.path.exists(self.__settings__['serial_dir']):
            pass
        else:
            os.mkdir(self.__settings__['serial_dir'])
        print('\ncPickle data stored in: ' + self.__settings__['serial_dir'])
        File = open(os.path.join(self.__settings__['serial_dir'], filename),
            'wb')
        try:
            pickle.dump(data, File, -1)
            print('Serialization complete')
        except Exception as E:
            print('Serialize error:\n', E)
        File.flush()
        File.close()
        del data

    def SerialDecode(self, filename):
        """
    Decode and return a serialised object saved with SerialEncode.

    This method reads a binary file created by the `SerialEncode` method and 
    reconstructs the original Python object. The method assumes that the filename 
    does not include the '.pscdat' extension, which it adds automatically. The 
    method will raise an error if the specified file does not exist in the 
    directory set in `self.__settings__['serial_dir']`.

    Parameters
    ----------
    filename : str
        The filename without the '.pscdat' extension for the serialized object.

    Returns
    -------
    object
        The deserialized Python object that was originally serialized and saved 
        using the `SerialEncode` method.

    Raises
    ------
    RuntimeError
        If the serialized data file does not exist in the specified directory.

    Examples
    --------
    >>> mod = PysMod()
    >>> obj = {'key': 'value'}
    >>> mod.SerialEncode(obj, 'sample_file')  # Assuming this method works
    >>> decoded_obj = mod.SerialDecode('sample_file')
    >>> print(decoded_obj)
    {'key': 'value'}
    """
        filename = str(filename) + '.pscdat'
        if not os.path.exists(os.path.join(self.__settings__['serial_dir'],
            filename)):
            raise RuntimeError('Serialized data ' + os.path.join(self.
                __settings__['serial_dir'], filename) + ' does not exist')
        data = None
        try:
            File = open(os.path.join(self.__settings__['serial_dir'],
                filename), 'rb')
            data = pickle.load(File)
            File.close()
            print('Serial decoding complete')
        except Exception as E:
            print('Serial decode error:\n', E)
        return data

    def InitialiseCompartments(self):
        """
    Initialize and rescale model compartments.

    This method initializes compartments according to the specified configurations
    and rescales their sizes based on a predefined fudge factor. If compartments
    are present and need adjustment, it updates their sizes and ensures species
    are placed in the correct compartment. Species concentrations are adjusted
    accordingly if they are specified in terms of concentration, and a scale
    factor is applied.

    The following attributes are modified or updated:
    - `__settings__['compartment_fudge_factor']`
    - Attributes corresponding to compartment names with their respective sizes.
    - Species initial values in `__sDict__` adjusted by the fudge factor.
    - `__CsizeAllIdx__` and `__CsizeAll__` are populated with the current index and sizes of compartments.

    Warnings:
    If a species defined in the model is not assigned to a compartment, and multiple compartments are available,
    an assertion error is raised, indicating that assignment to a compartment is necessary.

    Example
    -------
    >>> model = PysMod()  # Assume PysMod is initialized with compartments and species
    >>> model.InitialiseCompartments()
    INFO: Rescaling compartment with size 0.0001 to 1.0
    COMPARTMENT WARNING: this model has a compartment defined ['compartment1'] but "species1" is not in it ... I'll try it for now
    INFO: Rescaling species (species1) to size 5.0.
    INFO: Rescaling species (species2) to size 10.0.

    Notes
    -----
    - This method assumes that the `__HAS_COMPARTMENTS__` attribute specifies whether compartments are used.
    - The `startFudging` value is used to determine if rescaling is necessary.

    Raises
    ------
    AssertionError
        If a species is not assigned to a compartment when multiple compartments are defined.
    """
        self.__settings__['compartment_fudge_factor'] = 1.0
        startFudging = 1e-08
        I_AM_FUDGING = False
        if self.__HAS_COMPARTMENTS__:
            tmp = min([abs(self.__compartments__[c]['size']) for c in list(
                self.__compartments__.keys())])
            if tmp < startFudging:
                self.__settings__['compartment_fudge_factor'] = tmp
            for c in list(self.__compartments__.keys()):
                if self.__settings__['compartment_fudge_factor'
                    ] < startFudging:
                    newsize = self.__compartments__[c]['size'
                        ] / self.__settings__['compartment_fudge_factor']
                    print('INFO: Rescaling compartment with size {} to {}'.
                        format(self.__compartments__[c]['size'], newsize))
                    self.__compartments__[c]['size'] = newsize
                    self.__compartments__[c].update({'scale': self.
                        __settings__['compartment_fudge_factor']})
                    I_AM_FUDGING = True
                setattr(self, self.__compartments__[c]['name'], self.
                    __compartments__[c]['size'])
                setattr(self, '{}_init'.format(self.__compartments__[c][
                    'name']), self.__compartments__[c]['size'])
        if self.__HAS_COMPARTMENTS__:
            for sp in list(self.__sDict__.keys()):
                if self.__sDict__[sp]['compartment'] is None and len(list(
                    self.__compartments__.keys())) == 1:
                    self.__sDict__[sp]['compartment'] = self.__compartments__[
                        list(self.__compartments__.keys())[0]]['name']
                    print(
                        'COMPARTMENT WARNING: this model has a compartment defined {} but "{}" is not in it ... I\'ll try it for now'
                        .format([self.__compartments__[c]['name'] for c in
                        list(self.__compartments__.keys())], sp))
                elif self.__sDict__[sp]['compartment'] is None and len(list
                    (self.__compartments__.keys())) > 1:
                    assert self.__sDict__[sp]['compartment'
                        ] != None, """
COMPARTMENT ERROR: this model has multiple compartments defined {} but "{}" is not in one!""".format(
                        [self.__compartments__[c]['name'] for c in list(
                        self.__compartments__.keys())], sp)
                if not self.__KeyWords__['Species_In_Conc'] and I_AM_FUDGING:
                    self.__sDict__[sp]['initial'] = self.__sDict__[sp][
                        'initial'] / self.__settings__[
                        'compartment_fudge_factor']
                    setattr(self, sp, self.__sDict__[sp]['initial'])
                    setattr(self, '{}_init'.format(sp), self.__sDict__[sp][
                        'initial'])
                    print('INFO: Rescaling species ({}) to size {}.'.format
                        (sp, self.__sDict__[sp]['initial']))
        self.__CsizeAllIdx__ = []
        self.__null_compartment__ = 1.0
        uses_null_compartments = False
        for s in self.__species__:
            if self.__HAS_COMPARTMENTS__:
                self.__CsizeAllIdx__.append(self.__sDict__[s]['compartment'])
            else:
                self.__CsizeAllIdx__.append('__null_compartment__')
                uses_null_compartments = True
        if uses_null_compartments:
            setattr(self, '__null_compartment___init', 1.0)
        if self.__HAS_COMPARTMENTS__:
            self.__CsizeAll__ = numpy.array([self.__compartments__[c][
                'size'] for c in self.__CsizeAllIdx__], 'd')
        else:
            self.__CsizeAll__ = numpy.ones(len(self.__CsizeAllIdx__), 'd')

    def InitialiseFunctions(self):
        """
    Initialise and set up function objects.

    This method initializes function objects from predefined functions in the
    model and assigns them as attributes of the instance. Each function is
    initialized with arguments and a formula. If any of the arguments reference
    other functions within the model, these are also set as attributes on the
    function object.

    Examples
    --------
    Assuming `model` is an instance of the `PysMod` class:

    >>> model.InitialiseFunctions()
    >>> print(model.some_function)
    <Function object representing some_function>

    Notes
    -----
    - This method assumes that `self.__functions__` is a dictionary where each
      key is the name of a function, and the value is a dictionary with 'args'
      and 'formula' keys.
    - The `Function` class must have `setArg` and `addFormula` methods for setting
      arguments and the function formula respectively.
    """
        for func in self.__functions__:
            F = Function(func, self)
            for arg in self.__functions__[func]['args']:
                F.setArg(arg.strip())
            F.addFormula(self.__functions__[func]['formula'])
            setattr(self, func, F)
        for func in self.__functions__:
            fobj = getattr(self, func)
            for f in fobj.functions:
                if f in self.__functions__:
                    setattr(fobj, f, getattr(self, f))

    def InitialiseEvents(self):
        """
    Initialize events for the PysMod model instance.

    This method sets up the event-handling mechanisms by creating event 
    objects based on the model's event definitions. It configures the event's
    trigger condition, delay, priority, persistent flag, and assignments.

    The events are stored in the instance's `__events__` attribute, 
    and each event object is also set as an attribute of the instance with 
    the event's name.

    Additionally, the method updates the model settings to ensure that 
    CVODE, the solver, returns event time points and sets a flag 
    (`__HAS_EVENTS__`) to indicate whether any events have been defined.

    Attributes
    ----------
    __settings__ : dict
        Configuration settings for the model, updated to return event time points.
    __events__ : list
        A list of initialized event objects.
    __HAS_EVENTS__ : bool
        A flag indicating whether the model contains any events.

    Raises
    ------
    None

    Examples
    --------
    Assuming `model` is an instance of `PysMod` with events defined in
    `__eDict__`:

    >>> model.InitialiseEvents()
    >>> print(model.__events__)
    [<Event object at 0x...>, ...]

    Each event defined in `model.__eDict__` is initialized and can be accessed via:
    
    >>> print(model.some_event_name)
    <Event object at 0x...>

    """
        self.__settings__['cvode_return_event_timepoints'] = True
        self.__events__ = []
        for e in self.__eDict__:
            ev = Event(e, self)
            ev._time_symbol = self.__eDict__[e]['tsymb']
            ev.setTrigger(self.__eDict__[e]['trigger'], self.__eDict__[e][
                'delay'])
            ev.setPriority(self.__eDict__[e]['priority'])
            ev.setPersistent(self.__eDict__[e]['persistent'])
            for ass in self.__eDict__[e]['assignments']:
                ev.setAssignment(ass, self.__eDict__[e]['assignments'][ass])
            self.__events__.append(ev)
            setattr(self, ev.name, ev)
        if len(self.__events__) > 0:
            self.__HAS_EVENTS__ = True
            print('Event(s) detected.')
        else:
            self.__HAS_EVENTS__ = False

    def InitialiseModel(self):
        """
    Initializes and sets up dynamic model attributes and methods based on the model 
    defined in the associated PSC file.

    This method configures various settings, initializes model state variables, 
    precompiles necessary expressions, and prepares the integrator settings. 
    It ensures that the model is ready for simulation by checking and configuring 
    event handling, rate rules, and assignment rules.

    Parameters
    ----------
    None

    Attributes Initialized
    ----------------------
    __inspec__ : numpy.ndarray
        An array initialized to zeros, representing species concentrations.
    __vvec__ : numpy.ndarray
        An array initialized to zeros, representing vector representations.
    mode_solver : str
        Default solver set to 'HYBRD'.
    mode_integrator : str
        Default integrator mode set, switched to 'CVODE' if events, rate rules, 
        or forced functions are detected and if the Assimulo package is available.
    __settings__ : dict
        A dictionary containing various simulation settings, including tolerance 
        levels, integration methods, and output formatting.
    sim_time : numpy.ndarray
        An array representing the time points for simulation.

    Example
    -------
    >>> mod = PysMod()  # Instantiate your model class
    >>> mod.InitialiseModel()  # Initialize the model
    >>> print(mod.mode_integrator)
    'LSODA'  # Mode integrator, switched if conditions are met

    This setup enables PySCeS to efficiently simulate computational models by 
    ensuring that all parameters and arrays are correctly initialized and configured.
    """
        self.__inspec__ = numpy.zeros(len(self.__species__), 'd')
        self._D_s_Order, DvarUpString = self.__initvar__()
        self.__CDvarUpString__ = compile(DvarUpString, 'DvarUpString', 'exec')
        self.state_species = None
        self.state_flux = None
        self.elas_var_u = None
        self.elas_var = None
        self.__vvec__ = numpy.zeros(self.__Nshape__[1], 'd')
        self.__tvec_a__ = None
        self.__tvec_c__ = None
        self._CVODE_Vtemp = numpy.zeros(self.__Nshape__[1], 'd')
        par_map2store, par_remap = self.__initfixed__()
        self.__par_map2storeC = compile(par_map2store, 'par_map2store', 'exec')
        self.__par_remapC = compile(par_remap, 'par_remap', 'exec')
        mapString, mapString_R, vString = self.__initREQ__()
        self.__mapFunc__ = compile(mapString, 'mapString', 'exec')
        self.__mapFunc_R__ = compile(mapString_R, 'mapString_R', 'exec')
        self.__vFunc__ = compile(vString, 'vString', 'exec')
        self.__settings__['mach_floateps'] = mach_eps
        self.__settings__['mode_number_format'] = '%2.4e'
        self.__settings__['mode_sim_init'] = 0
        self.__settings__['mode_sim_max_iter'] = 3
        self.mode_state_init = 0
        self.mode_solver = 'HYBRD'
        self.mode_solver_fallback = 1
        self.STATE_extra_output = []
        self.__settings__['mode_state_init3_factor'] = 0.1
        self.__mode_state_init2_array__ = numpy.logspace(0, 5, 18)
        self.__settings__['mode_state_mesg'] = 1
        self.__settings__['mode_state_nan_on_fail'] = False
        self.__settings__['solver_switch_warning'] = True
        self.__settings__['mode_solver_fallback_integration'] = 1
        self.__settings__['mode_elas_deriv_order'] = 3
        self.__settings__['mode_mca_scaled'] = 1
        self.__settings__['mode_eigen_output'] = 0
        self.__settings__['mode_suppress_info'] = 0
        self.__settings__['write_array_header'] = 1
        self.__settings__['write_array_spacer'] = 1
        self.__settings__['write_array_html_header'] = 1
        self.__settings__['write_array_html_footer'] = 1
        self.__settings__['write_array_html_format'] = '%2.4f'
        self.__settings__['write_arr_lflush'] = 5
        self.__settings__['hybrd_xtol'] = 1e-12
        self.__settings__['hybrd_maxfev'] = 0
        self.__settings__['hybrd_epsfcn'] = copy.copy(self.__settings__[
            'mach_floateps'])
        self.__settings__['hybrd_factor'] = 100
        self.__settings__['hybrd_mesg'] = 1
        self.__settings__['fintslv_tol'] = 0.001
        self.__settings__['fintslv_step'] = 5
        self.__fintslv_range__ = numpy.array([1, 10, 100, 1000, 5000, 10000,
            50000, 50100, 50200, 50300, 50400, 50500, 50600, 50700, 50800, 
            50850, 50900, 50950, 51000], 'd')
        self.__settings__['fintslv_rmult'] = 1.0
        self.__settings__['nleq2_advanced_mode'] = False
        self.__settings__['nleq2_growth_factor'] = 10
        self.__settings__['nleq2_iter'] = 3
        self.__settings__['nleq2_iter_max'] = 10
        self.__settings__['nleq2_rtol'] = 1e-08
        self.__settings__['nleq2_jacgen'] = 2
        self.__settings__['nleq2_iscaln'] = 0
        self.__settings__['nleq2_nonlin'] = 4
        self.__settings__['nleq2_qrank1'] = 1
        self.__settings__['nleq2_qnscal'] = 0
        self.__settings__['nleq2_ibdamp'] = 0
        self.__settings__['nleq2_nitmax_growth_factor'] = 4
        self.__settings__['nleq2_iormon'] = 2
        self.__settings__['nleq2_mesg'] = 1
        self.nleq2_nitmax = 50
        self.__settings__['small_concentration'] = 1e-06
        self.__StateOK__ = True
        self.data_sstate = None
        self._sim = None
        self.data_sim = None
        self.data_stochsim = None
        self.__scan2d_results__ = None
        self.__scan2d_pars__ = None
        self.__settings__['lsoda_atol'] = 1e-12
        self.__settings__['lsoda_rtol'] = 1e-07
        self.__settings__['lsoda_mxstep'] = 0
        self.__settings__['lsoda_h0'] = 0.0
        self.__settings__['lsoda_hmax'] = 0.0
        self.__settings__['lsoda_hmin'] = 0.0
        self.__settings__['lsoda_mxordn'] = 12
        self.__settings__['lsoda_mxords'] = 5
        self.__settings__['lsoda_mesg'] = 1
        self.mode_integrator = 'LSODA'
        if self.__HAS_EVENTS__:
            if _HAVE_ASSIMULO:
                print(
                    """
INFO: events detected and we have Assimulo installed,
switching to CVODE (mod.mode_integrator='CVODE')."""
                    )
                self.mode_integrator = 'CVODE'
            else:
                print(
                    """
WARNING: PySCeS needs Assimulo installed for event handling,
PySCeS will continue with LSODA (NOTE: ALL EVENTS WILL BE IGNORED!)."""
                    )
                print(
                    """Assimulo may be installed from conda-forge or compiled from source.
See: https://jmodelica.org/assimulo"""
                    )
                self.__events__ = []
        if self.__HAS_RATE_RULES__:
            if _HAVE_ASSIMULO:
                print(
                    """INFO: RateRules detected and Assimulo installed,
switching to CVODE (mod.mode_integrator='CVODE').
"""
                    )
                self.mode_integrator = 'CVODE'
            else:
                print(
                    """
WARNING: RateRules detected! PySCeS prefers CVODE but will continue with LSODA
(NOTE: VARIABLE COMPARTMENTS ARE NOT SUPPORTED WITH LSODA!)"""
                    )
                print(
                    """Assimulo may be installed from conda-forge or compiled from source.
See: https://jmodelica.org/assimulo"""
                    )
        if self.__HAS_FORCED_FUNCS__:
            if _HAVE_ASSIMULO:
                print(
                    """INFO: Assignment Rules detected and Assimulo installed,
switching to CVODE (mod.mode_integrator='CVODE').
"""
                    )
                self.mode_integrator = 'CVODE'
            else:
                print(
                    """
WARNING: Assignment Rules detected! PySCeS prefers CVODE but will continue with LSODA
(NOTE: THE VALUES OF ASSIGNMENT RULES DURING THE SIMULATION CANNOT BE TRACKED WITH LSODA!)"""
                    )
                print(
                    """Assimulo may be installed from conda-forge or compiled from source.
See: https://jmodelica.org/assimulo"""
                    )
        self.mode_integrate_all_odes = False
        self.__settings__['cvode_abstol'] = 1e-09
        self.__settings__['cvode_reltol'] = 1e-09
        self.__settings__['cvode_mxstep'] = 5000
        self.__settings__['cvode_stats'] = False
        self.__settings__['cvode_h0'] = 0.0
        self.__settings__['cvode_hmax'] = 0.0
        self.__settings__['cvode_hmin'] = 0.0
        self.__settings__['cvode_mxord'] = 5
        self.__settings__['cvode_access_solver'] = True
        if self.__HAS_PIECEWISE__ and self.__settings__['cvode_reltol'
            ] <= 1e-09:
            self.__settings__['cvode_reltol'] = 1e-06
            print(
                """INFO: Piecewise functions detected increasing CVODE tolerance slightly
(mod.__settings__["cvode_reltol"] = 1.0e-6 )."""
                )
        self.sim_start = 0.0
        self.sim_end = 10.0
        self.sim_points = 20.0
        sim_steps = (self.sim_end - self.sim_start) / self.sim_points
        self.sim_time = numpy.arange(self.sim_start, self.sim_end +
            sim_steps, sim_steps)
        self._CVODE_xdata = None
        self.__SIMPLOT_OUT__ = []
        self.__settings__['elas_evar_upsymb'] = 1
        self.__settings__['elas_epar_upsymb'] = 1
        self.__settings__['elas_evar_remap'] = 1
        self.__settings__['mode_elas_deriv_factor'] = 0.0001
        self.__settings__['mode_elas_deriv_min'] = 1e-12
        self.__settings__['elas_zero_flux_fix'] = False
        self.__settings__['elas_zero_conc_fix'] = False
        self.__settings__['elas_scaling_div0_fix'] = False
        self.__settings__['mca_ccj_upsymb'] = 1
        self.__settings__['mca_ccs_upsymb'] = 1
        self.__settings__['mca_ccall_fluxout'] = 1
        self.__settings__['mca_ccall_concout'] = 1
        self.__settings__['mca_ccall_altout'] = 0
        self.__settings__['pitcon_fix_small'] = 0
        self.pitcon_par_space = numpy.logspace(-1, 3, 10)
        self.pitcon_iter = 10
        self.__settings__['pitcon_allow_badstate'] = 0
        self.__settings__['pitcon_flux_gen'] = 1
        self.__settings__['pitcon_filter_neg'] = 1
        self.__settings__['pitcon_filter_neg_res'] = 0
        self.__settings__['pitcon_max_step'] = 30.0
        self.pitcon_target_points = []
        self.pitcon_limit_points = []
        self.pitcon_targ_val = 0.0
        self.__settings__['pitcon_init_par'] = 1
        self.__settings__['pitcon_par_opt'] = 0
        self.__settings__['pitcon_jac_upd'] = 0
        self.__settings__['pitcon_targ_val_idx'] = 0
        self.__settings__['pitcon_limit_point_idx'] = 0
        self.__settings__['pitcon_output_lvl'] = 0
        self.__settings__['pitcon_jac_opt'] = 1
        self.__settings__['pitcon_abs_tol'] = 1e-05
        self.__settings__['pitcon_rel_tol'] = 1e-05
        self.__settings__['pitcon_min_step'] = 0.01
        self.__settings__['pitcon_start_step'] = 0.3
        self.__settings__['pitcon_start_dir'] = 1.0
        self.__settings__['pitcon_max_grow'] = 3.0
        self.InitialiseConservationArrays()
        self.scan_in = ''
        self.scan_out = []
        self._scan = None
        self.scan_res = None
        self.__settings__['scan1_mca_mode'] = 0
        self.__settings__['scan1_dropbad'] = 0
        self.__settings__['scan1_nan_on_bad'] = True
        self.__settings__['scan1_mesg'] = True
        self.__scan_errors_par__ = None
        self.__scan_errors_idx__ = None

    def InitialiseConservationArrays(self):
        """
    Initialise conservation related vectors/arrays.

    This method initializes vectors and arrays related to conservation laws in
    the model. It was originally part of the `InitialiseModel` method but has
    been separated so that it can be invoked when the stoichiometry of the model
    is reanalyzed. The method also adjusts specific settings related to the model
    configuration and moiety conservation.

    Attributes
    ----------
    __SI__ : numpy.ndarray
        A vector initialized to zeros, with a length equal to the number of columns
        in the zero matrix (`__lzeromatrix__`).
    __SALL__ : numpy.ndarray
        A vector initialized to zeros, with a length equal to the number of species
        in the model.
    __reordered_lcons : numpy.ndarray
        A reordered conservation matrix, only initialized if moiety conservation
        is detected and active in the model.
    
    Examples
    --------
    >>> mod = PysMod()
    >>> mod.InitialiseConservationArrays()
    """
        self.__SI__ = numpy.zeros(self.__lzeromatrix__.shape[1], 'd')
        self.__SALL__ = numpy.zeros(len(self.__species__), 'd')
        self.__settings__['pitcon_max_step'] = 10 * (len(self.__SI__) + 1)
        if self.__HAS_MOIETY_CONSERVATION__ == True:
            idx = [list(self.conservation_matrix_col).index(n) for n in
                range(len(self.conservation_matrix_col))]
            self.__reordered_lcons = self.conservation_matrix.take(idx, 1)
            if self.__HAS_COMPARTMENTS__ and self.__KeyWords__[
                'Species_In_Conc']:
                self.__Build_Tvec__(amounts=False)
            else:
                self.__Build_Tvec__(amounts=True)
            self.showConserved(screenwrite=0)

    def InitialiseOldFunctions(self):
        """
    InitialiseOldFunctions()

    Parses and initializes user-defined functions specified by the `!T` and `!U` 
    markers in the PSC input file. It compiles and executes the user-defined code 
    to prepare the model for execution.

    Notes
    -----
    This method compiles the user-defined functions from the provided PSC input 
    file and executes them. If there are any issues with non-instantiated attributes, 
    it prompts the user to manually resolve these using the `ReloadUserFunc()` method. 
    Additionally, it attempts to initialize functions critical to the model setup 
    via `ReloadInitFunc()`.

    Warnings
    --------
    - If there is an error in executing the user-defined functions, it may be due 
      to non-instantiated attributes like MCA/state attributes.
    - Ensure attributes used in your functions exist and are loaded correctly.

    Examples
    --------
    Initialize user-defined functions within a PysMod instance:

    >>> model = PysMod()
    >>> model.InitialiseOldFunctions()
    
    Ensure user functions are working without errors by using `ReloadUserFunc()` 
    if necessary:

    >>> model.ReloadUserFunc()  # Reload user functions if there were issues
    """
        self.__CODE_user = compile(self._Function_user, '_Function_user',
            'exec')
        try:
            exec(self.__CODE_user)
        except Exception as e:
            print('WARNING: User function error\n', e)
            print(
                'This might be due to non-instantiated (e.g. MCA/state) attributes'
                )
            print(
                'Make sure the attributes that are used in your function exist ...'
                )
            print('and manually load using the self.ReloadUserFunc() method')
        try:
            self.ReloadInitFunc()
        except Exception as e:
            print('WARNING: Init function error\n', e)
            print(
                'This function is meant to be used exclusively for the initialization'
                )
            print(
                'of PySCeS properties and expressions based on model attributes defined in the input file.'
                )
            print('Reinitialize with ReloadInitFunc()')

    def ReloadUserFunc(self):
        """
    Recompile and execute the user function (!U) from the input file.

    This method attempts to parse and re-execute user-defined functions specified
    in the PSC input file using the `!U` directive. It recompiles the function
    code and then executes it. If errors occur during execution, a warning
    message is printed to notify the user of potential issues related
    to non-instantiated attributes.

    Parameters
    ----------
    None

    Examples
    --------
    Assuming that you have a user-defined function specified with `!U`
    in your PSC input file, you would typically use this method if you
    have made changes to the function and need to recompile and reload it::

        model = PysMod()
        # After making changes to the user function in the input file
        model.ReloadUserFunc()

    In case the user function has dependencies or requires certain model 
    state attributes, ensure these attributes are initialized before 
    executing the function.
    """
        self.__CODE_user = compile(self._Function_user, '_Function_user',
            'exec')
        try:
            exec(self.__CODE_user)
        except Exception as e:
            print('WARNING: User function load error\n', e)

    def ReloadInitFunc(self):
        """
    Recompile and execute the user initialisations (!I) as defined in the PSC input file and in `mod.__InitFuncs__`.

    This method recompiles and executes user-defined initializations found in the PSC input file, allowing
    for dynamic setup based on the user's model specifications. As of 2015, it also supports defining
    InitialAssignments without requiring a `self.*` prefix in the input file. This facilitates the customization 
    of model properties and state variables upon initialization.

    The method processes and sorts the initialization tasks to respect any dependencies between them, ensuring 
    that cyclic dependencies are raised as errors. It distinguishes between simple assignments (e.g., constants) 
    and more complex expressions (e.g., evaluated at runtime).

    Examples
    --------
    Suppose your PSC input file defines initializations with placeholders for model-specific expressions 
    or state assignments. Use this method to enforce those initial setup routines:

    >>> model = PysMod()
    >>> model.ReloadInitFunc()

    This will execute the initializations as defined, providing the necessary setup for further simulations 
    or analyses.

    Raises
    ------
    RuntimeError
        If a cyclic dependency is detected among the initializations.
    
        Repeatedly go through all of the nodes in the graph, moving each of
        the nodes that has all its edges resolved, onto a sequence that
        forms our sorted graph. A node has all of its edges resolved and
        can be moved once all the nodes its edges point to, have been moved
        from the unsorted graph onto the sorted one.
        """

        def topolgical_sort(graph_unsorted):
            """
            Repeatedly go through all of the nodes in the graph, moving each of
            the nodes that has all its edges resolved, onto a sequence that
            forms our sorted graph. A node has all of its edges resolved and
            can be moved once all the nodes its edges point to, have been moved
            from the unsorted graph onto the sorted one.
            """
            graph_sorted = []
            graph_unsorted = dict(graph_unsorted)
            while graph_unsorted:
                acyclic = False
                for node, edges in list(graph_unsorted.items()):
                    for edge in edges:
                        if edge in graph_unsorted:
                            break
                    else:
                        acyclic = True
                        del graph_unsorted[node]
                        graph_sorted.append((node, edges))
                if not acyclic:
                    raise RuntimeError('A cyclic dependency occurred')
            return graph_sorted
        simpleAss = []
        exprsAss = []
        symbolsX = []
        for init in self.__InitFuncs__:
            if hasattr(self, init):
                if type(self.__InitFuncs__[init]) == str:
                    if self.__InitFuncs__[init].isdigit():
                        simpleAss.append((init, eval(self.__InitFuncs__[init]))
                            )
                    else:
                        InfixParser.setNameStr('self.', '')
                        InfixParser.parse(str(self.__InitFuncs__[init]))
                        piecewises = InfixParser.piecewises
                        symbols = InfixParser.names
                        for s_ in symbols:
                            if s_ not in symbolsX:
                                symbolsX.append(s_)
                        if init not in symbolsX:
                            symbolsX.append(init)
                        functions = InfixParser.functions
                        code_string = 'self.{} = {}'.format(init,
                            InfixParser.output)
                        xcode = compile(code_string, '_InitAss_', 'exec')
                        exprsAss.append((init, xcode, symbols))
                else:
                    setattr(self, init, self.__InitFuncs__[init])
            else:
                assert hasattr(self, init), '\nModelInit error'
        for init, val in simpleAss:
            setattr(self, init, val)
        execOrder = []
        for init, xcode, symbols2 in exprsAss:
            symbols2 = [symbolsX.index(s_) for s_ in symbols2]
            execOrder.append((symbolsX.index(init), symbols2))
            eval(xcode)

    def __initREQ__(self):
        """
    Compile and initialize the model rate equations and associated arrays. 
    
    This method constructs the rate equations as lambda functions, allowing 
    them to be callable as `mod.R1(mod)`. It processes species and reactions, 
    potentially using substitution for assignment rules if configured. It is also 
    responsible for warning if reactions are not associated with any compartments, 
    given certain settings.

    If `display_debug` setting is enabled, additional debugging information will 
    be printed to the console, including mapping strings for species and reactions.

    Returns
    -------
    tuple of str
        A tuple containing three strings:
        - `mapString`: Maps species to indices.
        - `mapString_R`: Initializes species from their initial states.
        - `vString`: Creates temporary variables for reaction rates.

    Examples
    --------
    >>> mod = PysMod()
    >>> map_str, map_str_r, v_str = mod.__initREQ__()
    
    In this example, `mod` is an instance of `PysMod`. The `__initREQ__` method 
    initializes the rate equations for the model. The method returns strings 
    representing the mappings and variable creations crucial for further simulation 
    processes. Debugging information will be included in the output if the 
    `display_debug` setting is enabled during the execution.
    """
        if self.__settings__['display_debug'] == 1:
            pass
        mapString = ''
        mapString_R = ''
        for x in range(0, len(self.__species__)):
            mapString += 'self.{} = s[{}]\n'.format(self.__species__[x], x)
            mapString_R += 'self.__inspec__[{}] = self.{}_init\n'.format(x,
                self.__species__[x])
        if self.__settings__['display_debug'] == 1:
            print('mapString')
            print(mapString)
            print('mapString_R')
            print(mapString_R)
        vString = ''
        EStat = []
        symbR = {}
        if len(list(self.__rules__.keys())) > 0 and self.__settings__[
            'mode_substitute_assignment_rules']:
            for ass in self.__rules__:
                symbR.update({self.__rules__[ass]['name']: self.__rules__[
                    ass]['code_string']})
        for x in range(0, len(self.__reactions__)):
            req1 = self.__nDict__[self.__reactions__[x]]['RateEq']
            vString += 'Vtemp[{}] = self.{}()\n'.format(x, self.
                __reactions__[x])
            if len(list(self.__rules__.keys())) > 0 and self.__settings__[
                'mode_substitute_assignment_rules']:
                InfixParser.setNameStr('self.', '')
                InfixParser.FunctionReplacements = symbR
                InfixParser.parse(req1.replace('self.', ''))
                req2 = InfixParser.output
                rObj = ReactionObj(self, self.__reactions__[x], req2, 'self.')
            else:
                rObj = ReactionObj(self, self.__reactions__[x], req1, 'self.')
            self.__ParsePiecewiseFunctions__(rObj.piecewises)
            rObj.compartment = self.__nDict__[rObj.name]['compartment']
            setattr(self, self.__reactions__[x], rObj)
        if self.__HAS_COMPARTMENTS__:
            cnames = [self.__compartments__[c]['name'] for c in self.
                __compartments__]
            warnings = ''
            for rr in self.__reactions__:
                rrobj = getattr(self, rr)
                warn = False
                if rrobj.compartment is None and self.__settings__[
                    'display_compartment_warnings']:
                    warnings += ('# {} is not located in a compartment.\n'.
                        format(rrobj.name))
                else:
                    for comp in cnames:
                        if comp in rrobj.symbols:
                            warn = False
                            break
                        else:
                            warn = True
                    if warn:
                        warnings += (
                            '# {}: {}\n#  assuming kinetic constants are flow constants.\n'
                            .format(rr, rrobj.formula))
            if warnings != '' and self.__settings__[
                'display_compartment_warnings']:
                print('\n# -- COMPARTMENT WARNINGS --')
                print(warnings)
        if self.__settings__['display_debug'] == 1:
            print('vString')
            print(vString)
        return mapString, mapString_R, vString

    def __initmodel__(self):
        """
    Generate the stoichiometric matrix N from the parsed model description.

    This method generates the stoichiometric matrix `N` based on the model's
    parsed description. The matrix represents the relationship between
    species and reactions within the model.

    Returns
    -------
    numpy.ndarray
        A stoichiometric matrix `N` where each row corresponds to a species
        and each column corresponds to a reaction.

    Examples
    --------
    >>> model = PysMod()
    >>> stoichiometric_matrix = model.__initmodel__()
    >>> print(stoichiometric_matrix)
    [[-1.  1.  0.]
     [ 1. -1.  0.]
     [ 0.  0.  1.]]

    Notes
    -----
    The stoichiometric matrix is constructed such that each entry in the matrix
    maps a species to a reaction. Negative values in the matrix indicate 
    consumption of a species, while positive values indicate production.

    See Also
    --------
    PysMod: Class representing the model.
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

    def __initvar__(self):
        """
    Compile and initialize the model variable species and derivatives, along with associated mapping arrays.

    This method sets up internal mappings for species variables and their derivatives, initializing
    them based on the species defined in the model. It generates strings for variable initialization
    and updates, which can be executed at runtime to update the state of species within the model.

    Parameters
    ----------
    None

    Returns
    -------
    tuple
        A tuple containing:
        - sOrder : list of str
            An ordered list of species variable names, prefixed with 'self.'.
        - DvarUpString : str
            A string representing Python code to update species variables with input values.

    Examples
    --------
    Suppose `mod` is an instance of the `PysMod` class, and the species in the model are defined as 'A', 'B', and 'C'.
    
    >>> sOrder, DvarUpString = mod.__initvar__()
    >>> print(sOrder)
    ['self.A', 'self.B', 'self.C']
    >>> print(DvarUpString)
    self.A = input[0]
    self.B = input[1]
    self.C = input[2]
    
    The `DvarUpString` can be executed within the context of `mod` to update the model's species variables.
    """
        mvarString = ''
        sOrder = []
        self.__remaps = ''
        DvarUpString = ''
        for x in range(0, len(self.__species__)):
            key = self.__species__[x]
            self.__inspec__[x] = getattr(self, key)
            mvarString += 'self.__inspec__[{}] = float(self.{})\n'.format(x,
                key)
            self.__remaps += 'self.{} = self.{}_ss\n'.format(key, key)
            sOrder.append('self.' + key)
            DvarUpString += 'self.{} = input[{}]\n'.format(self.__species__
                [x], x)
        if self.__settings__['display_debug'] == 1:
            print('s_initDeriv')
            print(s_initDeriv)
            print('\n__remaps')
            print(self.__remaps)
            print('\nDvarUpString')
            print(DvarUpString)
        mvarFunc = compile(mvarString, 'mvarString', 'exec')
        return sOrder, DvarUpString

    def __initfixed__(self):
        """
    Compile and initialise the fixed species and associated mapping arrays.

    This method generates the mapping and remapping strings for the fixed species,
    based on the defined parameters in the model. It returns two strings: one for
    mapping the parameter values to a temporary storage array, and another for
    remapping those values back to the class attributes.

    Parameters
    ----------
    None

    Returns
    -------
    par_map2store : str
        A string of Python code to map current parameter values to a temporary
        storage array `parVal_hold`.
    par_remap : str
        A string of Python code to remap values from `parVal_hold` back to the
        original class attributes.

    Example
    -------
    >>> model = PysMod()
    >>> par_map2store, par_remap = model.__initfixed__()
    >>> print(par_map2store)
    parVal_hold[0] = self.parameter1
    parVal_hold[1] = self.parameter2
    ...

    >>> print(par_remap)
    self.parameter1 = parVal_hold[0]
    self.parameter2 = parVal_hold[1]
    ...
    """
        runmapString = ''
        if self.__settings__['display_debug'] == 1:
            print('InitParams2')
            print(self.__parameters__)
            print('InitStrings2')
        par_map2store = ''
        par_remap = ''
        for x in range(len(self.__parameters__)):
            par_map2store += 'parVal_hold[{}] = self.{}\n'.format(x, self.
                __parameters__[x])
            par_remap += 'self.{} = parVal_hold[{}]\n'.format(self.
                __parameters__[x], x)
        if self.__settings__['display_debug'] == 1:
            print('par_map2store')
            print(par_map2store)
            print('par_remap')
            print(par_remap)
        return par_map2store, par_remap

    def _SpeciesAmountToConc(self, s):
        """
    Convert species amount to concentration.

    This method converts the given species amount `s` to concentration using
    the compartment sizes specified in the `__CsizeAll__` array of the class instance.
    It first updates `__CsizeAll__` with the current values from the attributes listed
    in `__CsizeAllIdx__`. The species amount is then divided by these compartment sizes
    to yield the concentrations.

    Parameters
    ----------
    s : numpy.ndarray or float
        The species amount to be converted into concentration.

    Returns
    -------
    numpy.ndarray or float
        The concentration of the species, calculated based on the compartment sizes.

    Example
    -------
    >>> model = PysMod()
    >>> model.__CsizeAllIdx__ = ['compartment_size_1', 'compartment_size_2']
    >>> model.compartment_size_1 = 2.0
    >>> model.compartment_size_2 = 3.0
    >>> amount = numpy.array([10, 15])
    >>> concentration = model._SpeciesAmountToConc(amount)
    >>> print(concentration)
    [5. 5.]

    Notes
    -----
    This method assumes that the compartment sizes have already been set to 
    the correct values and updated prior to calling the method.
    """
        for x in range(len(self.__CsizeAllIdx__)):
            self.__CsizeAll__[x] = getattr(self, self.__CsizeAllIdx__[x])
        return s / self.__CsizeAll__

    def _SpeciesConcToAmount(self, s):
        """
    Converts a concentration value to an amount by using internal size attributes.

    This method multiplies the input concentration by internal size scaling factors
    to convert it into an amount. The size scaling factors are dynamically fetched
    from the object instance.

    Parameters
    ----------
    s : float or array_like
        The concentration value(s) that need to be converted to amount(s).

    Returns
    -------
    float or array_like
        The amount(s) corresponding to the input concentration(s), scaled by
        the internal size attributes.

    Notes
    -----
    The method accesses and uses internal attributes `__CsizeAllIdx__` and
    `__CsizeAll__`. It assumes these are pre-defined and properly configured
    within the object.

    Examples
    --------
    >>> pys_mod_instance = PysMod()
    >>> concentration = 10.0
    >>> amount = pys_mod_instance._SpeciesConcToAmount(concentration)
    >>> print(amount)
    [10.0, 20.0, 30.0]  # Example output assuming `__CsizeAll__` = [1, 2, 3]
    """
        for x in range(len(self.__CsizeAllIdx__)):
            self.__CsizeAll__[x] = getattr(self, self.__CsizeAllIdx__[x])
        return s * self.__CsizeAll__

    def _FixedSpeciesAmountToConc(self):
        """
    Converts the amount of fixed species to concentration.

    This method iterates over the list of fixed species and adjusts their values 
    from amount to concentration based on the corresponding scaling factor. The 
    result replaces the original fixed species amount with its concentration in 
    the object state.

    Parameters
    ----------
    None

    Attributes
    ----------
    __fixed_species__ : list of str
        List of attribute names representing the fixed species whose amounts are 
        to be converted into concentrations.

    __CsizeFS__ : list of float
        List of scaling factors corresponding to each fixed species. Used for 
        converting amounts into concentrations.
        
    __CsizeFSIdx__ : list of str
        List of attribute names representing the scaling factor for each fixed 
        species.

    Example
    -------
    Assume `self` is an instance of `PysMod` with fixed species 'species_a', 
    'species_b' and corresponding scaling factors:
    
    >>> self.__fixed_species__ = ['species_a', 'species_b']
    >>> self.species_a = 10.0  # amount
    >>> self.species_b = 15.0  # amount
    >>> self.__CsizeFSIdx__ = ['scale_a', 'scale_b']
    >>> self.scale_a = 2.0
    >>> self.scale_b = 3.0

    After calling `_FixedSpeciesAmountToConc`, the object attributes will be:

    >>> self._FixedSpeciesAmountToConc()
    >>> self.species_a 
    5.0  # concentration
    >>> self.species_b
    5.0  # concentration

    Notes
    -----
    This method modifies the instance attributes in place.
    """
        for x in range(len(self.__fixed_species__)):
            am = getattr(self, self.__fixed_species__[x])
            self.__CsizeFS__[x] = getattr(self, self.__CsizeFSIdx__[x])
            setattr(self, self.__fixed_species__[x], am / self.__CsizeFS__[x])

    def _FixedSpeciesConcToAmount(self):
        """Convert fixed species concentrations to amounts.

This method updates the internal attributes corresponding to a predefined list
of fixed species by converting their concentrations to amounts. The conversion
uses the compartment sizes specified in `self.__CsizeFS__` and `self.__CsizeFSIdx__`.

The conversion is performed as follows:
    amount = concentration * compartment_size

Notes
-----
- It is assumed that `self.__fixed_species__` contains attribute names
  corresponding to the fixed species concentrations to be converted.
- `self.__CsizeFS__` and `self.__CsizeFSIdx__` should map correctly to the
  compartment sizes and their indices to perform the conversion properly.

Examples
--------
>>> class PysMod:
...     def __init__(self):
...         self.species_A_conc = 2.0
...         self.species_B_conc = 3.5
...         self.__fixed_species__ = ['species_A_conc', 'species_B_conc']
...         self.__CsizeFS__ = [0.0, 0.0]
...         self.__CsizeFSIdx__ = ['compartment_A_size', 'compartment_B_size']
...         self.compartment_A_size = 10.0
...         self.compartment_B_size = 20.0
...
...     def _FixedSpeciesConcToAmount(self):
...         for x in range(len(self.__fixed_species__)):
...             cn = getattr(self, self.__fixed_species__[x])
...             self.__CsizeFS__[x] = getattr(self, self.__CsizeFSIdx__[x])
...             setattr(self, self.__fixed_species__[x], cn * self.__CsizeFS__[x])
...
>>> pysmod_instance = PysMod()
>>> pysmod_instance._FixedSpeciesConcToAmount()
>>> pysmod_instance.species_A_conc
20.0
>>> pysmod_instance.species_B_conc
70.0
"""
        for x in range(len(self.__fixed_species__)):
            cn = getattr(self, self.__fixed_species__[x])
            self.__CsizeFS__[x] = getattr(self, self.__CsizeFSIdx__[x])
            setattr(self, self.__fixed_species__[x], cn * self.__CsizeFS__[x])

    def _EvalREq2_alt(self, s, Vtemp):
        """
    Evaluate the rate equations with mass conservation correction.

    This method adjusts the full species conservation vector `s` with
    respect to moiety awareness and then evaluates the rate equations.
    It modifies the species vector `s` based on internal moiety balances
    before proceeding with the rate evaluation. 

    Parameters
    ----------
    s : array_like
        A vector representing species conservation. This should be a list
        or array where each entry corresponds to a species involved in the
        reaction network.
    
    Vtemp : array_like
        The rate vector. This typically represents the temporary or transient
        rates associated with each step or component of the reaction network.

    Returns
    -------
    array_like
        The updated species vector after adjusting for moiety conservation 
        and evaluating the reaction rates using the provided rate vector.

    Example
    -------
    >>> model = PysMod()
    >>> s = [1.0, 2.0, 0.5, 3.0]  # Example species vector
    >>> Vtemp = [0.1, 0.3, 0.05, 0.2]  # Example rate vector
    >>> updated_s = model._EvalREq2_alt(s, Vtemp)
    >>> print(updated_s)
    [0.85, 1.9, 0.5, 2.7]  # Hypothetical updated values

    Notes
    -----
    This function requires that the instance of `PysMod` has correctly 
    initialized its internal moiety properties (like `__SI__`, `__tvec_a__`, 
    and others) before calling this method.

    """
        for x in range(0, len(self.lzeromatrix_col)):
            self.__SI__[x] = s[self.lzeromatrix_col[x]]
        for x in range(0, len(self.lzeromatrix_row)):
            s[self.lzeromatrix_row[x]] = self.__tvec_a__[x] + numpy.add.reduce(
                self.__lzeromatrix__[x, :] * self.__SI__)
        return self._EvalREq(s, Vtemp)

    def _EvalREq(self, s, Vtemp):
        """
    Evaluate the rate equations for the given species and rate vector.

    This method evaluates the rate equations by first checking for the presence of rate rules 
    and compartment-specific details. It dynamically executes mapped functions and attempts to 
    correct for any errors during the computation process. Any errors encountered during the execution
    of rate functions are caught and printed, and the rate vector (`Vtemp`) is updated accordingly.

    Parameters
    ----------
    s : array-like
        A vector representing the species concentrations or amounts depending on system requirements.
    Vtemp : array-like
        A pre-allocated vector for storing the computed rates.

    Returns
    -------
    Vtemp : array-like
        The updated rate vector after evaluating the rate equations.

    Raises
    ------
    Exception
        General exceptions are caught during the forcing function execution.
    ArithmeticError, AttributeError, NameError, ZeroDivisionError, ValueError
        Specific errors during evaluations that may occur with rate rules 
        or function mappings.

    Examples
    --------
    >>> pys_mod_instance = PysMod()
    >>> species_vector = [1.0, 0.5, 0.3]  # Example species concentrations or amounts
    >>> rate_vector = [0.0, 0.0, 0.0]     # Pre-allocated rate vector
    >>> result = pys_mod_instance._EvalREq(species_vector, rate_vector)
    >>> print(result)
    [<computed_rate_1>, <computed_rate_2>, <computed_rate_3>]
    """
        if self.__HAS_RATE_RULES__:
            exec(self._CODE_raterule_map)
        if self.__HAS_COMPARTMENTS__ and self.__KeyWords__['Species_In_Conc']:
            s = self._SpeciesAmountToConc(s)
        exec(self.__mapFunc__)
        try:
            self.Forcing_Function()
        except Exception as de:
            print('INFO: forcing function failure', de)
        try:
            exec(self.__vFunc__)
        except (ArithmeticError, AttributeError, NameError,
            ZeroDivisionError, ValueError) as detail:
            print('INFO: REq evaluation failure:', detail)
            Vtemp[:] = self.__settings__['mach_floateps']
        if self.__HAS_RATE_RULES__:
            try:
                exec(self._CODE_raterule)
            except (ArithmeticError, AttributeError, NameError,
                ZeroDivisionError, ValueError) as detail:
                print('INFO: RateRule evaluation failure:', detail)
        return Vtemp

    def _EvalREq2(self, s, Vtemp):
        """
    Evaluate the rate equations with mass conservation correction.

    This method evaluates the rate equations considering mass conservation 
    and regenerates the full species vector. PySCeS uses a reduced set of 
    ordinary differential equations (ODEs) for its core operations. This method 
    accepts a reduced species vector and returns the full species vector.

    Parameters
    ----------
    s : array-like
        The reduced species vector (concentrations or amounts).
    Vtemp : array-like
        The rate vector to be updated with the calculated values.

    Returns
    -------
    array-like
        The updated rate vector after evaluating the rate equations.

    Notes
    -----
    This method expands the reduced species vector into a full species vector 
    by ensuring mass balance through the use of pre-defined stoichiometric and 
    linkage matrix mappings.

    Examples
    --------
    >>> # Assuming 'pysmod_instance' is an instance of 'PysMod'
    >>> s = [0.1, 0.2, 0.3]  # Example reduced species vector
    >>> Vtemp = [0, 0, 0]    # Initial rate vector
    >>> updated_Vtemp = pysmod_instance._EvalREq2(s, Vtemp)
    >>> print(updated_Vtemp)
    [0.05, 0.04, 0.07]  # Example output; actual values depend on internal model setup
    """
        for x in range(len(s)):
            self.__SALL__[self.nrmatrix_row[x]] = s[x]
        for x in range(len(self.lzeromatrix_row)):
            self.__SALL__[self.lzeromatrix_row[x]] = self.__tvec_a__[x
                ] + numpy.add.reduce(self.__lzeromatrix__[x, :] * s)
        return self._EvalREq(self.__SALL__, Vtemp)

    def _EvalODE(self, s, Vtemp):
        """
    _EvalODE(s, Vtemp)

    Core ODE evaluation routine for evaluating the reduced set of ordinary 
    differential equations (ODEs). This method considers whether mass 
    conservation is present to decide between the use of either the full system 
    (N*v) or the reduced system (Nr*v) of equations. Depending on the situation, 
    either `_EvalREq` or `_EvalREq2` is invoked. The method returns the rate 
    of change of the species vector `sdot`.

    Parameters
    ----------
    s : array_like
        The species vector representing the concentrations or amounts of 
        each species in the system.
        
    Vtemp : array_like
        The temporary rate vector used during the evaluation of reaction rates.

    Returns
    -------
    sdot : ndarray
        The derivatives of the species concentrations, representing the rate 
        of change of each species with respect to time.

    Example
    -------
    Suppose you have a `PysMod` object initialized, and species and rate vectors:

    >>> model = PysMod()
    >>> species_vector = numpy.array([1.0, 2.0, 0.5])  # Example species concentrations
    >>> rate_vector = numpy.zeros(len(species_vector))  # Initialize rate vector
    >>> derivatives = model._EvalODE(species_vector, rate_vector)
    >>> print(derivatives)
    array([...])  # The result will depend on the model's equations and parameters

    Notes
    -----
    Ensure that the model has been properly initialized and that all the 
    required pre-computation steps are completed before calling this method. 
    This includes setting up the reaction rate rules and any mass conservation 
    constraints that apply to the model. The returned `sdot` can be used 
    in numerical integration methods to simulate the system dynamics over time.
    """
        if self.__HAS_RATE_RULES__:
            s, self.__rrule__ = numpy.split(s, [self.Nrmatrix.shape[0]])
        if self.__HAS_MOIETY_CONSERVATION__ == True:
            for x in range(len(s)):
                self.__SALL__[self.nrmatrix_row[x]] = s[x]
            for x in range(len(self.lzeromatrix_row)):
                self.__SALL__[self.lzeromatrix_row[x]] = self.__tvec_a__[x
                    ] + numpy.add.reduce(self.__lzeromatrix__[x, :] * s)
            self.__vvec__ = self._EvalREq(self.__SALL__, Vtemp)
            s = numpy.add.reduce(self.__nrmatrix__ * self.__vvec__, 1)
        else:
            self.__vvec__ = self._EvalREq(s, Vtemp)
            s = numpy.add.reduce(self.__nmatrix__ * self.__vvec__, 1)
        if self.__HAS_RATE_RULES__:
            s = numpy.concatenate([s, self.__rrule__]).copy()
        return s

    def _EvalODE_CVODE(self, s, Vtemp):
        """
    Evaluate the full set of ordinary differential equations (ODEs).

    This method evaluates the entire set of ODEs by utilizing either the reaction rate matrix `N*v`
    via `_EvalREq` or the reduced reaction rate matrix `Nr*v` via `_EvalREq2_alt` based on the presence
    of moiety conservation. If rate rules are applicable, the input vector `s` is modified to include 
    rate rules.

    Parameters
    ----------
    s : ndarray
        Species vector representing the concentrations of chemical species in the system.
    Vtemp : ndarray
        Rate vector representing the current rates of reaction.

    Returns
    -------
    ndarray
        The computed derivative of the species vector, `sdot`, representing the rates of change
        of the species concentrations.

    Example
    -------
    Create an instance of `PysMod` and evaluate ODEs:

    >>> model = PysMod()
    >>> species_vector = numpy.array([1.0, 0.5, 0.3])
    >>> rate_vector = numpy.array([0.0, 0.0])
    >>> sdot = model._EvalODE_CVODE(species_vector, rate_vector)
    >>> print(sdot)
    [rate_of_change_species_1, rate_of_change_species_2, rate_of_change_species_3]
    
    Notes
    -----
    Internally, the method handles the presence of rate rules by splitting the species vector and 
    concatenating it back after computations. The decision on whether to use the full or reduced 
    reaction matrices is based on the presence of moiety conservation.
    """
        if self.__HAS_RATE_RULES__:
            s, self.__rrule__ = numpy.split(s, [self.Nmatrix.shape[0]])
        self.__vvec__ = self._EvalREq(s, Vtemp)
        s = numpy.add.reduce(self.__nmatrix__ * self.__vvec__, 1)
        if self.__HAS_RATE_RULES__:
            s = numpy.concatenate([s, self.__rrule__]).copy()
        return s

    def Forcing_Function(self):
        """
    Forcing_Function()

    User-defined forcing function that is either specified in the PSC input file as !F or by overwriting this method.
    This method is evaluated before every rate equation evaluation, allowing dynamic adjustment of any state or rate 
    parameters as needed.

    Parameters
    ----------
    None

    Example
    -------
    Define a custom forcing function in a subclass:

    >>> class CustomPysMod(PysMod):
    ...     def Forcing_Function(self):
    ...         # Custom logic for adjusting parameters dynamically
    ...         self.some_parameter += 1
    ...         print("Forcing function executed.")
    ...
    >>> model = CustomPysMod()
    >>> model.Forcing_Function()
    Forcing function executed.

    Notes
    -----
    This method executes the self.__CODE_forced script, which can modify model state variables or parameters
    for computational experiments or simulations that require parameter adjustment at every evaluation step.
    """
        exec(self.__CODE_forced)

    def User_Function(self):
        """
    Execute a user-defined function specified within the class instance.

    **Deprecated**

    This method is deprecated and may be removed in future versions. It currently executes a block
    of user-defined code stored in the `self.__CODE_user` attribute. Users are encouraged to use
    alternative, more robust mechanisms for custom code execution within their models.

    Notes
    -----
    - The functionality provided by this method can be replaced by using the `Forcing_Function` or
      other similar appropriate methods.
    - Directly executing code via `exec` poses significant security risks. It is recommended to 
      validate or sanitize any input before use.

    Examples
    --------
    Define a block of custom Python code and execute it:

    >>> class MyPysMod(PysMod):
    ...     def __init__(self, custom_code):
    ...         self.__CODE_user = custom_code
    ...
    >>> model = MyPysMod("print('Executing user-defined function')")
    >>> model.User_Function()
    Executing user-defined function

    Warnings
    --------
    - As this method uses Python's `exec()` for code execution, it's important to ensure that the 
      code passed to `self.__CODE_user` is safe to execute.
    """
        exec(self.__CODE_user)

    def savePSC(self, filename=None, directory=None, iValues=True):
        """
    Saves the model object to a PSC file.

    This method writes the current state of the model to a PSC (Portable Storage Container) file, 
    which can be used for later retrieval or sharing. The file can be saved under a specified 
    filename and location, and the user can choose whether to record the model's initial values 
    or the current state values.

    Parameters
    ----------
    filename : str, optional
        The desired name of the PSC file without extension. If not provided, the model's name 
        will be used with the '.psc' extension.
    directory : str, optional
        The directory path where the PSC file will be stored. If not specified, the current 
        working directory is used.
    iValues : bool, default=True
        Determines which values are saved in the PSC file:
        - If True, the initial values of the model are recorded.
        - If False, the current values of the model are recorded.

    Examples
    --------
    Save the model using a specific filename in the current directory with initial values:

    >>> model = PysMod()
    >>> model.savePSC(filename='my_model')

    Save the model using its default name in a specified directory with current values:

    >>> model = PysMod()
    >>> model.savePSC(directory='/path/to/save', iValues=False)

    Save the model using the provided filename within a custom directory:

    >>> model = PysMod()
    >>> model.savePSC(filename='custom_model', directory='/path/to/directory')
    """
        interface.writeMod2PSC(self, filename, directory, iValues)

    def exportSBML(self, filename=None, directory=None, iValues=True):
        """
    Exports the model object to an SBML file.

    Parameters
    ----------
    filename : str, optional
        The desired name of the SBML file, without extension. If `None`, the 
        model's name will be used as the filename. The file will have the 
        `.xml` extension.
    directory : str, optional
        The directory where the SBML file will be saved. If not specified, 
        the current working directory is used.
    iValues : bool, optional
        Determines whether to export the model's initial values or the current 
        values. If `True`, initial values are used; if `False`, current values 
        are used. Default is `True`.

    Examples
    --------
    Exporting a model using the default parameters:
    
    >>> model = PysMod()
    >>> model.exportSBML()

    Exporting a model with a specific filename and to a specified directory:

    >>> model.exportSBML(filename='custom_model', directory='/path/to/export')

    Exporting a model using current values instead of initial values:

    >>> model.exportSBML(iValues=False)

    Notes
    -----
    This function requires an interface function `writeMod2SBML` that 
    performs the actual export operation.
    """
        interface.writeMod2SBML(self, filename, directory, iValues)

    def _EvalExtraData(self, xdata):
        """
    Evaluates a list of model attributes and returns an array of their values.

    Parameters
    ----------
    xdata : list of str
        A list of attribute names (as strings) that belong to the model object.
        Each attribute must be accessible via `getattr`.

    Returns
    -------
    numpy.ndarray
        An array containing the values of the specified attributes in the same order
        they appear in the `xdata` list.

    Raises
    ------
    AttributeError
        If any attribute in `xdata` does not exist in the model object.

    Examples
    --------
    Suppose you are working with an instance of `PysMod` class, and you want to
    collect the current values of attributes named `attr1`, `attr2`, and `attr3`.

    >>> model = PysMod()
    >>> model.attr1 = 10
    >>> model.attr2 = 20
    >>> model.attr3 = 30
    >>> xdata = ['attr1', 'attr2', 'attr3']
    >>> values = model._EvalExtraData(xdata)
    >>> print(values)
    array([10, 20, 30])

    """
        return numpy.array([getattr(self, d) for d in xdata])
    __CVODE_initialise__ = True
    CVODE_continuous_result = None
    __CVODE_initial_num__ = None

    def CVODE_continue(self, tvec):
        """
    **Experimental:** Continues a simulation over a new time vector. The
    CVODE memory object is reused and not reinitialized. Model parameters can be
    changed between calls to this method. The `mod.data_sim` objects from
    the initial simulation and all calls to this method are stored in the list
    `mod.CVODE_continuous_result`.

    Parameters
    ----------
    tvec : numpy.ndarray
        A numpy array of time points for the continuation of the simulation.

    Notes
    -----
    This method is experimental and requires CVODE to be used as the integration
    algorithm (i.e., `self.mode_integrator` must be 'CVODE').

    Upon simulation continuation, the method records simulation time and prints
    out the time taken for simulation points. It also updates the simulation
    data and checks for success, printing 'Simulation failure' if unsuccessful.

    Examples
    --------
    Continue a simulation with a new time vector:

    >>> import numpy as np
    >>> model = PysMod()  # Assume PysMod is properly initialized elsewhere
    >>> new_time_vector = np.linspace(10, 100, num=1000)
    >>> model.CVODE_continue(new_time_vector)

    If successful, the new simulation results are placed in `model.CVODE_continuous_result`.

    
    For what should be rather obvious reasons, this method requires CVODE to be used as the default integration algorithm.
    """
        assert self.mode_integrator == 'CVODE', """
For what should be rather obvious reasons, this method requires CVODE to be used as the default integration algorithm.
"""
        if self.CVODE_continuous_result is None:
            self.CVODE_continuous_result = [self.data_sim]
        self.__CVODE_initialise__ = False
        self.sim_time = tvec
        Tsim0 = time.time()
        sim_res, rates, simOK = self.CVODE(None)
        Tsim1 = time.time()
        print('{} time for {} points: {}'.format(self.mode_integrator, len(
            self.sim_time), Tsim1 - Tsim0))
        if self.__HAS_RATE_RULES__:
            sim_res, rrules = numpy.split(sim_res, [len(self.__species__)],
                axis=1)
            print('RateRules evaluated and added to mod.data_sim.')
        self.data_sim = IntegrationDataObj()
        self.IS_VALID = simOK
        self.data_sim.setTime(self.sim_time)
        self.data_sim.setSpecies(sim_res, self.__species__)
        self.data_sim.setRates(rates, self.__reactions__)
        if self.__HAS_RATE_RULES__:
            self.data_sim.setRules(rrules, self.__rate_rules__)
        if len(self._CVODE_extra_output) > 0:
            self.data_sim.setXData(self._CVODE_xdata, lbls=self.
                _CVODE_extra_output)
            self._CVODE_xdata = None
        if not simOK:
            print('Simulation failure')
        del sim_res
        self.CVODE_continuous_result.append(self.data_sim)
        self.__CVODE_initialise__ = True

    def CVODE(self, initial):
        """
    Interface to the CVODE integration algorithm.

    This method integrates the set of ordinary differential equations
    specified by the model using the CVODE solver from Assimulo. The
    integration is performed starting from the given initial conditions.

    Parameters
    ----------
    initial : array-like
        Vector containing the initial concentrations of the species.

    Returns
    -------
    sim_res : ndarray
        Array containing the simulation results, which includes the temporal
        evolution of the species concentrations.
    rates : ndarray
        Array containing the computed reaction rates at each time point.
    success : bool
        Boolean flag indicating whether the integration was successful.

    Raises
    ------
    AssertionError
        If Assimulo is not installed or did not import correctly.

    Notes
    -----
    - The CVODE solver allows for complex models involving conservation laws
      and rate rules.
    - The method can handle events and moiety conservation, and it returns
      detailed integration statistics if enabled.

    Examples
    --------
    >>> model = PysMod()  # Assuming a properly initialized PysMod object
    >>> initial_conditions = [1.0, 0.0, 0.0]  # Example initial concentrations
    >>> sim_res, rates, success = model.CVODE(initial_conditions)
    >>> if success:
    ...     print("Integration completed successfully.")
    ...     print("Simulation results:", sim_res)
    ...     print("Reaction rates:", rates)
    ... else:
    ...     print("Integration failed.")

    """
        assert _HAVE_ASSIMULO, '\nAssimulo is not installed or did not import correctly\n{}'.format(
            _ASSIMULO_LOAD_ERROR)
        Vtemp = numpy.zeros(self.__Nshape__[1], 'd')

        def findi(t, s):
            self._TIME_ = t
            return self._EvalODE(s, Vtemp)

        def ffull(t, s):
            self._TIME_ = t
            return self._EvalODE_CVODE(s, Vtemp)
        if self.mode_integrate_all_odes:
            rhs = ffull
        else:
            rhs = findi
        if self.__CVODE_initialise__:
            if self.__HAS_RATE_RULES__:
                initial, rrules = numpy.split(initial, [self.Nrmatrix.shape[0]]
                    )
            if (self.__HAS_MOIETY_CONSERVATION__ and self.
                mode_integrate_all_odes):
                initial = self.Fix_S_indinput(initial, amounts=True)
            if self.__HAS_RATE_RULES__:
                initial = numpy.concatenate([initial, rrules])
        if not self.__settings__['cvode_track_assignment_rules']:
            self._CVODE_extra_output = []
        self._CVODE_XOUT = False
        if len(self._CVODE_extra_output) > 0:
            out = []
            for d in self._CVODE_extra_output:
                if (hasattr(self, d) and d not in self.__species__ + self.
                    __reactions__ + self.__rate_rules__):
                    out.append(d)
                else:
                    print(
                        """
Warning: CVODE is ignoring extra data ({}), it either doesn't exist or it's a species or rate.
"""
                        .format(d))
            if len(out) > 0:
                self._CVODE_extra_output = out
                self._CVODE_XOUT = True
            del out
        problem = EventsProblem(self, rhs=rhs, y0=initial)
        self._problem = problem
        sim = CVode(problem)
        if self.__settings__['cvode_access_solver']:
            self._solver = sim
        if self.__settings__['cvode_stats']:
            sim.verbosity = 10
        else:
            sim.verbosity = 40
        sim.atol = self.__settings__['cvode_abstol']
        sim.rtol = self.__settings__['cvode_reltol']
        sim.maxsteps = self.__settings__['cvode_mxstep']
        sim.inith = self.__settings__['cvode_h0']
        sim.maxh = self.__settings__['cvode_hmax']
        sim.minh = self.__settings__['cvode_hmin']
        sim.maxord = self.__settings__['cvode_mxord']
        t, sim_res = sim.simulate(self.sim_end, ncp=0, ncp_list=self.sim_time)
        t = numpy.array(t)
        idx = [0] + [numpy.max(numpy.where(t == i)) for i in problem.
            event_times] + [len(t)]
        rates = numpy.zeros((sim_res.shape[0], len(self.__reactions__)))
        if (self.__HAS_MOIETY_CONSERVATION__ and not self.
            mode_integrate_all_odes):
            for i in range(len(idx) - 1):
                for j in range(len(self.parameters)):
                    setattr(self, self.parameters[j], problem.parvals[i][j])
                for r in range(idx[i], idx[i + 1]):
                    self._EvalODE(sim_res[r].copy(), self._CVODE_Vtemp)
                    rates[r] = self.__vvec__
            if self.__HAS_RATE_RULES__:
                sim_res, rrules = numpy.split(sim_res, [self.Nrmatrix.shape
                    [0]], axis=1)
            res = numpy.zeros((sim_res.shape[0], len(self.__species__)))
            for x in range(sim_res.shape[0]):
                res[x] = self.Fix_S_indinput(sim_res[x], amounts=True)
            sim_res = res
            del res
            if self.__HAS_COMPARTMENTS__ and self.__KeyWords__['Output_In_Conc'
                ]:
                for x in range(0, sim_res.shape[0]):
                    sim_res[x] = self._SpeciesAmountToConc(sim_res[x])
            if self.__HAS_RATE_RULES__:
                sim_res = numpy.concatenate([sim_res, rrules], axis=1)
        else:
            for i in range(len(idx) - 1):
                for j in range(len(self.parameters)):
                    setattr(self, self.parameters[j], problem.parvals[i][j])
                for r in range(idx[i], idx[i + 1]):
                    self._EvalODE_CVODE(sim_res[r].copy(), self._CVODE_Vtemp)
                    rates[r] = self.__vvec__
            if self.__HAS_RATE_RULES__:
                sim_res, rrules = numpy.split(sim_res, [self.Nmatrix.shape[
                    0]], axis=1)
            if self.__HAS_COMPARTMENTS__ and self.__KeyWords__['Output_In_Conc'
                ]:
                for x in range(sim_res.shape[0]):
                    sim_res[x] = self._SpeciesAmountToConc(sim_res[x])
            if self.__HAS_RATE_RULES__:
                sim_res = numpy.concatenate([sim_res, rrules], axis=1)
        if self.__settings__['cvode_return_event_timepoints']:
            self.sim_time = t
            self._ev_idx = idx
        else:
            tidx = [numpy.where(t == i)[0][0] for i in self.sim_time]
            sim_res = sim_res[tidx]
            rates = rates[tidx]
            self._ev_idx = [0] + [numpy.min(numpy.where(self.sim_time >= i)
                ) for i in problem.event_times] + [len(t)]
        return sim_res, rates, True

    def CVODE_VPYTHON(self, s):
        """
    Future VPython hook for CVODE.

    This method is intended as a placeholder for future implementation of VPython
    integration with the CVODE solver. Currently, it does not carry out any operations 
    and serves as a conceptual hook for potential future functionalities.

    Parameters
    ----------
    s : any
        Currently, this parameter serves no purpose but is reserved for future use. 
        It is intended that in a future implementation, `s` might be used to pass 
        relevant simulation state or parameters for visualization in VPython.

    Notes
    -----
    This method is a stub and has no functionality at the moment. Users should not 
    expect any output or behavior from calling this method until it is fully implemented.

    Examples
    --------
    Currently, this method does not do anything and is meant for future expansion. 
    When implemented, its usage might look something like this:

    >>> model = PysMod()
    >>> state = ...  # some state or data to use with VPython
    >>> model.CVODE_VPYTHON(state)
    """
        pass

    def LSODA(self, initial):
        """
    PySCeS interface to the LSODA integration algorithm.

    This method interfaces with the LSODA integration algorithm to simulate dynamic behavior
    of a biochemical system starting with a given set of initial conditions. It returns an
    array of species concentrations and a status flag indicating the success of the simulation.
    Integration parameters and LSODA controls are accessible as attributes of the PysMod instance
    with the prefix `lsoda_`.

    Parameters
    ----------
    initial : array_like
        A vector containing the initial concentrations of the species.

    Returns
    -------
    sim_res : ndarray
        The simulation results; an array of species concentrations for each time point in the
        simulation.
    rates : ndarray
        The reaction rates corresponding to each time point in the simulation.
    status : bool
        A boolean flag indicating whether the integration was successful (`True`) or not (`False`).

    Warnings
    --------
    Integration of all ODEs is not supported with LSODA if `mode_integrate_all_odes` is set.
    PySCeS will then switch to integrating a reduced set of ODEs, and compute the dependent
    conserved species from the L-matrix. Consider using CVODE if this functionality is required.

    Notes
    -----
    - If integration errors occur, the method may automatically adjust the `lsoda_mxstep` setting
      to ensure convergence. If persistent errors remain after maximum allowed retries, consider
      switching to CVODE by setting `mod.mode_integrator = 'CVODE'`.

    Example
    -------
    >>> mod = PysMod()
    >>> initial_conditions = [1.0, 0.5, 2.0]  # Example initial concentrations for the species
    >>> sim_res, rates, status = mod.LSODA(initial_conditions)
    >>> if status:
    ...     print("Integration was successful.")
    ... else:
    ...     print("Integration failed. Consider trying CVODE.")

    """
        if self.mode_integrate_all_odes:
            print(
                """
NOTE: Integration of all ODEs is not supported with LSODA. PySCeS will integrate
a reduced set of ODEs and the dependent conserved species will be calculated from the L-matrix.
The `mod.mode_integrate_all_odes` flag is ignored. To explicitly integrate all ODEs, switch
to CVODE (mod.mode_integrator='CVODE')"""
                )
        Vtemp = numpy.zeros(self.__Nshape__[1], 'd')

        def function_sim(s, t):
            self._TIME_ = t
            return self._EvalODE(s, Vtemp)
        iter = 0
        go = True
        status = 0
        lsoda_mxstep_bak = self.__settings__['lsoda_mxstep']
        while go:
            sim_res, infodict = scipy.integrate.odeint(function_sim,
                initial, self.sim_time, atol=self.__settings__['lsoda_atol'
                ], rtol=self.__settings__['lsoda_rtol'], mxstep=self.
                __settings__['lsoda_mxstep'], h0=self.__settings__[
                'lsoda_h0'], hmax=self.__settings__['lsoda_hmax'], hmin=
                self.__settings__['lsoda_hmin'], mxordn=self.__settings__[
                'lsoda_mxordn'], mxords=self.__settings__['lsoda_mxords'],
                printmessg=self.__settings__['lsoda_mesg'], full_output=1)
            if infodict['message'] == 'Integration successful.':
                status = 0
            else:
                status = 1
            if status > 0 and iter < self.__settings__['mode_sim_max_iter']:
                if self.__settings__['lsoda_mxstep'] == 0:
                    print(
                        """
Integration error

Setting self.__settings__["lsoda_mxstep"] = 1000 and reSimulating ..."""
                        )
                    self.__settings__['lsoda_mxstep'] = 1000
                else:
                    print(
                        'Integration error\n\nSetting self.__settings__["lsoda_mxstep"] = '
                         + repr(self.__settings__['lsoda_mxstep'] * 3) +
                        ' and reSimulating ...')
                    self.__settings__['lsoda_mxstep'] = self.__settings__[
                        'lsoda_mxstep'] * 3
                iter += 1
            elif status > 0 and iter == self.__settings__['mode_sim_max_iter']:
                print(
                    """
This simulation is going nowhere fast
Consider trying CVODE (mod.mode_integrator = 'CVODE')
"""
                    )
                print('self.__settings__["lsoda_mxstep"] = ' + repr(self.
                    __settings__['lsoda_mxstep']))
                print("__settings__['mode_sim_max_iter'] = " + repr(iter))
                go = False
            else:
                go = False
        self.__settings__['lsoda_mxstep'] = lsoda_mxstep_bak
        rates = numpy.zeros((sim_res.shape[0], len(self.__reactions__)))
        if status == 0:
            tmp = None
            for r in range(sim_res.shape[0]):
                tmp = self._EvalODE(sim_res[r].copy(), self._CVODE_Vtemp)
                rates[r] = self.__vvec__
            del tmp
            if self.__HAS_RATE_RULES__:
                sim_res, rrules = numpy.split(sim_res, [self.Nrmatrix.shape
                    [0]], axis=1)
            if self.__HAS_MOIETY_CONSERVATION__ == True:
                res = numpy.zeros((sim_res.shape[0], len(self.__species__)))
                for x in range(0, sim_res.shape[0]):
                    res[x, :] = self.Fix_S_indinput(sim_res[x, :], amounts=True
                        )
                sim_res = res
                del res
            if self.__HAS_COMPARTMENTS__ and self.__KeyWords__['Output_In_Conc'
                ]:
                for x in range(0, sim_res.shape[0]):
                    sim_res[x] = self._SpeciesAmountToConc(sim_res[x])
            if self.__HAS_RATE_RULES__:
                sim_res = numpy.concatenate([sim_res, rrules], axis=1)
            return sim_res, rates, True
        else:
            if (self.__HAS_MOIETY_CONSERVATION__ == True and self.
                __HAS_RATE_RULES__):
                sim_res = numpy.zeros((sim_res.shape[0], len(self.
                    __species__) + len(self.__rrule__)), 'd')
            elif self.__HAS_MOIETY_CONSERVATION__ == True:
                sim_res = numpy.zeros((sim_res.shape[0], len(self.
                    __species__)), 'd')
            sim_res[:] = numpy.nan
            return sim_res, rates, False

    def HYBRD(self, initial):
        """
    HYBRD(initial)

    PySCeS interface to the HYBRD solver. Returns a steady-state solution and
    error flag. This solver is a good general-purpose solver for finding the
    roots of a system of nonlinear equations.

    Algorithm controls are available as mod.hybrd_<control>.

    Parameters
    ----------
    initial : array_like
        A vector of initial species concentrations used as a starting point
        for the solver.

    Returns
    -------
    steady_state : ndarray
        An array representing the steady-state solution of the species concentrations.
    success : bool
        A boolean flag indicating whether the solver successfully found a steady state.
        Returns `True` if successful, otherwise `False`.

    Examples
    --------
    >>> import numpy as np
    >>> from pysmod import PysMod
    >>> model = PysMod()
    >>> initial_concentrations = np.array([0.5, 0.3, 0.2])
    >>> steady_state, success = model.HYBRD(initial_concentrations)
    >>> if success:
    ...     print("Steady state found:", steady_state)
    ... else:
    ...     print("Failed to find a steady state.")

    Notes
    -----
    - The available HYBRD controls such as tolerance and maximum function evaluations
      can be accessed and possibly adjusted through the `mod.hybrd_<control>` attributes.
    - Ensure that `initial` is of compatible dimensions and units as expected by the model.
    """
        Vtemp = numpy.zeros(self.__Nshape__[1], 'd')

        def function_state(s):
            return self._EvalODE(s, Vtemp)
        state_out = scipy.optimize.fsolve(function_state, initial, args=(),
            xtol=self.__settings__['hybrd_xtol'], maxfev=self.__settings__[
            'hybrd_maxfev'], epsfcn=self.__settings__['hybrd_epsfcn'],
            factor=self.__settings__['hybrd_factor'], col_deriv=0,
            full_output=1)
        if state_out[2] == 1:
            if self.__settings__['hybrd_mesg'] == 1:
                print('(hybrd)', state_out[3])
            return state_out[0], True
        else:
            if self.__settings__['hybrd_mesg']:
                print('INFO: (hybrd) Invalid steady state:')
            if self.__settings__['hybrd_mesg'] == 1:
                print('(hybrd)', state_out[3])
            return state_out[0], False

    def FINTSLV(self, initial):
        """
    FINTSLV(initial)

    Forward integration steady-state solver. This method finds a steady state when the
    maximum change in species concentration falls within a specified tolerance.
    It returns the steady-state solution and an error flag indicating success or failure.
    Algorithm controls are available as `mod.fintslv_<control>`.

    Parameters
    ----------
    initial : array_like
        Vector of initial species concentrations.

    Returns
    -------
    steady_state : array_like
        The species concentrations at steady state.
    success : bool
        True if the steady state was successfully found within the specified tolerance,
        False otherwise.

    Example
    -------
    >>> mod = PysMod()
    >>> initial_conditions = [1.0, 0.5, 0.2]
    >>> steady_state, success = mod.FINTSLV(initial_conditions)
    >>> if success:
    ...     print("Steady state found:", steady_state)
    ... else:
    ...     print("Failed to find steady state.")

    Notes
    -----
    The method integrates the system forward in time using the LSODA algorithm
    within `scipy.integrate.odeint`. It checks for convergence to a steady state
    when the maximum change in species concentration over successive steps is
    below a specified tolerance.
    Algorithm controls like `fintslv_tol` (tolerance) and `fintslv_step` (number of successive small changes)
    are configurable via the model settings.
    """
        sim_time = self.__fintslv_range__ * self.__settings__['fintslv_rmult']
        Vtemp = numpy.zeros(self.__Nshape__[1], 'd')

        def function_sim(s, t):
            return self._EvalODE(s, Vtemp)
        res, infodict = scipy.integrate.odeint(function_sim, initial.copy(),
            sim_time, atol=self.__settings__['lsoda_atol'], rtol=self.
            __settings__['lsoda_rtol'], mxstep=10000, h0=self.__settings__[
            'lsoda_h0'], hmax=self.__settings__['lsoda_hmax'], hmin=self.
            __settings__['lsoda_hmin'], mxordn=self.__settings__[
            'lsoda_mxordn'], mxords=self.__settings__['lsoda_mxords'],
            printmessg=self.__settings__['lsoda_mesg'], full_output=1)
        if infodict['message'] == 'Integration successful.':
            status = True
        else:
            status = False
        OK = 0
        if status:
            for x in range(len(res)):
                if OK >= self.__settings__['fintslv_step']:
                    break
                if x > 0 and OK < self.__settings__['fintslv_step']:
                    dif = abs(res[x] - res[x - 1])
                    if max(dif) < self.__settings__['fintslv_tol']:
                        OK += 1
                    else:
                        pass
        if OK >= 5:
            return res[-1], True
        else:
            return res[-1], False

    def NLEQ2(self, initial):
        """
    Solves for a steady state using the NLEQ2 algorithm via the PySCeS interface.

    The `NLEQ2` method is a powerful steady-state solver that is often effective 
    when the `HYBRD` method fails. It offers various algorithm controls, accessible 
    through settings with the prefix `mod.nleq2_<control>`. The function returns a 
    steady-state solution along with an error flag indicating the success of the 
    solution process.

    Parameters
    ----------
    initial : array_like
        A vector representing the initial species concentrations.

    Returns
    -------
    array_like
        The computed steady-state solution vector.
    bool
        Error flag where `True` indicates successful convergence of the solver, 
        and `False` indicates failure.

    Notes
    -----
    The algorithm controls for the NLEQ2 solver can be customized through a variety 
    of options, such as setting advanced modes or iteration limits, to increase the 
    robustness and efficiency of the solution process.

    Example
    -------
    >>> import PySCeS
    >>> model = PySCeS.PysMod()
    >>> initial_concentrations = [1.0, 0.5, 0.1]
    >>> solution, success = model.NLEQ2(initial_concentrations)
    >>> if success:
    ...     print("Steady-state solution found:", solution)
    ... else:
    ...     print("Solution did not converge.")
    """
        Vtemp = numpy.zeros(self.__Nshape__[1], 'd')
        initial0 = initial.copy()
        s_scale = initial.copy()
        N = len(initial)
        iwk = numpy.zeros(N + 52, 'i')
        rwk = numpy.zeros((N + max(N, 10) + 15) * N + 61, 'd')
        iopt = numpy.zeros(50, 'i')
        if self.__settings__['nleq2_jacgen'] == 1:
            print(
                "(nleq2)User supplied Jacobian not supported yet ... setting __settings__['nleq2_jacgen'] = 0"
                )
            self.__settings__['nleq2_jacgen'] = 0
        rtol = mach_eps * 10.0 * N
        iopt[2] = self.__settings__['nleq2_jacgen']
        iopt[8] = self.__settings__['nleq2_iscaln']
        iopt[10] = 0
        iopt[11] = 6
        iopt[12] = 0
        iopt[13] = 6
        iopt[14] = 0
        iopt[15] = 6
        iopt[30] = self.__settings__['nleq2_nonlin']
        iopt[31] = self.__settings__['nleq2_qrank1']
        iopt[34] = self.__settings__['nleq2_qnscal']
        iopt[37] = self.__settings__['nleq2_ibdamp']
        iopt[38] = self.__settings__['nleq2_iormon']
        iwk[30] = self.nleq2_nitmax

        def func(s, ifail):
            s = self._EvalODE(s, Vtemp)
            if numpy.isnan(s).any():
                ifail = -1
            elif (numpy.abs(s) > 1e+150).any():
                ifail = -1
            return s, ifail

        def jacfunc(s):
            return s
        ierr = 0
        BRETT_DEBUG_MODE = False
        ADVANCED_MODE = self.__settings__['nleq2_advanced_mode']
        if ADVANCED_MODE:
            max_iter = self.__settings__['nleq2_iter']
            max_iter_ceiling = self.__settings__['nleq2_iter_max']
        else:
            max_iter_ceiling = self.__settings__['nleq2_iter']
            max_iter = self.__settings__['nleq2_iter']
        iter = 1
        GO = True
        while GO:
            if BRETT_DEBUG_MODE:
                print('s_scale', s_scale)
                print('ierr({}) = {}'.format(iter, ierr))
                print('rtol({}) = {}'.format(iter, rtol))
                print('nitmax({}) = {}'.format(iter, iwk[30]))
                print('s_scale({}) = {}'.format(iter, s_scale))
                print('max_iter({}) = {}'.format(iter, max_iter))
            res, s_scale, rtol, iopt, ierr = nleq2.nleq2(func, jacfunc,
                initial, s_scale, rtol, iopt, iwk, rwk)
            if ierr == 0:
                GO = False
            elif ierr == 21:
                GO = False
            elif ierr == 2:
                iwk[30] = iwk[30] * self.__settings__['nleq2_growth_factor']
            if iter >= max_iter and iter >= max_iter_ceiling:
                GO = False
            iter += 1
        if BRETT_DEBUG_MODE and ierr > 0:
            print('ierr = {}'.format(ierr))
            print(res)
            time.sleep(5)
        if ierr > 0:
            if self.__settings__['nleq2_mesg']:
                print('(nleq2) exits with ierr = {}'.format(ierr))
        elif self.__settings__['nleq2_mesg']:
            print('(nleq2) The solution converged.')
        if ierr > 0:
            return res, False
        else:
            return res, True

    def PITCON(self, scanpar, scanpar3d=None):
        """
    PITCON(scanpar, scanpar3d=None)

    PySCeS interface to the PITCON continuation algorithm. This method performs
    a single parameter continuation, implemented as a "scan", starting in
    `mod.pitcon_par_space`. The second argument does not affect the continuation
    but can be used to include an additional parameter along the third axis in
    the results. It returns an array of results corresponding to the parameter
    continuation. Various algorithm controls are accessible via `mod.pitcon_<control>` settings.

    Parameters
    ----------
    scanpar : str
        The model parameter to perform the continuation scan on.
        Should be a string representing a parameter in the model.
        
    scanpar3d : float or None, optional
        An additional parameter that can be included for 3D plot outputs. 
        This does not affect the continuation process. Default is None.

    Returns
    -------
    numpy.ndarray
        An array containing the results of the parameter continuation scan.

    Raises
    ------
    AssertionError
        If PITCON is not installed or if `scanpar` is not a string.
        
    NameError
        If `scanpar` is not a parameter in the model.
        
    NotImplementedError
        If the model contains RateRules, as bifurcation analysis is not 
        currently supported for such models.

    Notes
    -----
    The PITCON algorithm allows for advanced continuation study of model parameters,
    and it provides an interface for external configurations through several 
    algorithm control settings.

    Examples
    --------
    To perform a parameter continuation scan on a model parameter 'x5' with 
    an additional plotting parameter for 3D visualization:

    >>> model = PysMod()
    >>> scan_results = model.PITCON('x5', scanpar3d=1.0)
    >>> print(scan_results)
    array([...])  # An array containing scan results.
    """
        assert pitcon_switch, 'PITCON is not installed! Continuation analysis not available.'
        if self.__HAS_RATE_RULES__:
            raise NotImplementedError(
                """
Bifurcation analysis not currently available for models containing RateRules"""
                )
        assert type(scanpar
            ) == str, '\nscanpar must be a <string> representing a model parameter'
        modpar = list(self.__parameters__)
        try:
            a = modpar.index(scanpar)
        except:
            raise NameError(repr(scanpar) + ' is not a parameter of this model'
                )
        if scanpar3d != None:
            if type(scanpar3d) == str:
                scanpar3d = float(scanpar3d)
                assert type(scanpar3d) == float, 'scanpar3d must be a <float>'
        par_hold = getattr(self, scanpar)
        if self.__settings__['pitcon_jac_opt'] < 1:
            self.__settings__['pitcon_jac_opt'] = 1
            print(
                """
INFO: .__settings__["pitcon_jac_opt"] set to 1 - user defined jacobian function not yet supported"""
                )

        def fx(s):
            setattr(self, scanpar, s[-1])
            try:
                sdot[:-1] = self._EvalODE(s[:-1], Vtemp)
            except Exception as ex:
                print('PITCON EXCEPTION 1', ex)
                sdot[:-1] = 0.0
            sdot[-1] = s[-1]
            return sdot
        parm = len(self.__SI__) + 1
        Vtemp = numpy.zeros(self.__Nshape__[1], 'd')
        xr2 = numpy.zeros(parm, 'd')
        sdot = numpy.zeros(parm, 'd')
        iwork = numpy.zeros(30 + parm, 'i')
        rwork = numpy.zeros(30 + 6 * parm * parm, 'd')
        ipar = numpy.zeros(parm, 'i')
        fpar = numpy.zeros(parm, 'd')
        res = []
        for xscan in self.pitcon_par_space:
            Vtemp[:] = 0.0
            xr2[:] = 0.0
            sdot[:] = 0.0
            ipar[:] = 0
            fpar[:] = 0.0
            iwork[0] = 0
            iwork[1] = self.__settings__['pitcon_init_par']
            iwork[2] = self.__settings__['pitcon_par_opt']
            iwork[3] = self.__settings__['pitcon_jac_upd']
            iwork[4] = self.__settings__['pitcon_targ_val_idx']
            iwork[5] = self.__settings__['pitcon_limit_point_idx']
            iwork[6] = self.__settings__['pitcon_output_lvl']
            iwork[7] = 6
            iwork[8] = self.__settings__['pitcon_jac_opt']
            iwork[9] = 0
            iwork[10] = 0
            iwork[11] = 0
            iwork[12] = 30
            iwork[13] = len(iwork)
            iwork[14] = 30 + 4 * parm
            iwork[15] = len(rwork)
            iwork[16] = self.__settings__['pitcon_max_step']
            iwork[17] = 0
            iwork[18] = 0
            iwork[19] = 0
            iwork[20] = 0
            iwork[21] = 0
            iwork[22] = 0
            iwork[23] = 0
            iwork[24] = 0
            iwork[25] = 0
            iwork[26] = 0
            iwork[27] = 0
            rwork[0] = self.__settings__['pitcon_abs_tol']
            rwork[1] = self.__settings__['pitcon_rel_tol']
            rwork[2] = self.__settings__['pitcon_min_step']
            rwork[3] = self.__settings__['pitcon_max_step']
            rwork[4] = self.__settings__['pitcon_start_step']
            rwork[5] = self.__settings__['pitcon_start_dir']
            rwork[6] = self.pitcon_targ_val
            rwork[7] = 0.0
            rwork[8] = 0.0
            rwork[9] = 0.0
            rwork[10] = 0.0
            rwork[11] = 0.0
            rwork[12] = 0.0
            rwork[13] = 0.0
            rwork[14] = 0.0
            rwork[15] = 0.0
            rwork[16] = 0.0
            rwork[17] = 0.0
            rwork[18] = 0.0
            rwork[19] = self.__settings__['pitcon_max_grow']
            rwork[20] = 0.0
            rwork[21] = 0.0
            rwork[22] = 0.0
            rwork[23] = 0.0
            rwork[24] = 0.0
            rwork[25] = 0.0
            rwork[26] = 0.0
            rwork[27] = 0.0
            rwork[28] = 0.0
            setattr(self, scanpar, xscan)
            self.State()
            if self.__HAS_COMPARTMENTS__ and self.__KeyWords__['Output_In_Conc'
                ]:
                self.__inspec__ = copy.copy(self._SpeciesConcToAmount(self.
                    state_species))
            else:
                self.__inspec__ = copy.copy(self.state_species)
            go = False
            if self.__StateOK__:
                go = True
            elif not self.__StateOK__ and self.__settings__[
                'pitcon_allow_badstate']:
                go = True
            else:
                go = False
            if go:
                if self.__HAS_MOIETY_CONSERVATION__ == True:
                    temp = numpy.zeros(len(self.__SI__), 'd')
                    for x in range(0, len(self.lzeromatrix_col)):
                        xr2[x] = self.__inspec__[self.nrmatrix_row[x]]
                else:
                    xr2[:-1] = copy.copy(self.state_species[:])
                xr2[-1] = copy.copy(xscan)
                for x in range(int(self.pitcon_iter)):
                    ierror, iwork, rwork, xr2 = pitcon.pitcon1(fx, fpar, fx,
                        ipar, iwork, rwork, xr2)
                    if iwork[0] == 2:
                        if self.__settings__['pitcon_flux_gen']:
                            if min(xr2[:-1]) < 0.0 and self.__settings__[
                                'pitcon_filter_neg']:
                                pass
                            else:
                                xout = xr2.tolist()
                                a = xout.pop(-1)
                                if self.__HAS_MOIETY_CONSERVATION__ == True:
                                    xout = self.Fix_S_indinput(xout,
                                        amounts=True)
                                else:
                                    xout = numpy.array(xout)
                                if (self.__HAS_COMPARTMENTS__ and self.
                                    __KeyWords__['Output_In_Conc']):
                                    xout = self._SpeciesAmountToConc(xout)
                                xout2 = self.__FluxGen__(xout).tolist()
                                xout = xout.tolist()
                                xout.insert(0, a)
                                if scanpar3d != None:
                                    xout.insert(0, scanpar3d)
                                xout2 = xout + xout2
                                res.append(xout2)
                        elif min(xr2[:-1]) < 0.0 and self.__settings__[
                            'pitcon_filter_neg']:
                            pass
                        else:
                            xout = xr2.tolist()
                            a = xout.pop(-1)
                            if self.__HAS_MOIETY_CONSERVATION__ == True:
                                xout = self.Fix_S_indinput(xout, amounts=True)
                            else:
                                xout = numpy.array(xout)
                            if self.__HAS_COMPARTMENTS__ and self.__KeyWords__[
                                'Output_In_Conc']:
                                xout = self._SpeciesAmountToConc(xout)
                            xout = xout.tolist()
                            xout.insert(0, a)
                            if scanpar3d != None:
                                xout.insert(0, scanpar3d)
                            res.append(xout)
                    elif iwork[0] == 3:
                        print('\nTarget point:')
                        xout = xr2.tolist()
                        a = xout.pop(-1)
                        if self.__HAS_MOIETY_CONSERVATION__ == True:
                            xout = self.Fix_S_indinput(xout, amounts=True)
                        else:
                            xout = numpy.array(xout)
                        if self.__HAS_COMPARTMENTS__ and self.__KeyWords__[
                            'Output_In_Conc']:
                            xout = self._SpeciesAmountToConc(xout)
                        xout = xout.tolist()
                        xout.insert(0, a)
                        print(xout)
                        self.pitcon_target_points.append(xout)
                    elif iwork[0] == 4:
                        print('\nLimit point')
                        xout = xr2.tolist()
                        a = xout.pop(-1)
                        if self.__HAS_MOIETY_CONSERVATION__ == True:
                            xout = self.Fix_S_indinput(xout, amounts=True)
                        else:
                            xout = numpy.array(xout)
                        if self.__HAS_COMPARTMENTS__ and self.__KeyWords__[
                            'Output_In_Conc']:
                            xout = self._SpeciesAmountToConc(xout)
                        xout = xout.tolist()
                        xout.insert(0, a)
                        print(xout)
                        self.pitcon_limit_points.append(xout)
                    elif iwork[0] == 1:
                        pass
                    else:
                        print(iwork[0])
            else:
                print('\nInvalid steady state, skipping ...')
        if self.__settings__['pitcon_filter_neg_res']:
            for result in range(len(res) - 1, -1, -1):
                if min(res[result][:len(self.species) + 1]) < 0.0:
                    res.pop(result)
        setattr(self, scanpar, par_hold)
        return numpy.array(res)

    def Simulate(self, userinit=0):
        """
    Simulate the system dynamics over time.

    This method is the PySCeS integration driver routine that evolves the system
    over the specified time interval. The resulting array of species concentrations
    is stored in the `mod.data_sim` object. Initial species concentrations can be
    determined based on the `mod.__settings__['mode_sim_init']` parameter.

    Parameters
    ----------
    userinit : int, optional
        Specifies how the initial species concentrations and time array should be set
        (default is 0):
        
        - 0: Initial concentrations are set from `mod.S_init`, and the time array 
          is calculated from `sim_start`, `sim_end`, and `sim_points`.
        - 1: Initial concentrations are set from `mod.S_init`, and the existing
          `mod.sim_time` is used.
        - 2: Initial concentrations are read from `mod.__inspec__`, and the existing
          `mod.sim_time` is used.

    `mod.__settings__['mode_sim_init']` allows different initialization strategies:

    - 0: Initializes with initial concentrations.
    - 1: Initializes with very small (close to zero) values.
    - 2: Initializes with the steady state from a previous calculation.
    - 3: Initializes with the final point of a previous simulation.

    Notes
    -----
    The method checks for the validity of species conservation laws and incorporates
    the effects of moiety conservation if present. It supports integrating the model
    using either LSODA or CVODE solvers, based on the `mode_integrator` setting.

    Raises
    ------
    NotImplementedError
        If `__StoichOK` is not set, the stoichiometry is reanalyzed.

    Examples
    --------
    Simulate the system using default initial conditions and time array:

    >>> model = PysMod()
    >>> model.Simulate()

    Use the results of a previously calculated steady state for initial conditions:

    >>> model.__settings__['mode_sim_init'] = 2
    >>> model.Simulate()

    Specify initial species concentrations from `mod.S_init` and use custom time array:

    >>> model.Simulate(userinit=1)

    After simulation, access the results from the `mod.data_sim` object:

    >>> species_data = model.data_sim.species
    >>> time_data = model.data_sim.time
    """
        if not self.__StoichOK:
            self.Stoichiometry_ReAnalyse()
        self._sim_time_bak = None
        if userinit != 0:
            if self.sim_time[0] != 0:
                self._sim_time_bak = copy.copy(self.sim_time)
                self.sim_time = [0.0] + list(self._sim_time_bak)
            self.sim_end = self.sim_time[-1]
        if userinit == 1:
            eval(self.__mapFunc_R__)
        elif userinit == 2:
            try:
                assert len(self.__inspec__) == len(self.__species__)
            except:
                print(
                    '\nINFO: mod.__inspec__ is the incorrect length, initialising with .sX_init'
                    )
                self.__inspec__ = numpy.zeros(len(self.__species__))
                eval(self.__mapFunc_R__)
        else:
            self.sim_start = float(self.sim_start)
            self.sim_end = float(self.sim_end)
            self.sim_points = int(self.sim_points)
            if self.sim_points == 1.0:
                print(
                    """*****
WARNING: simulations require a minimum of 2 points,setting sim_points = 2.0
*****"""
                    )
                self.sim_points = 2.0
            self.sim_time = numpy.linspace(self.sim_start, self.sim_end,
                self.sim_points, endpoint=True, retstep=False)
            eval(self.__mapFunc_R__)
        if self.__HAS_RATE_RULES__:
            for r in range(len(self.__rate_rules__)):
                self.__rrule__[r] = getattr(self, '{}_init'.format(self.
                    __rate_rules__[r]))
        for c in range(len(self.__CsizeAllIdx__)):
            cval = getattr(self, '{}_init'.format(self.__CsizeAllIdx__[c]))
            setattr(self, self.__CsizeAllIdx__[c], cval)
            self.__CsizeAll__[c] = cval
        if self.__HAS_COMPARTMENTS__ and self.__KeyWords__['Species_In_Conc']:
            self.__inspec__ = self._SpeciesConcToAmount(self.__inspec__)
        for s in range(len(self.__species__)):
            setattr(self, self.__species__[s], self.__inspec__[s])
        if self.__HAS_MOIETY_CONSERVATION__ == True:
            self.__Build_Tvec__(amounts=True)
            self.showConserved(screenwrite=0)
            if self.__settings__['display_debug'] == 1:
                print(self.conserved_sums)
        if self.__settings__['mode_sim_init'] == 0:
            s0_sim_init = copy.copy(self.__inspec__)
        elif self.__settings__['mode_sim_init'] == 1:
            s0_sim_init = copy.copy(self.__inspec__)
            s0_sim_init[:] = self.__settings__['small_concentration']
        elif self.__settings__['mode_sim_init'] == 2:
            if self.state_species != None:
                s0_sim_init = copy.copy(self.state_species) * 1.0
                for x in range(len(s0_sim_init)):
                    if s0_sim_init[x] < 0.0:
                        s0_sim_init[x] = self.__settings__[
                            'small_concentration']
                        print(
                            'Negative concentration detected in SimInit: s[' +
                            repr(x) + '] set to ' + repr(self.__settings__[
                            'small_concentration']))
            else:
                s0_sim_init = copy.copy(self.__inspec__)
        elif self.__settings__['mode_sim_init'] == 3:
            try:
                s0_sim_init = copy.copy(self.data_sim.species[-1])
            except:
                s0_sim_init = copy.copy(self.__inspec__)
        else:
            s0_sim_init = copy.copy(self.__inspec__)
        if self.__HAS_MOIETY_CONSERVATION__ == True:
            temp = numpy.zeros(len(self.__SI__), 'd')
            for x in range(0, len(self.lzeromatrix_col)):
                temp[x] = s0_sim_init[self.nrmatrix_row[x]]
            s0_sim_init = temp
            del temp
        if self.__HAS_RATE_RULES__:
            s0_sim_init = numpy.concatenate([s0_sim_init, self.__rrule__])
        self._sim = None
        Tsim0 = time.time()
        if self.mode_integrator == 'LSODA':
            sim_res, rates, simOK = self.LSODA(copy.copy(s0_sim_init))
        elif self.mode_integrator == 'CVODE':
            sim_res, rates, simOK = self.CVODE(copy.copy(s0_sim_init))
        if self._sim_time_bak is not None:
            self.sim_time = copy.copy(self._sim_time_bak)
            sim_res = sim_res[1:]
            rates = rates[1:]
        Tsim1 = time.time()
        if self.__settings__['lsoda_mesg']:
            print('{} time for {} points: {}'.format(self.mode_integrator,
                len(self.sim_time), Tsim1 - Tsim0))
        if self.__HAS_RATE_RULES__:
            sim_res, rrules = numpy.split(sim_res, [len(self.__species__)],
                axis=1)
            print('RateRules evaluated and added to mod.data_sim.')
        self.data_sim = IntegrationDataObj()
        self.IS_VALID = simOK
        self.data_sim.setTime(self.sim_time)
        self.data_sim.setSpecies(sim_res, self.__species__)
        self.data_sim.setRates(rates, self.__reactions__)
        if self.__HAS_RATE_RULES__:
            self.data_sim.setRules(rrules, self.__rate_rules__)
        if self._CVODE_XOUT:
            self._CVODE_xdata = numpy.zeros((len(self.sim_time), len(self.
                _CVODE_extra_output)))
            r = self.__rules__
            ars = {name: r[name] for name in r if r[name]['type'] ==
                'assignment'}
            for k in range(len(self._CVODE_extra_output)):
                name = self._CVODE_extra_output[k]
                self._update_assignment_rule_code(ars[name])
                for i in range(len(self._ev_idx) - 1):
                    for j in range(len(self.parameters)):
                        setattr(self, self.parameters[j], self._problem.
                            parvals[i][j])
                    p = self._ev_idx[i]
                    q = self._ev_idx[i + 1]
                    self._CVODE_xdata[p:q, k] = eval(ars[name][
                        'data_sim_string'])
                self.data_sim.setXData(self._CVODE_xdata, lbls=self.
                    _CVODE_extra_output)
            self._CVODE_xdata = None
        if not simOK:
            print('Simulation failure')
        del sim_res

    def _update_assignment_rule_code(self, rule):
        """
    Updates the data simulation string for a given assignment rule.

    This method transforms the `code_string` in the `rule` to a `data_sim_string`
    suitable for numerical simulations. It replaces symbols in the rule's code with
    references to simulation data stored in `self.data_sim`. This is necessary for
    post-simulation analysis where symbols can represent model species, reactions, or
    other rules.

    Parameters
    ----------
    rule : dict
        A dictionary representing an assignment rule. It must contain the keys:
        - `symbols`: A list of symbols present in `code_string`.
        - `code_string`: The original string containing the rule's computation logic.
        - `data_sim_string`: The output string where replacements have been made to 
          reference simulation data.

    Notes
    -----
    By replacing symbols with calls to `self.data_sim.getSimData("{symbol}")[p:q, 1]`,
    the method aligns rule evaluation with the integrated simulation timeline `p:q`.
    This is especially useful in cases involving symbolic overlaps or substrings between
    different rule components.

    Examples
    --------
    Suppose you have a rule defined as follows:

    >>> rule = {
    ...     'symbols': ['S1', 'J1'],
    ...     'code_string': 'S1 + J1'
    ... }
    >>> model._update_assignment_rule_code(rule)
    >>> print(rule['data_sim_string'])
    'self.data_sim.getSimData("S1")[p:q,1] + self.data_sim.getSimData("J1")[p:q,1]'

    In this example, the method adjusts the `code_string` to refer to simulation data
    for species `S1` and reaction `J1`, allowing for time-course evaluation post 
    simulation.
    """
        replacements = []
        rule['data_sim_string'] = rule['code_string']
        for s in sorted(rule['symbols'], key=len, reverse=True):
            if (s in self.__reactions__ or s in self.__rules__ or s in self
                .__species__):
                partialmatches = []
                for sym in rule['symbols']:
                    if not re.fullmatch(s, sym) and re.search(s, sym):
                        partialmatches.append(sym)
                for e in enumerate(partialmatches):
                    replacements.append((e[1], f'_zzzz_{e[0]}_z'))
                replacements.append((f'self.{s}',
                    f'self.data_sim.getSimData("{s}")[p:q,1]'))
                for e in enumerate(partialmatches):
                    replacements.append((f'_zzzz_{e[0]}_z', e[1]))
        for old, new in replacements:
            rule['data_sim_string'] = rule['data_sim_string'].replace(old, new)

    @property
    def sim(self):
        """
    Retrieve the simulation data in the specified custom data type format.

    This property fetches the simulation data from the `data_sim` object and
    formats it based on the user's settings. If the custom data type specified
    in the settings is `pandas`, the data is returned as a Pandas DataFrame;
    otherwise, it is returned as a NumPy structured array.

    Returns
    -------
    sim_data : pandas.DataFrame or numpy.recarray
        The simulation data formatted according to the custom data type
        specified in the settings.

    Examples
    --------
    Assuming `model` is an instance of `PysMod` and a simulation has been
    completed, access the simulation data as follows:

    >>> model.__settings__['custom_datatype'] = 'pandas'
    >>> sim_data = model.sim
    >>> type(sim_data)
    <class 'pandas.core.frame.DataFrame'>

    Or for a NumPy structured array:

    >>> model.__settings__['custom_datatype'] = 'numpy'
    >>> sim_data = model.sim
    >>> type(sim_data)
    <class 'numpy.recarray'>
    """
        if self._sim is None and self.data_sim is not None:
            data = self.data_sim.getAllSimData(lbls=True)
            if self.__settings__['custom_datatype'] == 'pandas':
                self._sim = self.__pandas.DataFrame(data[0], columns=data[1])
            else:
                self._sim = numpy.rec.fromrecords(data[0], names=data[1])
        return self._sim

    def State(self):
        """
    Solve for a steady state using non-linear solver algorithms.

    This method drives the PySCeS non-linear solver routine to find a steady 
    state solution for a model using one of the following algorithms: HYBRD, 
    NLEQ2, or FINTSLV. The results are stored in `mod.state_species` and 
    `mod.state_flux`. These results can be viewed using the `mod.showState()` 
    method.

    The solver can be initialized in three different ways using the 
    `mode_state_init` attribute:
    
    - `mod.mode_state_init = 0`: Initialize with species initial values.
    - `mod.mode_state_init = 1`: Initialize with small values.
    - `mod.mode_state_init = 2`: Initialize with the final value of a 
      10-logstep simulation (`numpy.logspace(0, 5, 18)`).

    Parameters
    ----------
    None

    Attributes
    ----------
    state_species : ndarray
        Steady-state species concentrations.
    state_flux : ndarray
        Steady-state flux values.
    data_sstate : StateDataObj
        Object containing the steady-state species and flux data.

    Warnings
    --------
    This method will print warnings if negative concentrations are detected or 
    if extremely small fluxes are present. It will also print an error message 
    if an invalid steady state solution is found.

    Example
    -------
    To compute a steady state solution for a model, ensure the model is 
    properly configured and call the State method:
    
    >>> model = PysMod()  # Assume PysMod is properly initialized.
    >>> model.mode_state_init = 0
    >>> model.State()
    >>> model.showState()

    This will calculate the steady state and display the results using 
    `showState`.
    """
        if not self.__StoichOK:
            self.Stoichiometry_ReAnalyse()
        self.__StateOK__ = True
        Vtemp = numpy.zeros(self.__Nshape__[1], 'd')

        def function_sim(s, t):
            return self._EvalODE(s, Vtemp)
        eval(self.__mapFunc_R__)
        if self.__HAS_RATE_RULES__:
            for r in range(len(self.__rate_rules__)):
                self.__rrule__[r] = getattr(self, '{}_init'.format(self.
                    __rate_rules__[r]))
        for c in range(len(self.__CsizeAllIdx__)):
            cval = getattr(self, '{}_init'.format(self.__CsizeAllIdx__[c]))
            setattr(self, self.__CsizeAllIdx__[c], cval)
            self.__CsizeAll__[c] = cval
        if self.__HAS_COMPARTMENTS__ and self.__KeyWords__['Species_In_Conc']:
            self.__inspec__ = self._SpeciesConcToAmount(self.__inspec__)
        for s in range(len(self.__species__)):
            setattr(self, self.__species__[s], self.__inspec__[s])
        if self.__HAS_MOIETY_CONSERVATION__ == True:
            self.__Build_Tvec__(amounts=True)
            self.showConserved(screenwrite=0)
            if self.__settings__['display_debug'] == 1:
                print(self.conserved_sums)
        s0_ss_init = None
        hybrd_factor_temp = copy.copy(self.__settings__['hybrd_factor'])
        if self.mode_state_init == 0:
            s0_ss_init = self.__inspec__.copy()
        elif self.mode_state_init == 1:
            s0_ss_init = self.__inspec__.copy()
            s0_ss_init[:] = self.__settings__['small_concentration']
        elif self.mode_state_init == 2:
            print(
                'This initialisation mode has been disabled, using initial values.'
                )
            s0_ss_init = self.__inspec__.copy()
        elif self.mode_state_init == 3:
            if (self.state_species != None and self.__HAS_COMPARTMENTS__ and
                self.__KeyWords__['Output_In_Conc']):
                s0_ss_init = self._SpeciesConcToAmount(abs(copy.copy(self.
                    state_species)) * float(self.__settings__[
                    'mode_state_init3_factor']))
            else:
                s0_ss_init = copy.copy(self.__inspec__)
        else:
            s0_ss_init = copy.copy(self.__inspec__)
        self.__settings__['hybrd_factor'] = hybrd_factor_temp
        if self.__HAS_MOIETY_CONSERVATION__ == True:
            temp = numpy.zeros(len(self.__SI__), 'd')
            for x in range(0, len(self.lzeromatrix_col)):
                temp[x] = s0_ss_init[self.nrmatrix_row[x]]
            if self.__HAS_RATE_RULES__:
                temp = numpy.concatenate([temp, self.__rrule__])
            s0_ss_init = temp
            del temp
        available_solvers = ['HYBRD']
        if self.mode_solver == 0:
            self.mode_solver = 'HYBRD'
        elif self.mode_solver == 1:
            self.mode_solver = 'FINTSLV'
        elif self.mode_solver == 2:
            self.mode_solver = 'NLEQ2'
        if nleq2_switch:
            available_solvers.append('NLEQ2')
        else:
            if self.mode_solver == 'NLEQ2':
                self.mode_solver = 'HYBRD'
            print(
                """INFO: switching to HYBRD.
Nleq2 solver not available see /nleq/readme.txt for details"""
                )
        if self.__settings__['mode_solver_fallback_integration']:
            available_solvers.append('FINTSLV')
        if not self.mode_solver_fallback == 1:
            assert self.mode_solver in available_solvers, '\nERROR: {} is not a valid ({}) solver!'.format(
                solver, str(available_solvers))
            available_solvers = [self.mode_solver]
        if self.mode_solver != 'HYBRD':
            available_solvers.insert(0, available_solvers.pop(
                available_solvers.index(self.mode_solver)))
        STATE_XOUT = False
        STATE_xdata = None
        if len(self.STATE_extra_output) > 0:
            out = []
            for d in self.STATE_extra_output:
                if (hasattr(self, d) and d not in self.__species__ + self.
                    __reactions__ + self.__rate_rules__):
                    out.append(d)
                else:
                    print(
                        """
WARNING: STATE is ignoring extra data ({}), it either doesn't exist or it's a species, rate or rule.
"""
                        .format(d))
            if len(out) > 0:
                self.STATE_extra_output = out
                STATE_XOUT = True
            del out
        self.__StateOK__ = True
        state_species = None
        rrules = None
        state_flux = None
        for solver in available_solvers:
            if solver == 'HYBRD':
                state_species, self.__StateOK__ = self.HYBRD(s0_ss_init.copy())
            elif solver == 'NLEQ2':
                state_species, self.__StateOK__ = self.NLEQ2(s0_ss_init.copy())
            elif solver == 'FINTSLV':
                state_species, self.__StateOK__ = self.FINTSLV(s0_ss_init.
                    copy())
            if numpy.isscalar(state_species):
                state_species = numpy.array([state_species], 'd')
            tmp = self._EvalODE(state_species.copy(), Vtemp)
            state_flux = self.__vvec__
            if self.__HAS_RATE_RULES__:
                state_species, rrules = numpy.split(state_species, [self.
                    Nrmatrix.shape[0]])
            if (state_species < 0.0).any():
                self.__StateOK__ = False
                if self.__settings__['mode_state_mesg']:
                    print('WARNING!! Negative concentrations detected.')
            if self.__StateOK__:
                break
            elif self.mode_solver_fallback and self.__settings__[
                'solver_switch_warning']:
                slv_idx = available_solvers.index(solver)
                if slv_idx != len(available_solvers) - 1:
                    print('INFO: STATE is switching to {} solver.'.format(
                        available_solvers[slv_idx + 1]))
                else:
                    print('INFO: STATE calculation failed!')
        if STATE_XOUT:
            STATE_xdata = self._EvalExtraData(self.STATE_extra_output)
        if self.__HAS_MOIETY_CONSERVATION__ == True:
            state_species = self.Fix_S_indinput(state_species, amounts=True)
        if self.__HAS_COMPARTMENTS__ and self.__KeyWords__['Output_In_Conc']:
            self.state_species = self._SpeciesAmountToConc(state_species.copy()
                )
        else:
            self.state_species = state_species
        self.state_flux = state_flux
        del state_species, state_flux
        if (self.state_flux == self.__settings__['mach_floateps']).any():
            print(
                '\nWARNING: extremely small flux detected! proceed with caution:'
                )
            print(self.state_flux)
        if not self.__StateOK__:
            print(
                """
***
WARNING: invalid steady state solution (species concentrations and fluxes)
***
"""
                )
            if self.__settings__['mode_state_nan_on_fail']:
                self.state_species[:] = numpy.nan
                self.state_flux[:] = numpy.nan
        self.SetStateSymb(self.state_flux, self.state_species)
        self.data_sstate = StateDataObj()
        self.data_sstate.setSpecies(self.state_species, self.__species__)
        self.data_sstate.setFluxes(self.state_flux, self.__reactions__)
        if self.__HAS_RATE_RULES__:
            self.data_sstate.setRules(rrules, self.__rate_rules__)
        if STATE_XOUT:
            self.data_sstate.setXData(STATE_xdata, lbls=self.STATE_extra_output
                )
            del STATE_xdata
        self.data_sstate.IS_VALID = self.__StateOK__
        if self.__settings__['display_debug'] == 1:
            print('self.state_species')
            print(self.state_species)
            print('self.state_flux')
            print(self.state_flux)

    def SetStateSymb(self, flux, metab):
        """
    Sets the individual steady-state flux and concentration attributes.

    This method assigns the steady-state flux values to class attributes named 
    `mod.J_<reaction>` and the steady-state concentrations to `mod.<species>_ss`.
    
    Parameters
    ----------
    flux : array-like
        The steady-state flux array. Each element corresponds to the flux of a specific reaction.
        
    metab : array-like
        The steady-state concentration array. Each element represents the concentration 
        of a specific species in the steady state.

    Examples
    --------
    Suppose you have a `PysMod` instance `mod` with predefined species and reactions:
    
    >>> mod = PysMod()
    
    After running a steady-state analysis, you obtain steady-state fluxes and concentrations:
    
    >>> steady_state_flux = [0.5, 1.2, 0.8]  # Corresponding to reactions A, B, C
    >>> steady_state_conc = [1.0, 0.75]      # Corresponding to species X, Y
    
    You can set these steady-state attributes using the `SetStateSymb` method:
    
    >>> mod.SetStateSymb(steady_state_flux, steady_state_conc)
    
    This will set the following attributes in `mod`:
    - `mod.J_A = 0.5`
    - `mod.J_B = 1.2`
    - `mod.J_C = 0.8`
    - `mod.X_ss = 1.0`
    - `mod.Y_ss = 0.75`
    """
        for x in range(0, len(self.state_species)):
            setattr(self, self.__species__[x] + '_ss', metab[x])
        for x in range(0, len(self.state_flux)):
            setattr(self, 'J_' + self.__reactions__[x], flux[x])

    def __Build_Tvec__(self, amounts=True):
        """
    Constructs vectors of conserved moiety totals from the `__inspec__` attribute.
    
    This method computes and assigns vectors representing the totals of conserved moieties 
    based on the boolean argument `amounts`. If `amounts` is `True`, the vectors represent 
    amounts; otherwise, they represent concentrations. The resulting vectors are stored in 
    `self.__tvec_a__` for amounts and `self.__tvec_c__` for concentrations.

    Parameters
    ----------
    amounts : bool, optional
        If `True`, constructs vectors based on moiety amounts. If `False`, constructs vectors
        based on moiety concentrations. The default is `True`.

    Notes
    -----
    This function should only be used if `self.__HAS_MOIETY_CONSERVATION__` is `True`, 
    indicating that the model involves moiety conservation laws.

    Example
    -------
    >>> model = PysMod()
    >>> model.__HAS_MOIETY_CONSERVATION__ = True
    >>> model.__inspec__ = numpy.array([10.0, 5.0, 1.0])
    >>> model.__reordered_lcons = numpy.identity(3)
    >>> model.__Build_Tvec__(amounts=True)
    >>> print("Amount vector:", model.__tvec_a__)
    Amount vector: [10.0, 5.0, 1.0]
    >>> print("Concentration vector:", model.__tvec_c__)
    Concentration vector: [10.0, 5.0, 1.0]
    """
        if self.__HAS_MOIETY_CONSERVATION__ == True:
            if amounts:
                self.__tvec_a__ = numpy.add.reduce(copy.copy(self.
                    __reordered_lcons) * self.__inspec__, 1)
                self.__tvec_c__ = numpy.add.reduce(copy.copy(self.
                    __reordered_lcons) * self._SpeciesAmountToConc(self.
                    __inspec__), 1)
            else:
                self.__tvec_c__ = numpy.add.reduce(copy.copy(self.
                    __reordered_lcons) * self.__inspec__, 1)
                self.__tvec_a__ = numpy.add.reduce(copy.copy(self.
                    __reordered_lcons) * self._SpeciesConcToAmount(self.
                    __inspec__), 1)

    def showConserved(self, File=None, screenwrite=1, fmt='%2.3f'):
        """
    Print the moiety conserved cycles present in the system.

    This method displays or writes to a file the summary of moiety conservation
    cycles that occur within the biochemical system modeled by the class. It
    processes the conservation matrix to enumerate these cycles and formats them
    for output.

    Parameters
    ----------
    File : file object, optional
        An open writable Python file object. If provided, the results are written
        to the file instead of or in addition to being printed to the console.
        Default is `None`.
    screenwrite : int, optional
        If set to 1, the results are printed to the console. If set to 0, no
        console output is performed. Default is `1`.
    fmt : str, optional
        A format string to format the numerical coefficients in the conservation
        relationships. Default is `'%2.3f'`.

    Notes
    -----
    The method relies on the class having a properly configured conservation matrix
    and is dependent on the presence of moiety conservation within the system.

    Examples
    --------
    Consider `mod` to be an instance of `PysMod` with moiety conservation computed:

    >>> mod.showConserved()
    Conserved relationships
    + {1.000}A - {2.000}B = 0.500

    Writing the output to a file:

    >>> with open('conserved_relationships.txt', 'w') as f:
    ...     mod.showConserved(File=f, screenwrite=0)
    """
        if self.__HAS_MOIETY_CONSERVATION__ == True:
            Tlist = list(range(0, len(self.__tvec_a__)))
            if Tlist != []:
                ConSumPstr = ''
                for x in range(0, len(Tlist)):
                    for y in range(0, len(self.conservation_matrix_col)):
                        if self.conservation_matrix[Tlist[x], y] > 0.0:
                            coeff = fmt % abs(self.conservation_matrix[
                                Tlist[x], y])
                            met = self._D_s_Order[self.
                                conservation_matrix_col[y]]
                            ConSumPstr += ' + {' + coeff + '}' + met.replace(
                                'self.', '')
                        elif self.conservation_matrix[Tlist[x], y] < 0.0:
                            coeff = fmt % abs(self.conservation_matrix[
                                Tlist[x], y])
                            met = self._D_s_Order[self.
                                conservation_matrix_col[y]]
                            ConSumPstr += ' - {' + coeff + '}' + met.replace(
                                'self.', '')
                    if self.__HAS_COMPARTMENTS__ and self.__KeyWords__[
                        'Output_In_Conc']:
                        ConSumPstr += ' = ' + self.__settings__[
                            'mode_number_format'] % self.__tvec_c__[Tlist[x]
                            ] + '\n'
                    else:
                        ConSumPstr += ' = ' + self.__settings__[
                            'mode_number_format'] % self.__tvec_a__[Tlist[x]
                            ] + '\n'
            self.conserved_sums = ConSumPstr
        else:
            self.conserved_sums = 'No moiety conservation'
        if File != None:
            print('\nConserved relationships')
            File.write('\n## Conserved relationships\n')
            File.write(self.conserved_sums)
        elif screenwrite:
            print('\nConserved relationships')
            print(self.conserved_sums)

    def showFluxRelationships(self, File=None):
        """
    showFluxRelationships(File=None)

    Print the flux relationships present in the system.

    Parameters
    ----------
    File : file-like object, optional
        An open writable Python file object where the flux relationships
        will be written. If `None`, the results are printed to the console.
        The default is `None`.

    Notes
    -----
    This method computes and prints the flux relationships derived from the
    system's matrix of zero-flux relationships (`__kzeromatrix__`). These
    relationships are presented by iterating over each row and displaying
    the respective relationships with their coefficients from the matrix.

    Examples
    --------
    Assume an instance `model` of the `PysMod` class with properly initialized
    matrices. To display the flux relationships:

    >>> model.showFluxRelationships()

    To write the flux relationships to a file:

    >>> with open('flux_relationships.txt', 'w') as f:
    ...     model.showFluxRelationships(File=f)
    """
        Ostr = ''
        for row in range(self.__kzeromatrix__.shape[0]):
            Ostr += '{} ='.format(self.reactions[self.kzeromatrix_row[row]])
            for col in range(self.__kzeromatrix__.shape[1]):
                if self.__kzeromatrix__[row, col] != 0.0:
                    if self.__kzeromatrix__[row, col] > 0.0:
                        Ostr += ' + {%2.2f}%s' % (abs(self.__kzeromatrix__[
                            row, col]), self.reactions[self.kzeromatrix_col
                            [col]])
                    else:
                        Ostr += ' - {%2.2f}%s' % (abs(self.__kzeromatrix__[
                            row, col]), self.reactions[self.kzeromatrix_col
                            [col]])
            Ostr += '\n'
        if File != None:
            print('\nFlux relationships')
            File.write('\n## Flux relationships\n')
            File.write(Ostr)
        else:
            print('\nFlux relationships')
            print(Ostr)

    def Fix_S_fullinput(self, s_vec, amounts=True):
        """
    Evaluate dependent species using the full concentration vector.

    This method updates the provided concentration vector `s_vec` by
    calculating the values of dependent species based on a pre-defined
    stoichiometric matrix. It supports two different modes based on
    the `amounts` parameter.

    Parameters
    ----------
    s_vec : array-like
        A full-length concentration vector that includes both
        independent and dependent species.
    amounts : bool, optional
        A flag indicating whether the calculation should be done in
        terms of amounts (`True`) or concentrations (`False`). Default
        is `True`.

    Returns
    -------
    array-like
        Updated concentration vector with evaluated dependent species.

    Examples
    --------
    Suppose you have an instance of `PysMod` class called `model` and a
    concentration vector `concentration_vector` available:

    >>> concentration_vector = [1.0, 2.0, 0.0, 0.0, 0.0]
    >>> updated_vector = model.Fix_S_fullinput(concentration_vector, amounts=True)
    >>> print(updated_vector)
    [1.0, 2.0, <calculated_value_1>, <calculated_value_2>, ...]

    Notes
    -----
    This method assumes that `self.__SI__`, `self.__lzeromatrix_col`, 
    `self.__lzeromatrix_row`, `self.__tvec_a__`, and `self.__tvec_c__` 
    are all pre-defined and initialized before calling this method.

    """
        for x in range(0, len(self.lzeromatrix_col)):
            self.__SI__[x] = s_vec[self.lzeromatrix_col[x]]
        for x in range(0, len(self.lzeromatrix_row)):
            if amounts:
                s_vec[self.lzeromatrix_row[x]] = self.__tvec_a__[x
                    ] + numpy.add.reduce(self.__lzeromatrix__[x, :] * self.
                    __SI__)
            else:
                s_vec[self.lzeromatrix_row[x]] = self.__tvec_c__[x
                    ] + numpy.add.reduce(self.__lzeromatrix__[x, :] * self.
                    __SI__)
        return s_vec

    def Fix_S_indinput(self, s_vec, amounts=True):
        """
    Evaluate and return a full concentration vector from independent species.

    Given a vector of independent species, this method evaluates and returns
    a full concentration vector by calculating the dependent species concentrations.
    The method uses either `self.__tvec_a__` or `self.__tvec_c__` depending on
    the `amounts` parameter.

    Parameters
    ----------
    s_vec : array-like
        A vector of concentrations for independent species.

    amounts : bool, optional
        A flag indicating whether to use `self.__tvec_a__` (default) or `self.__tvec_c__`
        for calculations. If `True`, `self.__tvec_a__` is used; otherwise, `self.__tvec_c__`
        is used.

    Returns
    -------
    ndarray
        A full concentration vector including both independent and dependent species.

    Examples
    --------
    >>> model = PysMod()
    >>> independent_species = [1.0, 2.0, 3.0] 
    >>> full_concentration_vector = model.Fix_S_indinput(independent_species, amounts=True)
    >>> print(full_concentration_vector)
    [1.0, 2.0, 3.0, 4.0, 5.0, ...]  # an example of the resulting full concentration vector
    """
        for x in range(len(s_vec)):
            self.__SALL__[self.nrmatrix_row[x]] = s_vec[x]
        for x in range(len(self.lzeromatrix_row)):
            if amounts:
                self.__SALL__[self.lzeromatrix_row[x]] = self.__tvec_a__[x
                    ] + numpy.add.reduce(self.__lzeromatrix__[x, :] * s_vec)
            else:
                self.__SALL__[self.lzeromatrix_row[x]] = self.__tvec_c__[x
                    ] + numpy.add.reduce(self.__lzeromatrix__[x, :] * s_vec)
        return copy.copy(self.__SALL__)

    def __FluxGen__(self, s):
        """
    Calculate the flux vector for given species concentrations.

    **Deprecated**: This method is primarily used by **PITCON**.

    Parameters
    ----------
    s : array-like
        A species vector/array that represents the concentrations of species 
        for which the flux needs to be calculated. It can be a 1D array for a 
        single set of concentrations or a 2D array for multiple sets.

    Returns
    -------
    numpy.ndarray
        An array representing the calculated flux values. The shape of the 
        output depends on the input `s`. If `s` is a 2D array, the output 
        will be a 2D array with the same number of rows.

    Notes
    -----
    This method is used to evaluate the fluxes in a system of reactions 
    based on the species concentrations provided. If moiety conservation 
    is present, a different evaluation path is taken.

    Examples
    --------
    >>> pys_mod_instance = PysMod()
    >>> species_array = numpy.array([0.1, 0.2, 0.3])
    >>> flux = pys_mod_instance.__FluxGen__(species_array)
    >>> print(flux)
    array([...])  # The output will be a flux vector

    >>> species_matrix = numpy.array([[0.1, 0.2, 0.3],
    ...                              [0.4, 0.5, 0.6]])
    >>> flux_matrix = pys_mod_instance.__FluxGen__(species_matrix)
    >>> print(flux_matrix)
    array([...])  # The output will be a matrix of flux vectors
    """
        Vtemp = numpy.zeros(self.__Nshape__[1], 'd')
        try:
            s.shape[0]
            Vout = numpy.zeros((len(s), len(self.__vvec__)), 'd')
            for x in range(len(s)):
                if self.__HAS_MOIETY_CONSERVATION__ == True:
                    Vout[x, :] = self._EvalREq2_alt(s[x, :], Vtemp)
                elif self.__HAS_MOIETY_CONSERVATION__ == False:
                    Vout[x, :] = self._EvalREq(s[x, :], Vtemp)
        except:
            if self.__HAS_MOIETY_CONSERVATION__ == True:
                Vout = self._EvalREq2_alt(s, Vtemp)
            elif self.__HAS_MOIETY_CONSERVATION__ == False:
                Vout = self._EvalREq(s, Vtemp)
        return Vout

    def FluxGenSim(self, s):
        """
    **Deprecated**

    This method is deprecated and may be removed in future versions. 
    Please refer to the documentation for updated methods or alternatives.

    Parameters
    ----------
    s : array-like
        A vector or array of species.

    Examples
    --------
    Deprecated usage; for demonstration purposes only:

    >>> model = PysMod()
    >>> species_vector = [0.1, 0.2, 0.3]
    >>> model.FluxGenSim(species_vector)  # No output, method is deprecated

    Notes
    -----
    This method does not perform any operations and is included only for compatibility reasons. 
    Users are encouraged to avoid using deprecated functionality in new code.
    """
        pass

    def ParGenSim(self):
        """
    **Deprecated**

    This method is deprecated and may be removed in a future version. 
    It currently does not perform any operations and exists for legacy code compatibility.

    Examples
    --------
    >>> obj = PysMod()
    >>> obj.ParGenSim()

    Notes
    -----
    Use of this method is discouraged as it is deprecated and does not perform any functionality. 
    Consider using updated methods or contacting the library maintainers for alternatives.
    """
        pass

    def Fix_Sim(self, metab, flux=0, par=0):
        """
    **Deprecated**

    This method is deprecated and should not be used in new code. It was previously used for simulating 
    fixes in metabolic networks by modifying both fluxes and parameters based on the provided 
    metabolites.

    Parameters
    ----------
    metab : str
        The identifier of the metabolite to be used in the simulation.
    flux : int, optional
        The flux value to be set in the simulation, by default 0.
    par : int, optional
        The parameter value to be set in the simulation, by default 0.

    Notes
    -----
    This method is deprecated and may be removed in future versions. Users are encouraged to find 
    alternative methods to achieve equivalent functionality.

    Examples
    --------
    Deprecated and should not be used. For historical reference only.

    >>> model = PysMod()
    >>> model.Fix_Sim(metab='ATP', flux=5, par=2)
    
    In the above example, 'ATP' is the metabolite whose simulation is being fixed with a flux and 
    parameter value of 5 and 2, respectively. However, this method is deprecated.
    """
        pass

    def EvalEvar(self, input=None, input2=None):
        """
    Calculate reaction elasticities towards the variable species.

    This function evaluates the elasticities of reactions with respect to
    variable species concentrations. The elasticities are calculated within
    the framework of Metabolic Control Analysis (MCA). Both `input` and
    `input2` must be valid steady-state solutions for the calculations to
    be accurate. If either is missing, the function will use the last
    known state values automatically. The elasticities are then scaled
    using the provided vectors.

    Parameters
    ----------
    input : array-like, optional
        The species concentration vector. Defaults to `None`, in which case
        the last state values are used.
    input2 : array-like, optional
        The reaction rate vector. Defaults to `None`, in which case the last
        state fluxes are used.

    Settings (mod.__settings__)
    ---------------------------
    elas_evar_upsymb : bool, default=True
        If True, attach individual elasticity symbols to the model instance.
    elas_zero_conc_fix : bool, default=False
        If True, zero concentrations detected in a steady-state solution are replaced
        by a very small number.
    elas_zero_flux_fix : bool, default=False
        If True, zero fluxes detected in a steady-state solution are replaced
        by a very small number.
    elas_scaling_div0_fix : bool, default=False
        If True, infinite values detected after scaling are set to zero.
    display_debug : int, default=0
        If set to 1, detailed debug information is printed during execution.

    Raises
    ------
    AssertionError
        If the length of `input` does not match the number of species, or if the
        length of `input2` does not match the number of reactions.

    Examples
    --------
    >>> model = PysMod()
    >>> species_conc = [0.02, 0.03, 0.04]
    >>> reaction_rates = [1.5, 1.2, 1.3]
    >>> model.EvalEvar(input=species_conc, input2=reaction_rates)
    >>> print(model.elas_var)
    # Outputs a matrix of elasticities scaled by the given inputs.
    """
        if input is None or input2 is None:
            input = self.state_species
            input2 = self.state_flux
        else:
            assert len(input) == len(self.species
                ), 'length error this array must have ' + str(len(self.species)
                ) + ' elements'
            assert len(input2) == len(self.reactions
                ), 'length error this array must have ' + str(len(self.
                reactions)) + ' elements'
            exec(self.__CDvarUpString__)
        if self.__settings__['elas_zero_flux_fix']:
            val = self.__settings__['mach_floateps'] ** 2
            for x in range(len(input2)):
                if abs(input2[x]) <= val:
                    print('Info: zero flux detected: J_{} set to {}'.format
                        (self.reactions[x], val))
                    input2[x] = val
        if self.__settings__['elas_zero_conc_fix']:
            val = self.__settings__['mach_floateps'] ** 2
            for x in range(len(input)):
                if abs(input[x]) <= val:
                    print('Info: zero concentration detected: {} set to {}'
                        .format(self.species[x], val))
                    input[x] = val
        if self.__settings__['display_debug'] == 1:
            print('\nVarinput = ' + repr(input) + '\n')
        self.__evmatrix__ = numpy.zeros((len(self.__reactions__), len(self.
            __species__)), 'd')
        self.ScaleKL(input, input2)
        self.elas_var_row = tuple(copy.copy(self.__reactions__))
        col = copy.copy(self._D_s_Order)
        for x in range(len(self._D_s_Order)):
            col[x] = col[x].replace('self.', '')
        self.elas_var_col = tuple(col)
        del col
        for react in range(len(self.__reactions__)):
            if self.__settings__['display_debug'] == 1:
                print('\nReaction: ' + self.__reactions__[react])
            for met in range(len(self.__species__)):
                countV = 0
                countMet = 0
                try:
                    """
                    I borrowed the idea of scaling the stepsize to So from Herbert Sauro's Jarnac TModel.uEEOp function
                    brett - 2004-08-18
                    """
                    hstep = input[met] * self.__settings__[
                        'mode_elas_deriv_factor']
                    if abs(hstep) < self.__settings__['mode_elas_deriv_min']:
                        hstep = self.__settings__['mode_elas_deriv_min']
                    a = ScipyDerivative(self.__num_deriv_function__, input[
                        met], order=self.__settings__[
                        'mode_elas_deriv_order'], dx=hstep, n=1, args=(self
                        .__reactions__[react], self.__species__[met]))
                except Exception as ex:
                    print(ex)
                    print('\nINFO: Elasticity evaluation failure in ', self
                        .__reactions__[react], self.__species__[met])
                    print('Elasticity has been set to zero')
                    print(
                        'A stepsize that is too large might cause this ... try decreasing the factor and or min stepsize'
                        )
                    print('Keep in mind machine precision, is', self.
                        __settings__['mach_floateps'], ' and if min stepsize')
                    print(
                        'becomes too small numeric error can become significant'
                        )
                    a = 0.0
                if abs(a) >= self.__settings__['mach_floateps']:
                    if self.__settings__['display_debug'] == 1:
                        print('species: ' + self.__species__[met])
                        print('--> d(' + self.__reactions__[react] + ')d(' +
                            self.__species__[met] + ') = ' + str(a))
                    self.__evmatrix__[react, met] = a
        if self.__settings__['elas_evar_remap']:
            eval(compile(self.__remaps, '__remaps', 'exec'))
        self.elas_var_u = self.__evmatrix__
        if self.__settings__['mode_mca_scaled'] == 1:
            Ds = numpy.zeros((len(input), len(input)), 'd')
            Dj = numpy.zeros((len(input2), len(input2)), 'd')
            for x in range(0, len(input)):
                Ds[x, x] = input[x]
            for x in range(0, len(input2)):
                Dj[x, x] = 1.0 / input2[x]
                if self.__settings__['elas_scaling_div0_fix'] and numpy.isinf(
                    Dj[x, x]):
                    print(
                        'Infinite elasticity detected during scaling setting to zero ({})'
                        .format(self.reactions[x]))
                    Dj[x, x] = 1e-16
            Dj_e = numpy.dot(Dj, self.__evmatrix__)
            self.__evmatrix__ = numpy.dot(Dj_e, Ds)
            self.elas_var = self.__evmatrix__
        else:
            self.elas_var = None
        if self.__settings__['elas_evar_upsymb'] == 1:
            r, c = self.__evmatrix__.shape
            output2 = ''
            for x in range(0, r):
                react = self.__reactions__[x]
                for y in range(0, c):
                    met = self._D_s_Order[y]
                    if self.__settings__['mode_mca_scaled']:
                        setattr(self, 'ec' + react + '_' + met.replace(
                            'self.', ''), self.__evmatrix__[x, y])
                    else:
                        setattr(self, 'uec' + react + '_' + met.replace(
                            'self.', ''), self.__evmatrix__[x, y])
            self.ec = BagOfStuff(self.__evmatrix__, self.elas_var_row, self
                .elas_var_col)
            self.ec.load()
            if self.__settings__['mode_mca_scaled'] == 1:
                self.ec.scaled = True
            else:
                self.ec.scaled = False
            if self.__settings__['display_debug'] == 1:
                print('\n\n********************************\n')
                print(output2)
                print('\n********************************\n')
        else:
            pass
        if self.__settings__['display_debug'] == 1:
            print('\ne_vmatrix')
            print(repr(self._D_s_Order).replace('self.', ''))
            print(self.__reactions__)
            print(self.__evmatrix__)

    def CleanNaNsFromArray(self, arr, replace_val=0.0):
        """
    Scan a matrix for NaN values and replace them with a specified value.

    This method iterates over a given 2D array, checks for NaN values, and replaces each occurrence with the specified `replace_val`, which defaults to 0.0.

    Parameters
    ----------
    arr : numpy.ndarray
        The array to be cleaned. It should be a 2D numpy array where NaNs need to be replaced.
    replace_val : scalar, optional
        The value to replace NaNs with. Default is 0.0.

    Examples
    --------
    >>> import numpy as np
    >>> from pysmod import PysMod  # Assuming the class is imported from a module named pysmod
    >>> obj = PysMod()
    >>> data = np.array([[1, 2, np.nan], [4, np.nan, 6]])
    >>> print("Before cleaning:", data)
    Before cleaning: [[ 1.  2. nan]
                      [ 4. nan  6.]]
    >>> obj.CleanNaNsFromArray(data, replace_val=0)
    >>> print("After cleaning:", data)
    After cleaning: [[1. 2. 0.]
                     [4. 0. 6.]]
    """
        nantest = numpy.isnan(arr)
        if nantest.any():
            for r in range(nantest.shape[0]):
                if nantest[r].any():
                    for c in range(nantest.shape[1]):
                        if numpy.isnan(arr[r, c]):
                            arr[r, c] = replace_val

    def EvalEpar(self, input=None, input2=None):
        """
    Calculate reaction elasticities towards the parameters.

    Both inputs (`input`, `input2`) should be valid steady-state
    solutions for Metabolic Control Analysis (MCA) and provided in the
    correct order for them to be used. If either or both are missing,
    the last state values are used automatically.
    Elasticities are scaled using `input` and `input2`.

    Parameters
    ----------
    input : array-like, optional
        Species concentration vector. If not provided, the last recorded
        state species are used.
    input2 : array-like, optional
        Reaction rate vector. If not provided, the last recorded state
        flux is used.

    Settings
    --------
    - elas_epar_upsymb : int, optional, default is 1
        Attaches individual elasticity symbols to the model instance.
    - elas_scaling_div0_fix : bool, optional, default is False
        If NaN's are detected in the variable and parameter elasticity
        matrix, they are replaced with zero.

    Raises
    ------
    AssertionError
        If the lengths of `input` or `input2` do not match the expected
        sizes corresponding to the species or reactions.

    Example
    -------
    >>> mod = PysMod()
    >>> species_concentrations = [0.5, 0.3, 0.2]
    >>> reaction_rates = [1.0, 0.5]
    >>> mod.EvalEpar(input=species_concentrations, input2=reaction_rates)
    """
        if input is None or input2 is None:
            input = self.state_species
            input2 = self.state_flux
        else:
            assert len(input) == len(self.species
                ), 'length error this array must have ' + str(len(self.species)
                ) + ' elements'
            assert len(input2) == len(self.reactions
                ), 'length error this array must have ' + str(len(self.
                reactions)) + ' elements'
            exec(self.__CDvarUpString__)
        parVal_hold = numpy.zeros(len(self.__parameters__), 'd')
        if self.__settings__['display_debug'] == 1:
            print('\nParinput = ' + repr(input) + '\n')
            print('\nparVal_hold1')
            print(parVal_hold)
        exec(self.__par_map2storeC)
        parVal2 = copy.copy(parVal_hold)
        self.__epmatrix__ = numpy.zeros((len(self.__reactions__), len(self.
            __parameters__)), 'd')
        self.elas_par_row = tuple(copy.copy(self.__reactions__))
        self.elas_par_col = tuple(self.__parameters__)
        for react in range(len(self.__reactions__)):
            if self.__settings__['display_debug'] == 1:
                print('\nReaction: ' + self.__reactions__[react])
            for par in range(len(self.__parameters__)):
                countV = 0
                countPar = 0
                try:
                    """ I got the idea of scaling the stepsize to So from Herbert Sauro's Jarnac TModel.uEEOp function

                    brett - 20040818"""
                    hstep = getattr(self, self.__parameters__[par]
                        ) * self.__settings__['mode_elas_deriv_factor']
                    if abs(hstep) < self.__settings__['mode_elas_deriv_min']:
                        hstep = self.__settings__['mode_elas_deriv_min']
                    a = ScipyDerivative(self.__num_deriv_function__,
                        getattr(self, self.__parameters__[par]), order=self
                        .__settings__['mode_elas_deriv_order'], dx=hstep, n
                        =1, args=(self.__reactions__[react], self.
                        __parameters__[par]))
                except Exception as ex:
                    print('\nNumeric derivative evaluation failure in ',
                        self.__reactions__[react], self.__parameters__[par])
                    print('Elasticity has been set to NaN')
                    print(
                        'A stepsize that is too large might cause this ... try decreasing the factor and or min stepsize'
                        )
                    print('Keep in mind machine precision, is', self.
                        __settings__['mach_floateps'], ' and if min stepsize')
                    print(
                        'becomes too small numeric error can become significant'
                        )
                    print(ex)
                    a = numpy.nan
                if numpy.isnan(a) or abs(a) > self.__settings__['mach_floateps'
                    ]:
                    if self.__settings__['display_debug'] == 1:
                        print('parameter: ' + self.__parameters__[par])
                        print('--> d(' + self.__reactions__[react] + ')d(' +
                            self.__parameters__[par] + ') = ' + repr(a))
                    self.__epmatrix__[react, par] = a
        self.elas_par_u = self.__epmatrix__
        exec(self.__par_remapC)
        if self.__settings__['display_debug'] == 1:
            print('\nparVal_hold2')
            print(parVal_hold)
        if self.__settings__['mode_mca_scaled'] == 1:
            Dp = numpy.zeros((len(self.__parameters__), len(self.
                __parameters__)), 'd')
            Dj = numpy.zeros((len(input2), len(input2)), 'd')
            for x in range(0, len(self.__parameters__)):
                Dp[x, x] = parVal_hold[x]
            for x in range(0, len(input2)):
                Dj[x, x] = 1.0 / input2[x]
                if self.__settings__['elas_scaling_div0_fix'] and numpy.isinf(
                    Dj[x, x]):
                    print(
                        'Infinite elasticity detected during scaling setting to zero ({})'
                        .format(self.reactions[x]))
                    Dj[x, x] = 1e-16
            Dj_e = numpy.dot(Dj, self.__epmatrix__)
            self.__epmatrix__ = numpy.dot(Dj_e, Dp)
            self.elas_par = self.__epmatrix__
        else:
            self.elas_par = None
        del parVal_hold
        if self.__settings__['elas_epar_upsymb'] == 1:
            r, c = self.__epmatrix__.shape
            output2 = ''
            for x in range(0, r):
                react = self.__reactions__[x]
                for y in range(0, c):
                    met = self.__parameters__[y]
                    if self.__settings__['mode_mca_scaled']:
                        setattr(self, 'ec' + react + '_' + met.replace(
                            'self.', ''), self.__epmatrix__[x, y])
                    else:
                        setattr(self, 'uec' + react + '_' + met.replace(
                            'self.', ''), self.__epmatrix__[x, y])
            self.ecp = BagOfStuff(self.__epmatrix__, self.elas_par_row,
                self.elas_par_col)
            self.ecp.load()
            if self.__settings__['mode_mca_scaled'] == 1:
                self.ecp.scaled = True
            else:
                self.ecp.scaled = False
            if self.__settings__['display_debug'] == 1:
                print('\n\n********************************\n')
                print(output2)
                print('\n********************************\n')
        if self.__settings__['display_debug'] == 1:
            print('\ne_pmatrix')
            print(self.__parameters__)
            print(self.__reactions__)
            print(self.__epmatrix__)

    def showEvar(self, File=None):
        """
    Write out all variable elasticities as 'LaTeX' formatted strings.

    This method writes out all variable elasticities as 'LaTeX' formatted strings
    and prints them to the console. Alternatively, if a writable file object 
    is provided, it writes the results to that file.

    Parameters
    ----------
    File : file-like object, optional
        An open writable Python file object. If provided, the results are 
        written to this file instead of being printed to the console. 
        Default is None.

    Notes
    -----
    The elasticities are presented as either scaled or unscaled, depending 
    on the settings in `self.__settings__`. The `mode_number_format` setting 
    determines the numerical display format of the elasticities.

    If no variable elasticities are calculated (i.e., `EvalEvar()` 
    has not been run), an informative message will be displayed.

    Examples
    --------
    Display the variable elasticities on the console:

    >>> model = PysMod()
    >>> model.showEvar()

    Write the variable elasticities to a file:

    >>> with open('elasticities.txt', 'w') as f:
    ...     model.showEvar(File=f)

    """
        try:
            r, c = self.__evmatrix__.shape
            evar_output = ''
            for x in range(0, r):
                react = self.__reactions__[x]
                evar_output += '\n' + repr(self.__reactions__[x]) + '\n'
                for y in range(0, c):
                    rtemp = self.__settings__['mode_number_format'
                        ] % self.__evmatrix__[x, y]
                    met = self._D_s_Order[y]
                    if self.__evmatrix__[x, y] != 0.0:
                        if self.__settings__['mode_mca_scaled']:
                            elas = ('\\ec{' + react + '}{' + met + '} = ' +
                                rtemp)
                        else:
                            elas = ('\\uec{' + react + '}{' + met + '} = ' +
                                rtemp)
                        evar_output += elas + '\n'
            evar_output = evar_output.replace('self.', '')
        except:
            evar_output = (
                'No variable elasticities - run EvalEvar() to calculate')
        if File != None:
            print('\nspecies elasticities')
            File.write('\n## species elasticities\n')
            File.write(evar_output)
        else:
            print('\nspecies elasticities')
            print(evar_output)

    def showEpar(self, File=None):
        """
    Write out all nonzero parameter elasticities as 'LaTeX' formatted strings, 
    alternatively write them to a file.

    Parameters
    ----------
    File : file-like object, optional
        An open writable Python file object. If provided, the output will be 
        written to this file. If None (the default), the output will be printed 
        to the console.

    Notes
    -----
    Parameter elasticities are calculated using the `EvalEpar` method, and the 
    results are presented in LaTeX format. This method outputs only nonzero 
    elasticities for better readability.

    This method relies on a matrix of elasticities (`__epmatrix__`) that has been 
    populated through prior computation. If `__epmatrix__` is uninitialized, 
    the method will output a message prompting the user to run `EvalEpar()`.

    The LaTeX output differentiates between scaled elasticities (`\\\\ec{reaction}{parameter}`) 
    and unscaled elasticities (`\\\\uec{reaction}{parameter}`), based on the module's 
    settings configuration (`mode_mca_scaled`).

    Examples
    --------
    >>> model = PysMod()
    >>> model.EvalEpar(input=species_vector, input2=rate_vector)
    >>> model.showEpar()
    
    Alternatively, write to a file:
    
    >>> with open('elasticities_output.txt', 'w') as file:
    >>>     model.showEpar(File=file)
    """
        try:
            r, c = self.__epmatrix__.shape
            epar_output = ''
            for x in range(0, r):
                react = self.__reactions__[x]
                epar_output += '\n' + repr(self.__reactions__[x]) + '\n'
                for y in range(0, c):
                    rtemp = self.__settings__['mode_number_format'
                        ] % self.__epmatrix__[x, y]
                    met = self.__parameters__[y]
                    if self.__settings__['mode_mca_scaled']:
                        elas = '\\ec{' + react + '}{' + met + '} = ' + rtemp
                    else:
                        elas = '\\uec{' + react + '}{' + met + '} = ' + rtemp
                    epar_output += elas + '\n'
            epar_output = epar_output.replace('self.', '')
        except:
            epar_output = (
                'No parameter elasticities - run EvalEpar() to calculate')
        if File != None:
            print('\nParameter elasticities')
            File.write('\n## Parameter elasticities\n')
            File.write(epar_output)
        else:
            print('\nParameter elasticities')
            print(epar_output)

    def showElas(self, File=None):
        """
    Print all elasticities to screen or file as 'LaTeX' compatible strings.
    
    This method writes all variable and parameter elasticities in 'LaTeX'
    formatted strings to the console or an optional file. It internally
    calls the `showEvar` and `showEpar` methods to obtain and format
    the elasticities.

    Parameters
    ----------
    File : file object, optional
        An open writable Python file object. If provided, the elasticities
        are written to this file. If `None`, the elasticities are printed
        to the console. The default is `None`.

    Examples
    --------
    To print all elasticities to the console:

    >>> model = PysMod()
    >>> model.showElas()

    To write all elasticities to a file:

    >>> with open('elasticities.txt', 'w') as f:
    >>>     model.showElas(f)

    Notes
    -----
    The method assumes that the elasticities have already been calculated 
    using the appropriate evaluation methods such as `EvalEvar` and 
    `EvalEpar`.
    """
        if File != None:
            self.showEvar(File)
            self.showEpar(File)
        else:
            self.showEvar()
            self.showEpar()

    def ScaleKL(self, input, input2):
        """
    Scale the K and L matrices with the current steady state or the user-provided input.

    This method scales the K and L matrices based on vectors of species concentrations and
    reaction rates. If either `input` or `input2` is `None`, the method uses the current 
    steady state values instead. The method also checks whether the provided inputs match
    the required dimensions of the species and reaction vectors.

    Parameters
    ----------
    input : array-like or None
        A vector of species concentrations. If `None`, the current steady state values are used.
    input2 : array-like or None
        A vector of reaction rates. If `None`, the current steady state values are used.

    Raises
    ------
    AssertionError
        If the length of `input` does not match the number of species.
    AssertionError
        If the length of `input2` does not match the number of reactions.

    Notes
    -----
    The method generates warning messages if zero species or null fluxes are detected.

    Examples
    --------
    >>> model = PysMod()
    >>> species_concentrations = [1.0, 2.0, 1.5, 0.5]
    >>> reaction_rates = [0.1, 0.2, 0.3]
    >>> model.ScaleKL(species_concentrations, reaction_rates)
    
    Use current steady state values:
    
    >>> model.ScaleKL(None, None)

    """
        lmat = copy.copy(self.__lmatrix__)
        if type(input) == type(None) or type(input2) == type(None):
            input = self.state_species
            input2 = self.state_flux
        else:
            assert len(input) == len(self.species
                ), 'length error this array must have ' + str(len(self.species)
                ) + ' elements'
            assert len(input2) == len(self.reactions
                ), 'length error this array must have ' + str(len(self.
                reactions)) + ' elements'
        s_1_scale = numpy.zeros((len(input), len(input)), 'd')
        for x in range(0, len(input)):
            try:
                s_1_scale[x, x] = 1.0 / input[self.lmatrix_row[x]]
            except:
                print('Zero species detected: ' + self.__species__[self.
                    lmatrix_row[x]] + '_ss = ' + repr(input[self.
                    lmatrix_row[x]]))
                s_1_scale[x, x] = 0.0
        Si_order = numpy.zeros(len(self.lmatrix_col))
        Si_scale = numpy.zeros((len(self.lmatrix_col), len(self.lmatrix_col
            )), 'd')
        for x in range(0, len(self.lmatrix_col)):
            Si_order[x] = self.lmatrix_col[x]
            Si_scale[x, x] = input[self.lmatrix_col[x]]
        L_Si = numpy.dot(lmat, Si_scale)
        L_scaled = numpy.dot(s_1_scale, L_Si)
        self.lmatrix_scaled = L_scaled
        kmat = copy.copy(self.__kmatrix__)
        j_1_scale = numpy.zeros((len(input2), len(input2)), 'd')
        for x in range(0, len(input2)):
            try:
                j_1_scale[x, x] = 1.0 / input2[self.kmatrix_row[x]]
            except:
                print('\nNull flux detected: ' + self.__reactions__[self.
                    kmatrix_row[x]] + '_ss = ' + repr(input2[self.
                    kmatrix_row[x]]))
                j_1_scale[x, x] = 0.0
        Ji_order = numpy.zeros(len(self.kmatrix_col))
        Ji_scale = numpy.zeros((len(self.kmatrix_col), len(self.kmatrix_col
            )), 'd')
        for x in range(0, len(self.kmatrix_col)):
            Ji_order[x] = self.kmatrix_col[x]
            Ji_scale[x, x] = input2[self.kmatrix_col[x]]
        K_Ji = numpy.dot(kmat, Ji_scale)
        K_scaled = numpy.dot(j_1_scale, K_Ji)
        self.kmatrix_scaled = K_scaled

    def __num_deriv_function__(self, x, react, met):
        """
    Evaluate the rate equation for a given reaction and metabolite.

    This system function is utilized by numeric perturbation methods to derive
    elasticities. It assigns the value `x` to the metabolite `met`, evaluates
    the reaction `react` for the rate `v`, and is designed to work with the
    `ScipyDerivative()` function using precompiled function strings.

    Parameters
    ----------
    x : float
        The value to assign to the metabolite `met`.
    react : str
        The name of the reaction whose rate equation is to be evaluated.
    met : str
        The metabolite, which will be temporarily set to `x` during the
        evaluation.

    Returns
    -------
    vout : float
        The evaluated rate of the reaction `react` after assigning `x` to `met`.

    Examples
    --------
    Consider a system with a reaction `R1` and a metabolite `A`.

    >>> model = PysMod()
    >>> reaction_rate = model.__num_deriv_function__(0.5, 'R1', 'A')
    >>> print(f"The rate of reaction R1 when metabolite A is set to 0.5 is {reaction_rate}")

    The above example demonstrates the use of `__num_deriv_function__` to
    compute the reaction rate for `R1` when the concentration of `A` is set
    to 0.5.
    
    Notes
    -----
    This function is not intended for direct use but rather as a utility within
    numeric differentiation routines.
    """
        bk = getattr(self, met)
        setattr(self, met, x)
        self.Forcing_Function()
        v = getattr(self, react)
        vout = v()
        setattr(self, met, bk)
        return vout

    def EvalCC(self):
        """
    Calculate the Metabolic Control Analysis (MCA) control coefficients using the current steady-state solution.

    This method computes control coefficients based on the steady-state configuration of the system's
    metabolic network. It generates both flux and concentration control coefficients and assigns them
    to the model instance depending on the settings.

    By setting `mod.__settings__["mca_ccj_upsymb"] = 1`, the flux control coefficients are attached to
    the model instance, while setting `mod.__settings__["mca_ccs_upsymb"] = 1` attaches the concentration
    control coefficients.

    The calculation might fail if there are NaN values present, which could be due to zero fluxes or 
    concentrations. The settings `elas_scaling_div0_fix` and `elas_zero_flux_fix` can help handle such cases.

    Parameters
    ----------
    None

    Returns
    -------
    None : NoneType
        This method performs calculations and updates internal attributes of the instance but does not 
        return a value.

    Example
    -------
    >>> # Assuming `mod` is an instance of PysMod
    >>> mod.EvalCC()
    >>> # Accessing computed control coefficients
    >>> print(mod.cc_flux)
    >>> print(mod.cc_conc)

    Notes
    -----
    Ensure the model's steady-state is correctly set before invoking this method. This method updates 
    the instance with scaled or unscaled control coefficients, as well as symbolic references to them, 
    based on the configuration set in `mod.__settings__`.
    """
        e2 = numpy.zeros((len(self.kmatrix_row), len(self.lmatrix_row)), 'd')
        for x in range(0, len(self.kmatrix_row)):
            for y in range(0, len(self.lmatrix_row)):
                e2[x, y] = self.elas_var_u[self.kmatrix_row[x], self.
                    lmatrix_row[y]]
        if self.__settings__['display_debug'] == 1:
            print(self.lmatrix_row)
            print(self.kmatrix_row)
            print('---')
            print(self.__nmatrix__.shape)
            print(self.nmatrix_row)
            print(self.nmatrix_col)
            print(self.nmatrix)
            print('---')
            print(self.__lmatrix__.shape)
            print(self.__lmatrix__)
            print('---')
            print(self.__kmatrix__.shape)
            print(self.__kmatrix__)
            print('---')
            print(e2.shape)
            print(e2)
        EL = -1.0 * numpy.dot(e2, self.__lmatrix__)
        KEL = numpy.concatenate((self.__kmatrix__, EL), 1)
        self.kel_unscaled = KEL
        go = 0
        try:
            Ci_u = scipy.linalg.inv(KEL)
            self.__structural__.MatrixFloatFix(Ci_u, val=1e-15)
            go = 1
        except Exception as ex:
            print(ex)
            print(
                """
INFO: K-EL matrix inversion failed this is possibly due to NaN values in the Elasticity matrix"""
                )
            print(
                'NaN elasticities can be caused by zero fluxes or concentrations look at the settings (in mod.__settings__)'
                )
            print("'elas_scaling_div0_fix' and 'elas_zero_flux_fix'")
            go = 0
        if go == 1 and not self.__settings__['mode_mca_scaled']:
            self.mca_ci = Ci_u
            self.__mca_CCi_unscaled = Ci_u
            del Ci_u
        elif go == 1 and self.__settings__['mode_mca_scaled']:
            Vi_l = self.__kmatrix__.shape[1]
            Vd_l = abs(self.__kmatrix__.shape[0] - self.__kmatrix__.shape[1])
            Si_l = self.__lmatrix__.shape[1]
            Sd_l = abs(self.__lmatrix__.shape[0] - self.__lmatrix__.shape[1])
            scale1_l = Vi_l + Si_l
            scale1 = numpy.zeros((scale1_l, scale1_l), 'd')
            for x in range(0, Vi_l):
                scale1[x, x] = 1.0 / self.state_flux[self.kmatrix_row[x]]
            for x in range(Vi_l, scale1_l):
                scale1[x, x] = 1.0 / self.state_species[self.lmatrix_row[x -
                    Vi_l]]
            scale2 = numpy.zeros((Vi_l + Vd_l, Vi_l + Vd_l), 'd')
            for x in range(0, len(self.state_flux)):
                scale2[x, x] = self.state_flux[self.kmatrix_row[x]]
            left = numpy.dot(scale1, Ci_u)
            Ci = numpy.dot(left, scale2)
            self.mca_ci = Ci
            del scale1_l, scale1, scale2, left, Ci, Ci_u
        Cirow = []
        Cicol = []
        for x in range(0, len(self.kmatrix_row)):
            Cicol.append(self.__reactions__[self.kmatrix_row[x]])
        for x in range(0, len(self.kmatrix_col)):
            Cirow.append(self.__reactions__[self.kmatrix_col[x]])
        for x in range(0, len(self.lmatrix_col)):
            Cirow.append(self._D_s_Order[self.lmatrix_col[x]].replace(
                'self.', ''))
        self.mca_ci_row = Cirow
        self.mca_ci_col = Cicol
        del e2, EL, KEL
        if self.__settings__['display_debug'] == 1:
            print('print self.mca_ci_row')
            print(self.mca_ci_row)
            print('print self.mca_ci_col')
            print(self.mca_ci_col)
            print('print self.mca_ci')
            print(self.mca_ci)
        CJi = numpy.zeros((len(self.kmatrix_col), self.mca_ci.shape[1]), 'd')
        CSi = numpy.zeros((self.mca_ci.shape[0] - len(self.kmatrix_col),
            self.mca_ci.shape[1]), 'd')
        for x in range(0, self.mca_ci.shape[0]):
            for y in range(0, self.mca_ci.shape[1]):
                if x < len(self.kmatrix_col):
                    CJi[x, y] = self.mca_ci[x, y]
                else:
                    CSi[x - len(self.kmatrix_col), y] = self.mca_ci[x, y]
        Ko_row = []
        Ko_col = []
        sKo = numpy.zeros(self.__kzeromatrix__.shape, 'd')
        xFactor = self.__kmatrix__.shape[0] - self.__kzeromatrix__.shape[0]
        for x in range(self.__kzeromatrix__.shape[0] - 1, -1, -1):
            Ko_row.append(self.__reactions__[self.kzeromatrix_row[x]])
            for y in range(self.__kzeromatrix__.shape[1] - 1, -1, -1):
                sKo[x, y] = self.kmatrix_scaled[x + xFactor, y]
                if x == 0:
                    Ko_col.append(self.__reactions__[self.kzeromatrix_col[y]])
        Ko_row.reverse()
        Ko_col.reverse()
        if self.__settings__['display_debug'] == 1:
            print('CJi')
            print(CJi)
            print('sKo')
            print(sKo)
            print('self.kzeromatrix')
            print(self.__kzeromatrix__)
        if self.__settings__['mode_mca_scaled']:
            self.mca_cjd = numpy.dot(sKo, CJi)
        else:
            self.mca_cjd = numpy.dot(self.__kzeromatrix__, CJi)
        self.mca_cjd_row = Ko_row
        self.mca_cjd_col = copy.copy(self.mca_ci_col)
        del sKo, CJi
        if self.__settings__['display_debug'] == 1:
            print('self.mca_cjd_row')
            print(self.mca_cjd_row)
            print('self.mca_cjd_col')
            print(self.mca_cjd_col)
            print('self.mca_cjd')
            print(self.mca_cjd)
        Lo_row = []
        Lo_col = []
        sLo = numpy.zeros(self.__lzeromatrix__.shape, 'd')
        xFactor = self.__lmatrix__.shape[0] - self.__lzeromatrix__.shape[0]
        for x in range(self.__lzeromatrix__.shape[0] - 1, -1, -1):
            Lo_row.append(self._D_s_Order[self.lzeromatrix_row[x]].replace(
                'self.', ''))
            for y in range(self.__lzeromatrix__.shape[1] - 1, -1, -1):
                sLo[x, y] = self.lmatrix_scaled[x + xFactor, y]
                if x == 0:
                    Lo_col.append(self._D_s_Order[self.lzeromatrix_col[y]])
        Lo_row.reverse()
        Lo_col.reverse()
        if self.__HAS_MOIETY_CONSERVATION__ == True:
            if self.__settings__['mode_mca_scaled']:
                self.mca_csd = numpy.dot(sLo, CSi)
            else:
                self.mca_csd = numpy.dot(self.__lzeromatrix__, CSi)
            self.mca_csd_row = Lo_row
            self.mca_csd_col = copy.copy(self.mca_ci_col)
        else:
            self.mca_csd = None
            self.mca_csd_row = None
            self.mca_csd_col = None
        del CSi, sLo
        if self.__settings__['display_debug'] == 1:
            print('self.mca_csd')
            print(self.mca_csd)
            print('self.mca_csd_row')
            print(self.mca_csd_row)
            print('self.mca_csd_col')
            print(self.mca_csd_col)
        if self.__HAS_FLUX_CONSERVATION__ and self.__HAS_MOIETY_CONSERVATION__:
            self.cc_flux = numpy.concatenate((self.mca_ci[:self.__kmatrix__
                .shape[1], :], self.mca_cjd))
            self.cc_flux_row = copy.copy(self.mca_ci_row[:self.__kmatrix__.
                shape[1]] + self.mca_cjd_row)
            self.cc_flux_col = copy.copy(self.mca_ci_col)
            self.cc_conc = numpy.concatenate((self.mca_ci[self.__kmatrix__.
                shape[1]:, :], self.mca_csd))
            self.cc_conc_row = copy.copy(self.mca_ci_row[self.__kmatrix__.
                shape[1]:] + self.mca_csd_row)
            self.cc_conc_col = copy.copy(self.mca_ci_col)
        elif self.__HAS_FLUX_CONSERVATION__:
            self.cc_flux = numpy.concatenate((self.mca_ci[:self.__kmatrix__
                .shape[1], :], self.mca_cjd))
            self.cc_flux_row = copy.copy(self.mca_ci_row[:self.__kmatrix__.
                shape[1]] + self.mca_cjd_row)
            self.cc_flux_col = copy.copy(self.mca_ci_col)
            self.cc_conc = copy.copy(self.mca_ci[self.__kmatrix__.shape[1]:, :]
                )
            self.cc_conc_row = copy.copy(self.mca_ci_row[self.__kmatrix__.
                shape[1]:])
            self.cc_conc_col = copy.copy(self.mca_ci_col)
        elif self.__HAS_MOIETY_CONSERVATION__:
            print(
                "INFO: this is interesting no dependent flux cc's only dependent conc cc's!"
                )
            self.cc_flux = copy.copy(self.mca_ci[:self.__kmatrix__.shape[1], :]
                )
            self.cc_flux_row = copy.copy(self.mca_ci_row[:self.__kmatrix__.
                shape[1]])
            self.cc_flux_col = copy.copy(self.mca_ci_col)
            self.cc_conc = numpy.concatenate((self.mca_ci[self.__kmatrix__.
                shape[1]:, :], self.mca_csd))
            self.cc_conc_row = copy.copy(self.mca_ci_row[self.__kmatrix__.
                shape[1]:] + self.mca_csd_row)
            self.cc_conc_col = copy.copy(self.mca_ci_col)
        else:
            print(
                'INFO: this is interesting no dependent flux/conc coefficients!'
                )
            self.cc_flux = copy.copy(self.mca_ci[:self.__kmatrix__.shape[1], :]
                )
            self.cc_flux_row = copy.copy(self.mca_ci_row[:self.__kmatrix__.
                shape[1]])
            self.cc_flux_col = copy.copy(self.mca_ci_col)
            self.cc_conc = copy.copy(self.mca_ci[self.__kmatrix__.shape[1]:, :]
                )
            self.cc_conc_row = copy.copy(self.mca_ci_row[self.__kmatrix__.
                shape[1]:])
            self.cc_conc_col = copy.copy(self.mca_ci_col)
        self.cc_all = numpy.concatenate((self.cc_flux, self.cc_conc))
        self.cc_all_row = self.cc_flux_row + self.cc_conc_row
        self.cc_all_col = copy.copy(self.mca_ci_col)
        self.cc = BagOfStuff(self.cc_all, self.cc_all_row, self.cc_all_col)
        self.cc.load()
        if self.__settings__['mode_mca_scaled'] == 1:
            self.cc.scaled = True
        else:
            self.cc.scaled = False
        if self.__settings__['display_debug'] == 1:
            print('self.cc_flux')
            print(self.cc_flux_row)
            print(self.cc_flux_col)
            print(self.cc_flux.shape)
            print('self.cc_conc')
            print(self.cc_conc_row)
            print(self.cc_conc_col)
            print(self.cc_conc.shape)
            print('self.cc_all')
            print(self.cc_all_row)
            print(self.cc_all_col)
            print(self.cc_all.shape)
        if self.__settings__['mca_ccj_upsymb']:
            r, c = self.cc_all.shape
            CJoutput = ''
            for x in range(0, len(self.__reactions__)):
                for y in range(0, c):
                    if self.__settings__['mode_mca_scaled']:
                        setattr(self, 'ccJ' + self.cc_all_row[x] + '_' +
                            self.cc_all_col[y], self.cc_all[x, y])
                    else:
                        setattr(self, 'uccJ' + self.cc_all_row[x] + '_' +
                            self.cc_all_col[y], self.cc_all[x, y])
        if self.__settings__['mca_ccs_upsymb']:
            CSoutput = ''
            for x in range(len(self.__reactions__), r):
                for y in range(0, c):
                    if self.__settings__['mode_mca_scaled']:
                        setattr(self, 'cc' + self.cc_all_row[x] + '_' +
                            self.cc_all_col[y], self.cc_all[x, y])
                    else:
                        setattr(self, 'ucc' + self.cc_all_row[x] + '_' +
                            self.cc_all_col[y], self.cc_all[x, y])
        if self.__settings__['display_debug'] == 1:
            print('CJoutput')
            print(CJoutput)
            print('CSoutput')
            print(CSoutput)

    def showCC(self, File=None):
        """
    Print all control coefficients as 'LaTeX' formatted strings to the screen or file.

    This method outputs the control coefficients of the system in a LaTeX format.
    The output can either be displayed on the screen or written to a specified
    file. The coefficients are optionally grouped by reaction or separated into
    flux and concentration coefficients based on the model's settings.

    Parameters
    ----------
    File : file-like object, optional
        An open, writable Python file object where the LaTeX formatted control
        coefficients will be written. If not provided, the output will be printed
        to the screen.

    Examples
    --------
    To print the control coefficients to the screen:

    >>> model_instance = PysMod()
    >>> model_instance.showCC()

    To write the control coefficients to a file:

    >>> with open('control_coefficients.tex', 'w') as f:
    >>>     model_instance.showCC(File=f)
    """
        r, c = self.cc_all.shape
        if self.__settings__['mca_ccall_altout']:
            CAltoutput = ''
            for x in range(0, c):
                col = self.cc_all_col[x]
                CAltoutput += '\n`Reaction ' + col + "'\n"
                for y in range(0, r):
                    if y == 0:
                        CAltoutput += '"Flux control coefficients"\n'
                    if y == len(self.__reactions__):
                        CAltoutput += '"Concentration control coefficients"\n'
                    row = self.cc_all_row[y]
                    rtemp = self.__settings__['mode_number_format'
                        ] % self.cc_all[y, x]
                    if y < len(self.__reactions__):
                        if self.__settings__['mode_mca_scaled']:
                            cc = '\\cc{J' + row + '}{' + col + '} = ' + rtemp
                        else:
                            cc = '\\ucc{J' + row + '}{' + col + '} = ' + rtemp
                    elif self.__settings__['mode_mca_scaled']:
                        cc = '\\cc{' + row + '}{' + col + '} = ' + rtemp
                    else:
                        cc = '\\ucc{' + row + '}{' + col + '} = ' + rtemp
                    CAltoutput += cc + '\n'
        else:
            if self.__settings__['mca_ccall_fluxout']:
                CJoutput = ''
                for x in range(0, len(self.__reactions__)):
                    row = self.cc_all_row[x]
                    CJoutput += '\n`J' + row + "'\n"
                    for y in range(0, c):
                        col = self.cc_all_col[y]
                        rtemp = self.__settings__['mode_number_format'
                            ] % self.cc_all[x, y]
                        if self.__settings__['mode_mca_scaled']:
                            cc = '\\cc{J' + row + '}{' + col + '} = ' + rtemp
                        else:
                            cc = '\\ucc{J' + row + '}{' + col + '} = ' + rtemp
                        CJoutput += cc + '\n'
            if self.__settings__['mca_ccall_concout']:
                CSoutput = ''
                for x in range(len(self.__reactions__), r):
                    row = self.cc_all_row[x]
                    CSoutput += '\n`' + row + "'\n"
                    for y in range(0, c):
                        col = self.cc_all_col[y]
                        rtemp = self.__settings__['mode_number_format'
                            ] % self.cc_all[x, y]
                        if self.__settings__['mode_mca_scaled']:
                            cc = '\\cc{' + row + '}{' + col + '} = ' + rtemp
                        else:
                            cc = '\\ucc{' + row + '}{' + col + '} = ' + rtemp
                        CSoutput += cc + '\n'
        if File != None:
            if self.__settings__['mca_ccall_altout']:
                print('\nControl coefficients grouped by reaction')
                File.write('\n## Control coefficients grouped by reaction\n')
                File.write(CAltoutput)
            else:
                if self.__settings__['mca_ccall_fluxout']:
                    print('\nFlux control coefficients')
                    File.write('\n## Flux control coefficients\n')
                    File.write(CJoutput)
                if self.__settings__['mca_ccall_concout']:
                    print('\nConcentration control coefficients')
                    File.write('\n## Concentration control coefficients\n')
                    File.write(CSoutput)
        elif self.__settings__['mca_ccall_altout']:
            print('\nControl coefficients grouped by reaction')
            print(CAltoutput)
        else:
            if self.__settings__['mca_ccall_fluxout']:
                print('\nFlux control coefficients')
                print(CJoutput)
            if self.__settings__['mca_ccall_concout']:
                print('\nConcentration control coefficients')
                print(CSoutput)

    def EvalRC(self):
        """
    EvalRC()

    Calculate the MCA response coefficients using the current steady-state solution.

    This method computes the response coefficients in Metabolic Control Analysis (MCA) using previously
    calculated control coefficients and elasticity matrices. These response coefficients provide insight
    into how changes in parameters affect the steady-state flux and concentration levels.

    Returns
    -------
    None

    Side Effects
    ------------
    - Modifies `self.mca_rc_par`, storing the calculated response coefficients.
    - Sets `self.mca_rc_par_row` and `self.mca_rc_par_col` for row and column labels, respectively.
    - Updates `self.rc`, a `BagOfStuff` instance containing the response coefficients matrix.
    - Adjusts the `scaled` attribute of `self.rc` based on current settings.

    Examples
    --------
    >>> pys_mod_instance = PysMod()
    >>> pys_mod_instance.EvalRC()
    >>> print(pys_mod_instance.mca_rc_par)
    [[...], [...], ...]  # Example display of calculated matrix of response coefficients.

    Notes
    -----
    - Ensure that the `cc_all` matrix and elasticity matrices have been appropriately initialized
      before calling this method.
    - The method relies on several class attributes related to elasticity and control coefficients
      matrices.
    """
        reordered_cc_all = numpy.zeros(self.cc_all.shape, 'd')
        for reac in range(len(self.elas_par_row)):
            reordered_cc_all[:, reac] = self.cc_all[:, list(self.cc_all_col
                ).index(self.elas_par_row[reac])]
        self.mca_rc_par = numpy.dot(reordered_cc_all, self.elas_par)
        self.mca_rc_par_row = self.cc_all_row
        self.mca_rc_par_col = self.elas_par_col
        del reordered_cc_all
        self.__structural__.MatrixFloatFix(self.mca_rc_par)
        self.rc = BagOfStuff(self.mca_rc_par, self.mca_rc_par_row, self.
            mca_rc_par_col)
        self.rc.load()
        if self.__settings__['mode_mca_scaled'] == 1:
            self.rc.scaled = True
        else:
            self.rc.scaled = False

    def EvalRCT(self):
        """
    Calculate the MCA response coefficients.

    This method calculates the Metabolic Control Analysis (MCA) response coefficients
    using the current steady-state solution. It also calculates the responses to changes
    in the sums of moiety conserved cycles.

    Parameters
    ----------
    None

    Attributes
    ----------
    mca_rct_conc : numpy.ndarray
        The calculated MCA response coefficients for concentrations.
    mca_rct_flux : numpy.ndarray
        The calculated MCA response coefficients for fluxes.
    mca_rct_conc_row : list
        Row labels for the concentration MCA response coefficients.
    mca_rct_conc_col : list
        Column labels for the concentration MCA response coefficients, corresponding to moiety sums.
    mca_rct_flux_row : list
        Row labels for the flux MCA response coefficients.
    mca_rct_flux_col : list
        Column labels for the flux MCA response coefficients, corresponding to moiety sums.
    rc : BagOfStuff
        An object containing the full MCA response coefficients matrix, with rows and columns
        labeled appropriately.

    Notes
    -----
    The calculations in this method rely heavily on the steady-state condition of the system,
    as well as the specified scaling settings in `__settings__`.

    Examples
    --------
    Calculate the MCA response coefficients and print the results:

    >>> pysmod_instance = PysMod()
    >>> pysmod_instance.EvalRCT()
    >>> print(pysmod_instance.mca_rct_conc)
    >>> print(pysmod_instance.mca_rct_flux)
    >>> print(pysmod_instance.rc)
    """
        reaction_reordered_cc_conc = numpy.zeros(self.cc_conc.shape, 'd')
        for reac in range(len(self.elas_var_row)):
            reaction_reordered_cc_conc[:, reac] = self.cc_conc[:, list(self
                .cc_conc_col).index(self.elas_var_row[reac])]
        reaction_reordered_cc_flux = numpy.zeros(self.cc_flux.shape, 'd')
        for reac in range(len(self.elas_var_row)):
            reaction_reordered_cc_flux[:, reac] = self.cc_flux[:, list(self
                .cc_flux_col).index(self.elas_var_row[reac])]
        reordered_cc_conc = numpy.zeros(self.cc_conc.shape, 'd')
        for metab in range(len(self.elas_var_col)):
            reordered_cc_conc[metab, :] = reaction_reordered_cc_conc[list(
                self.cc_conc_row).index(self.elas_var_col[metab]), :]
        reordered_cc_flux = numpy.zeros(self.cc_flux.shape, 'd')
        for flux in range(len(self.elas_var_row)):
            reordered_cc_flux[flux, :] = reaction_reordered_cc_flux[list(
                self.cc_flux_row).index(self.elas_var_row[flux]), :]
        m = len(self.species)
        r = len(self.Lmatrix.col)
        Im = numpy.matrix(numpy.eye(m))
        reordered_zero_I = numpy.zeros((m, m - r), 'd')
        for species in self.elas_var_col:
            if species in self.Consmatrix.row:
                row = self.elas_var_col.index(species)
                col = self.Consmatrix.row.index(species)
                reordered_zero_I[row, col] = 1.0
        if self.__settings__['mode_mca_scaled'] == 1:
            diag_T = None
            if self.__KeyWords__['Species_In_Conc']:
                diag_T = numpy.matrix(numpy.diag(self.__tvec_c__))
            else:
                diag_T = numpy.matrix(numpy.diag(self.__tvec_a__))
            diag_S = numpy.matrix(numpy.diag(self.data_sstate.getStateData(
                *self.elas_var_col)))
            self.mca_rct_conc = numpy.dot(numpy.dot(numpy.dot(numpy.dot(
                reordered_cc_conc, self.elas_var) + Im, numpy.linalg.inv(
                diag_S)), reordered_zero_I), diag_T)
            self.mca_rct_flux = numpy.dot(numpy.dot(numpy.dot(numpy.dot(
                reordered_cc_flux, self.elas_var), numpy.linalg.inv(diag_S)
                ), reordered_zero_I), diag_T)
        else:
            self.mca_rct_conc = numpy.dot(numpy.dot(reordered_cc_conc, self
                .elas_var_u) + Im, reordered_zero_I)
            self.mca_rct_flux = numpy.dot(numpy.dot(reordered_cc_flux, self
                .elas_var_u), reordered_zero_I)
        del reordered_cc_conc
        del reordered_cc_flux
        moiety_sum_names = ['T_{0}'.format(dependent_species) for
            dependent_species in self.Consmatrix.row]
        self.__structural__.MatrixFloatFix(self.mca_rct_conc)
        self.mca_rct_conc_row = self.elas_var_col
        self.mca_rct_conc_col = moiety_sum_names
        self.__structural__.MatrixFloatFix(self.mca_rct_flux)
        self.mca_rct_flux_row = self.elas_var_row
        self.mca_rct_flux_col = moiety_sum_names
        vstacked_col = list(self.elas_var_row) + list(self.elas_var_col)
        vstacked = numpy.vstack((self.mca_rct_flux, self.mca_rct_conc))
        reordered_vstack = numpy.zeros(vstacked.shape, 'd')
        for row_index in range(len(self.rc.row)):
            reordered_vstack[row_index, :] = vstacked[list(vstacked_col).
                index(self.rc.row[row_index]), :]
        hstacked = numpy.hstack((self.rc.matrix, vstacked))
        self.rc = BagOfStuff(hstacked, self.rc.row, self.rc.col +
            moiety_sum_names)
        self.rc.load()
        if self.__settings__['mode_mca_scaled'] == 1:
            self.rc.scaled = True
        else:
            self.rc.scaled = False

    def EvalEigen(self):
        """
    Calculate the eigenvalues or vectors of the unscaled Jacobian matrix and 
    thereby analyze the stability of a system.

    This method computes the eigenvalues and eigenvectors of the Jacobian matrix 
    based on the current settings. The Jacobian matrix is constructed from the 
    stoichiometry and elasticity matrices. The method is used to investigate 
    system stability and to understand the system's dynamic behavior near 
    equilibrium points.

    Parameters
    ----------
    None

    Attributes
    ----------
    jacobian : tuple
        Tuple representation of the Jacobian matrix.
    jacobian_row : tuple
        Row names of the Jacobian matrix indicating system variables.
    jacobian_col : tuple
        Column names of the Jacobian matrix indicating system variables.
    eigen_values : ndarray
        Array of eigenvalues calculated from the Jacobian matrix.
    eigen_vecleft : ndarray, optional
        Left eigenvectors, available if `mode_eigen_output` is True.
    eigen_vecright : ndarray, optional
        Right eigenvectors, available if `mode_eigen_output` is True.

    Notes
    -----
    Eigenvalues are crucial in determining the stability. If any eigenvalue has 
    a positive real part, the system is unstable at the steady state.

    Examples
    --------
    Calculate the eigenvalues and eigenvectors of the system:

    >>> pysmod_instance = PysMod()
    >>> pysmod_instance.EvalEigen()
    >>> print("Eigenvalues: ", pysmod_instance.eigen_values)
    >>> if pysmod_instance.__settings__['mode_eigen_output']:
    >>>     print("Left Eigenvectors: ", pysmod_instance.eigen_vecleft)
    >>>     print("Right Eigenvectors: ", pysmod_instance.eigen_vecright)

    Displaying debug information if enabled:

    >>> pysmod_instance.__settings__['display_debug'] = 1
    >>> pysmod_instance.EvalEigen()
    """
        Nr = copy.copy(self.__nrmatrix__)
        Es = numpy.zeros(self.elas_var_u.shape, 'd')
        for x in range(0, len(self.nmatrix_col)):
            for y in range(0, len(self.lmatrix_row)):
                Es[x, y] = self.elas_var_u[self.nmatrix_col[x], self.
                    lmatrix_row[y]]
        if self.__HAS_MOIETY_CONSERVATION__ == True:
            jac1 = numpy.dot(Nr, Es)
            jacobian = numpy.dot(jac1, copy.copy(self.__lmatrix__))
        else:
            jacobian = numpy.dot(Nr, Es)
        Si = []
        for x in range(0, len(self.lmatrix_col)):
            Si.append(self.__species__[self.lmatrix_col[x]])
            self.jacobian = tuple(jacobian)
            self.jacobian_row = tuple(Si)
            self.jacobian_col = tuple(Si)
            if self.__settings__['mode_eigen_output']:
                eigenval, self.eigen_vecleft, self.eigen_vecright = (scipy.
                    linalg.eig(jacobian, numpy.identity(jacobian.shape[0],
                    'd'), left=1, right=1))
                if self.__settings__['display_debug'] == 1:
                    print('\nEigenvalues')
                    print(eigenval)
                    print('\nLeft Eigenvector')
                    print(self.eigen_vecleft)
                    print('\nRight Eigenvector')
                    print(self.eigen_vecright)
            else:
                eigenval = scipy.linalg.eigvals(jacobian)
                if self.__settings__['display_debug'] == 1:
                    print('\nEigenvalues')
                    print(eigenval)
        self.eigen_values = eigenval
        for x in range(len(self.eigen_values)):
            setattr(self, 'lambda' + str(x + 1), self.eigen_values[x])
        if self.__settings__['display_debug'] == 1:
            print('\nEigenvalues attached as lambda1 ... lambda' + repr(len
                (eigenorder)) + '\n')

    def showEigen(self, File=None):
        """
    Print the eigenvalues and stability analysis of a system generated with EvalEigen()
    to the screen or file.

    Parameters
    ----------
    File : file-like object, optional
        An open, writable Python file object to which the output is written.
        If None, the output is printed to the screen.

    Examples
    --------
    Analyze and display the eigenvalues and stability of the system:

    >>> model = PysMod()
    >>> model.EvalEigen()  # Calculates the eigenvalues and updates the system
    >>> model.showEigen()  # Prints the eigenvalues and stability to the screen

    Save the eigenvalues and stability analysis to a file:

    >>> with open('eigen_analysis.txt', 'w') as file:
    ...     model.showEigen(file)  # Writes the output to 'eigen_analysis.txt'
    
    Notes
    -----
    This method depends on the prior execution of `EvalEigen()` to ensure
    that eigenvalues are available for analysis.
    """
        eigenstats = ''
        eigenstats += 'Max real part: ' + self.__settings__[
            'mode_number_format'] % max(self.eigen_values.real) + '\n'
        eigenstats += 'Min real part: ' + self.__settings__[
            'mode_number_format'] % min(self.eigen_values.real) + '\n'
        eigenstats += 'Max absolute imaginary part: ' + self.__settings__[
            'mode_number_format'] % max(abs(self.eigen_values.imag)) + '\n'
        eigenstats += 'Min absolute imaginary part: ' + self.__settings__[
            'mode_number_format'] % min(abs(self.eigen_values.imag)) + '\n'
        eigenstats += 'Stiffness: ' + self.__settings__['mode_number_format'
            ] % (max(abs(self.eigen_values.real)) / min(abs(self.
            eigen_values.real))) + '\n'
        pure_real = 0
        pure_imag = 0
        negcplx = 0
        pure_cplx = 0
        pure_zero = 0
        pos_real = 0
        neg_real = 0
        for x in self.eigen_values:
            if x.real != 0.0 and x.imag == 0.0:
                pure_real += 1
            if x.real == 0.0 and x.imag != 0.0:
                pure_imag += 1
            if x.real == 0.0 and x.imag == 0.0:
                pure_zero += 1
            if x.real > 0.0:
                pos_real += 1
            if x.real < 0.0:
                neg_real += 1
            if x.real < 0.0 and x.imag != 0.0:
                negcplx += 1
        eigenstats += '\n## Stability'
        eiglen = len(self.eigen_values)
        if neg_real == eiglen:
            eigenstats += ' --> Stable state'
        elif pure_zero > 0:
            eigenstats += ' --> Undetermined'
        elif pos_real == eiglen:
            eigenstats += ' --> Unstable state'
        eigenstats += '\nPurely real: ' + repr(pure_real)
        eigenstats += '\nPurely imaginary: ' + repr(pure_imag)
        eigenstats += '\nZero: ' + repr(pure_zero)
        eigenstats += '\nPositive real part: ' + repr(pos_real)
        eigenstats += '\nNegative real part: ' + repr(neg_real)
        if File != None:
            print('\nEigen values')
            File.write('\n## Eigen values\n')
            scipy.io.write_array(File, self.eigen_values, precision=2,
                keep_open=1)
            print('\nEigen statistics')
            File.write('\n## Eigen statistics\n')
            File.write(eigenstats)
            if self.__settings__['mode_eigen_output']:
                print('\nLeft eigen vector')
                File.write('\n## Left eigen vector\n')
                scipy.io.write_array(File, self.eigen_vecleft, precision=2,
                    keep_open=1)
                print('\nRight eigen vector')
                File.write('\n## Right eigen vector\n')
                scipy.io.write_array(File, self.eigen_vecright, precision=2,
                    keep_open=1)
        else:
            print('\nEigen values')
            print(self.eigen_values)
            print('\nEigen statistics')
            print(eigenstats)
            if self.__settings__['mode_eigen_output']:
                print('\nLeft eigen vector')
                print(self.eigen_vecleft)
                print('\nRight eigen vector')
                print(self.eigen_vecright)

    def doSim(self, end=10.0, points=21):
        """
    Run a time simulation from t=0 to t=sim_end with sim_points.

    This method sets up a time simulation by defining the simulation end time
    and the number of data points to be computed. After setting up, it 
    invokes the `Simulate` method to execute the simulation.

    Parameters
    ----------
    end : float, optional
        The simulation end time, default is 10.0.
    points : int, optional
        The number of points in the simulation, default is 21.

    Calls
    -----
    Simulate
        Executes the simulation using the specified end time and number of points.

    Example
    -------
    >>> model = PysMod()
    >>> model.doSim(end=5.0, points=100)
    This will run a simulation from t=0 to t=5.0 with 100 time points.
    """
        self.sim_end = end
        self.sim_points = points
        self.Simulate()

    def doSimPlot(self, end=10.0, points=21, plot='species', fmt='lines',
        filename=None):
        """
    Run a time simulation from t=0 to `end` with `points` and plot the results.
    
    Parameters
    ----------
    end : float, optional
        The end time of the simulation, default is 10.0.
    points : int, optional
        The number of points in the simulation, default is 21.
    plot : str, optional
        Select the output data to plot, default is 'species'. Options are:
        
        - 'species': Plot the species data.
        - 'rates': Plot the rates data.
        - 'all': Plot both species and rates data.
    fmt : str, optional
        The format of the plot, depends on the UPI backend, default is 'lines'.
        Options can be:
        
        - '': Default format.
        - 'lines': Use line style in plot.
        - 'points': Use points only in plot.
    filename : str or None, optional
        The name of the file to which the plot is exported as a PNG image.
        If `None`, the plot is not exported. Default is `None`.

    Calls
    -----
    Simulate()
        Starts the simulation process.
    SimPlot(plot, format, filename)
        Handles the plotting of simulation results according to specified parameters.

    Notes
    -----
    This method requires the simulation process to be complete before plotting.
    Ensure that the simulation settings are correctly configured according to your needs.

    Examples
    --------
    To run a simulation with the default settings and plot species:

    >>> p = PysMod()
    >>> p.doSimPlot()

    To run a simulation up to time 20.0 with 50 points, plot both species and rates, 
    and export the plot as 'simulation_output.png':

    >>> p.doSimPlot(end=20.0, points=50, plot='all', filename='simulation_output')
    """
        self.sim_end = end
        self.sim_points = points
        self.Simulate()
        self.SimPlot(plot=plot, format=fmt, filename=filename)

    def exportSimAsSedML(self, output='files', return_sed=False, vc_given=
        'PySCeS', vc_family='Software', vc_email='', vc_org=
        'pysces.sourceforge.net'):
        """
    Exports the current simulation as SED-ML in various formats, creating and storing the 
    SED-ML files in a folder derived from the model name.

    Parameters
    ----------
    output : str, optional
        The type of SED-ML export format. Options are one or more comma-separated values from: 
        'files', 'archive', 'combine'. Default is 'files'.
        - 'files': Exports the plain SBML and SEDML XML files.
        - 'archive': Exports as a SED-ML archive `<file>.sedx` containing SBML and SEDML XML files.
        - 'combine': Exports as a COMBINE archive `<file>.omex` containing the SBML, SEDML, 
          manifest (XML), and metadata (RDF).
          
    return_sed : bool, optional
        If True, returns the SED object after export. Default is False.

    vc_given : str, optional
        Version control given name. Default is 'PySCeS'.

    vc_family : str, optional
        Version control family name. Default is 'Software'.

    vc_email : str, optional
        Version control email associated with the software. Default is 'bgoli@users.sourceforge.net'.

    vc_org : str, optional
        Organization for the version control. Default is '<pysces.sourceforge.net>'.

    Returns
    -------
    SED.SED, optional
        Returns the SED object if `return_sed` is True; otherwise, returns None.

    Calls
    -----
    SED.SED.writeSedScript : Writes the SED script to a file when 'files' is in output.
    SED.SED.writeSedXML : Writes the SED XML to a file when 'files' is in output.
    SED.SED.writeSedXArchive : Writes the SED-ML as a `.sedx` archive when 'archive' is in output.
    SED.SED.writeCOMBINEArchive : Writes the SED-ML as a `.omex` COMBINE archive when 'combine' is in output. 

    Example
    -------
    >>> pysmod_instance = PysMod()
    >>> pysmod_instance.exportSimAsSedML(output='files,combine', return_sed=False, 
                                          vc_given='JohnDoe', vc_family='BioSoftware', 
                                          vc_email='johndoe@example.com', vc_org='example.org')

    Notes
    -----
    The SED-ML files are saved in a directory named after the model, ensuring that all 
    exports get a distinct and organized file structure.
    """
        sedname = self.ModelFile.replace('.psc', '')
        sedout = os.path.join(OUTPUT_DIR, 'sedout')
        if not os.path.exists(sedout):
            os.makedirs(sedout)
        S = SED.SED(sedname + '_sed', sedout)
        S.addModel(sedname, self)
        S.addSimulation('sim0', self.sim_start, self.sim_end, self.
            sim_points, self.__SIMPLOT_OUT__, initial=None, algorithm=
            'KISAO:0000019')
        S.addTask('task0', 'sim0', sedname)
        S.addTaskDataGenerators('task0')
        S.addTaskPlot('task0')
        output = [s_.strip() for s_ in output.split(',')]
        for s_ in output:
            if s_ == 'files':
                S.writeSedScript()
                S.writeSedXML()
            elif s_ == 'archive':
                S.writeSedXArchive()
            elif s_ == 'combine':
                S.writeCOMBINEArchive(vc_given=vc_given, vc_family=
                    vc_family, vc_email=vc_email, vc_org=vc_org)
        S._SED_CURRENT_ = False
        if return_sed:
            return S
        else:
            del S

    def doSimPerturb(self, pl, end):
        """
    Run a simulation with perturbations.

    .. deprecated:: 1.0
       This method is deprecated and will be removed in future versions.
       Please use the `events` method instead for handling perturbations
       during simulations.

    Parameters
    ----------
    pl : object
        The perturbation list defining the changes to be applied during 
        the simulation.
    end : float
        The end time of the simulation.

    Examples
    --------
    >>> model = PysMod()
    >>> perturbations = [('parameter1', 1.5), ('parameter2', 0.8)]
    >>> model.doSimPerturb(perturbations, 10.0)
    """
        pass

    def doState(self):
        """
    Calculate the steady-state solution of the system.

    This method computes the steady-state of the current model by invoking the
    `State()` method. The steady-state solution is essential for understanding
    the equilibrium behavior of the system's variables under constant conditions.

    Calls
    -----
    State
        Computes the steady-state solution for the model.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Examples
    --------
    >>> model = PysMod()
    >>> model.doState()
    The system's steady-state solution is calculated and stored within the model.
    """
        self.State()

    def doStateShow(self):
        """
    doStateShow()

    Calculate the steady-state solution of a system and show the results.

    This method computes the steady-state solution using the `State` method
    and then displays the results using the `showState` method. The method
    asserts that the steady state is valid by checking the `__StateOK__`
    attribute.

    Calls
    -----
    State()
        Calculates the steady-state solution of the system.

    showState()
        Displays the steady-state results.

    Arguments
    ---------
    None

    See Also
    --------
    doState : Calculates the steady-state without showing the results.

    Examples
    --------
    >>> pymod = PysMod()
    >>> pymod.doStateShow()
    INFO: All systems functional. Steady-state results displayed.

    Notes
    -----
    It's essential to ensure `__StateOK__` is set to `True` after the state
    calculation to confirm the results are valid before proceeding to display them. 
    If `__StateOK__` is `False`, the system will raise an `AssertionError`.
    
INFO: Invalid steady state: run mod.doState() and check for errors"""
        self.State()
        assert self.__StateOK__ == True, """

INFO: Invalid steady state: run mod.doState() and check for errors"""
        self.showState()

    def doElas(self):
        """
    Calculate the model elasticities and automatically computes the steady state.

    This method executes a series of calculations needed to determine the model 
    elasticities. It begins by calculating the steady-state solution of the system 
    and proceeds to evaluate both variable and parameter elasticities. The method assumes 
    that the state has been correctly computed, and checks for its validity before 
    proceeding further.

    Calls
    -----
    State()
        Computes the steady-state solution of the system.
    EvalEvar()
        Evaluates variable elasticities.
    EvalEpar()
        Evaluates parameter elasticities.

    Raises
    ------
    AssertionError
        If the steady-state solution is not valid, indicating that `doState()` 
        should be run and errors should be checked.

    Examples
    --------
    >>> mod = PysMod()
    >>> mod.doElas()
    INFO: Steady state and elasticities computed successfully.

    Notes
    -----
    Make sure that the object `mod` is an initialized instance of the `PysMod` class 
    with all necessary parameters configured before calling `doElas()`.
    
    
INFO: Invalid steady state: run mod.doState() and check for errors"""
        self.State()
        assert self.__StateOK__ == True, """

INFO: Invalid steady state: run mod.doState() and check for errors"""
        self.EvalEvar()
        self.EvalEpar()

    def doEigen(self):
        """
    Calculate the eigenvalues, and automatically perform steady state
    and elasticity analysis.

    This method calculates the eigenvalues of the system and ensures that
    a steady-state and elasticity analysis is performed prior to the eigenvalue
    calculation. It is automatically assumed that the steady state is valid.

    Calls
    -----
    State()
        Calculates the steady state of the system.
    EvalEvar()
        Evaluates elasticity variables.
    Evaleigen()
        Computes the eigenvalues of the system.

    Raises
    ------
    AssertionError
        If the steady-state validation fails.

    Example
    -------
    >>> mod = PysMod()
    >>> mod.doEigen()
    Calculating steady state...
    Steady state valid: True
    Evaluating elasticities...
    Elasticities evaluated successfully.
    Calculating eigenvalues...
    Eigenvalues calculation complete.

    Notes
    -----
    It is important to ensure that the system is properly initialized before
    calling this method. The steady state result must be valid for the eigenvalues
    calculation to proceed without errors.
    
    
INFO: Invalid steady state: run mod.doState() and check for errors"""
        self.State()
        assert self.__StateOK__ == True, """

INFO: Invalid steady state: run mod.doState() and check for errors"""
        self.__settings__['elas_evar_upsymb'] = 1
        self.EvalEvar()
        self.__settings__['elas_evar_upsymb'] = 1
        self.EvalEigen()

    def doEigenShow(self):
        """
    doEigenShow()

    Calculate the eigenvalues, automatically performs a steady state and elasticity analysis,
    and displays the results.

    This method calls the `doEigen()` function to perform the calculations necessary for determining 
    the eigenvalues, including ensuring a steady state and conducting an elasticity analysis. 
    It then calls the `showEigen()` function to display the computed eigenvalues and analysis results.

    Calls
    -----
    doEigen() : None
        Performs the calculation of eigenvalues and the necessary preliminary analyses.
    showEigen() : None
        Displays the results of the eigenvalue calculations and analysis.

    Arguments
    ---------
    None

    Example
    -------
    >>> model = PysMod()
    >>> model.doEigenShow()
    INFO: Steady state calculated successfully.
    Eigenvalues and elasticity analysis results
    -------------------------------------------
    # Output of the showEigen() method

    Notes
    -----
    - Ensure that the model's state is correctly initialized before invoking `doEigenShow()` 
      by running appropriate state validation methods if necessary.
    - The `doEigen()` method assumes a valid steady state; any assertion error indicates 
      a need to manually run `doState()` and resolve potential issues.
    """
        self.doEigen()
        self.showEigen()

    def doEigenMca(self):
        """
    Calculate a full control analysis and eigenvalues, automatically performing 
    steady state, elasticity, and control analysis.

    This method asserts that the system is in a valid steady state before proceeding 
    to evaluate elasticity, control coefficients, and eigenvalues. It is a comprehensive 
    analysis method in the `PysMod` class.

    Calls
    -----
    State : method
        Checks and sets the steady state of the system.
    EvalEvar : method
        Evaluates the system's elasticity variables.
    EvalCC : method
        Evaluates the control coefficients.
    Evaleigen : method
        Calculates the eigenvalues of the system for analysis.

    Raises
    ------
    AssertionError
        If the system is not in a valid steady state, raises an assertion error 
        with information to run `doState()` for error checking.

    Examples
    --------
    >>> mod = PysMod()
    >>> mod.doEigenMca()
    INFO: Process successfully completed with valid steady state, elasticity, and control analysis.

    Notes
    -----
    Ensure that `doState()` has been successfully executed and that the system is 
    in a valid state before invoking this method to avoid AssertionErrors.

    

INFO: Invalid steady state: run mod.doState() and check for errors"""
        self.State()
        assert self.__StateOK__ == True, """

INFO: Invalid steady state: run mod.doState() and check for errors"""
        self.EvalEvar()
        self.EvalCC()
        self.EvalEigen()

    def doMca(self):
        """
    Perform a complete Metabolic Control Analysis on the model, automatically calculating a steady state.

    This method conducts a comprehensive Metabolic Control Analysis (MCA), which includes computations 
    for a steady state, elasticity coefficients, and control coefficients. The method automates the 
    process of ensuring that the model is in a valid steady state before proceeding with further 
    analyses.

    Calls
    -----
    State()
        Computes the steady state of the model.
    EvalEvar()
        Evaluates the elasticity variables.
    EvalEpar()
        Evaluates the elasticity parameters.
    EvalCC()
        Evaluates the control coefficients.

    Raises
    ------
    AssertionError
        If the steady state is not valid, prompting the user to run `doState()`.

    Examples
    --------
    >>> mod = PysMod()
    >>> mod.doState()  # Ensure the model is set up correctly, may require additional setup.
    >>> mod.doMca()    # Perform a full MCA on the model.
    

INFO: Invalid steady state run: mod.doState() and check for errors"""
        self.State()
        assert self.__StateOK__ == True, """

INFO: Invalid steady state run: mod.doState() and check for errors"""
        self.EvalEvar()
        self.EvalEpar()
        self.EvalCC()

    def doMcaRC(self):
        """
    Perform a comprehensive Metabolic Control Analysis (MCA) with Response Coefficients on the model, 
    automatically calculates a steady state.

    This method executes the necessary steps to perform a full MCA by computing the steady state, 
    elasticities, control coefficients, and response coefficients. It ensures that a valid steady state 
    is achieved before proceeding with the analysis. The method calls a series of internal analysis functions.

    Calls
    -----
    State() :
        Calculates the steady state of the model.
    EvalEvar() :
        Evaluates the elasticities of the model variables.
    EvalEpar() :
        Evaluates the elasticities of the model parameters.
    EvalCC() :
        Computes the control coefficients.
    EvalRC() :
        Computes the response coefficients.

    Returns
    -------
    None

    Raises
    ------
    AssertionError
        If the steady state is not valid (`self.__StateOK__` is False).

    Examples
    --------
    >>> mod = PysMod()
    >>> mod.doMcaRC()
    # Ensure that the model's steady state is valid and that all associated analyses are completed.

    Notes
    -----
    It is recommended to ensure that the model parameters are set appropriately before calling this method.
    
    
INFO: Invalid steady state: run mod.doState() and check for errors"""
        self.State()
        assert self.__StateOK__ == True, """

INFO: Invalid steady state run: mod.doState() and check for errors"""
        self.EvalEvar()
        self.EvalEpar()
        self.EvalCC()
        self.EvalRC()

    def doMcaRCT(self):
        """
    Perform a complete Metabolic Control Analysis on the model, automatically calculating a steady state.
    
    Additionally, response coefficients to the sums of moiety-conserved cycles are calculated.
    
    This method calls several other methods to evaluate various components of the metabolic model. It first checks 
    if the steady state is valid and then proceeds to evaluate elasticities, control coefficients, and response 
    coefficients related to moiety conservation if applicable.
    
    Calls
    -----
    State()
        Initialize and validate the steady state of the model.
    EvalEvar()
        Evaluate elasticities of variables.
    EvalEpar()
        Evaluate elasticities of parameters.
    EvalCC()
        Evaluate control coefficients.
    EvalRC()
        Evaluate response coefficients.
    EvalRCT()
        Evaluate response coefficients to moiety-conserved cycles (if applicable).
    
    Raises
    ------
    AssertionError
        If the steady state is invalid, an assertion error is raised with a message to check the model state.
    
    Notes
    -----
    If moiety conservation is not detected in the model, the method will revert to a simpler analysis using 
    `doMcaRC()` instead.

    Examples
    --------
    >>> model = PysMod()
    >>> model.doMcaRCT()
    Perform a complete Metabolic Control Analysis with additional moiety-conserved cycles response coefficients.
    
    This execution checks and uses the internal state of `model` to evaluate all components.
    
    
INFO: Invalid steady state run: mod.doState() and check for errors"""
        self.State()
        assert self.__StateOK__ == True, """

INFO: Invalid steady state run: mod.doState() and check for errors"""
        self.EvalEvar()
        self.EvalEpar()
        self.EvalCC()
        self.EvalRC()
        if self.__HAS_MOIETY_CONSERVATION__:
            self.EvalRCT()
        else:
            print(
                'No moiety conservation detected, reverting to simple doMcaRC().'
                )

    def showSpecies(self, File=None):
        """
    Prints the current value of the model's variable species (`mod.X`) to the screen or a specified file.

    Parameters
    ----------
    File : file-like object, optional
        An open, writable Python file object. If provided, the species values
        are written to this file. If `None` (default), the values are printed
        to the screen.

    Example
    -------
    To display species values on the screen:

    >>> model = PysMod()
    >>> model.showSpecies()

    To write species values to a file:

    >>> with open('species_values.txt', 'w') as f:
    ...     model.showSpecies(File=f)

    Notes
    -----
    The function iterates over the list of species in `self.__species__`.
    It formats these values according to the 'mode_number_format' specified
    in `self.__settings__`.
    """
        out_list = []
        print('\nSpecies values')
        out_list.append('\n## species values\n')
        for x in range(len(self.__species__)):
            if File is None:
                print(self.__species__[x] + ' = ' + self.__settings__[
                    'mode_number_format'] % getattr(self, self.__species__[x]))
            else:
                out_list.append(self.__species__[x] + ' = ' + self.
                    __settings__['mode_number_format'] % getattr(self, self
                    .__species__[x]) + '\n')
        if File != None:
            for x in out_list:
                File.write(x)

    def showSpeciesI(self, File=None):
        """
    showSpeciesI(File=None)

    Prints the current value of the model's variable species initial values (`mod.X_init`) to screen or to a file.

    Parameters
    ----------
    File : writable file object, optional
        An open, writable Python file object where the species initial values will be written.
        If `None`, the values are printed to the screen. Default is `None`.

    Examples
    --------
    >>> # Create an instance of the PysMod class
    >>> mod = PysMod()
    >>> # Display species initial values to the screen
    >>> mod.showSpeciesI()
    Species initial values
    species1_init = 1.0
    species2_init = 2.5
    ...

    >>> # Write species initial values to a file
    >>> with open('species_initial_values.txt', 'w') as file:
    ...     mod.showSpeciesI(File=file)

    Notes
    -----
    This method utilizes the class attribute `__species__`, which should contain the species names,
    and `__settings__`, which is expected to include a 'mode_number_format' for formatting numbers.
    """
        out_list = []
        print('\nSpecies initial values')
        out_list.append('\n## species initial values\n')
        for x in range(len(self.__species__)):
            if File is None:
                print(self.__species__[x] + '_init = ' + self.__settings__[
                    'mode_number_format'] % getattr(self, self.__species__[
                    x] + '_init'))
            else:
                out_list.append(self.__species__[x] + ' = ' + self.
                    __settings__['mode_number_format'] % getattr(self, self
                    .__species__[x] + '_init') + '\n')
        if File != None:
            for x in out_list:
                File.write(x)

    def showSpeciesFixed(self, File=None):
        """
    Prints the current value of the model's fixed species values to screen or file.

    This method outputs the current values of the fixed species in the model,
    either to the standard output or to a specified file. If no file is provided, 
    the values are printed to the console. Otherwise, they are written to the 
    provided file object.

    Parameters
    ----------
    File : file object, optional
        An open, writable Python file object where the fixed species values 
        will be written. If None, the values are printed to the standard output.

    Examples
    --------
    Printing to console:

    >>> model = PysMod()
    >>> model.showSpeciesFixed()

    Writing to a file:

    >>> with open('fixed_species.txt', 'w') as f:
    >>>     model.showSpeciesFixed(File=f)

    Notes
    -----
    The method relies on the `__fixed_species__` attribute, which should be a list 
    of species names, and on the `__settings__['mode_number_format']` string for 
    formatting the output values. Each species value is accessed as an attribute 
    of the instance, corresponding to the name in `__fixed_species__`.
    
    """
        out_list = []
        print('\nFixed species')
        out_list.append('\n## fixed species\n')
        for x in range(len(self.__fixed_species__)):
            if File is None:
                print(self.__fixed_species__[x] + ' = ' + self.__settings__
                    ['mode_number_format'] % getattr(self, self.
                    __fixed_species__[x]))
            else:
                out_list.append(self.__fixed_species__[x] + ' = ' + self.
                    __settings__['mode_number_format'] % getattr(self, self
                    .__fixed_species__[x]) + '\n')
        if File != None:
            for x in out_list:
                File.write(x)

    def showPar(self, File=None):
        """
    Prints the current values of the model's parameter values (`mod.P`) to the screen or a specified file.

    Parameters
    ----------
    File : file-like object, optional
        An open, writable file object where the parameter values will be written. If `None` (default),
        the output will be printed to the screen.

    Example
    -------
    >>> # Assuming `model` is an instance of PysMod
    >>> model.showPar()
    Parameters
    parameter1 = 12.345
    parameter2 = 67.890
    ...
    
    >>> with open('parameters.txt', 'w') as file:
    ...     model.showPar(File=file)
    """
        out_list = []
        print('\nParameters')
        out_list.append('\n## parameters\n')
        for x in range(len(self.__parameters__)):
            if File is None:
                print(self.__parameters__[x] + ' = ' + self.__settings__[
                    'mode_number_format'] % getattr(self, self.
                    __parameters__[x]))
            else:
                out_list.append(self.__parameters__[x] + ' = ' + self.
                    __settings__['mode_number_format'] % getattr(self, self
                    .__parameters__[x]) + '\n')
        if File != None:
            for x in out_list:
                File.write(x)

    def showModifiers(self, File=None):
        """
    Prints the current value of the model's modifiers per reaction to screen or file.

    Parameters
    ----------
    File : file-like object, optional
        An open, writable Python file object. If provided, the modifier information
        will be written to this file; otherwise, it will be printed to the screen.

    Notes
    -----
    This method accesses `self.__modifiers__`, which is expected to be a list
    of tuples. Each tuple represents a reaction and its corresponding modifiers.

    Example
    -------
    >>> model = PysMod()
    >>> model.showModifiers()
    Modifiers per reaction:
    Reaction1 has modifier: ModifierA 
    Reaction2 has modifiers: ModifierB ModifierC 
    Reactions with no modifiers:
    Reaction3 Reaction4 

    # To write to a file instead of printing to the screen:
    >>> with open('modifiers.txt', 'w') as file:
    >>>     model.showModifiers(file)

    """
        noMod = []
        out_list = []
        print('\nModifiers per reaction:')
        out_list.append('\n## modifiers per reaction\n')
        for reac in self.__modifiers__:
            rstr = ''
            if len(reac[1]) == 1:
                if File is None:
                    print(reac[0] + ' has modifier:', end=' ')
                rstr = rstr + reac[0] + ' has modifier: '
                for x in reac[1]:
                    if File is None:
                        print(x, end=' ')
                    rstr = rstr + x + ' '
                if File is None:
                    print(' ')
                rstr += '\n'
            elif len(reac[1]) > 1:
                if File is None:
                    print(reac[0] + ' has modifiers: ', end=' ')
                rstr = rstr + reac[0] + ' has modifiers: '
                for x in reac[1]:
                    if File is None:
                        print(x, end=' ')
                    rstr = rstr + x + ' '
                if File is None:
                    print(' ')
                rstr += '\n'
            else:
                noMod.append(reac[0])
            out_list.append(rstr)
        if len(noMod) > 0:
            print('\nReactions with no modifiers:')
            out_list.append('\n## reactions with no modifiers\n')
            cntr = 0
            rstr = ''
            for n in range(len(noMod)):
                cntr += 1
                if cntr > 6:
                    if File is None:
                        print(' ')
                    rstr += '\n'
                    cntr = 0
                if File is None:
                    print(noMod[n], end=' ')
                rstr = rstr + noMod[n] + ' '
            if File is None:
                print(' ')
            rstr += '\n'
            out_list.append(rstr)
            if File != None:
                for x in out_list:
                    File.write(x)

    def showState(self, File=None):
        """
    Displays the results of the last steady-state analysis, including both 
    steady-state species concentrations and fluxes. The output can be directed 
    to the console or a file.

    Parameters
    ----------
    File : file-like object, optional
        An open, writable Python file object where the results will be 
        written. If `None`, the results will be printed to the console.
        Default is `None`.

    Notes
    -----
    If a valid steady state is not found, a message indicating this will be 
    printed or written to the file.

    Examples
    --------
    >>> model = PysMod()
    >>> # Assuming model has been set up and steady-state analysis is performed
    >>> model.showState()
    Steady-state species concentrations
    A_ss = 0.123
    B_ss = 0.456
    ...
      
    Steady-state fluxes
    J_R1 = 0.789
    J_R2 = 0.012
    ...

    >>> with open('steady_state_results.txt', 'w') as f:
    ...     model.showState(File=f)

    """
        out_list = []
        print('\nSteady-state species concentrations')
        out_list.append('\n## Current steady-state species concentrations\n')
        if self.__StateOK__:
            for x in range(len(self.state_species)):
                if File is None:
                    print(self.__species__[x] + '_ss = ' + self.
                        __settings__['mode_number_format'] % self.
                        state_species[x])
                else:
                    out_list.append(self.__species__[x] + '_ss = ' + self.
                        __settings__['mode_number_format'] % self.
                        state_species[x] + '\n')
        else:
            print('No valid steady state found')
            out_list.append('No valid steady state found.\n')
        print('\nSteady-state fluxes')
        out_list.append('\n## Steady-state fluxes\n')
        if self.__StateOK__:
            for x in range(len(self.state_flux)):
                if File is None:
                    print('J_' + self.__reactions__[x] + ' = ' + self.
                        __settings__['mode_number_format'] % self.state_flux[x]
                        )
                else:
                    out_list.append('J_' + self.__reactions__[x] + ' = ' + 
                        self.__settings__['mode_number_format'] % self.
                        state_flux[x] + '\n')
        else:
            print('No valid steady state found')
            out_list.append('No valid steady state found.\n')
        if File != None:
            for x in out_list:
                File.write(x)

    def showRate(self, File=None):
        """
    Prints the current rates of all the reactions using the current parameter values and species concentrations.

    This method calculates and prints the reaction rates based on the current parameter values and species
    concentrations. If a file object is provided, the output is written to the file; otherwise, it is printed
    to the screen.

    Parameters
    ----------
    File : file-like object, optional
        An open, writable Python file object. If provided, the reaction rates are written to this file.
        Default is None, which results in printing the output to the console.

    Examples
    --------
    >>> my_model = PysMod()
    >>> my_model.showRate()
    Reaction rates:
    Reaction1 = 1.23e-3
    Reaction2 = 4.56e-4
    ...

    # Writing to a file
    >>> with open("reaction_rates.txt", "w") as f:
    >>>     my_model.showRate(File=f)
    
    Notes
    -----
    Ensure that the model parameters and species concentrations are set before calling this method
    to obtain meaningful results. The format used for the output is determined by the model's settings.
    """
        Vtemp = [r() for r in getattr(self, self.__reactions__[r])]
        outrate = ''
        for x in range(len(self.__reactions__)):
            outrate += self.__reactions__[x] + ' = ' + self.__settings__[
                'mode_number_format'] % Vtemp[x] + '\n'
        del s, Vtemp
        if File != None:
            print('\nReaction rates')
            File.write('\n##Reaction rates\n')
            File.write(outrate)
        else:
            print('\nReaction rates')
            print(outrate)

    def showRateEq(self, File=None):
        """
    showRateEq(File=None)

    Prints the reaction stoichiometry and rate equations to screen or to a file.

    Parameters
    ----------
    File : file-like object, optional
        An open, writable Python file object. If provided, the output is written to this file.
        Defaults to None, in which case the output is printed to the console.

    Notes
    -----
    The method retrieves reaction stoichiometry and rate equations from internal structures
    and formats them for output. The stoichiometry is displayed in a typical chemical equation
    format, with reactants and products separated by an arrow ('>' for irreversible reactions
    and '=' for reversible reactions).

    Any assignment rules that are part of the model will also be output if they exist.

    Examples
    --------
    Below is a basic example of how to call the `showRateEq` method in a `PysMod` object:

    >>> model = PysMod()
    >>> model.showRateEq()  # prints to console

    To write the output to a file instead:

    >>> with open('reaction_equations.txt', 'w') as f:
    >>>     model.showRateEq(File=f)

    This will write the reaction stoichiometry and rate equations to 'reaction_equations.txt'.
    """
        out_list = []
        tnform = self.__settings__['mode_number_format']
        self.__settings__['mode_number_format'] = '%2.5f'
        print('\nReaction stoichiometry and rate equations')
        out_list.append('\n## Reaction stoichiometry and rate equations\n')
        for key in self.Kmatrix.row:
            if File is None:
                print(key + ':')
            else:
                out_list.append(key + ':\n')
            reagL = []
            reagR = []
            for reagent in self.__nDict__[key]['Reagents']:
                if self.__nDict__[key]['Reagents'][reagent] > 0:
                    if self.__nDict__[key]['Reagents'][reagent] == 1.0:
                        reagR.append(reagent.replace('self.', ''))
                    else:
                        reagR.append('{' + self.__settings__[
                            'mode_number_format'] % abs(self.__nDict__[key]
                            ['Reagents'][reagent]) + '}' + reagent.replace(
                            'self.', ''))
                elif self.__nDict__[key]['Reagents'][reagent] < 0:
                    if self.__nDict__[key]['Reagents'][reagent] == -1.0:
                        reagL.append(reagent.replace('self.', ''))
                    else:
                        reagL.append('{' + self.__settings__[
                            'mode_number_format'] % abs(self.__nDict__[key]
                            ['Reagents'][reagent]) + '}' + reagent.replace(
                            'self.', ''))
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
            if self.__nDict__[key]['Type'] == 'Rever':
                symbol = ' = '
            else:
                symbol = ' > '
            if File is None:
                print('\t' + substring + symbol + prodstring)
                print('\t' + self.__nDict__[key]['RateEq'].replace('self.', '')
                    )
            else:
                out_list.append('\t' + substring + symbol + prodstring + '\n')
                out_list.append('\t' + self.__nDict__[key]['RateEq'].
                    replace('self.', '') + '\n\n')
        if len(list(self.__rules__.keys())) > 0:
            out_list.append('\n# Assignment rules\n')
            for ass in self.__rules__:
                out_list.append('!F {} = {}\n'.format(self.__rules__[ass][
                    'name'], self.__rules__[ass]['formula']))
            out_list.append('\n')
        if File != None:
            for x in out_list:
                File.write(x)
        self.__settings__['mode_number_format'] = tnform

    def showODE(self, File=None, fmt='%2.3f'):
        """
    Print the full set of ordinary differential equations (ODEs) generated by PySCeS to either the screen or a file.

    Parameters
    ----------
    File : file-like object, optional
        An open, writable Python file object where the ODE representation will be written.
        If `None`, the ODEs will be printed to the screen. Default is `None`.
    fmt : str, optional
        A string representing the format for number output. Default is `'%2.3f'`.

    Returns
    -------
    None

    Notes
    -----
    This function constructs and prints a formatted representation of the unreduced ODEs,
    based on the species and reactions available in the model. If rate rules are present,
    they will also be included in the output.

    Examples
    --------
    To print the ODEs to the screen using the default number format:

    >>> model = PysMod()
    >>> model.showODE()

    To write the ODEs to a file with a specific number format:

    >>> with open('odes.txt', 'w') as file:
    >>>     model.showODE(File=file, fmt='%1.2f')
    """
        maxmetlen = 0
        for x in self.__species__:
            if len(x) > maxmetlen:
                maxmetlen = len(x)
        maxreaclen = 0
        for x in self.__reactions__:
            if len(x) > maxreaclen:
                maxreaclen = len(x)
        odes = "\n## ODE's (unreduced)\n"
        for x in range(self.__nmatrix__.shape[0]):
            odes += 'd' + self.__species__[x] + (maxmetlen - len(self.
                __species__[x])) * ' ' + '|'
            beginline = 0
            for y in range(self.__nmatrix__.shape[1]):
                if abs(self.__nmatrix__[x, y]) > 0.0:
                    if self.__nmatrix__[x, y] > 0.0:
                        if beginline == 0:
                            odes += '  ' + fmt % abs(self.__nmatrix__[x, y]
                                ) + '*' + self.__reactions__[y] + (maxreaclen -
                                len(self.__reactions__[y])) * ' '
                            beginline = 1
                        else:
                            odes += ' + ' + fmt % abs(self.__nmatrix__[x, y]
                                ) + '*' + self.__reactions__[y] + (maxreaclen -
                                len(self.__reactions__[y])) * ' '
                    else:
                        if beginline == 0:
                            odes += ' -' + fmt % abs(self.__nmatrix__[x, y]
                                ) + '*' + self.__reactions__[y] + (maxreaclen -
                                len(self.__reactions__[y])) * ' '
                        else:
                            odes += ' - ' + fmt % abs(self.__nmatrix__[x, y]
                                ) + '*' + self.__reactions__[y] + (maxreaclen -
                                len(self.__reactions__[y])) * ' '
                        beginline = 1
            odes += '\n'
        if self.__HAS_RATE_RULES__:
            odes += '\n## Rate Rules:\n'
            for rule in self.__rate_rules__:
                odes += 'd{} | {}\n'.format(rule, self.__rules__[rule][
                    'formula'].replace('()', ''))
            odes += '\n'
        if File != None:
            print("\nODE's (unreduced)\n")
            File.write(odes)
        else:
            print(odes)

    def showODEr(self, File=None, fmt='%2.3f'):
        """
    Print a representation of the reduced set of ODE's generated by PySCeS to screen or file.

    This method outputs a representation of the reduced set of Ordinary Differential Equations (ODEs),
    derived using the PySCeS library. The ODEs can be printed to the screen or written to a specified file.
    This function is particularly useful for analyzing the dynamics of a reaction system in a reduced form.

    Parameters
    ----------
    File : file-like object, optional
        An open, writable Python file object where the ODEs will be printed. If not provided or set to 
        None, the ODEs are printed to the standard output (screen).
    fmt : str, optional
        The format string for the output numbers in the ODEs. Defaults to '%2.3f', which specifies the 
        floating-point precision.

    Examples
    --------
    >>> model = PysMod()  # Assume PysMod is initialized with necessary data
    >>> model.showODEr(fmt='%2.2f')

    This will print the reduced ODEs to the console with numbers formatted to two decimal places.

    >>> with open('odes_output.txt', 'w') as file:
    ...     model.showODEr(File=file)

    This will write the reduced ODEs to the file 'odes_output.txt' using the default number format.
    """
        maxmetlen = 0
        for x in self.__species__:
            if len(x) > maxmetlen:
                maxmetlen = len(x)
        maxreaclen = 0
        for x in self.__reactions__:
            if len(x) > maxreaclen:
                maxreaclen = len(x)
        odes = "\n## ODE's (reduced)\n"
        for x in range(self.nrmatrix.shape[0]):
            odes += 'd' + self.__species__[self.nrmatrix_row[x]] + (maxmetlen -
                len(self.__species__[self.nrmatrix_row[x]])) * ' ' + '|'
            beginline = 0
            for y in range(self.__nrmatrix__.shape[1]):
                if abs(self.__nrmatrix__[x, y]) > 0.0:
                    if self.__nrmatrix__[x, y] > 0.0:
                        if beginline == 0:
                            odes += '  ' + fmt % abs(self.__nrmatrix__[x, y]
                                ) + '*' + self.__reactions__[self.
                                nrmatrix_col[y]] + (maxreaclen - len(self.
                                __reactions__[self.nrmatrix_col[y]])) * ' '
                            beginline = 1
                        else:
                            odes += ' + ' + fmt % abs(self.__nrmatrix__[x, y]
                                ) + '*' + self.__reactions__[self.
                                nrmatrix_col[y]] + (maxreaclen - len(self.
                                __reactions__[self.nrmatrix_col[y]])) * ' '
                    else:
                        if beginline == 0:
                            odes += ' -' + fmt % abs(self.__nrmatrix__[x, y]
                                ) + '*' + self.__reactions__[self.
                                nrmatrix_col[y]] + (maxreaclen - len(self.
                                __reactions__[self.nrmatrix_col[y]])) * ' '
                        else:
                            odes += ' - ' + fmt % abs(self.__nrmatrix__[x, y]
                                ) + '*' + self.__reactions__[self.
                                nrmatrix_col[y]] + (maxreaclen - len(self.
                                __reactions__[self.nrmatrix_col[y]])) * ' '
                        beginline = 1
            odes += '\n'
        if self.__HAS_RATE_RULES__:
            odes += '\n## Rate Rules:\n'
            for rule in self.__rate_rules__:
                odes += 'd{} | {}\n'.format(rule, self.__rules__[rule][
                    'formula'].replace('()', ''))
            odes += '\n'
        if File != None:
            print("\nODE's (reduced)\n")
            File.write(odes)
            self.showConserved(File)
        else:
            print(odes)
            self.showConserved()

    def showModel(self, filename=None, filepath=None, skipcheck=0):
        """
    showModel(filename=None, filepath=None, skipcheck=0)

    The PySCeS 'save' command, prints the entire model to the screen or to a file in a 
    PSC (PySCeS Scripted Computing) format. This includes only basic model attributes 
    and does not save any additional functions.

    Parameters
    ----------
    filename : str or None, optional
        The name of the output PSC file. If `None`, the model is printed to the console. 
        Default is `None`.
    filepath : str or None, optional
        The directory path where the output PSC file should be saved. If `None`, 
        it uses the default model directory. Default is `None`.
    skipcheck : int, optional
        If set to a non-zero value, the method skips the check to see if the file 
        already exists and overwrites it automatically. Default is `0`.

    Notes
    -----
    This function only saves model attributes and does not save additional functions 
    defined in the model.

    Example
    -------
    To print the model to the console without saving to a file:

    >>> model = PysMod()
    >>> model.showModel()

    To save the model to a file named 'mymodel.psc' in the current directory:

    >>> model.showModel(filename='mymodel.psc')

    To save the model to a file named 'mymodel.psc' in a specific directory and skip 
    the file existence check:

    >>> model.showModel(filename='mymodel.psc', filepath='/path/to/directory', skipcheck=1)
    """
        if filepath is None:
            filepath = self.ModelDir
        if filename is None:
            print('\nFixed species')
            if len(self.__fixed_species__) == 0:
                print('<none>')
            else:
                for x in self.__fixed_species__:
                    print(x, end=' ')
            print(' ')
            self.showRateEq()
            self.showSpeciesI()
            self.showPar()
        else:
            FBuf = 0
            if type(filename) == str:
                filename = chkpsc(filename)
                go = 0
                loop = 0
                if skipcheck != 0:
                    loop = 1
                    go = 1
                while loop == 0:
                    try:
                        filex = os.path.join(filepath, filename)
                        if os.path.isfile(filex):
                            inp = input('\nfile ' + filex +
                                ' exists.\nOverwrite? ([y]/n) ')
                        else:
                            raise IOError
                        if inp == 'y' or inp == '':
                            go = 1
                            loop = 1
                        elif inp == 'n':
                            filename = input('\nfile "' + filename +
                                '" exists. Enter a new filename: ')
                            go = 1
                            filex = os.path.join(filepath, filename)
                            filename = chkpsc(filename)
                        else:
                            print('\nInvalid input')
                    except IOError:
                        print('\nfile "' + filex +
                            '" does not exist, proceeding')
                        loop = 1
                        go = 1
            else:
                print('\nI hope we have a filebuffer')
                if type(filename).__name__ == 'file' or type(filename
                    ).__name__ == 'StringIO' or type(filename
                    ).__name__ == 'TextIOWrapper':
                    go = 1
                    FBuf = 1
                    print('Seems like it')
                else:
                    go = 0
                    print('Are you sure you know what ur doing')
            if go == 1:
                if not FBuf:
                    if filex[-4:] == '.psc':
                        pass
                    else:
                        print('Assuming extension is .psc')
                        filex += '.psc'
                    outFile = open(filex, 'w')
                else:
                    outFile = filename
                header = '# PySCeS (' + __version__ + ') model input file \n'
                header += '# http://pysces.sourceforge.net/ \n'
                header += '# \n'
                header += '# Original input file: ' + self.ModelFile + '\n'
                header += '# This file generated: ' + time.strftime(
                    '%a, %d %b %Y %H:%M:%S') + '\n\n'
                outFile.write(header)
                outFile.write('## Fixed species\n')
                if len(self.__fixed_species__) == 0:
                    outFile.write('#  <none>')
                else:
                    outFile.write('FIX: ')
                    for x in self.__fixed_species__:
                        outFile.write(x + ' ')
                outFile.write('\n')
                self.showRateEq(File=outFile)
                self.showSpeciesI(File=outFile)
                self.showPar(File=outFile)
                if type(filename) == str or type(filename
                    ).__name__ == 'file' or type(filename
                    ).__name__ == 'TextIOWrapper':
                    outFile.close()

    def showN(self, File=None, fmt='%2.3f'):
        """
    Print the stoichiometric matrix (N), including row and column labels to the screen or a file.

    Parameters
    ----------
    File : file object, optional
        An open, writable Python file object where the matrix will be printed. 
        If None, the matrix is printed to the screen. Default is None.
    fmt : str, optional
        The format specification for numbers in the output. 
        Default is '%2.3f' which formats the numbers to three decimal places.

    Example
    -------
    >>> model = PysMod()
    >>> model.showN()
    Stoichiometric matrix (N)
    
        reac1    reac2   reac3   
    sp1    1.000  0.000 -1.000
    sp2    0.000  1.000  1.000
    sp3   -1.000 -1.000  0.000

    The above example assumes a model with three reactions (`reac1`, `reac2`, `reac3`) 
    and three species (`sp1`, `sp2`, `sp3`), illustrating how the stoichiometric matrix 
    is presented with the default formatting.

    Notes
    -----
    If any reaction or species name is longer than the space allocated to its label, 
    it will be truncated and an indication will be provided at the bottom of the output.
    """
        km_trunc = 0
        km_trunc2 = 0
        km = ''
        rtr = []
        ctr = []
        for a in range(len(self.nmatrix_col)):
            if a == 0:
                km += 9 * ' '
            if len(self.__reactions__[self.nmatrix_col[a]]) > 7:
                km += self.__reactions__[self.nmatrix_col[a]][:6] + '* '
                km_trunc = 1
                ctr.append(a)
            else:
                km += self.__reactions__[a] + (9 - len(self.__reactions__[a
                    ]) - 1) * ' '
        for x in range(len(self.nmatrix_row)):
            if len(self.__species__[x]) > 6:
                km += '\n' + self.__species__[x][:5] + '*' + 1 * ' '
                km_trunc2 = 1
                rtr.append(x)
            else:
                km += '\n' + self.__species__[x] + (8 - len(self.
                    __species__[x]) - 1) * ' '
            for y in range(len(self.nmatrix_col)):
                if y != 0:
                    km += 3 * ' '
                if self.__nmatrix__[x, y] >= 10.0:
                    km += ' ' + fmt % self.__nmatrix__[x, y]
                elif self.__nmatrix__[x, y] >= 0.0:
                    km += '  ' + fmt % self.__nmatrix__[x, y]
                elif abs(self.__nmatrix__[x, y]) >= 10.0:
                    km += fmt % self.__nmatrix__[x, y]
                else:
                    km += ' ' + fmt % self.__nmatrix__[x, y]
        if km_trunc:
            km += '\n\n' + '* indicates truncated (col) reaction\n('
            for x in range(len(self.nmatrix_col)):
                try:
                    ctr.index(x)
                    km += self.__reactions__[x] + ', '
                except:
                    pass
            km += ')'
        if km_trunc2:
            km += '\n\n' + '* indicates truncated (row) species:\n('
            for x in range(len(self.nmatrix_row)):
                try:
                    rtr.index(x)
                    km += self.__species__[x] + ', '
                except:
                    pass
            km += ')'
        if File != None:
            print('\nStoichiometric matrix (N)')
            File.write('\n\n## Stoichiometric matrix (N)\n')
            File.write(km)
            File.write('\n')
        else:
            print('\nStoichiometric matrix (N)')
            print(km)

    def showNr(self, File=None, fmt='%2.3f'):
        """
    Prints the reduced stoichiometric matrix (Nr), including row and column labels
    to either the screen or a specified file.

    Parameters
    ----------
    File : file object, optional
        An open, writable Python file object where the output should be written.
        If not provided, the output is printed to the screen.
    fmt : str, optional
        A format string for specifying the number format in the output. 
        Default is '%2.3f', which represents a floating point number with at
        least two characters and three digits after the decimal point.

    Examples
    --------
    >>> model = PysMod()
    >>> model.showNr()
    Prints the reduced stoichiometric matrix to the screen.

    >>> with open('output.txt', 'w') as file:
    ...     model.showNr(File=file)
    Writes the reduced stoichiometric matrix to 'output.txt'.

    Notes
    -----
    The method handles the truncation of long reaction and species names by 
    appending an asterisk (*) and provides an indicator for truncated entries 
    at the end of the output.

    """
        km_trunc = 0
        km_trunc2 = 0
        km = ''
        ctr = []
        rtr = []
        for a in range(len(self.nrmatrix_col)):
            if a == 0:
                km += 9 * ' '
            if len(self.__reactions__[self.nrmatrix_col[a]]) > 7:
                km += self.__reactions__[self.nrmatrix_col[a]][:6] + '* '
                km_trunc = 1
                ctr.append(a)
            else:
                km += self.__reactions__[a] + (9 - len(self.__reactions__[a
                    ]) - 1) * ' '
        for x in range(len(self.nrmatrix_row)):
            if len(self.__species__[self.nrmatrix_row[x]]) > 6:
                km += '\n' + self.__species__[self.nrmatrix_row[x]][:5
                    ] + '*' + 1 * ' '
                km_trunc2 = 1
                rtr.append(x)
            else:
                km += '\n' + self.__species__[self.nrmatrix_row[x]] + (8 -
                    len(self.__species__[self.nrmatrix_row[x]]) - 1) * ' '
            for y in range(len(self.nrmatrix_col)):
                if y != 0:
                    km += 3 * ' '
                if self.__nrmatrix__[x, y] >= 10.0:
                    km += ' ' + fmt % self.__nrmatrix__[x, y]
                elif self.__nrmatrix__[x, y] >= 0.0:
                    km += '  ' + fmt % self.__nrmatrix__[x, y]
                elif abs(self.__nrmatrix__[x, y]) >= 10.0:
                    km += fmt % self.__nrmatrix__[x, y]
                else:
                    km += ' ' + fmt % self.__nrmatrix__[x, y]
        if km_trunc:
            km += '\n\n' + '* indicates truncated (col) reaction:\n('
            for x in range(len(self.nrmatrix_col)):
                try:
                    ctr.index(x)
                    km += self.__reactions__[x] + ', '
                except:
                    pass
            km += ')'
        if km_trunc2:
            km += '\n\n' + '* indicates truncated (row) species:\n('
            for x in range(len(self.nrmatrix_row)):
                try:
                    rtr.index(x)
                    km += self.__species__[self.nrmatrix_row[x]] + ', '
                except:
                    pass
            km += ')'
        if File != None:
            print('\nReduced stoichiometric matrix (Nr)')
            File.write('\n\n## Reduced stoichiometric matrix (Nr)\n')
            File.write(km)
            File.write('\n')
        else:
            print('\nReduced stoichiometric matrix (Nr)')
            print(km)

    def showK(self, File=None, fmt='%2.3f'):
        """
    Print the Kernel matrix (K), including row and column labels. This can be 
    directed to the screen or an open, writable file object.

    Parameters
    ----------
    File : file-like object, optional
        An open, writable Python file object where the output will be written. 
        If `None`, the output is printed to the screen by default.
    fmt : str, optional
        The format string to control the output number format. Default is `'%2.3f'`.

    Notes
    -----
    If any reaction or row names are longer than 7 characters, they will be 
    truncated and marked with an asterisk (*). A section at the bottom will list 
    truncated names and their full form.

    If there is no flux conservation, a message will be printed stating "No flux 
    conservation" instead of the matrix. The presence or absence of flux 
    conservation is dependent on the attribute `__HAS_FLUX_CONSERVATION__`.

    Examples
    --------
    >>> model = PysMod()
    >>> model.showK()
    Kernel matrix (K)
    Reaction1 Reaction2* 
    ReactionA        0.000  1.000
    ReactionB*       1.000  0.000

    >>> with open('output.txt', 'w') as f:
    ...     model.showK(File=f)
    >>> # The matrix will be written to 'output.txt' with the specified format.
    """
        km_trunc = 0
        km_trunc2 = 0
        km = ''
        ctr = []
        rtr = []
        for a in range(len(self.kmatrix_col)):
            if a == 0:
                km += 8 * ' '
            if len(self.__reactions__[self.kmatrix_col[a]]) > 7:
                km += self.__reactions__[self.kmatrix_col[a]][:6] + '* '
                km_trunc = 1
                ctr.append(a)
            else:
                km += self.__reactions__[self.kmatrix_col[a]] + (9 - len(
                    self.__reactions__[self.kmatrix_col[a]]) - 1) * ' '
        for x in range(len(self.kmatrix_row)):
            if len(self.__reactions__[self.kmatrix_row[x]]) > 6:
                km += '\n' + self.__reactions__[self.kmatrix_row[x]][:5
                    ] + '*' + 1 * ' '
                km_trunc2 = 1
                rtr.append(x)
            else:
                km += '\n' + self.__reactions__[self.kmatrix_row[x]] + (8 -
                    len(self.__reactions__[self.kmatrix_row[x]]) - 1) * ' '
            for y in range(len(self.kmatrix_col)):
                if y != 0:
                    km += 4 * ' '
                if abs(self.__kmatrix__[x, y]) == 0.0:
                    km += ' ' + fmt % abs(self.__kmatrix__[x, y])
                elif self.__kmatrix__[x, y] < 0.0:
                    km += fmt % self.__kmatrix__[x, y]
                else:
                    km += ' ' + fmt % self.__kmatrix__[x, y]
        if km_trunc:
            km += '\n\n' + '* indicates truncated (col) reaction:\n('
            for x in range(len(self.kmatrix_col)):
                try:
                    ctr.index(x)
                    km += self.__reactions__[self.kmatrix_col[x]] + ', '
                except:
                    pass
            km += ')'
        if km_trunc2:
            km += '\n\n' + '* indicates truncated (row) reaction:\n('
            for x in range(len(self.kmatrix_row)):
                try:
                    rtr.index(x)
                    km += self.__reactions__[self.kmatrix_row[x]] + ', '
                except:
                    pass
            km += ')'
        if File != None:
            print('\nKernel matrix (K)')
            File.write('\n\n## Kernel matrix (K)\n')
            if not self.__HAS_FLUX_CONSERVATION__:
                File.write('No flux conservation\n')
            else:
                File.write(km)
                File.write('\n')
        else:
            print('\nKernel matrix (K)')
            if not self.__HAS_FLUX_CONSERVATION__:
                print('No flux conservation')
            else:
                print(km)

    def showL(self, File=None, fmt='%2.3f'):
        """
    Print the Link matrix (L), including row and column labels to screen or File.

    This method displays the Link matrix which provides insight into the relationship between different species in the model. 
    The matrix is printed on the screen by default, or can be written to a file if a file object is provided. 
    The elements of the matrix can be formatted using a specified format.

    Parameters
    ----------
    File : file-like object, optional
        An open, writable Python file object where the Link matrix will be written. 
        If None (default), the matrix is printed to the console.
    fmt : str, optional
        A string representing the number format to be used when displaying the matrix elements. 
        Defaults to '%2.3f', which formats the numbers with 3 decimal places.

    Notes
    -----
    If any species' name is longer than 7 characters, it is truncated in the display, and an asterisk (*) is appended to indicate truncation.
    If there is no moiety conservation within the model, this will be noted in the output.

    Examples
    --------
    >>> model = PysMod()
    >>> model.showL()
    Link matrix (L)
    # Output displayed on screen

    >>> with open('link_matrix.txt', 'w') as f:
    ...     model.showL(File=f)
    # Link matrix written to 'link_matrix.txt'
    """
        km_trunc = 0
        km_trunc2 = 0
        km = ''
        ctr = []
        rtr = []
        for a in range(len(self.lmatrix_col)):
            if a == 0:
                km += 8 * ' '
            if len(self.__species__[self.lmatrix_col[a]]) > 7:
                km += self.__species__[self.lmatrix_col[a]][:6] + '* '
                km_trunc = 1
                ctr.append(a)
            else:
                km += self.__species__[self.lmatrix_col[a]] + (9 - len(self
                    .__species__[self.lmatrix_col[a]]) - 1) * ' '
        for x in range(len(self.lmatrix_row)):
            if len(self.__species__[self.lmatrix_row[x]]) > 6:
                km += '\n' + self.__species__[self.lmatrix_row[x]][:5] + '* '
                km_trunc2 = 1
                rtr.append(x)
            else:
                km += '\n' + self.__species__[self.lmatrix_row[x]] + (8 -
                    len(self.__species__[self.lmatrix_row[x]]) - 1) * ' '
            for y in range(len(self.lmatrix_col)):
                if y != 0:
                    km += 4 * ' '
                if abs(self.__lmatrix__[x, y]) == 0.0:
                    km += ' ' + fmt % abs(self.__lmatrix__[x, y])
                elif self.__lmatrix__[x, y] < 0.0:
                    km += fmt % self.__lmatrix__[x, y]
                else:
                    km += ' ' + fmt % self.__lmatrix__[x, y]
        if km_trunc:
            km += '\n\n' + '* indicates truncated (col) species:\n('
            for x in range(len(self.lmatrix_col)):
                try:
                    ctr.index(x)
                    km += self.__species__[self.lmatrix_col[x]] + ', '
                except:
                    pass
            km += ')'
        if km_trunc2:
            km += '\n\n' + '* indicates truncated (row) species:\n('
            for x in range(len(self.lmatrix_row)):
                try:
                    rtr.index(x)
                    km += self.__species__[self.lmatrix_row[x]] + ', '
                except:
                    pass
            km += ')'
        if File != None:
            print('\nLink matrix (L)')
            File.write('\n\n## Link matrix (L)\n')
            if not self.__HAS_MOIETY_CONSERVATION__:
                File.write('No moiety conservation\n')
                File.write(km)
                File.write('\n')
            else:
                File.write(km)
                File.write('\n')
        else:
            print('\nLink matrix (L)')
            if not self.__HAS_MOIETY_CONSERVATION__:
                print('"No moiety conservation"\n')
                print(km)
            else:
                print(km)

    def TestSimState(self, endTime=10000, points=101, diff=1e-05):
        """
    **Deprecated**

    This method is deprecated and may be removed in future versions. 
    It was intended to test the simulation state over a specified period,
    but it no longer performs any operation.

    Parameters
    ----------
    endTime : int, optional
        The end time of the simulation in arbitrary units. Default is 10000.
    points : int, optional
        The number of data points in the simulation. Default is 101.
    diff : float, optional
        The difference threshold for some unspecified simulation state comparison. Default is 1e-05.

    Notes
    -----
    Since this method is deprecated, it is advised not to use it in new code.
    There may not be any functionality linked with the specified parameters.

    Examples
    --------
    This method does not perform any operation; thus, calling it has no effect:

    >>> model = PysMod()
    >>> model.TestSimState(endTime=5000, points=50, diff=1e-4)
    """
        pass

    def ResetNumberFormat(self):
        """
    Reset the PySCeS default number format.

    This method resets the default number format used in PySCeS, stored as 
    `mod.mode_number_format`, to the scientific notation format `%2.4e`.

    Parameters
    ----------
    None

    Examples
    --------
    >>> pysmod_instance = PysMod()
    >>> pysmod_instance.ResetNumberFormat()
    >>> print(pysmod_instance._PysMod__settings__['mode_number_format'])
    '%2.4e'

    Notes
    -----
    This method modifies the internal settings of the PysMod class to ensure that numerical
    outputs are formatted consistently in scientific notation with four decimal places of precision.
    """
        self.__settings__['mode_number_format'] = '%2.4e'

    def Write_array(self, inputd, File=None, Row=None, Col=None, close_file
        =0, separator='  '):
        """
    Write_array(input, File=None, Row=None, Col=None, close_file=0, separator='  ')

    Writes an array to a file with optional row and column labels. You can specify a  
    comma `,` separator to generate a CSV-style file.

    The following settings influence the behavior of this method:
    
    - `mod.__settings__['write_array_header']`: Adds the filename as a header line (1 = yes, 0 = no).
    - `mod.__settings__['write_array_spacer']`: Adds a space after the header line (1 = yes, 0 = no).
    - `mod.__settings__['write_arr_lflush']`: Sets the flush rate for large file writes.

    Parameters
    ----------
    input : array-like
        The data array to be written to the file.
    File : file object, optional
        An open, writable Python file object. If `None`, an error will be logged.
    Row : list, optional
        A list of row labels.
    Col : list, optional
        A list of column labels. Must match the number of columns in `input`.
    close_file : int, optional
        Whether to close the file after writing (1) or leave it open (0). Default is 0.
    separator : str, optional
        The column separator to use. Default is `'  '` (double space).

    Examples
    --------
    Write an array to a file with headers and separators:
    
    >>> import numpy as np
    >>> array = np.array([[1.5, 2.3], [3.0, 4.8]])
    >>> with open('output.txt', 'w') as f:
    ...     mod = PysMod()
    ...     mod.__settings__ = {'write_array_header': 1, 'write_array_spacer': 1, 'write_arr_lflush': 3}
    ...     mod.Write_array(array, File=f, Row=['Row1', 'Row2'], Col=['Col1', 'Col2'], separator=', ')

    The output in 'output.txt' will be formatted based on the given settings. Be sure to update `PysMod` 
    instance settings before calling `Write_array` for the appropriate file format.
    """
        fname = 'Write_array_' + self.ModelFile.replace('.psc', ''
            ) + '_' + time.strftime('%H:%M:%S')
        if type(inputd) == tuple or type(inputd) == list:
            inputd = numpy.array(inputd)
        if Row != None:
            assert len(Row) == inputd.shape[0
                ], 'len(Row) must be equal to len(input_row)'
        if Col != None:
            assert len(Col) == inputd.shape[1
                ], 'len(Col) must be equal to len(input_col)'
        if File != None:
            if self.__settings__['write_array_header']:
                File.write('\n## ' + fname + '\n')
            if self.__settings__['write_array_spacer']:
                File.write('\n')
            if Col != None:
                print('Writing column')
                File.write('#')
                try:
                    input_width = len(self.__settings__[
                        'mode_number_format'] % inputd[0, 0]) + 1 + len(str
                        (separator))
                except:
                    input_width = len(self.__settings__[
                        'mode_number_format'] % inputd[0]) + 1 + len(str(
                        separator))
                for x in Col:
                    if len(str(x)) <= len(str(inputd[0, 0])):
                        spacer = (input_width - len(str(x))) * ' '
                        File.write(str(x) + spacer)
                    else:
                        spacer = len(str(separator)) * ' '
                        File.write(str(x[:input_width]) + spacer)
                File.write('\n')
            try:
                print('\nWriting array (normal) to file')
                flush_count = 0
                if self.__settings__['write_arr_lflush'] < 2:
                    print('INFO: LineFlush must be >= 2')
                    self.__settings__['write_arr_lflush'] = 2
                for x in range(inputd.shape[0]):
                    flush_count += 1
                    for y in range(inputd.shape[1]):
                        if inputd[x, y] < 0.0:
                            File.write(self.__settings__[
                                'mode_number_format'] % inputd[x, y])
                        else:
                            File.write(' ' + self.__settings__[
                                'mode_number_format'] % inputd[x, y])
                        if y < inputd.shape[1] - 1:
                            File.write(separator)
                    File.write('\n')
                    if flush_count == self.__settings__['write_arr_lflush']:
                        File.flush()
                        flush_count = 0
            except IndexError as e:
                print('\nWriting vector (normal) to file')
                for x in range(len(inputd)):
                    File.write(self.__settings__['mode_number_format'] %
                        inputd[x])
                    if x < len(inputd) - 1:
                        File.write(separator)
                File.write('\n')
            if Row != None:
                print('Writing row')
                File.write('# Row: ')
                for x in Row:
                    File.write(str(x) + separator)
                File.write('\n')
            if close_file:
                File.close()
        else:
            print(
                'INFO: You need to supply an open writable file as the 2nd argument to this function'
                )

    def Write_array_latex(self, inputa, File=None, Row=None, Col=None,
        close_file=0):
        """
    Write an array to an open file as a 'LaTeX' {array}.

    This method writes a given array to a specified file in the LaTeX array format. It provides
    options to include row and column labels, and allows for the file to be closed after writing.

    Parameters
    ----------
    inputa : array-like
        The array to be written.
    File : file object, optional
        An open, writable Python file object. If not provided, an error message will be printed.
    Row : list, optional
        A list of row labels. The length of the list must be equal to the number of rows in `inputa`.
    Col : list, optional
        A list of column labels. The length of the list must be equal to the number of columns in `inputa`.
    close_file : int, optional, default: 0
        Specifies whether to close the file after writing (1 to close, 0 to leave open).

    Raises
    ------
    AssertionError
        If the length of `Row` does not match the number of rows in `inputa` or if the length
        of `Col` does not match the number of columns in `inputa`.

    Notes
    -----
    - Ensure that the file passed to the `File` parameter is open and writable.
    - If the `Row` or `Col` options are provided, the arrays must be two-dimensional.

    Example
    -------
    Here is an example of how to use `Write_array_latex`:

    >>> import numpy as np
    >>> from io import StringIO
    >>> mod = PysMod()
    >>> array_data = np.array([[1.0, 2.0], [3.0, 4.0]])
    >>> file_buffer = StringIO()
    >>> mod.Write_array_latex(array_data, File=file_buffer, Row=['Row1', 'Row2'], Col=['Col1', 'Col2'], close_file=1)
    >>> print(file_buffer.getvalue())
    %% Write_array_latex_example_model_HH:MM:SS
    \\[
    \\begin{array}{r|rr}
     & \\small{Col1} & \\small{Col2} \\\\
    \\hline
    \\small{Row1} &  1.0000 &  2.0000 \\\\
    \\small{Row2} &  3.0000 &  4.0000 \\\\
    \\end{array}
    \\]
    """
        if type(inputa) == tuple or type(inputa) == list:
            inputa = numpy.array(inputa)
        try:
            a = inputa.shape[1]
        except:
            inputa.shape = 1, inputa.shape[0]
        if Row != None:
            assert len(Row) == inputa.shape[0
                ], 'len(Row) must be equal to len(input_row)'
            print('Writing row')
        if Col != None:
            assert len(Col) == inputa.shape[1
                ], 'len(Col) must be equal to len(input_col)'
            print('Writing column')
        if self.__settings__['write_arr_lflush'] < 2:
            print('INFO: LineFlush must be >= 2')
            self.__settings__['write_arr_lflush'] = 2
        fname = 'Write_array_latex_' + self.ModelFile.replace('.psc', ''
            ) + '_' + time.strftime('%H:%M:%S')
        if File != None:
            File.write('\n%% ' + fname + '\n')
            try:
                a = inputa.shape
                print('\nWriting array (LaTeX) to file')
                File.write('\\[\n')
                File.write('\\begin{array}{')
                if Row != None:
                    File.write('r|')
                for x in range(inputa.shape[1]):
                    File.write('r')
                File.write('}\n ')
                flush_count = 0
                for x in range(inputa.shape[0]):
                    if Col != None and x == 0:
                        for el in range(len(Col)):
                            elx = str(Col[el]).replace('_', '\\_')
                            if Row != None:
                                File.write(' & $\\small{' + elx + '}$')
                            elif el == 0:
                                File.write(' $\\small{' + elx + '}$')
                            else:
                                File.write(' & $\\small{' + elx + '}$')
                        File.write(' \\\\ \\hline\n ')
                    flush_count += 1
                    if Row != None:
                        el2 = str(Row[x]).replace('_', '\\_')
                        File.write('$\\small{' + el2 + '}$ &')
                    for y in range(inputa.shape[1]):
                        if inputa[x, y] < 0.0:
                            val = '%2.4f' % inputa[x, y]
                        else:
                            val = ' ' + '%2.4f' % inputa[x, y]
                        File.write(val)
                        if y < inputa.shape[1] - 1:
                            File.write(' &')
                    File.write(' \\\\\n ')
                    if flush_count == self.__settings__['write_arr_lflush']:
                        File.flush()
                        flush_count = 0
                File.write('\\end{array}\n')
                File.write('\\]\n\n')
            except:
                print(
                    '\nINFO: Only arrays can currently be processed with this method.\n'
                    )
            if close_file:
                File.close()
        else:
            print(
                'INFO: You need to supply an open writable file as the 2nd argument to this method'
                )

    def Write_array_html(self, inputa, File=None, Row=None, Col=None, name=
        None, close_file=0):
        """
    Write an array as an HTML table or complete document, including optional row and column labels.
    
    Tables are formatted with colored columns if they exceed a specified size. The method can include
    or exclude an HTML header and footer based on module settings.

    Parameters
    ----------
    input : ndarray
        The array to be written as an HTML table.
    File : file-like object, optional
        An open, writable Python file object where the HTML content will be written. Defaults to `None`.
    Row : list, optional
        A list of row labels corresponding to the array's rows. Defaults to `None`.
    Col : list, optional
        A list of column labels corresponding to the array's columns. Defaults to `None`.
    name : str, optional
        An HTML table description line to be included above the table. Defaults to `None`.
    close_file : int, optional
        Specify whether to close the file after writing (1) or keep it open (0). Defaults to `0`.

    Settings
    --------
    mod.__settings__['write_array_html_header'] : bool
        Determines whether to write the HTML document header.
    mod.__settings__['write_array_html_footer'] : bool
        Determines whether to write the HTML document footer.
    
    Notes
    -----
    If the `File` argument is not provided, this method will print an informational message
    reminding you to supply a writable file object.

    Example
    -------
    >>> import numpy as np
    >>> data = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    >>> with open('output.html', 'w') as f:
    ...     mod.Write_array_html(data, File=f, Row=['Row1', 'Row2', 'Row3'], Col=['Col1', 'Col2', 'Col3'], name='My Table', close_file=1)
    This will create an HTML file named 'output.html' with a table representation of `data`, 
    including specified row and column labels, and a table description.

    """
        if type(inputa) == tuple or type(inputa) == list:
            inputa = numpy.array(inputa)
        try:
            a = inputa.shape[1]
        except:
            inputa.shape = 1, inputa.shape[0]
        if Row != None:
            assert len(Row) == inputa.shape[0
                ], 'len(Row) must be equal to len(input_row)'
            print('Writing row')
        if Col != None:
            assert len(Col) == inputa.shape[1
                ], 'len(Col) must be equal to len(input_col)'
            print('Writing column')
        if self.__settings__['write_arr_lflush'] < 2:
            print('INFO: LineFlush must be >= 2')
            self.__settings__['write_arr_lflush'] = 2
        if name != None:
            fname = ('PySCeS data "' + name +
                '" generated from model file: ' + self.ModelFile + ' ' +
                time.strftime('%H:%M:%S'))
        else:
            fname = ('PySCeS data generated from model file: ' + self.
                ModelFile + ' ' + time.strftime('%H:%M:%S'))
        if File != None:
            header = '\n'
            header += (
                '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">\n'
                )
            header += '<html>\n'
            header += '<head>\n'
            header += ('<title>PySCeS data generated from: ' + self.
                ModelFile + ' - ' + time.strftime('%H:%M:%S (%Z)') +
                '</title>\n')
            header += (
                '<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">\n'
                )
            header += '</head>\n'
            header += '<body>\n\n'
            header += (
                '<h4><a href="http://pysces.sourceforge.net">PySCeS</a> data generated from: '
                 + self.ModelFile + '</h4>\n\n')
            if self.__settings__['write_array_html_header']:
                File.write(header)
            File.write('<!-- ' + fname + '-->\n')
            if name != None:
                File.write(
                    '\n<p>\n<font face="Arial, Helvetica, sans-serif"><strong>'
                     + name + '</strong></font>')
            else:
                File.write('\n<p>\n')
            File.write(
                '\n<table border="1" cellpadding="2" cellspacing="2" bgcolor="#FFFFFF">'
                )
            try:
                a = inputa.shape
                print('Writing array (HTML) to file')
                double_index = 15
                if Col != None:
                    File.write('\n<tr>\n')
                    if Row != None:
                        File.write('  <td>&nbsp;</td>\n')
                    for col in Col:
                        File.write(
                            '  <td bgcolor="#CCCCCC"><div align="center"><b>' +
                            str(col) + '</b></div></td>\n')
                    if Row != None and inputa.shape[1] + 1 >= double_index:
                        File.write('  <td>&nbsp;</td>\n')
                    File.write('</tr>')
                flush_count = 0
                colour_count_row = 6
                colour_count_col = 6
                rowcntr = 0
                html_colour = '#FFFFCC'
                for x in range(inputa.shape[0]):
                    rowcntr += 1
                    if rowcntr == colour_count_row and inputa.shape[0
                        ] > 3 * colour_count_row:
                        File.write('\n<tr bgcolor="' + html_colour + '">\n')
                        rowcntr = 0
                    else:
                        File.write('\n<tr>\n')
                    if Row != None:
                        File.write(
                            '  <td bgcolor="#CCCCCC"><div align="center"><b>' +
                            str(Row[x]) + '</b></div></td>\n')
                    colcntr = 0
                    for y in range(inputa.shape[1]):
                        colcntr += 1
                        if colcntr == colour_count_col and inputa.shape[1
                            ] > 3 * colour_count_row:
                            File.write('  <td nowrap bgcolor="' +
                                html_colour + '">' + self.__settings__[
                                'write_array_html_format'] % inputa[x, y] +
                                '</td>\n')
                            colcntr = 0
                        else:
                            File.write('  <td nowrap>' + self.__settings__[
                                'write_array_html_format'] % inputa[x, y] +
                                '</td>\n')
                    if Row != None and inputa.shape[1] + 1 >= double_index:
                        File.write(
                            '  <td bgcolor="#CCCCCC"><div align="center"><b>' +
                            Row[x] + '</b></div></td>\n')
                    File.write('</tr>')
                    if flush_count == self.__settings__['write_arr_lflush']:
                        File.flush()
                        flush_count = 0
                if Col != None and inputa.shape[0] + 1 >= double_index:
                    File.write('\n<tr>\n')
                    if Row != None:
                        File.write('  <td>&nbsp;</td>\n')
                    for col in Col:
                        File.write(
                            '  <td bgcolor="#CCCCCC"><div align="center"><b>' +
                            col + '</b></div></td>\n')
                    if Row != None and inputa.shape[1] + 1 >= double_index:
                        File.write('  <td>&nbsp;</td>\n')
                    File.write('</tr>\n')
            except:
                print(
                    '\nINFO: Only arrays can currently be processed with this method.\n'
                    )
            File.write('\n</table>\n')
            File.write('</p>\n\n')
            if self.__settings__['write_array_html_footer']:
                try:
                    File.write(
                        '<p><a href="http://pysces.sourceforge.net"><font size="3">PySCeS '
                         + __version__ +
                        '</font></a><font size="2"> HTML output (model <i>' +
                        self.ModelFile + '</i> analysed at ' + time.
                        strftime('%H:%M:%S') + ' by <i>' + getuser() +
                        """</i>)</font></p>
""")
                except:
                    File.write(
                        '<p><a href="http://pysces.sourceforge.net"><font size="3">PySCeS '
                         + __version__ +
                        '</font></a><font size="2"> HTML output (model <i>' +
                        self.ModelFile + '</i> analysed at ' + time.
                        strftime('%H:%M:%S - %Z') + ')</font></p>\n')
                File.write('</body>\n')
                File.write('</html>\n')
            else:
                File.write('\n')
            if close_file:
                File.close()
        else:
            print(
                'INFO: You need to supply an open writable file as the 2nd argument to this method'
                )

    def SimPlot(self, plot='species', filename=None, title=None, log=None,
        format='lines'):
        """
    Plot the simulation results using the new UPI pysces.plt interface.

    Parameters
    ----------
    plot : str or list, optional
        The output to plot. Default is 'species'. Options are:
        - 'all': Plot both rates and species.
        - 'species': Plot species.
        - 'rates': Plot reaction rates.
        - list: A list of model attributes (species, rates), e.g., `['S1', 'R1']`.
    filename : str, optional
        If not None, the file is exported to this filename. Default is None.
    title : str, optional
        The plot title. Default is None.
    log : str, optional
        Use log axis for 'x', 'y', 'xy'. Default is None.
    format : str, optional
        Line format. It is backend dependent. Default is 'lines'.

    Raises
    ------
    RuntimeError
        If the `plot` argument is not one of the allowed options.

    Notes
    -----
    This method will plot the simulation results. Depending on the chosen `plot` argument, 
    it can plot all data, only species, only rates, or selected attributes. The matplotlib 
    equivalent of axes labels and title are set based on whether 'Output_In_Conc' is enabled.

    Examples
    --------
    >>> model = PysMod()
    >>> # Plot all data (species and rates) and save to 'output.png'
    >>> model.SimPlot(plot='all', filename='output.png')

    >>> # Plot specific species
    >>> model.SimPlot(plot=['S1', 'S2'], title='Selected Species Plot')

    >>> # Plot rates with a log scale for x-axis
    >>> model.SimPlot(plot='rates', log='x', format='dashed')

    """
        data = None
        labels = None
        allowedplots = ['all', 'species', 'rates']
        if type(plot) != list and plot not in allowedplots:
            raise RuntimeError('\nPlot must be one of {} not "{}"'.format(
                str(allowedplots), plot))
        if plot == 'all':
            data, labels = self.data_sim.getAllSimData(lbls=True)
        elif plot == 'species':
            data, labels = self.data_sim.getSpecies(lbls=True)
        elif plot == 'rates':
            data, labels = self.data_sim.getRates(lbls=True)
        else:
            plot = [at for at in plot if at in self.__species__ + self.
                __reactions__ + [self.data_sim.time_label]]
            kwargs = {'lbls': True}
            if len(plot) > 0:
                data, labels = self.data_sim.getSimData(*plot, **kwargs)
        del allowedplots
        self.__SIMPLOT_OUT__ = labels
        plt.plotLines(data, 0, list(range(1, data.shape[1])), titles=labels,
            formats=[format])
        RngTime = self.data_sim.getTime()
        end = RngTime[-1] + 0.2 * RngTime[-1]
        plt.setRange('x', RngTime[0], end)
        del RngTime
        if self.__KeyWords__['Output_In_Conc']:
            M = 'Concentration'
        else:
            M = (
                'Amount (%(multiplier)s x %(kind)s x 10**%(scale)s)**%(exponent)s'
                 % self.__uDict__['substance'])
        xu = (
            'Time (%(multiplier)s x %(kind)s x 10**%(scale)s)**%(exponent)s' %
            self.__uDict__['time'])
        plt.setAxisLabel('x', xu)
        if plot == 'all':
            yl = 'Rates, {}'.format(M)
        elif plot == 'rates':
            yl = 'Rate'
        elif plot == 'species':
            yl = '{}'.format(M)
        else:
            yl = 'Rates, {}'.format(M)
        plt.setAxisLabel('y', yl)
        if log != None:
            plt.setLogScale(log)
        if title is None:
            plt.setGraphTitle('PySCeS Simulation (' + self.ModelFile + ') ' +
                time.strftime('%a, %d %b %Y %H:%M:%S'))
        else:
            plt.setGraphTitle(title)
        plt.replot()
        if filename != None:
            plt.export(filename, directory=self.ModelOutput, outtype='png')

    def Scan1(self, range1=[], runUF=0):
        """
    Perform a single dimension parameter scan using the steady-state solvers.

    The parameter to be scanned is defined as a model attribute in `mod.scan_in`,
    while the required output is specified in the list `mod.scan_out`. This method
    allows exploration of model behavior across parameter variations and can be
    easily visualized using `Scan1Plot()`.

    Attributes
    ----------
    mod.scan_in : str
        A model attribute that defines the parameter to be scanned (e.g., 'P', 'Vmax1').
    mod.scan_out : list of str
        A list specifying the desired output attributes (e.g., ['A','T2', 'ecR1_s1', 'ccJR1_R1', 'rcJR1_s1']).
    mod.scan_res : ndarray
        Contains the results of the parameter scan.
    mod.scan : numpy.record array
        The scan results including `scan_in` and `scan_out`. Accessed as `mod.scan.Vmax`, `mod.scan.A_ss`, `mod.scan.J_R1`, etc.
    mod.__settings__['scan1_mca_mode'] : int
        Modifies the scan algorithm to evaluate elasticities (1) or control coefficients (2).

    Parameters
    ----------
    range1 : list, optional
        A predefined range over which to scan. Default is an empty list `[]`.
    runUF : int, optional
        Indicates whether to run the user-defined function `mod.User_Function` (1) before evaluating the steady state. Default is `0`.

    Examples
    --------
    >>> model = PysMod()
    >>> model.scan_in = 'P'
    >>> model.scan_out = ['A', 'T2', 'ccJR1_R1']
    >>> model.Scan1(range1=[0.1, 0.2, 0.3, 0.4], runUF=1)
    >>> print(model.scan_res)
    [[0.1, 0.023, 0.458]
     [0.2, 0.027, 0.541]
     [0.3, 0.031, 0.632]
     [0.4, 0.035, 0.714]]

    Notes
    -----
    - The `mod.scan_in` parameter must be set to a valid model attribute before running this method.
    - This method also sets specific settings for evaluating elasticities and control coefficients, which can be auto-detected or explicitly forced using the `mod.__settings__['scan1_mca_mode']`.
    - Results of the scan are stored in `mod.scan_res` and can be visualized with the `Scan1Plot()` method if desired.
    """
        if self.__settings__['scan1_dropbad'] != 0:
            print(
                """
****
INFO: Dropping invalid steady states can mask interesting behaviour
****
"""
                )
        run = 1
        parameters = list(self.__parameters__) + [(a + '_init') for a in
            self.__species__]
        scanpar = []
        scanIn = [self.scan_in]
        for x in scanIn:
            try:
                a = parameters.index(x)
                scanpar.append(x)
            except:
                print(x, ' is not a parameter')
        if len(scanpar) < 1:
            print('mod.scan_in should contain something')
            run = 0
        elif len(scanpar) > 1:
            print('Only 1D scans possible - first element will be used')
            scanpar = scanpar[0]
            print('self.scan_in = ', scanpar)
        else:
            scanpar = scanpar[0]
        self.scan_in = scanpar
        wrong = []
        for x in self.scan_out:
            if x[:2] == 'ec' and self.__settings__['scan1_mca_mode'] < 2:
                self.__settings__['scan1_mca_mode'] = 1
            elif x[:2] == 'cc':
                self.__settings__['scan1_mca_mode'] = 2
            elif x[:2] == 'rc':
                self.__settings__['scan1_mca_mode'] = 3
        if self.__settings__['scan1_mca_mode'] == 1:
            self.doElas()
        elif self.__settings__['scan1_mca_mode'] == 2:
            self.doMca()
        elif self.__settings__['scan1_mca_mode'] == 3:
            self.doMcaRC()
        else:
            self.doState()
        for x in self.scan_out:
            try:
                if x.startswith('ec'):
                    try:
                        getattr(self.ec, x[2:])
                    except AttributeError:
                        getattr(self.ecp, x[2:])
                elif x.startswith('cc'):
                    if x.startswith('ccJ'):
                        getattr(self.cc, x[3:])
                    else:
                        getattr(self.cc, x[2:])
                elif x.startswith('rc'):
                    if x.startswith('rcJ'):
                        getattr(self.rc, x[3:])
                    else:
                        getattr(self.rc, x[2:])
                elif x in self.reactions:
                    getattr(self.data_sstate, x)
                    print(
                        'INFO: using steady-state flux for reaction ({} --> J_{})'
                        .format(x, x))
                    self.scan_out[self.scan_out.index(x)] = 'J_{}'.format(x)
                elif x.startswith('J_'):
                    getattr(self.data_sstate, x[2:])
                elif x in self.species:
                    getattr(self.data_sstate, x)
                    print(
                        'INFO: using steady-state concentration for species ({} --> {}_ss)'
                        .format(x, x))
                    self.scan_out[self.scan_out.index(x)] = '{}_ss'.format(x)
                elif x.endswith('_ss'):
                    getattr(self.data_sstate, x[:-3])
                else:
                    getattr(self, x)
            except:
                print(x + ' is not a valid attribute')
                wrong.append(x)
            if x == self.scan_in:
                wrong.append(x)
                print(x, ' is mod.scan_in ')
        if len(wrong) != 0:
            try:
                for x in wrong:
                    print(self.scan_out)
                    self.scan_out.remove(x)
                    print(self.scan_out)
            except:
                print('No valid output')
                run = 0
        assert len(self.scan_out
            ) != 0, 'Output parameter list (mod.scan_out) empty - do the model attributes exist ... see manual for details.'
        self._scan = None
        result = []
        cntr = 0
        cntr2 = 1
        cntr3 = 0
        self.__scan_errors_par__ = None
        self.__scan_errors_idx__ = None
        if self.__settings__['scan1_mesg']:
            print('\nScanning ...')
        if len(self.scan_in) > 0 and run == 1:
            badList = []
            badList_idx = []
            if self.__settings__['scan1_mesg']:
                print(len(range1) - cntr * cntr2, end=' ')
            for xi in range(len(range1)):
                x = range1[xi]
                setattr(self, self.scan_in, x)
                if self.__settings__['scan1_mca_mode'] == 1:
                    self.doElas()
                elif self.__settings__['scan1_mca_mode'] == 2:
                    self.doMca()
                elif self.__settings__['scan1_mca_mode'] == 3:
                    self.doMcaRC()
                else:
                    self.State()
                rawres = [x]
                if runUF:
                    self.User_Function()
                for res in self.scan_out:
                    if res.startswith('ec'):
                        try:
                            rawres.append(getattr(self.ec, res[2:]))
                        except AttributeError:
                            rawres.append(getattr(self.ecp, res[2:]))
                    elif res.startswith('cc'):
                        if res.startswith('ccJ'):
                            rawres.append(getattr(self.cc, res[3:]))
                        else:
                            rawres.append(getattr(self.cc, res[2:]))
                    elif res.startswith('rc'):
                        if res.startswith('rcJ'):
                            rawres.append(getattr(self.rc, res[3:]))
                        else:
                            rawres.append(getattr(self.rc, res[2:]))
                    elif res in self.reactions:
                        rawres.append(getattr(self.data_sstate, res))
                    elif res.startswith('J_'):
                        rawres.append(getattr(self.data_sstate, res[2:]))
                    elif res in self.species:
                        rawres.append(getattr(self.data_sstate, res))
                    elif res.endswith('_ss'):
                        rawres.append(getattr(self.data_sstate, res[:-3]))
                    else:
                        rawres.append(getattr(self, res))
                if not self.__StateOK__ and self.__settings__['scan1_dropbad'
                    ] != 0:
                    pass
                elif not self.__StateOK__:
                    if self.__settings__['scan1_nan_on_bad']:
                        result.append([numpy.nan] * len(rawres))
                    else:
                        result.append(rawres)
                    badList.append(x)
                    badList_idx.append(xi)
                else:
                    result.append(rawres)
                cntr += 1
                cntr3 += 1
                if cntr == 20:
                    if self.__settings__['scan1_mesg']:
                        print(len(range1) - cntr * cntr2, end=' ')
                    cntr = 0
                    cntr2 += 1
                if cntr3 == 101:
                    if self.__settings__['scan1_mesg']:
                        print(' ')
                    cntr3 = 0
        if self.__settings__['scan1_mesg']:
            print('\ndone.\n')
        if len(badList) != 0:
            self.__scan_errors_par__ = badList
            self.__scan_errors_idx__ = badList_idx
            print('\nINFO: ' + str(len(badList)) +
                ' invalid steady states detected at ' + self.scan_in +
                ' values:')
            print(badList)
        if len(result) == 0:
            self.scan_res = numpy.zeros((len(range1), len(self.scan_out) + 1))
        else:
            self.scan_res = numpy.array(result)

    @property
    def scan(self):
        """
    Retrieve the results of a parameter scan as a data structure.

    This method returns the results of the parameter scan performed
    by the `Scan1` method. The results can be formatted as a
    Pandas DataFrame or a NumPy record array, depending on the
    configuration setting `custom_datatype`.

    Returns
    -------
    pandas.DataFrame or numpy.recarray
        The parameter scan results, formatted according to the `custom_datatype`
        setting. The columns include both the scanned parameter and the
        specified outputs.

    Notes
    -----
    - If `scan_res` is not available, the method returns None.
    - Before generating the results, ensure that a parameter scan has been
      executed using the `Scan1` method to populate `scan_res`.

    Example
    -------
    >>> mod = PysMod()
    >>> mod.scan_in = 'Vmax1'
    >>> mod.scan_out = ['A', 'T2']
    >>> mod.Scan1(range1=[0.1, 0.2, 0.3])
    >>> scan_results = mod.scan
    >>> print(scan_results)
          Vmax1    A_ss    T2_ss
    0     0.1     1.25    0.85
    1     0.2     1.30    0.88
    2     0.3     1.35    0.90

    """
        if self._scan is None and self.scan_res is not None:
            if self.__settings__['custom_datatype'] == 'pandas':
                self._scan = self.__pandas.DataFrame(self.scan_res, columns
                    =[self.scan_in] + self.scan_out)
            else:
                self._scan = numpy.rec.fromrecords(self.scan_res, names=[
                    self.scan_in] + self.scan_out)
        return self._scan

    def Scan1Plot(self, plot=[], title=None, log=None, format='lines',
        filename=None):
        """
    Plots the results of a parameter scan generated with `Scan1()`.

    This method visualizes the data obtained from a parameter scan,
    allowing the user to specify various plotting options such as the
    title, log scale, and file export options.

    Parameters
    ----------
    plot : list of str, optional
        Specifies which outputs to plot. If empty, all outputs specified in
        `mod.scan_out` are used. Must be a subset of `mod.scan_out` (default=[]).
    title : str, optional
        The plot title. If None, a default title is used based on the model file
        and current date and time (default=None).
    log : {'x', 'xy', 'xyz'}, optional
        Specifies the axis to use a logarithmic scale. If None, linear axes are
        used (default=None).
    format : str, optional
        The line format for the plot, which can be backend dependent. Common
        styles are 'lines' or 'points' (default='lines').
    filename : str, optional
        The filename for exporting the plot as a PNG file. If None, the plot
        is not exported (default=None).

    Examples
    --------
    Assuming `mod` is an instance of `PysMod` with a completed scan:

    >>> mod.Scan1(range1=[0, 1, 2], runUF=1)
    >>> mod.Scan1Plot(plot=['A', 'B'], title='Scan Results', log='x', format='points', filename='scan_plot.png')

    This will plot the variables 'A' and 'B' from the scan results with a custom
    title and x-axis logarithmic scale. The plot is exported to 'scan_plot.png'.

    """
        data = self.scan_res
        labels = None
        if type(plot) == str:
            plot = [plot]
        if len(plot) == 0:
            plot = self.scan_out
        else:
            plot = [at for at in plot if hasattr(self, at)]
        if len(plot) == 0:
            print('No plottable output specified using self.scan_out')
            plot = self.scan_out
        plotidx = [(self.scan_out.index(c) + 1) for c in plot]
        labels = [self.scan_in] + self.scan_out
        plt.plotLines(data, 0, plotidx, titles=labels, formats=[format])
        end = data[-1, 0] + 0.2 * data[-1, 0]
        plt.setRange('x', data[0, 0], end)
        plt.setAxisLabel('x', self.scan_in)
        yl = 'Steady-state variable'
        plt.setAxisLabel('y', yl)
        if log != None:
            plt.setLogScale(log)
        if title is None:
            plt.setGraphTitle('PySCeS Scan1 (' + self.ModelFile + ') ' +
                time.strftime('%a, %d %b %Y %H:%M:%S'))
        else:
            plt.setGraphTitle(title)
        plt.replot()
        if filename != None:
            plt.export(filename, directory=self.ModelOutput, outtype='png')

    def Scan2D(self, p1, p2, output, log=False):
        """
    Generate a 2-dimensional parameter scan using the steady-state solvers.

    This method conducts a two-dimensional parameter scan by varying two parameters
    over specified ranges and evaluating a given steady-state variable or property.
    The results can be used to analyze the impact of parameter variations on system behavior.
    
    Parameters
    ----------
    p1 : list
        A list specifying the first parameter scan settings in the format 
        [parameter1, start, end, points], where 'parameter1' is the name of the 
        parameter, 'start' is the starting value, 'end' is the ending value, and 
        'points' is the number of points for the scan.
    p2 : list
        A list specifying the second parameter scan settings in the format 
        [parameter2, start, end, points], with components similar to `p1`.
    output : str
        The steady-state variable or property to be evaluated during the scan, 
        such as 'J_R1', 'A_ss', or 'ecR1_s1'.
    log : bool, optional
        If True, the scan will employ a logarithmic scale for both axes. 
        Default is False, which uses a linear scale.
        
    Raises
    ------
    RuntimeError
        If `p1[0]` or `p2[0]` are not valid model attributes.

    Examples
    --------
    >>> model = PysMod()
    >>> p1 = ['parameter1', 0.1, 1.0, 50]
    >>> p2 = ['parameter2', 0.2, 2.0, 50]
    >>> output = 'A_ss'
    >>> model.Scan2D(p1, p2, output)
    >>> print(model.__scan2d_results__)

    Here, a 2D parameter scan is conducted by varying 'parameter1' from 0.1 to 1.0
    and 'parameter2' from 0.2 to 2.0, assessing the steady-state variable 'A_ss'.
    The results are stored in `model.__scan2d_results__`.
    """
        for p in [p1, p2]:
            if not hasattr(self, p[0]):
                raise RuntimeError('"{}" is not a valid model attribute'.
                    format(p[0]))
        p1 = list(p1) + [log, False]
        p2 = list(p2) + [log, False]
        sc1 = Scanner(self)
        sc1.quietRun = True
        sc1.addScanParameter(*p1)
        sc1.addScanParameter(*p2)
        sc1.addUserOutput(output)
        sc1.Run()
        self.__scan2d_pars__ = [p1[0], p2[0], output]
        self.__scan2d_results__ = sc1.getResultMatrix()
        del sc1

    def Scan2DPlot(self, title=None, log=None, format='lines', filename=None):
        """
    Plot the results of a 2D scan generated with `Scan2D`.

    This method visualizes the results of a two-dimensional parameter scan on
    the model, plotting the specified steady-state variable or property against
    two parameters that have been varied over defined ranges.

    Parameters
    ----------
    title : str, optional
        The plot title. If `None`, a default title including the model filename
        and the current date and time will be used. Default is `None`.
    log : {'x', 'xy', 'xyz'}, optional
        Determines which axis, if any, should be displayed in logarithmic scale.
        If `None`, a linear axis is assumed. Default is `None`.
    format : str, optional
        Specifies the backend dependent line format, which could either be
        'lines' or 'points'. Default is 'lines'.
    filename : str, optional
        The filename of the PNG file for exporting the plot. If `None`, the
        plot is not exported. Default is `None`.

    Examples
    --------
    Create a 2D scan and plot the results:

    >>> model = PysMod()
    >>> model.Scan2D(['param1', 0, 10, 100], ['param2', 0, 5, 50], 'A_ss')
    >>> model.Scan2DPlot(title='Steady-state Analysis', log='xy', format='points', filename='scan_results.png')

    Notes
    -----
    Ensure that the method `Scan2D` has been called prior to invoking this
    method, as it relies on the results generated from that scan.

    Raises
    ------
    AttributeError
        If `Scan2D` has not been successfully executed before this method.
    """
        plt.splot(self.__scan2d_results__, 0, 1, 2, titles=self.
            __scan2d_pars__, format=format)
        plt.setRange('xyz')
        plt.setAxisLabel('x', self.__scan2d_pars__[0])
        plt.setAxisLabel('y', self.__scan2d_pars__[1])
        plt.setAxisLabel('z', 'Steady-state variable')
        if log != None:
            plt.setLogScale(log)
        if title is None:
            plt.setGraphTitle('PySCeS Scan2D (' + self.ModelFile + ') ' +
                time.strftime('4a, %d %b %Y %H:%M:%S'))
        else:
            plt.setGraphTitle(title)
        plt.replot()
        if filename != None:
            plt.export(filename, directory=self.ModelOutput, outtype='png')

    def SetQuiet(self):
        """
    Disable solver reporting noise.

    The `SetQuiet` method configures the internal settings of the PysMod instance
    to minimize the verbosity of solver-related outputs. This is useful for
    streamlining results and focusing on the essential data produced by the solvers.

    Specifically, it sets:
    - `hybrd_mesg` to 0: Disables hybrid solver messages.
    - `nleq2_mesg` to 0: Disables non-linear equation solver messages.
    - `lsoda_mesg` to 0: Disables LSODA solver messages.
    - `mode_state_mesg` to 0: Disables mode state messages.
    - `scan1_mesg` to 0: Disables scan-related messages.
    - `solver_switch_warning` to False: Disables warnings about solver switching.

    Parameters
    ----------
    None

    Example
    -------
    To reduce the verbosity of solver outputs in a PysMod instance, use:

    >>> model = PysMod()
    >>> model.SetQuiet()
    # Model solver messages will now be minimal during execution.
    """
        self.__settings__['hybrd_mesg'] = 0
        self.__settings__['nleq2_mesg'] = 0
        self.__settings__['lsoda_mesg'] = 0
        self.__settings__['mode_state_mesg'] = 0
        self.__settings__['scan1_mesg'] = 0
        self.__settings__['solver_switch_warning'] = False

    def SetLoud(self):
        """
    SetLoud()

    Enable detailed solver reporting by adjusting various settings to increase verbosity
    during the solver operations. This function sets the messaging for different solvers
    and reporting components to their highest level, providing more detailed feedback.

    The following settings are adjusted:
    
    - `hybrd_mesg`: Set to 1 to enable messaging for the HYBRD solver.
    - `nleq2_mesg`: Set to 1 to enable messaging for the NLEQ2 solver.
    - `lsoda_mesg`: Set to 1 to enable messaging for the LSODA solver.
    - `mode_state_mesg`: Set to 1 to enable state mode messaging.
    - `scan1_mesg`: Set to 1 to enable messaging during Scan1D operations.
    - `solver_switch_warning`: Set to True to enable warnings related to solver switching.

    Parameters
    ----------
    None

    Examples
    --------
    >>> model = PysMod()
    >>> model.SetLoud()
    >>> # Now, run a solver operation and detailed reporting will be enabled.
    
    Notes
    -----
    Use this method when you need to troubleshoot solver behavior or when detailed operational
    information is required for analysis.
    """
        self.__settings__['hybrd_mesg'] = 1
        self.__settings__['nleq2_mesg'] = 1
        self.__settings__['lsoda_mesg'] = 1
        self.__settings__['mode_state_mesg'] = 1
        self.__settings__['scan1_mesg'] = 1
        self.__settings__['solver_switch_warning'] = True

    def clone(self):
        """
    Returns a deep copy of this model object (experimental!).

    This method creates and returns a deep copy of the current `PysMod` object. 
    It is labeled as experimental, indicating that the functionality may not 
    be fully stable and should be used with caution. A deep copy duplicates the 
    object and all objects it references, ensuring that modifications to the 
    copied object do not affect the original one.

    Returns
    -------
    PysMod
        A deep copy of the current `PysMod` instance.

    Notes
    -----
    When using this method, be aware that it is experimental, and unexpected 
    behavior may occur. Extensive testing is recommended before deploying in 
    a production environment.

    Examples
    --------
    >>> original_mod = PysMod()
    >>> cloned_mod = original_mod.clone()
    INFO: Cloning function is currently experimental ... use with care.
    >>> # Now you have a deep copy of `original_mod` in `cloned_mod`.
    >>> # Changes to `cloned_mod` will not affect `original_mod`.
    >>> cloned_mod.some_attribute = 'new value'
    >>> assert original_mod.some_attribute != cloned_mod.some_attribute
    """
        print(
            '\nINFO: Cloning function is currently experimental ... use with care.'
            )
        return copy.deepcopy(self)


class WasteManagement(object):
    """
    WasteManagement

    A class responsible for handling waste management operations and related processes. 
    This class is a foundational component for managing various waste management tasks 
    within a system, relying on adjustable settings and processes for efficiency.

    Attributes
    ----------
    _gc_delay_hack : gplt
        An internal variable used to control the delay in garbage collection processes.
        This is implemented as a hack for optimizing system performance when managing 
        large amounts of data and ensuring smooth operation of waste management tasks.
        
    Notes
    -----
    The `WasteManagement` class may be expanded in the future to include more 
    attributes and methods that enhance and streamline waste-related processes. As it 
    stands, it focuses on essential operations with the capability to adjust system 
    behaviors through its internal variables.
    """
    _gc_delay_hack = gplt


__waste_manager = WasteManagement()
if __name__ == '__main__':
    print('\nTo use PySCeS import it from a Python Shell: \n\timport pysces\n')
