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
TODO: Parameter elasticities wrt the compartments
"""

from pysces.version import __version__

__doc__ = '''
            PyscesModel
            -----------

            This module contains the core PySCeS classes which
            create the model and associated data objects

            '''
import os, copy, gc, time
import math, operator, re
import pprint, pickle, io

try:
    input = raw_input  # Py2 compatibility
except NameError:
    pass
import numpy
import scipy
import scipy.linalg
import scipy.integrate

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

from getpass import getuser

from . import model_dir as MODEL_DIR
from . import output_dir as OUTPUT_DIR
from . import install_dir as INSTALL_DIR
from . import PyscesRandom as random
from .PyscesScan import Scanner
from .core2.InfixParser import MyInfixParser

from . import (
    nleq2,
    nleq2_switch,
    pitcon,
    plt,
    gplt,
    PyscesStoich,
    PyscesParse,
    __SILENT_START__,
    __CHGDIR_ON_START__,
    SED,
)

if __CHGDIR_ON_START__:
    CWD = OUTPUT_DIR
else:
    CWD = os.getcwd()

interface = None

# this is incredibly crude but effectively masks off unsupported random functions
del (
    random.setstate,
    random.getstate,
)  # random.division,
del random.randrange, random.Random, random.choice
del random.sample, random.shuffle, random.jumpahead
del random.SystemRandom, random.WichmannHill, random.triangular
# used by functions random.NV_MAGICCONST, random.SG_MAGICCONST, random.BPF, random.RECIP_BPF


# Scipy version check
if (
    int(scipy.__version__.split('.')[0]) < 1
    and int(scipy.__version__.split('.')[1]) < 6
):
    print(
        '\nINFO: Your version of SciPy ('
        + scipy.version.version
        + ') might be too old\n\tVersion 0.6.x or newer is strongly recommended\n'
    )
else:
    if not __SILENT_START__:
        print(
            'You are using NumPy ({}) with SciPy ({})'.format(
                numpy.__version__, scipy.__version__
            )
        )

_HAVE_PYSUNDIALS = False
_PYSUNDIALS_LOAD_ERROR = ''
try:
    from pysundials import cvode

    if not __SILENT_START__:
        print('PySundials available')
    _HAVE_PYSUNDIALS = True
except Exception as ex:
    _PYSUNDIALS_LOAD_ERROR = '{}'.format(ex)
    _HAVE_PYSUNDIALS = False

# We now welcome a prototype StomPy interface to PySCeS which will provide expanding stochastic simulation support.

# for future fun
_HAVE_VPYTHON = False
_VPYTHON_LOAD_ERROR = ''
__psyco_active__ = 0
'''
try:
    import visual
    _HAVE_VPYTHON = True
    vpyscene = visual.display(x=150, y=50, title='PySCeS '+__version__+' - C Brett G. Olivier, Stellenbosch 2006', width=640, height=480,\
    center=(0,0,0), background=(0,0,0))
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
'''
# machine specs
mach_spec = numpy.MachAr()
# grab a parser
pscParser = PyscesParse.PySCeSParser(debug=0)


def chkpsc(File):
    """
    chkpsc(File)

    Chekc whether the filename "File" has a '.psc' extension and adds one if not.

    Arguments:

    File: filename string

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
    chkmdir()

    Import and grab pysces.model_dir

    Arguments:
    None

    """
    ##  from pysces import model_dir as MODEL_DIR
    pass


class BagOfStuff(object):
    """
    A collection of attributes defined by row and column lists
    used by Response coefficients etc matrix is an array of values
    while row/col are lists of row colummn name strings
    """

    matrix = None
    row = None
    col = None

    def __init__(self, matrix, row, col):
        "Load attributes from arrays"
        assert len(row) == matrix.shape[0], "\nRow dimension mismatch"
        assert len(col) == matrix.shape[1], "\nCol dimension mismatch"
        self.matrix = matrix
        self.row = list(row)
        self.col = list(col)

    def load(self):
        for r in range(self.matrix.shape[0]):
            for c in range(self.matrix.shape[1]):
                setattr(
                    self, str(self.row[r]) + '_' + str(self.col[c]), self.matrix[r, c]
                )

    def get(self, attr1, attr2):
        "Returns a single attribute \"attr1_attr2\" or None"
        try:
            return getattr(self, attr1 + '_' + attr2)
        except Exception as ex:
            print(ex)
            return None

    def select(self, attr, search='a'):
        """Return a dictionary of <attr>_<name>, <name>_<attr> : val or {} if none
        If attr exists as an index for both left and right attr then:
        search='a' : both left and right attributes (default)
        search='l' : left attributes only
        search='r' : right attributes"""
        output = {}

        if attr in self.row and attr in self.col:
            if search in ['a', 'l']:
                row = self.row.index(attr)
                for col in range(self.matrix.shape[1]):
                    output.setdefault(attr + '_' + self.col[col], self.matrix[row, col])
            if search in ['a', 'r']:
                col = self.col.index(attr)
                for row in range(self.matrix.shape[0]):
                    output.setdefault(self.row[row] + '_' + attr, self.matrix[row, col])
        elif attr in self.row:
            row = self.row.index(attr)
            for col in range(self.matrix.shape[1]):
                output.setdefault(attr + '_' + self.col[col], self.matrix[row, col])
        elif attr in self.col:
            col = self.col.index(attr)
            for row in range(self.matrix.shape[0]):
                output.setdefault(self.row[row] + '_' + attr, self.matrix[row, col])
        return output

    def list(self):
        """
        Return all attributes as a attr:val dictionary

        """
        output = {}
        for row in range(self.matrix.shape[0]):
            for col in range(self.matrix.shape[1]):
                output.setdefault(
                    self.row[row] + '_' + self.col[col], self.matrix[row, col]
                )
        return output

    def listAllOrdered(self, order='descending', absolute=True):
        """
        Return an ordered list of (attr, value) tuples

         - *order* [default='descending'] the order to return as: 'descending' or 'ascending'
         - *absolute* [default=True] use the absolute value

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
    New class used to store parameter scan data (uses StateDataObj)
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
        """takes a list/array of input parameter values and the associated ssdata object"""
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
        if not self.HAS_MOD_DATA:
            self.HAS_MOD_DATA = True
            self.mod_data_labels = [attr for attr in args if hasattr(mod, attr)]
        self.mod_data.append(
            tuple([getattr(mod, attr) for attr in self.mod_data_labels])
        )

    def closeScan(self):
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
        """getScanData(\*args) feed this method species/flux/rule/mod labels and it
        will return an array of [parameter(s), sp1, f1, ....]
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
                output = numpy.hstack(
                    (
                        output,
                        self.species.take([self.species_labels.index(roc)], axis=-1),
                    )
                )
            if self.HAS_FLUXES and roc in self.flux_labels:
                lout.append(roc)
                output = numpy.hstack(
                    (output, self.fluxes.take([self.flux_labels.index(roc)], axis=-1))
                )
            if self.HAS_RULES and roc in self.rules_labels:
                lout.append(roc)
                output = numpy.hstack(
                    (output, self.rules.take([self.rules_labels.index(roc)], axis=-1))
                )
            if self.HAS_XDATA and roc in self.xdata_labels:
                lout.append(roc)
                output = numpy.hstack(
                    (output, self.xdata.take([self.xdata_labels.index(roc)], axis=-1))
                )
            if self.HAS_MOD_DATA and roc in self.mod_data_labels:
                lout.append(roc)
                output = numpy.hstack(
                    (
                        output,
                        self.mod_data.take([self.mod_data_labels.index(roc)], axis=-1),
                    )
                )
        if not lbls:
            return output
        else:
            return output, lout


class StateDataObj(object):
    """
    New class used to store steady-state data.
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

    ##  def setLabels(self, species=None, fluxes=None, rules=None):
    ##  """set the species, rate and rule label lists"""
    ##  if species != None:
    ##  self.species_labels = species
    ##  if fluxes != None:
    ##  self.flux_labels = fluxes
    ##  if rules != None:
    ##  self.rules_labels = rules

    def setSpecies(self, species, lbls=None):
        """Set the species array"""
        self.species = species
        self.HAS_SPECIES = True
        if lbls != None:
            self.species_labels = lbls
            for s in range(len(self.species_labels)):
                setattr(self, self.species_labels[s], self.species[s])

    def setFluxes(self, fluxes, lbls=None):
        """set the flux array"""
        self.fluxes = fluxes
        self.HAS_FLUXES = True
        if lbls != None:
            self.flux_labels = lbls
            for f in range(len(self.flux_labels)):
                setattr(self, self.flux_labels[f], self.fluxes[f])

    def setRules(self, rules, lbls=None):
        """Set the results of rate rules"""
        self.rules = rules
        self.HAS_RULES = True
        if lbls != None:
            self.rules_labels = lbls
            for r in range(len(self.rules_labels)):
                setattr(self, self.rules_labels[r], self.rules[r])

    def setXData(self, xdata, lbls=None):
        """Sets extra simulation data"""
        self.xdata = xdata
        self.HAS_XDATA = True
        if lbls != None:
            self.xdata_labels = lbls
            for x in range(len(self.xdata_labels)):
                setattr(self, self.xdata_labels[x], self.xdata[x])

    def getSpecies(self, lbls=False):
        """return species array"""
        output = None
        if self.HAS_SPECIES:
            output = self.species
        if not lbls:
            return output
        else:
            return output, self.species_labels

    def getFluxes(self, lbls=False):
        """return flux array"""
        output = None
        if self.HAS_FLUXES:
            output = self.fluxes
        if not lbls:
            return output
        else:
            return output, self.flux_labels

    def getRules(self, lbls=False):
        """Return rule array"""
        output = None
        if self.HAS_RULES:
            output = self.rules
        if not lbls:
            return output
        else:
            return output, self.rules_labels

    def getXData(self, lbls=False):
        """Return xdata array"""
        output = None
        if self.HAS_XDATA:
            output = self.xdata
        if not lbls:
            return output
        else:
            return output, self.xdata_labels

    def getAllStateData(self, lbls=False):
        """
        Return all available data as species+fluxes+rules
        if lbls=True returns (array,labels) else just array
        """
        labels = []
        output = None
        if self.HAS_SPECIES:
            output = self.species
            labels += self.species_labels
        if self.HAS_FLUXES:
            if output == None:
                output = self.fluxes
            else:
                output = numpy.hstack((output, self.fluxes))
            labels += self.flux_labels
        if self.HAS_RULES:
            if output == None:
                output = self.rules
            else:
                output = numpy.hstack((output, self.rules))
            labels += self.rules_labels
        if self.HAS_XDATA:
            if output == None:
                output = self.xdata
            else:
                output = numpy.hstack((output, self.xdata))
            labels += self.xdata_labels
        if not lbls:
            return output
        else:
            return output, labels

    def getStateData(self, *args, **kwargs):
        """getSimData(\*args) feed this method species/rate labels and it
        will return an array of [time, sp1, r1, ....]
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
                print('I don\'t have an attribute {} ... ignoring.'.format(roc))
        if not lbls:
            return output
        else:
            return numpy.array(output), lout


class IntegrationDataObj(object):
    """
    This class is specifically designed to store the results of a time simulation
    It has methods for setting the Time, Labels, Species and Rate data and
    getting Time, Species and Rate (including time) arrays. However, of more use:

    - getOutput(\*args) feed this method species/rate labels and it will return
      an array of [time, sp1, r1, ....]
    - getDataAtTime(time) the data generated at time point "time".
    - getDataInTimeInterval(time, bounds=None) more intelligent version of the above
      returns an array of all data points where: time-bounds <= time <= time+bounds
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
        """set the species, rate and rule label lists"""
        if species != None:
            self.species_labels = species
        if rates != None:
            self.rate_labels = rates
        if rules != None:
            self.rules_labels = rules

    def setTime(self, time, lbl=None):
        """Set the time vector"""
        self.time = time.reshape(len(time), 1)
        self.HAS_TIME = True
        if lbl != None:
            self.time_label = lbl

    def setSpecies(self, species, lbls=None):
        """Set the species array"""
        self.species = species
        self.HAS_SPECIES = True
        if lbls != None:
            self.species_labels = lbls

    def setRates(self, rates, lbls=None):
        """set the rate array"""
        self.rates = rates
        self.HAS_RATES = True
        if lbls != None:
            self.rate_labels = lbls

    def setRules(self, rules, lbls=None):
        """Set the results of rate rules"""
        self.rules = rules
        self.HAS_RULES = True
        if lbls != None:
            self.rules_labels = lbls

    def setXData(self, xdata, lbls=None):
        """Sets extra simulation data"""
        self.xdata = xdata
        self.HAS_XDATA = True
        if lbls != None:
            self.xdata_labels = lbls

    def getTime(self, lbls=False):
        """return the time vector"""
        output = None
        if self.HAS_TIME:
            output = self.time.reshape(len(self.time),)
        if not lbls:
            return output
        else:
            return output, [self.time_label]

    def getSpecies(self, lbls=False):
        """return time+species array"""
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
        """return time+rate array"""
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
        """Return time+rule array"""
        ##  assert self.rules != None, "\nNo rules"
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
        """Return time+xdata array"""
        ##  assert self.rules != None, "\nNo rules"
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
        """Return all data generated at "time" """
        # TODO add rate rule data
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
         getDataInTimeInterval(time, bounds=None) returns an array of all
         data points where: time-bounds <= time <= time+bounds
         where bound defaults to stepsize
        """
        # TODO add rate rule data
        temp_t = self.time.reshape(len(self.time),)
        if bounds == None:
            bounds = temp_t[1] - temp_t[0]
        c1 = temp_t >= time - bounds
        c2 = temp_t <= time + bounds
        print('Searching ({}:{}:{})'.format(time - bounds, time, time + bounds))

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
        Old alias for getSimData()
        getOutput(\*args) feed this method species/rate labels and it
        will return an array of [time, sp1, r1, ....]
        """
        return self.getSimData(*args, **kwargs)

    def getAllSimData(self, lbls=False):
        """
        Return all available data as time+species+rates+rules
        if lbls=True returns (array,lables) else just array
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
        """getSimData(\*args) feed this method species/rate labels and it
        will return an array of [time, sp1, r1, ....]
        """
        output = self.time
        ##  print argimrgs
        if 'lbls' in kwargs:
            lbls = kwargs['lbls']
        else:
            lbls = False
        lout = [self.time_label]
        for roc in args:
            if self.HAS_SPECIES and roc in self.species_labels:
                lout.append(roc)
                output = numpy.hstack(
                    (
                        output,
                        self.species.take([self.species_labels.index(roc)], axis=-1),
                    )
                )
            if self.HAS_RATES and roc in self.rate_labels:
                lout.append(roc)
                output = numpy.hstack(
                    (output, self.rates.take([self.rate_labels.index(roc)], axis=-1))
                )
            if self.HAS_RULES and roc in self.rules_labels:
                lout.append(roc)
                output = numpy.hstack(
                    (output, self.rules.take([self.rules_labels.index(roc)], axis=-1))
                )
            if self.HAS_XDATA and roc in self.xdata_labels:
                lout.append(roc)
                output = numpy.hstack(
                    (output, self.xdata.take([self.xdata_labels.index(roc)], axis=-1))
                )
        if not lbls:
            return output
        else:
            return output, lout


# this must stay in sync with core2
class NewCoreBase(object):
    """
    Core2 base class, needed here as we use Core2 derived classes
    in PySCes
    """

    name = None
    __DEBUG__ = False

    def getName(self):
        return self.name

    def setName(self, name):
        self.name = name

    def get(self, attr):
        """Return an attribute whose name is str(attr)"""
        return self.__getattribute__(attr)


# this must stay in sync with core2
class NumberBase(NewCoreBase):
    """
    Derived Core2 number class.
    """

    value = None
    value_initial = None

    def __call__(self):
        return self.value

    def getValue(self):
        return self.value

    def setValue(self, v):
        self.value = v


# Finally killed my lambda functions - brett07
class ReactionObj(NewCoreBase):
    """
    Defines a reaction with a KineticLaw *kl8, *formula* and *name* bound
    to a model instance, *mod*.
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
        mod : model instance
        name = reaction name
        kl = kinetic law
        klrepl = replacement string for KineticLaw
        """
        self.mod = mod
        self.name = name
        self.setKineticLaw(kl, klrepl='self.')

    def __call__(self, *args):
        exec(self.xcode)
        return self.rate

    def setKineticLaw(self, kl, klrepl='self.'):
        InfixParser.setNameStr('self.mod.', '')
        InfixParser.parse(kl.replace(klrepl, ''))
        self.symbols = InfixParser.names
        self.piecewises = InfixParser.piecewises
        formula = InfixParser.output
        # this code has to stay together #
        for pw in InfixParser.piecewises:
            formula = formula.replace(
                'self.mod.{}'.format(pw), 'self.mod.{}()'.format(pw)
            )
        # this code has to stay together #
        self.code_string = 'self.rate={}'.format(formula)
        self.xcode = compile(self.code_string, 'Req: {}'.format(self.name), 'exec')
        self.formula = (
            kl.replace(klrepl, '')
            .replace('numpy.', '')
            .replace('math.', '')
            .replace('operator.', '')
        )


InfixParser = MyInfixParser()
InfixParser.buildlexer()
InfixParser.buildparser(
    debug=0, debugfile='infix.dbg', tabmodule='infix_tabmodule', outputdir=OUTPUT_DIR
)
InfixParser.setNameStr('self.', '')
os.chdir(CWD)

# adapted from core2
class Function(NewCoreBase):
    """
    Function class ported from Core2 to enable the use of functions in PySCeS.
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
        self.name = name
        self.argsl = []
        self.functions = []
        self.mod = mod

    def __call__(self, *args):
        for ar in range(len(args)):
            self.__setattr__(self.argsl[ar], args[ar])
        exec(self.xcode)
        return self.value

    def setArg(self, var, value=None):
        self.__setattr__(var, value)
        self.argsl.append(var)

    def addFormula(self, formula):
        self.formula = formula
        InfixParser.setNameStr('self.', '')
        InfixParser.SymbolReplacements = {'_TIME_': 'mod._TIME_'}
        InfixParser.parse(formula)
        self.piecewises = InfixParser.piecewises
        self.symbols = InfixParser.names
        self.functions = InfixParser.functions
        self.code_string = 'self.value={}'.format(InfixParser.output)
        self.xcode = compile(self.code_string, 'Func: {}'.format(self.name), 'exec')


# adapted from core2
class EventAssignment(NumberBase):
    """
    Event assignments are actions that are triggered by an event.
    Ported from Core2 to build an event handling framework fro PySCeS
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
        setattr(self.mod, self.variable, self.value)
        if self.__DEBUG__:
            print('\tAssigning {} = {}'.format(self.variable, self.value))
        return True

    def __init__(self, name, mod):
        self.setName(name)
        self.mod = mod

    def setVariable(self, var):
        self.variable = var

    def setFormula(self, formula):
        self.formula = formula
        InfixParser.setNameStr('self.mod.', '')
        InfixParser.SymbolReplacements = {'_TIME_': 'self.mod._TIME_'}
        InfixParser.parse(formula)
        self.piecewises = InfixParser.piecewises
        self.symbols = InfixParser.names
        self.code_string = 'self.value={}'.format(InfixParser.output)
        self.xcode = compile(self.code_string, 'EvAs: {}'.format(self.name), 'exec')
        ##  print '\t', self.name, self.code_string

    def evaluateAssignment(self):
        exec(self.xcode)


# adapted from core2
class Event(NewCoreBase):
    """
    Event's have triggers and fire EventAssignments when required.
    Ported from Core2.
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
    symbols = None
    _time_symbol = None
    piecewises = None
    mod = None
    __DEBUG__ = True

    def __init__(self, name, mod):
        self.setName(name)
        self.assignments = []
        self.mod = mod

    def __call__(self, time):
        self._TIME_ = time
        exec(self.xcode)
        ret = False
        if self.state0 and not self.state:
            self.state0 = self.state
        if not self.state0 and self.state:
            for ass in self.assignments:
                ass.evaluateAssignment()
            self.state0 = self.state
            self._need_action = True
            self._ASS_TIME_ = time + self.delay
            if self.__DEBUG__:
                print('\nevent {} is evaluating at {}'.format(self.name, time))
            ret = True
        if self._need_action and self._TIME_ >= self._ASS_TIME_:
            for ass in self.assignments:
                ass()
            if self.__DEBUG__:
                print(
                    'event {} is assigning at {} (delay={})'.format(
                        self.name, time, self.delay
                    )
                )
            self._need_action = False
            ret = True
        return ret

    def setTrigger(self, formula, delay=0.0):
        self.formula = formula
        self.delay = delay
        InfixParser.setNameStr('self.mod.', '')
        ##  print formula, delay
        if self._time_symbol != None:
            InfixParser.SymbolReplacements = {self._time_symbol: '_TIME_'}
        else:
            InfixParser.SymbolReplacements = {'_TIME_': '_TIME_'}
        InfixParser.parse(formula)
        self.piecewises = InfixParser.piecewises
        self.symbols = InfixParser.names
        self.code_string = 'self.state={}'.format(InfixParser.output)
        self.xcode = compile(self.code_string, 'Ev: {}'.format(self.name), 'exec')
        ##  print self.name, self.code_string

    def setAssignment(self, var, formula):
        ass = EventAssignment(var, mod=self.mod)
        ass.setVariable(var)
        ass.setFormula(formula)
        self.assignments.append(ass)
        self.__setattr__('_' + var, ass)

    def reset(self):
        self.state0 = False
        self.state = False
        self._TIME_ = 0.0
        self._ASS_TIME_ = 0.0


class PieceWise(NewCoreBase):
    """
    Generic piecewise class adapted from Core2 that generates a compiled
    Python code block that allows evaluation of arbitrary length piecewise
    functions. Piecewise statements should be defined in assignment rules
    as `piecewise(<Piece>, <Conditional>, <OtherValue>)` where there can
    be an arbitrary number of `<Piece>, <Conditional>` pairs.

    - *args* a dictionary of piecewise information generated by the InfixParser as InfixParser.piecewises

    """

    name = None
    value = None
    formula = None
    code_string = None
    xcode = None
    _names = None
    _TIME_ = None

    def __init__(self, pwd, mod):
        pwd = pwd.copy()
        self.mod = mod
        if pwd['other'] != None:
            other = 'self.value = {}'.format(
                pwd.pop('other').replace('self.', 'self.mod.')
            )
        else:
            other = 'pass'
            pwd.pop('other')

        InfixParser.setNameStr('self.mod.', '')
        self._names = []
        if len(list(pwd.keys())) == 1:
            formula = pwd[0][0]
            InfixParser.parse(formula)
            for n in InfixParser.names:
                if n not in self._names and n != '_TIME_':
                    self._names.append(n)
            formula = InfixParser.output
            thenStat = pwd[0][1].replace('self.', 'self.mod.')
            ##  thenStat = pwd[0][1]
            ##  InfixParser.setNameStr('self.mod.', '')
            ##  InfixParser.parse(thenStat)
            ##  thenStat = InfixParser.output
            self.code_string = 'if {}:\n    self.value = {}\nelse:\n    {}'.format(
                formula, thenStat, other,
            )
            self.formula = self.code_string.replace('self.', '')
        else:
            formula = pwd[0][0]
            InfixParser.parse(formula)
            for n in InfixParser.names:
                if n not in self._names and n != '_TIME_':
                    self._names.append(n)
            formula = InfixParser.output
            thenStat = pwd[0][1].replace('self.', 'self.mod.')
            ##  thenStat = pwd[0][1]
            ##  InfixParser.setNameStr('self.mod.', '')
            ##  InfixParser.parse(thenStat)
            ##  thenStat = InfixParser.output
            self.code_string = 'if {}:\n    self.value = {}\n'.format(formula, thenStat)
            pwd.pop(0)
            for p in pwd:
                formula = pwd[p][0]
                InfixParser.parse(formula)
                for n in InfixParser.names:
                    if n not in self._names and n != '_TIME_':
                        self._names.append(n)

                formula = InfixParser.output
                thenStat = pwd[p][1].replace('self.', 'self.mod.')
                ##  thenStat = pwd[p][1]
                ##  InfixParser.setNameStr('self.mod.', '')
                ##  InfixParser.parse(thenStat)
                ##  thenStat = InfixParser.output
                self.code_string += 'elif {}:\n    self.value = {}\n'.format(
                    formula, thenStat,
                )
            self.code_string += 'else:\n    .format'.format(other)
            self.formula = self.code_string.replace('self.', '')
        self.xcode = compile(self.code_string, 'PieceWise', 'exec')

    def __call__(self):
        exec(self.xcode)
        return self.value


class PysMod(object):
    """
    Create a model object and instantiate a PySCeS model so that it can be used for further analyses. PySCeS
    model descriptions can be loaded from files or strings (see the *loader* argument for details).

    - *File* the name of the PySCeS input file if not explicit a \*.psc extension is assumed.
    - *dir* if specified, the path to the input file otherwise the default PyscesModel directory (defined in the pys_config.ini file) is assumed.
    - *autoload* autoload the model, pre 0.7.1 call mod.doLoad(). (default=True) **new**
    - *loader* the default behaviour is to load PSC file, however, if this argument is set to 'string' an input file can be supplied as the *fString* argument (default='file')
    - *fString* a string containing a PySCeS model file (use with *loader='string'*) the *File* argument now sepcifies the new input file name.

    """

    __version__ = __version__
    __pysces_directory__ = INSTALL_DIR
    __settings__ = None
    random = random
    # __STOMPY__ = None

    def __init__(self, File=None, dir=None, loader='file', fString=None, autoload=True):
        """
        Create a model object and instantiate a PySCeS model so that it can be used for further analyses. PySCeS
        model descriptions can be loaded from files or strings (see the *loader* argument for details).

        - *File* the name of the PySCeS input file if not explicit a \*.psc extension is assumed.
        - *dir* if specified, the path to the input file otherwise the default PyscesModel directory (defined in the pys_config.ini file) is assumed.
        - *autoload* autoload the model, pre 0.7.1 call mod.doLoad(). (default=True) **new**
        - *loader* the default behaviour is to load PSC file, however, if this argument is set to 'string' an input file can be supplied as the *fString* argument (default='file')
        - *fString* a string containing a PySCeS model file (use with *loader='string'*) the *File* argument now specifies the new input file name.

        """

        # if _HAVE_STOMPY:
        # self.__STOMPY__ = StomPyInterface(MODEL_DIR, OUTPUT_DIR)
        # print 'PySCeS/StochPy interface active'
        # else:
        # self.__STOMPY__ = None

        self.__settings__ = {}
        self.__settings__.update({'enable_deprecated_attr': True})
        self.WorkDir = CWD
        if loader == 'file':
            self.LoadFromFile(File, dir)
        elif loader == 'string':
            self.LoadFromString(File, fString)
        else:
            self.LoadFromFile(File, dir)
        # stuff that needs to be done before initmodel
        self.__settings__['mode_substitute_assignment_rules'] = False
        self.__settings__['display_compartment_warnings'] = False
        self._TIME_ = 0.0  # this will be the built-in time
        self.piecewise_functions = []
        self.__piecewises__ = {}
        self.__HAS_PIECEWISE__ = False
        if autoload:
            self.ModelLoad()
            self.__PSC_auto_load = True
        else:
            self.__PSC_auto_load = False

    def ModelLoad(self, stoich_load=0):
        """
        Load and instantiate a PySCeS model so that it can be used for further analyses. This function
        replaces the pre-0.7.1 doLoad() method.

        - *stoich_load* try to load a structural analysis saved with Stoichiometry_Save_Serial() (default=0)

        """
        self.InitialiseInputFile()
        assert self.__parseOK, '\nError in input file, parsing could not complete'
        self.Stoichiometry_Analyse(override=0, load=stoich_load)
        # add new Style functions to model
        self.InitialiseFunctions()
        self.InitialiseCompartments()
        self.InitialiseRules()
        self.InitialiseEvents()
        self.InitialiseModel()
        self.InitialiseRuleChecks()
        self.InitialiseOldFunctions()  # TODO replace this with initialisation functions

    def doLoad(self, stoich_load=0):
        """
        Load and instantiate a PySCeS model so that it can be used for further analyses. This function is
        being replaced by the ModelLoad() method.

        - *stoich_load* try to load a structural analysis saved with Stoichiometry_Save_Serial() (default=0)

        """
        if not self.__PSC_auto_load:
            self.ModelLoad(stoich_load=stoich_load)
        else:
            print(
                'PySCeS now automatically loads the model on model object instantiation. If you do not want this behaviour pass the autoload=False argument to the constructor, if you really want to reload the model, run reLoad().'
            )

    def reLoad(self, stoich_load=0):
        """
        Re-load and instantiate a PySCeS model so that it can be used for further analyses. This is just a convenience call to the ModelLoad() method.

        - *stoich_load* try to load a structural analysis saved with Stoichiometry_Save_Serial() (default=0)

        """
        self.ModelLoad(stoich_load=stoich_load)

    def LoadFromString(self, File=None, fString=None):
        """
        Docstring required
        """

        # grab model directory
        chkmdir()
        mdir = MODEL_DIR

        # check for .psc extension
        File = chkpsc(File)

        print('Using model directory: ' + mdir)

        if not os.path.isdir(os.path.join(MODEL_DIR, "orca")):
            os.mkdir(os.path.join(MODEL_DIR, "orca"))

        mdir = os.path.join(MODEL_DIR, "orca")

        # write string to file
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
            print(
                os.path.join(mdir, File)
                + ' does not exist please re-instantiate model .....',
                end=' ',
            )
        self.__settings__['display_debug'] = 0
        self.ModelOutput = OUTPUT_DIR

        # Initialize serial directory
        self.__settings__['serial_dir'] = os.path.join(self.ModelOutput, 'pscdat')
        # Initialize stoichiometric precision
        self.__settings__['stoichiometric_analysis_fp_zero'] = mach_spec.eps * 2.0e4
        self.__settings__['stoichiometric_analysis_lu_precision'] = self.__settings__[
            'stoichiometric_analysis_fp_zero'
        ]
        self.__settings__['stoichiometric_analysis_gj_precision'] = (
            self.__settings__['stoichiometric_analysis_lu_precision'] * 10.0
        )

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
        __init__(File=None,dir=None)

        Initialise a PySCeS model object with PSC file that can be found in optional directory.
        If a a filename is not supplied the pysces.model_dir directory contents is displayed and
        the model name can be entered at the promp (<ctrl>+C exits the loading process).

        Arguments:

        File [default=None]: the name of the PySCeS input file
        dir [default=pysces.model_dir]: the optional directory where the PSC file can be found

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
        if File == None:
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
                    'No models available in model directory, please set using pys_usercfg.ini or call\
                with pysces.model(\'file.psc\',dir=\'path\\to\\models\')'
                )
                dir = input('\nPlease enter full path to models <CTRL+C> exits: ')

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
            print(
                os.path.join(dir, File)
                + ' does not exist please re-instantiate model .....',
                end=' ',
            )
        self.__settings__['display_debug'] = 0
        self.ModelOutput = OUTPUT_DIR

        ##  # Initialise elementary modes
        ##  if os.sys.platform == 'win32':
        ##  self.eModeExe_int = os.path.join(self.__metatool,'meta43_int.exe')
        ##  self.eModeExe_dbl = os.path.join(self.__metatool,'meta43_double.exe')
        ##  else:
        ##  self.eModeExe_int = os.path.join(self.__metatool,'meta43_int')
        ##  self.eModeExe_dbl = os.path.join(self.__metatool,'meta43_double')
        ##  print 'Done.'

        # Initialize serial directory
        self.__settings__['serial_dir'] = os.path.join(self.ModelOutput, 'pscdat')
        # Initialize stoichiometric precision
        self.__settings__['stoichiometric_analysis_fp_zero'] = mach_spec.eps * 2.0e4
        self.__settings__['stoichiometric_analysis_lu_precision'] = self.__settings__[
            'stoichiometric_analysis_fp_zero'
        ]
        self.__settings__['stoichiometric_analysis_gj_precision'] = (
            self.__settings__['stoichiometric_analysis_lu_precision'] * 10.0
        )

    # def __LoadStomPyInterface__(self):
    # """
    # Load the StomPy Stochastic simulation interface
    # """
    # if _HAVE_STOMPY and self.__STOMPY__ == None:
    # self.__STOMPY__ = StomPyInterface(MODEL_DIR, OUTPUT_DIR)
    # print 'PySCeS/StomPy interface active'

    def __ParsePiecewiseFunctions__(self, piecewises):
        """
        THIS IS HIGHLY EXPERIMENTAL! Takes the a piecewises dictionary
        and creates a piecewise object while substituting the object name
        in the formula string.

        - *piecewises* a piecewise dictionary created by the InfixParser

        """
        if piecewises != None and len(list(piecewises.keys())) > 0:
            self.__HAS_PIECEWISE__ = True
            for p in piecewises:
                print('Info: adding piecewise object: {}'.format(p))
                # if len(piecewises[p].keys()) == 2:
                # piecewises[p][0].reverse() # only for libsbml generated infix
                self.__piecewises__.update({p: piecewises[p]})
                P = PieceWise(piecewises[p], self)
                P.setName(p)
                self.piecewise_functions.append(P)
                setattr(self, p, P)

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

        print('\nParsing file: {}'.format(os.path.join(self.ModelDir, self.ModelFile)))

        pscParser.ParsePSC(self.ModelFile, self.ModelDir, self.WorkDir)
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
            InfixParser.__pwcntr__ = 0
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
            if self.__KeyWords__['ModelType'] == None:
                self.__KeyWords__['ModelType'] = ['Deterministic']
            else:
                self.__KeyWords__['ModelType'] = tuple(
                    [t.strip() for t in self.__KeyWords__['ModelType'].split(',')]
                )

            # we now check for modeltype, if it specified as stochastic check if stochpy is available
            # if self.__KeyWords__['ModelType'] == ['Stochastic']:
            # if _HAVE_STOMPY:
            # print 'INFO: This model suggests that it requires discrete simulation and StochPy is installed ... mod.doStochSim() and mod.doStochSimPlot() are available for use.'
            # else:
            # print 'INFO: This model suggests that it requires stochastic simulation and StochPy (stompy.sf.net) is not installed ... PySCeS will treat this model as continuous.'

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
            self.__cbm_fluxbounds__ = pscParser.cbm_FluxBounds.copy()
            self.__cbm_objfuncs__ = pscParser.cbm_ObjectiveFunctions.copy()
            self.__cbm_userfluxconstraints__ = pscParser.cbm_UserFluxConstraints.copy()
            ##  if pscParser.ModelUsesNumpyFuncs:
            ##  print 'Numpy functions detected in kinetic laws.\n'
        else:
            print('\nERROR: model parsing error, please check input file.\n')
        # added in a check for model correctness and human error reporting (1=ok, 0=error)
        if len(pscParser.SymbolErrors) != 0:
            print('\nUndefined symbols:\n{}'.format(pscParser.SymbolErrors))
        if not pscParser.ParseOK:
            print(
                '\n\n*****\nModel parsing errors detected in input file '
                + self.ModelFile
                + '\n*****'
            )
            print('\nInput file errors')
            for error in pscParser.LexErrors:
                print(error)
                ##  try:
                ##  print error[0] + 'in line:\t' + str(error[1]) + ' ('+ error[2][:20] +' ...)'
                ##  except IndexError:
                ##  print 'Illegal character:', error.__repr__(), error.__str__()
            print('\nParser errors')
            for error in pscParser.ParseErrors:
                print(error)
                ##  try:
                ##  print error[0] + '- ' + error[2][:20]
                ##  except IndexError:
                ##  print 'Illegal character:', error
            assert pscParser.ParseOK == 1, 'Input File Error'

    def InitialiseRules(self):
        # we need to detect different types of rules etc
        # defmod
        self.__HAS_FORCED_FUNCS__ = False
        self.__HAS_RATE_RULES__ = False
        rate_rules = {}
        assignment_rules = {}
        for ar in self.__rules__:
            if self.__rules__[ar]['type'] == 'assignment':
                self.__HAS_FORCED_FUNCS__ = True
                assignment_rules.update({ar: self.__rules__[ar]})
            elif self.__rules__[ar]['type'] == 'rate':
                self.__HAS_RATE_RULES__ = True
                rate_rules.update({ar: self.__rules__[ar]})

        # THE NEW WAY (adds formula parsing for numpy functions etc)
        InfixParser.setNameStr('self.', '')
        self._NewRuleXCode = {}
        if (
            self.__HAS_FORCED_FUNCS__
            and not self.__settings__['mode_substitute_assignment_rules']
        ):
            code_string = ''
            all_names = []
            for ar in assignment_rules:
                name = assignment_rules[ar]['name']
                InfixParser.setNameStr('self.', '')
                InfixParser.parse(assignment_rules[ar]['formula'])
                assignment_rules[ar]['symbols'] = InfixParser.names
                formula = InfixParser.output
                # this code has to stay together #
                for pw in InfixParser.piecewises:
                    formula = formula.replace(
                        'self.{}'.format(pw), 'self.{}()'.format(pw)
                    )
                self.__ParsePiecewiseFunctions__(InfixParser.piecewises)
                # this code has to stay together #
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
            for ar in indep + keep:
                evalCode = 'self.{} = {}\n'.format(
                    assignment_rules[ar]['name'], assignment_rules[ar]['code_string'],
                )
                self._NewRuleXCode.update({assignment_rules[ar]['name']: evalCode})
                code_string += evalCode

            print('\nAssignment rule(s) detected.')
            self._Function_forced = code_string
        elif (
            self.__HAS_FORCED_FUNCS__
            and self.__settings__['mode_substitute_assignment_rules']
        ):
            # here we substitute nested assignment rules
            InfixParser.setNameStr('self.', '')
            symbR = {}
            for ass in assignment_rules:
                InfixParser.setNameStr('self.', '')
                InfixParser.parse(assignment_rules[ass]['formula'])
                formula = InfixParser.output
                # this code has to stay together #
                for pw in InfixParser.piecewises:
                    formula = formula.replace(
                        'self.{}'.format(pw), 'self.{}()'.format(pw)
                    )
                self.__ParsePiecewiseFunctions__(InfixParser.piecewises)
                # this code has to stay together #
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

        self.__CODE_forced = compile(self._Function_forced, 'AssignRules', 'exec')
        # tested in InitialiseRuleChecks()

        # defmod
        self.__rate_rules__ = []
        self.__rrule__ = None
        rr_code_block = ''
        rr_map_block = ''
        self.__CODE_raterule = None
        self._NewRateRuleXCode = {}
        if self.__HAS_RATE_RULES__:
            # create a rr vector
            self.__rrule__ = numpy.ones(len(list(rate_rules.keys())), 'd')
            # brett2008 debug stuff doesn't do any harm
            cntr = 0
            rr_keys = list(rate_rules.keys())
            rr_keys.sort()
            for ar in rr_keys:
                name = rate_rules[ar]['name']
                InfixParser.setNameStr('self.', '')
                InfixParser.parse(rate_rules[ar]['formula'])
                formula = InfixParser.output
                rate_rules[ar]['symbols'] = InfixParser.names
                # this code has to stay together #
                for pw in InfixParser.piecewises:
                    formula = formula.replace(
                        'self.{}'.format(pw), 'self.{}()'.format(pw)
                    )
                self.__ParsePiecewiseFunctions__(InfixParser.piecewises)
                # this code has to stay together #
                rate_rules[ar]['code_string'] = formula
                self.__rate_rules__.append(name)
                rr_code_block += 'self.__rrule__[{}] = {}\n'.format(cntr, formula)
                ##  rr_code_block += 'self.%s = self.__rrule__[%s]\n' % (name, cntr)
                rr_map_block += 'self.{} = self.__rrule__[{}]\n'.format(name, cntr)
                cntr += 1
                # create mod.<rule name>_init attributes
                setattr(self, '{}_init'.format(name), getattr(self, name))
            print('Rate rule(s) detected.')
        else:
            rr_code_block = 'pass\n'
            rr_map_block = 'pass\n'
        # TODO consider putting this in self.__HAS_RATE_RULES__
        self.__CODE_raterule = compile(rr_code_block, 'RateRules', 'exec')
        self.__CODE_raterule_map = compile(rr_map_block, 'RateRuleMap', 'exec')

        del rate_rules, assignment_rules

    def InitialiseRuleChecks(self):
        try:
            exec(self.__CODE_raterule)
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
                    exec(compile(self._NewRuleXCode[kk], 'NewRuleXCode', 'exec'))
                    if self.__HAS_RATE_RULES__:
                        exec(self.__CODE_raterule)
                        exec(self.__CODE_raterule_map)
                except Exception as ex:
                    zeroDivErr2.append(kk)
            exit -= 1
            zeroDivErr = zeroDivErr2
            if exit < 0:
                print('WARNING: ZeroDivision elimination failed')

        try:
            exec(self.__CODE_forced)
            # print 'done.'
        except ZeroDivisionError:
            print('WARNING: Assignment rule ZeroDivision on intialisation')
        except Exception as ex:
            print('WARNING: Assignment rule error (disabling all rules)\n', ex)
            self.__CODE_forced = compile('pass\n', 'AssignRules', 'exec')

    def Stoichiometry_Init(self, nmatrix, load=0):
        """
        Stoichiometry_Init(nmatrix,load=0)

        Initialize the model stoichiometry. Given a stoichiometric matrix N, this method will return an instantiated PyscesStoich instance and status flag. Alternatively, if load is enabled, PySCeS will
        attempt to load a previously saved stoichiometric analysis (saved with Stoichiometry_Save_Serial)
        and test it's correctness. The status flag indicates 0 = reanalyse stoichiometry or
        1 = complete structural analysis preloaded.

        Arguments:

        nmatrix: The input stoichiometric matrix, N
        load [default=0]: try to load a saved stoichiometry (1)

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
                        badList.append((stc.__species__[x], self.__species__[x]))
                        go = 0
                    for y in range(col):
                        if stc.__reactions__[y] != self.__reactions__[y]:
                            badList.append(
                                (stc.__reactions__[y], self.__reactions__[y])
                            )
                            go = 0
                        if (
                            abs(nmatrix[x, y] - stc.nmatrix[x, y])
                            > stc.stoichiometric_analysis_fp_zero
                        ):
                            badList.append(((x, y), nmatrix[x, y], stc.nmatrix[x, y]))
                            go = 0
            if not go:
                ##  print 'Stoichiometry mismatch'
                ##  for x in badList:
                ##  print x,
                print('\nProblem loading stoichiometry, reanalysing ...')
                stc = PyscesStoich.Stoich(nmatrix)
                status = 0
            else:
                print('Stoichiometry verified ... we have liftoff')
                status = 1
        else:
            # print 'Instantiating new stoichiometry ...'
            stc = PyscesStoich.Stoich(nmatrix)
            status = 0
        return stc, status

    def Stoichiometry_Save_Serial(self):
        """Serialize and save a Stoichiometric instance to binary pickle""" """
        Stoichiometry_Save_Serial()

        Serilaise and save the current model stoichiometry to a file with name <model>_stoichiometry.pscdat
        in the mod.__settings__['serial_dir'] directory (default: mod.model_output/pscdat)

        Arguments:
        None

        """
        # new plan, I introduce species and reaction arrays into the stoich instance for verification purposes
        # brett - 20050831
        self.__structural__.__species__ = copy.copy(self.__species__)
        self.__structural__.__reactions__ = copy.copy(self.__reactions__)
        self.SerialEncode(self.STOICH, self.ModelFile[:-4] + '_stoichiometry')

    def Stoichiometry_Load_Serial(self):
        """
        Stoichiometry_Load_Serial()

        Load a saved stoichiometry saved with mod.Stoichiometry_Save_Serial() and return
        a stoichiometry instance.

        Arguments:
        None

        """
        stc = self.SerialDecode(self.ModelFile[:-4] + '_stoichiometry')
        return stc

    def Stoichiometry_Analyse(self, override=0, load=0):
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
            self.__nmatrix__ = self.__initmodel__()  # Creates the model N
            # print '\nintializing N\n'
        else:
            print('\nStoichiometric override active\n')

        assert (
            len(self.__nmatrix__) > 0
        ), '\nUnable to generate Stoichiometric Matrix! model has:\n{} reactions\n{} species\nwhat did you have in mind?\n'.format(
            len(self.__reactions__), len(self.__species__)
        )

        self.__Nshape__ = self.__nmatrix__.shape  # Get the shape of N
        ##  self.__Vtemp__ = numpy.zeros((self.__Nshape__[1])) # going going ....

        # get stoich instance and whether it was analysed or loaded - brett 20050830
        self.__structural__, stc_load = self.Stoichiometry_Init(
            self.__nmatrix__, load=load
        )

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
            ##  self.__nmatrix__ = copy.copy(self.nmatrix)
            self.nmatrix = self.__nmatrix__  # done with caution brett2008
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

    def Stoichiometry_ReAnalyse(self):
        """
        Stoichiometry_ReAnalyse()

        Reanalyse the stoichiometry using the current N matrix ie override=1
        (for use with mod.Stoich_matrix_SetValue)

        Arguments:
        None

        """
        self.Stoichiometry_Analyse(override=1)
        self.InitialiseConservationArrays()

    def Stoich_nmatrix_SetValue(self, species, reaction, value):
        """
        Stoich_nmatrix_SetValue(species,reaction,value)

        Change a stoichiometric coefficient's value in the N matrix. Only a coefficients magnitude may be set, in other words a
        a coefficient's value must remain negative, positive or zero. After changing a coefficient
        it is necessary to Reanalyse the stoichiometry.

        Arguments:

        species: species name (s0)
        reaction: reaction name (R4)
        value: new coefficient value

        """
        index = self.__Stoich_nmatrix_FindIndex__(species, reaction)
        if self.__Stoich_nmatrix_CheckValue__(index, value):
            self.__Stoich_nmatrix_UpdateValue__(index, value)
            self.__StoichOK = 0

    def __Stoich_nmatrix_FindIndex__(self, species, reaction):
        """
        __Stoich_nmatrix_FindIndex__(species,reaction)

        Return the N matrix co-ordinates of the coefficient referenced by species and reaction name.

        Arguments:

        species: species name
        reaction: reaction name

        """
        rval = None
        cval = None
        for x in enumerate(self.species):
            if species == x[1]:
                rval = x[0]
        for y in enumerate(self.reactions):
            if reaction == y[1]:
                cval = y[0]
        return (rval, cval)

    def __Stoich_nmatrix_UpdateValue__(self, xxx_todo_changeme, val):
        """
        __Stoich_nmatrix_UpdateValue__((x,y), val)

        Update N matrix co-ordinates (x,y) with val

        Arguments:
        (x,y): row, column coordinates
        val: value

        """
        (x, y) = xxx_todo_changeme
        self.nmatrix[x, y] = val
        self.__nmatrix__[x, y] = val
        self.Nmatrix.array[x, y] = val

    def __Stoich_nmatrix_CheckValue__(self, xxx_todo_changeme1, val):
        """
        __Stoich_nmatrix_CheckValue__((x,y),val)

        Check validity of new coefficient value against existing N matrix coefficient (x,y) for existance and/or sign.
        Returns 1 (true) or 0.

        Arguments:
        (x,y): N matrix coordinates
        val: new value

        """
        (x, y) = xxx_todo_changeme1
        go = 1
        if x == None or y == None:
            print('\nSpecies (nmatrix) index  =', x)
            print('Reaction (nmatrix) index =', y)
            print(
                '\nI\'m confused, perhaps you entered an incorrect species or reaction name'
            )
            print('or they are in the wrong order?')
            go = 0
        elif abs(self.__nmatrix__[x, y]) == 0.0:
            if val != 0.0:
                go = 0
                print('\nZero coefficient violation')
                print(
                    '  nmatrix['
                    + repr(x)
                    + ','
                    + repr(y)
                    + '] can only be = 0.0 (input '
                    + str(val)
                    + ')'
                )
        elif self.__nmatrix__[x, y] > 0.0:
            if val <= 0.0:
                go = 0
                print('\nPositive coefficient violation')
                print(
                    '  nmatrix['
                    + repr(x)
                    + ','
                    + repr(y)
                    + '] can only be > 0.0 (input '
                    + str(val)
                    + ')'
                )
        elif self.__nmatrix__[x, y] < 0.0:
            if val >= 0.0:
                go = 0
                print('\nNegative coefficient violation')
                print(
                    '  nmatrix['
                    + repr(x)
                    + ','
                    + repr(y)
                    + '] can only be < 0.0 (input '
                    + str(val)
                    + ')'
                )
        if go:
            return 1
        else:
            return 0

    def SerialEncode(self, data, filename):
        """
        SerialEncode(data,filename)

        Serialise and save a Python object using a binary pickle to file. The serialised object
        is saved as <filename>.pscdat in the directory defined by mod.model_serial.

        Arguments:

        data: pickleable Python object
        filename: the ouput filename

        """
        filename = str(filename) + '.pscdat'
        if os.path.exists(self.__settings__['serial_dir']):
            pass
        else:
            os.mkdir(self.__settings__['serial_dir'])
        print('\ncPickle data stored in: ' + self.__settings__['serial_dir'])

        File = open(os.path.join(self.__settings__['serial_dir'], filename), 'wb')
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
        SerialDecode(filename)

        Decode and return a serialised object saved with SerialEncode.

        Arguments:

        filename: the filename (.pscdat is assumed)

        """
        filename = str(filename) + '.pscdat'
        if not os.path.exists(os.path.join(self.__settings__['serial_dir'], filename)):
            raise RuntimeError(
                'Serialized data '
                + os.path.join(self.__settings__['serial_dir'], filename)
                + ' does not exist'
            )
        data = None
        try:
            File = open(os.path.join(self.__settings__['serial_dir'], filename), 'rb')
            data = pickle.load(File)
            File.close()
            print('Serial decoding complete')
        except Exception as E:
            print('Serial decode error:\n', E)

        return data

    def InitialiseCompartments(self):
        # the name should say it all
        self.__settings__['compartment_fudge_factor'] = 1.0
        startFudging = 1.0e-8
        I_AM_FUDGING = False
        if self.__HAS_COMPARTMENTS__:
            tmp = min(
                [
                    abs(self.__compartments__[c]['size'])
                    for c in list(self.__compartments__.keys())
                ]
            )
            if tmp < startFudging:
                ##  self.__settings__['compartment_fudge_factor'] = 10.0**round(numpy.log10(tmp))
                self.__settings__[
                    'compartment_fudge_factor'
                ] = tmp  # sneaky b*))_)$^&*@rds
                ##  print 'compartment_fudge_factor', self.__settings__['compartment_fudge_factor']
            for c in list(self.__compartments__.keys()):
                if self.__settings__['compartment_fudge_factor'] < startFudging:
                    newsize = (
                        self.__compartments__[c]['size']
                        / self.__settings__['compartment_fudge_factor']
                    )
                    print(
                        'INFO: Rescaling compartment with size {} to {}'.format(
                            self.__compartments__[c]['size'], newsize
                        )
                    )
                    self.__compartments__[c]['size'] = newsize
                    self.__compartments__[c].update(
                        {'scale': self.__settings__['compartment_fudge_factor']}
                    )
                    I_AM_FUDGING = True
                setattr(
                    self,
                    self.__compartments__[c]['name'],
                    self.__compartments__[c]['size'],
                )
                setattr(
                    self,
                    '{}_init'.format(self.__compartments__[c]['name']),
                    self.__compartments__[c]['size'],
                )
        if self.__HAS_COMPARTMENTS__:
            for sp in list(self.__sDict__.keys()):
                if (
                    self.__sDict__[sp]['compartment'] == None
                    and len(list(self.__compartments__.keys())) == 1
                ):
                    self.__sDict__[sp]['compartment'] = self.__compartments__[
                        list(self.__compartments__.keys())[0]
                    ]['name']
                    print(
                        'COMPARTMENT WARNING: this model has a compartment defined {} but \"{}\" is not in it ... I\'ll try it for now'.format(
                            [
                                self.__compartments__[c]['name']
                                for c in list(self.__compartments__.keys())
                            ],
                            sp,
                        )
                    )
                    ##  print self.__sDict__[sp]
                elif (
                    self.__sDict__[sp]['compartment'] == None
                    and len(list(self.__compartments__.keys())) > 1
                ):
                    assert (
                        self.__sDict__[sp]['compartment'] != None
                    ), '\nCOMPARTMENT ERROR: this model has multiple compartments defined {} but \"{}\" is not in one!'.format(
                        [
                            self.__compartments__[c]['name']
                            for c in list(self.__compartments__.keys())
                        ],
                        sp,
                    )
                # brett 2008 fudge this!
                if not self.__KeyWords__['Species_In_Conc'] and I_AM_FUDGING:
                    self.__sDict__[sp]['initial'] = (
                        self.__sDict__[sp]['initial']
                        / self.__settings__['compartment_fudge_factor']
                    )
                    setattr(self, sp, self.__sDict__[sp]['initial'])
                    setattr(self, '{}_init'.format(sp), self.__sDict__[sp]['initial'])
                    print(
                        'INFO: Rescaling species ({}) to size {}.'.format(
                            sp, self.__sDict__[sp]['initial']
                        )
                    )

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
        ##  self.__CsizeFSIdx__ = []
        ##  for s in self.__fixed_species__:
        ##  if self.__HAS_COMPARTMENTS__:
        ##  self.__CsizeFSIdx__.append(self.__sDict__[s]['compartment'])
        ##  else:
        ##  self.__CsizeFSIdx__.append('__null_compartment__')

        # Init compartment divisions vector - brett 2008
        ##  self.__CsizeAll__ = numpy.ones(len(self.__CsizeAllIdx__), 'd')
        if self.__HAS_COMPARTMENTS__:
            self.__CsizeAll__ = numpy.array(
                [self.__compartments__[c]['size'] for c in self.__CsizeAllIdx__], 'd'
            )
        else:
            self.__CsizeAll__ = numpy.ones(len(self.__CsizeAllIdx__), 'd')
        ##  self.__CsizeFS__ = numpy.array([self.__compartments__[c]['size'] for c in self.__CsizeFSIdx__], 'd')

    # add new function types to model
    def InitialiseFunctions(self):
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

    # add event types to model
    def InitialiseEvents(self):
        self.__events__ = []
        # for each event
        for e in self.__eDict__:
            ev = Event(e, self)
            ev._time_symbol = self.__eDict__[e]['tsymb']
            ev.setTrigger(self.__eDict__[e]['trigger'], self.__eDict__[e]['delay'])
            # for each assignment
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
        InitialiseModel()

        Initialise and set up dynamic model attributes and methods using the model defined in
        the associated PSC file

        Arguments:
        None

        """
        # TODO killkillkill
        self.__inspec__ = numpy.zeros((len(self.__species__)), 'd')
        self._D_s_Order, DvarUpString = self.__initvar__()  # Initialise s[] array

        self.__CDvarUpString__ = compile(DvarUpString, 'DvarUpString', 'exec')

        self.state_species = None
        self.state_flux = None
        ##  self.sim_res = None
        self.elas_var_u = None
        self.elas_var = None
        self.__vvec__ = numpy.zeros((self.__Nshape__[1]), 'd')
        self.__tvec_a__ = None
        self.__tvec_c__ = None

        self._CVODE_Vtemp = numpy.zeros((self.__Nshape__[1]), 'd')

        #  #debug
        #  self.display_debugcount = 0
        #  self.cntr = 50

        par_map2store, par_remap = self.__initfixed__()  # Initialise fixed species

        # self.__par_map2storeC =  par_map2store
        # self.__par_remapC =  par_remap
        self.__par_map2storeC = compile(par_map2store, 'par_map2store', 'exec')
        self.__par_remapC = compile(par_remap, 'par_remap', 'exec')

        ##  self.__xMake = compile(xString,'xString','exec')
        ##  exec(self.__xMake)

        # initialise rate equations
        mapString, mapString_R, vString = self.__initREQ__()

        self.__mapFunc__ = compile(mapString, 'mapString', 'exec')
        self.__mapFunc_R__ = compile(mapString_R, 'mapString_R', 'exec')
        self.__vFunc__ = compile(vString, 'vString', 'exec')
        ##  exec(lambdaFuncsC)

        # Machine specific values (for IEEE compliant FP machines this is around 2.e-16)
        # FutureNote: investigate arbitrary precision FP in Python, probably GNU-GMP based - brett 20040326
        self.__settings__['mach_floateps'] = mach_spec.eps

        # PySCeS mode switches
        # this affects number output formatting in PySCeS)
        self.__settings__['mode_number_format'] = '%2.4e'
        # 0:initval, 1:zeroval, 2:lastss 3:sim_res[-1]
        self.__settings__['mode_sim_init'] = 0
        # maximum number of auto-stepsize adjustments
        self.__settings__['mode_sim_max_iter'] = 3

        # TODO UPGRADE
        self.mode_state_init = 0  # 0:initval, 1:zeroval, 2:sim, 3:%state
        self.mode_solver = 'HYBRD'  # 0:HYBRD, 1:FINTSLV, 2:NLEQ2
        self.mode_solver_fallback = 1  # 0:Solver failure fails,
        # 1:solver failure falls back to NLEQ2 if present then FINTSLV (see next)
        self.STATE_extra_output = []  # extra data added to data_sstate object
        # factor to multiply the previous ss used as initialiser
        self.__settings__['mode_state_init3_factor'] = 0.1
        # the simulation array used to initialise the steady state
        self.__mode_state_init2_array__ = numpy.logspace(0, 5, 18)
        self.__settings__['mode_state_mesg'] = 1  # happy exit message from State()
        # give last best solutions or NaN as solution
        self.__settings__['mode_state_nan_on_fail'] = False
        # the irritating switching to message
        self.__settings__['solver_switch_warning'] = True
        # 0:no fallback to forward integration 1:fallback to integration
        self.__settings__['mode_solver_fallback_integration'] = 1
        # this is now initialised in the model parsing sectiomn ParseModel
        # the order of ScipyDerivative - 3 seems good for most situations
        self.__settings__['mode_elas_deriv_order'] = 3
        # MCA scaling option 0:unscaled E+C+Eig 1:scaled E+C+Eig
        self.__settings__['mode_mca_scaled'] = 1
        # 0:normal, 1:extended eigen value + left/right vectors
        self.__settings__['mode_eigen_output'] = 0
        # under consideration
        # 0:warnings displayed, 1:warnings suppressed
        self.__settings__['mode_suppress_info'] = 0
        # new compartment stuff (using dictionary directly) - brett 2008

        ##  self.Species_In_Conc = False

        # misc settings
        # write an id header before the array
        self.__settings__['write_array_header'] = 1
        # write a empty line before the array
        self.__settings__['write_array_spacer'] = 1
        self.__settings__['write_array_html_header'] = 1  # write html page header
        self.__settings__['write_array_html_footer'] = 1  # write html page footer
        self.__settings__['write_array_html_format'] = '%2.4f'  # html number format
        # lines to write before flushing to disk
        self.__settings__['write_arr_lflush'] = 5

        # Hybrd options (zero value means routine decides)
        self.__settings__['hybrd_xtol'] = 1.0e-12  # relative error tolerance
        # Maximum number of calls, zero means then 100*(len(species)+1)
        self.__settings__['hybrd_maxfev'] = 0
        # A suitable step length for the forward-difference approximation of the Jacobian
        self.__settings__['hybrd_epsfcn'] = copy.copy(
            self.__settings__['mach_floateps']
        )
        # A parameter determining the initial step bound in interval (0.1,100)
        self.__settings__['hybrd_factor'] = 100
        self.__settings__['hybrd_mesg'] = 1  # print the exit status message

        # Finstslv options (forward integration solver) - brett (20030331)
        # max allowed deviation between max(sim_res) to be a steady state
        self.__settings__['fintslv_tol'] = 1.0e-3
        # threshold number of steps where deviation < atol to be declared a steady state
        self.__settings__['fintslv_step'] = 5
        # a "saturation" type range - should be ok for most systems
        self.__fintslv_range__ = numpy.array(
            [
                1,
                10,
                100,
                1000,
                5000,
                10000,
                50000,
                50100,
                50200,
                50300,
                50400,
                50500,
                50600,
                50700,
                50800,
                50850,
                50900,
                50950,
                51000,
            ],
            'd',
        )
        # self.__fintslv_range__ = numpy.array([1,10,100,1000,5000,10000,50000,100000,500000,1000000,1000100, 1000200,1000300,1000400,1000500,1000600,1000700,1000800],'d')
        self.__settings__['fintslv_rmult'] = 1.0

        # NLEQ2 options (optional non-linear solver) - brett (20030331)
        # IOPT
        # use advanced NLEQ2 features ... may speedup some parameter scans
        self.__settings__['nleq2_advanced_mode'] = False
        # nitmax growth factor per nleq2 iteration
        self.__settings__['nleq2_growth_factor'] = 10
        # number of itertions to loop the solver through (2 should almost always be sufficient)
        self.__settings__['nleq2_iter'] = 3
        # maximum numbe rof iterations in advanced mode
        self.__settings__['nleq2_iter_max'] = 10

        self.__settings__['nleq2_rtol'] = 1.0e-8  # automatically calculated by nleq12
        self.__settings__['nleq2_jacgen'] = 2  # 2:numdiff, 3:numdiff+feedback
        # 0:xscal lower thresholdof scaling vector, 1:always scaling vector
        self.__settings__['nleq2_iscaln'] = 0
        # Dangerous anything not zero seqfaults
        ## self.__settings__['nleq2_mprerr'] = 1       # 0:no output, 1:error, 2:+warning, 3:+info
        # 1:linear, 2:mildly non-lin, 3:highly non-lin, 4:extremely non-lin
        self.__settings__['nleq2_nonlin'] = 4
        # 0:no Broyden approx. rank-1 updates, 1:Broyden approx. rank-1 updates
        self.__settings__['nleq2_qrank1'] = 1
        # 0:Automatic row scaling is active, 1:inactive
        self.__settings__['nleq2_qnscal'] = 0
        # 0:auto damping strategy, 1:damping on, 2:damping off
        self.__settings__['nleq2_ibdamp'] = 0
        # on non-superlinear convergance nitmax *= this (ierr=5)
        self.__settings__['nleq2_nitmax_growth_factor'] = 4
        # 0:default(2) 1:convergance not checked, 2:+'weak stop', 3:+'hard stop' criterion
        self.__settings__['nleq2_iormon'] = 2
        self.__settings__['nleq2_mesg'] = 1  # print the exit status message
        # TODO:
        self.nleq2_nitmax = 50  # optimized and self adjusting :-)

        # Other SteadyState options
        # the initial values of 1:zeroval smaller than 1.0e-6 sometimes give problems
        self.__settings__['small_concentration'] = 1.0e-6
        ##  self.__state_set_conserve__ = 1
        self.__StateOK__ = True

        self.data_sstate = None
        self._sim = None  # new object for recarray with simulation results
        self.data_sim = None  # the new integration results data object
        self.data_stochsim = None  # the new stochastic integration results data object
        self.__scan2d_results__ = None  # yaaaay gnuplot is back
        self.__scan2d_pars__ = None

        # Simulate/lsoda options (zero value means routine decides)
        # The input parameters rtol and atol determine the error
        self.__settings__['lsoda_atol'] = 1.0e-12
        # control performed by the solver. -- johann 20050217 changed from 1e-5
        self.__settings__['lsoda_rtol'] = 1.0e-7
        # maximum number (internally defined) steps allowed per point. 0: x <= 500
        self.__settings__["lsoda_mxstep"] = 0
        # the step size to be attempted on the first step.
        self.__settings__["lsoda_h0"] = 0.0
        self.__settings__["lsoda_hmax"] = 0.0  # the maximum absolute step size allowed.
        self.__settings__["lsoda_hmin"] = 0.0  # the minimum absolute step size allowed.
        # maximum order to be allowed for the nonstiff (Adams) method.
        self.__settings__["lsoda_mxordn"] = 12
        # maximum order to be allowed for the stiff (BDF) method.
        self.__settings__["lsoda_mxords"] = 5
        self.__settings__["lsoda_mesg"] = 1  # print the exit status message

        # try to select the best integration algorithm
        self.mode_integrator = 'LSODA'  # LSODA/CVODE set the intgration algorithm
        if self.__HAS_EVENTS__:
            if _HAVE_PYSUNDIALS:
                print(
                    '\nINFO: events detected and we have PySundials installed, switching to CVODE (mod.mode_integrator=\'CVODE\').'
                )
                self.mode_integrator = 'CVODE'
            else:
                print(
                    '\nWARNING: PySCeS needs PySundials installed for event handling, PySCeS will continue with LSODA (NOTE: ALL EVENTS WILL BE IGNORED!).'
                )
                print(
                    'Windows users may install the "pysces_pysundials" package from pysces.sf.net'
                )
                print(
                    'For other OS\'s see pysundials.sf.net for installation details\n'
                )
                self.__events__ = []
        # if self.__HAS_RATE_RULES__ or self.__HAS_COMPARTMENTS__:
        if self.__HAS_RATE_RULES__:
            if _HAVE_PYSUNDIALS:
                print(
                    'INFO: RateRules detected and PySundials installed, switching to CVODE (mod.mode_integrator=\'CVODE\').\n'
                )
                self.mode_integrator = 'CVODE'
            else:
                print(
                    '\nWARNING: RateRules detected! PySCeS prefers CVODE but will continue with LSODA (NOTE: VARIABLE COMPARTMENTS ARE NOT SUPPORTED WITH LSODA!)'
                )
                print(
                    'Windows users may install the "pysces_pysundials" package from pysces.sf.net'
                )
                print(
                    'For other OS\'s see pysundials.sf.net for installation details\n'
                )

        # pysundials
        self.mode_integrate_all_odes = False  # only available with CVODE
        self.__settings__["cvode_abstol"] = 1.0e-15  # absolute tolerance
        ##  self.__settings__["cvode_abstol_max"]  = 1.0e-3 # not used anymore
        self.__settings__["cvode_abstol_factor"] = 1.0e-6
        self.__settings__["cvode_reltol"] = 1.0e-9  # relative tolerance
        self.__settings__["cvode_auto_tol_adjust"] = True  # auto reduce tolerances
        self.__settings__["cvode_mxstep"] = 1000  # max step default
        # print some pretty stuff after a simulation
        self.__settings__["cvode_stats"] = False

        if self.__HAS_PIECEWISE__ and self.__settings__["cvode_reltol"] <= 1.0e-9:
            self.__settings__["cvode_reltol"] = 1.0e-6
            print(
                'INFO: Piecewise functions detected increasing CVODE tolerance slightly (mod.__settings__[\"cvode_reltol\"] = 1.0e-9 ).'
            )

        # Normal simulation options
        self.sim_start = 0.0
        self.sim_end = 10.0
        self.sim_points = 20.0
        sim_steps = (self.sim_end - self.sim_start) / self.sim_points
        self.sim_time = numpy.arange(
            self.sim_start, self.sim_end + sim_steps, sim_steps
        )
        self.CVODE_extra_output = []
        self.CVODE_xdata = None
        self.__SIMPLOT_OUT__ = []

        # elasticity options
        self.__settings__["elas_evar_upsymb"] = 1  # attach elasticities
        self.__settings__["elas_epar_upsymb"] = 1  # attach parameter elasticities
        # remap steady state values ... no if SimElas
        self.__settings__["elas_evar_remap"] = 1
        # used to determine a stepsize dx=So*factor
        self.__settings__['mode_elas_deriv_factor'] = 0.0001
        # minimum value dx is allowed to have
        self.__settings__['mode_elas_deriv_min'] = 1.0e-12
        # replaces zero fluxes with a very small number
        self.__settings__['elas_zero_flux_fix'] = False
        # replaces zero concentrations with a very small number
        self.__settings__['elas_zero_conc_fix'] = False
        # if Infinite values are created when scaling set to zero
        self.__settings__['elas_scaling_div0_fix'] = False

        # MCA options
        self.__settings__["mca_ccj_upsymb"] = 1  # attach the flux control coefficients
        # attach the concentration control coefficients
        self.__settings__["mca_ccs_upsymb"] = 1
        self.__settings__["mca_ccall_fluxout"] = 1  # in .showCC() output flux cc's
        self.__settings__["mca_ccall_concout"] = 1  # in .showCC() output conc cc's
        # in .showCC() all CC's group by reaction
        self.__settings__["mca_ccall_altout"] = 0

        # gone in 60 seconds brett2008 (Refactored to PyscesLink)
        ##  # Elementary mode options
        ##  self.emode_userout = 0     # write metatool output to file 0:no, 1:yes
        ##  self.emode_intmode = 0  # 0:float metatool, 1:integer metatool
        ##  self.emode_file = self.ModelFile[:-4] + '_emodes'

        # pitcon (pitcon 6.1) continuation options - brett 20040429
        # my interface controls
        # in the REq evaluation values <1.e-15 = 1.e-15
        self.__settings__["pitcon_fix_small"] = 0
        # parameter space that pitcon must search in
        self.pitcon_par_space = numpy.logspace(-1, 3, 10)
        # number of iterations to search for every point in par_space
        self.pitcon_iter = 10
        # initialize with non steady-state values 0:no,1:yes
        self.__settings__["pitcon_allow_badstate"] = 0
        # generate fluxes as output
        self.__settings__["pitcon_flux_gen"] = 1
        # drop pitcon results containing negative concentrations 0:no,1:yes
        self.__settings__["pitcon_filter_neg"] = 1
        # drop output results containing negative concentrations 0:no,1:yes
        self.__settings__["pitcon_filter_neg_res"] = 0
        self.__settings__["pitcon_max_step"] = 30.0  # Maximum stepsize
        # TODO:
        self.pitcon_target_points = []  # list of calculated target points
        self.pitcon_limit_points = []  # list of calculated limit points
        self.pitcon_targ_val = 0.0  # Target value (Seek solution with iwork[4]=)

        # moved to InitialiseConservationArrays
        # self.__settings__["pitcon_max_step"] = 10*(len(self.__SI__)+1) #max corrector steps

        # pitcon integer options iwork in pitcon/dpcon61.f
        self.__settings__["pitcon_init_par"] = 1  # Use X(1) for initial parameter
        # Parameterization option 0:allows program
        self.__settings__["pitcon_par_opt"] = 0
        self.__settings__["pitcon_jac_upd"] = 0  # Update jacobian every newton step
        self.__settings__["pitcon_targ_val_idx"] = 0  # Seek target values for X(n)
        self.__settings__["pitcon_limit_point_idx"] = 0  # Seek limit points in X(n)
        self.__settings__["pitcon_output_lvl"] = 0  # Control amount of output.
        # Jacobian choice. 0:supply jacobian,1:use forward difference,2:central difference
        self.__settings__["pitcon_jac_opt"] = 1
        # pitcon float options rwork in pitcon/dpcon61.f
        self.__settings__["pitcon_abs_tol"] = 0.00001  # Absolute error tolerance
        self.__settings__["pitcon_rel_tol"] = 0.00001  # Relative error tolerance
        self.__settings__["pitcon_min_step"] = 0.01  # Minimum stepsize
        self.__settings__["pitcon_start_step"] = 0.3  # Starting stepsize
        self.__settings__["pitcon_start_dir"] = 1.0  # Starting direction +1.0/-1.0
        self.__settings__["pitcon_max_grow"] = 3.0  # maximum growth factor

        # setup conservation matrices (L_switch is set by stoich to 1 if it "detects" conservation)
        self.InitialiseConservationArrays()

        # scan option - removed from initsimscan
        self.scan_in = ''
        self.scan_out = []
        self._scan = None
        self.scan_res = None
        # should scan1 run mca analysis 0:no,1:elas,2:cc
        self.__settings__["scan1_mca_mode"] = 0
        # should scan1 drop invalid steady states (rare)?
        self.__settings__["scan1_dropbad"] = 0
        # invalid steady states are returned as NaN
        self.__settings__["scan1_nan_on_bad"] = True
        # print out progress messages for large (>20) scans
        self.__settings__["scan1_mesg"] = True
        self.__scan_errors_par__ = None  # collect errors, parameters values
        self.__scan_errors_idx__ = None  # collect errors, indexes

    def InitialiseConservationArrays(self):
        """
        Initialise conservation related vectors/array was in InitialiseModel but has been
        moved out so is can be called by when the stoichiometry is reanalysed

        """
        # Init and Sx vectors # these vectors are multiplying all by themselves - brett
        self.__SI__ = numpy.zeros((self.__lzeromatrix__.shape[1]), 'd')
        self.__SALL__ = numpy.zeros((len(self.__species__)), 'd')
        self.__settings__["pitcon_max_step"] = 10 * (
            len(self.__SI__) + 1
        )  # max corrector steps
        if self.__HAS_MOIETY_CONSERVATION__ == True:
            # reorder the conservation matrix so that it is compatible with L
            idx = [
                list(self.conservation_matrix_col).index(n)
                for n in range(len(self.conservation_matrix_col))
            ]
            self.__reordered_lcons = self.conservation_matrix.take(idx, 1)
            if self.__HAS_COMPARTMENTS__ and self.__KeyWords__['Species_In_Conc']:
                self.__Build_Tvec__(amounts=False)
            else:
                self.__Build_Tvec__(amounts=True)
            self.showConserved(screenwrite=0)

    def InitialiseOldFunctions(self):
        """
        InitialiseOldFunctions()

        Parse and initialise user defined functions specified by !T !U in the PSC input file

        Arguments:
        None

        """
        self.__CODE_user = compile(self._Function_user, '_Function_user', 'exec')
        try:
            exec(self.__CODE_user)
            # print 'done.'
        except Exception as e:
            print('WARNING: User function error\n', e)
            print('This might be due to non-instantiated (e.g. MCA/state) attributes')
            print('Make sure the attributes that are used in your function exist ...')
            print('and manually load using the self.ReloadUserFunc() method')

        # print '\nInitializing init function ...'
        # print self._Function_init

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
        ReloadUserFunc()

        Recompile and execute the user function (!U) from the input file.

        Arguments:
        None

        """
        ##  print "User functions disabled ..."
        self.__CODE_user = compile(self._Function_user, '_Function_user', 'exec')
        try:
            exec(self.__CODE_user)
        except Exception as e:
            print('WARNING: User function load error\n', e)

    def ReloadInitFunc(self):
        """
        ReloadInitFunc()

        Recompile and execute the user initialisations (!I) as defined in the PSC input file.
        and in mod.__InitFuncs__.

        UPDATE 2015: can now be used to define InitialAssignments (no need for self.* prefix in input file)

        Arguments:
        None

        """

        def topolgical_sort(graph_unsorted):
            """
            Repeatedly go through all of the nodes in the graph, moving each of
            the nodes that has all its edges resolved, onto a sequence that
            forms our sorted graph. A node has all of its edges resolved and
            can be moved once all the nodes its edges point to, have been moved
            from the unsorted graph onto the sorted one.
            """

            # This is the list we'll return, that stores each node/edges pair
            # in topological order.
            graph_sorted = []

            # Convert the unsorted graph into a hash table. This gives us
            # constant-time lookup for checking if edges are unresolved, and
            # for removing nodes from the unsorted graph.
            graph_unsorted = dict(graph_unsorted)

            # Run until the unsorted graph is empty.
            while graph_unsorted:

                # Go through each of the node/edges pairs in the unsorted
                # graph. If a set of edges doesn't contain any nodes that
                # haven't been resolved, that is, that are still in the
                # unsorted graph, remove the pair from the unsorted graph,
                # and append it to the sorted graph. Note here that by using
                # using the items() method for iterating, a copy of the
                # unsorted graph is used, allowing us to modify the unsorted
                # graph as we move through it. We also keep a flag for
                # checking that that graph is acyclic, which is true if any
                # nodes are resolved during each pass through the graph. If
                # not, we need to bail out as the graph therefore can't be
                # sorted.
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
                    # Uh oh, we've passed through all the unsorted nodes and
                    # weren't able to resolve any of them, which means there
                    # are nodes with cyclic edges that will never be resolved,
                    # so we bail out with an error.
                    raise RuntimeError("A cyclic dependency occurred")

            return graph_sorted

        simpleAss = []
        exprsAss = []
        symbolsX = []
        for init in self.__InitFuncs__:
            if hasattr(self, init):
                # TODO: to be continued by brett
                # print(init, self.__InitFuncs__[init])
                if type(self.__InitFuncs__[init]) == str:
                    if self.__InitFuncs__[init].isdigit():
                        # self.__InitFuncs__[init] = eval(self.__InitFuncs__[init])
                        # setattr(self, init, self.__InitFuncs__[init])
                        # self.__InitFuncs__[init] = compile(self.__InitFuncs__[init], '_InitFunc_', 'single')
                        simpleAss.append((init, eval(self.__InitFuncs__[init])))
                        # setattr(self, init, eval(self.__InitFuncs__[init]))
                    else:
                        InfixParser.setNameStr('self.', '')
                        # InfixParser.SymbolReplacements = {'_TIME_':'mod._TIME_'}
                        InfixParser.parse(str(self.__InitFuncs__[init]))
                        piecewises = InfixParser.piecewises
                        symbols = InfixParser.names
                        for s_ in symbols:
                            if s_ not in symbolsX:
                                symbolsX.append(s_)
                        if init not in symbolsX:
                            symbolsX.append(init)
                        functions = InfixParser.functions
                        code_string = 'self.{} = {}'.format(init, InfixParser.output)
                        xcode = compile(code_string, '_InitAss_', 'exec')
                        exprsAss.append((init, xcode, symbols))
                        # setattr(self, init, eval(xcode))
                    # else:
                    ##setattr(self, init, self.__InitFuncs__[init])
                    # self.__InitFuncs__[init] = compile(self.__InitFuncs__[init], '_InitFunc_', 'single')
                    # setattr(self, init, eval(self.__InitFuncs__[init]))
                else:
                    setattr(self, init, self.__InitFuncs__[init])
            else:
                assert hasattr(self, init), '\nModelInit error'
        # TODO: to be continued by brett
        # print('symbolsX', symbolsX)
        for init, val in simpleAss:
            # print('s', init,val)
            setattr(self, init, val)
        execOrder = []
        for init, xcode, symbols2 in exprsAss:
            symbols2 = [symbolsX.index(s_) for s_ in symbols2]
            execOrder.append((symbolsX.index(init), symbols2))
            # TODO: to be continued by brett
            # print('x', symbols2)
            # print('e', init,code_string)
            eval(xcode)
        # TODO: to be continued by brett
        # print(execOrder)
        # print(topolgical_sort(execOrder))
        # print(symbolsX)

    def __initREQ__(self):
        """
        __initREQ__()

        Compile and initialise the model rate equations and associated arrays. Rate equations are
        generated as lambda functions callable as mod.R1(mod)

        Arguments:
        None

        """
        # create work arrays for Reactions and variable species
        # s = numpy.zeros(len(self.__species__),'d')
        # self.Vtemp = numpy.zeros(len(self.__reactions__),'d')
        # This appears to be redundant leftover code - brett 20031215

        if self.__settings__['display_debug'] == 1:
            pass
            # print 'Vtemp'

        # create the map string by taking the VarReagent[] and assigning it to an s[index]
        # a reverse map function mapping self.sXi back to self.init[x] for simulation initiation
        mapString = ''
        mapString_R = ''
        for x in range(0, len(self.__species__)):
            mapString += 'self.{} = s[{}]\n'.format(self.__species__[x], x)
            mapString_R += 'self.__inspec__[{}] = self.{}_init\n'.format(
                x, self.__species__[x],
            )
        ##  print mapString_R
        # this creates the mapping string for the derivative functions

        if self.__settings__['display_debug'] == 1:
            print('mapString')
            print(mapString)
            print('mapString_R')
            print(mapString_R)

        # create the REq string in ReactionIDs order as a single multiline definition
        # using indexes: vString (Vtemp[x] =)
        # or a list of individual compile definitions: vArray (v + ReactionID = )
        vString = ''
        ##  vString2 = ''
        ##  DvOrder = []
        EStat = []
        ##  dispREq = ''

        symbR = {}
        if (
            len(list(self.__rules__.keys())) > 0
            and self.__settings__['mode_substitute_assignment_rules']
        ):
            # create substitution dictionary
            for ass in self.__rules__:
                symbR.update(
                    {self.__rules__[ass]['name']: self.__rules__[ass]['code_string']}
                )

        for x in range(0, len(self.__reactions__)):
            req1 = self.__nDict__[self.__reactions__[x]]['RateEq']
            ##  DvOrder.append(self.__reactions__[x])
            vString += 'Vtemp[{}] = self.{}()\n'.format(x, self.__reactions__[x])

            # Core update inspired by Core2, lambda functions replaced by Reaction instances

            if (
                len(list(self.__rules__.keys())) > 0
                and self.__settings__['mode_substitute_assignment_rules']
            ):
                # substitute assignment rules
                InfixParser.setNameStr('self.', '')
                InfixParser.FunctionReplacements = symbR
                InfixParser.parse(req1.replace('self.', ''))
                req2 = InfixParser.output
                ##  req1_names = InfixParser.names
                rObj = ReactionObj(self, self.__reactions__[x], req2, 'self.')
            else:
                rObj = ReactionObj(self, self.__reactions__[x], req1, 'self.')
            self.__ParsePiecewiseFunctions__(rObj.piecewises)
            rObj.compartment = self.__nDict__[rObj.name]['compartment']

            setattr(self, self.__reactions__[x], rObj)

        # this is to check that if a model defines compartments then REq's MUST be in amounts brett2008
        # brett2008
        if self.__HAS_COMPARTMENTS__:
            cnames = [self.__compartments__[c]['name'] for c in self.__compartments__]
            warnings = ''
            for rr in self.__reactions__:
                rrobj = getattr(self, rr)
                warn = False
                if (
                    rrobj.compartment == None
                    and self.__settings__['display_compartment_warnings']
                ):
                    warnings += "# {} is not located in a compartment.\n".format(
                        rrobj.name
                    )
                else:
                    for comp in cnames:
                        if comp in rrobj.symbols:
                            warn = False
                            break
                        else:
                            warn = True
                    if warn:
                        warnings += "# {}: {}\n#  assuming kinetic constants are flow constants.\n".format(
                            rr, rrobj.formula
                        )
            if warnings != '' and self.__settings__['display_compartment_warnings']:
                print('\n# -- COMPARTMENT WARNINGS --')
                print(warnings)
        if self.__settings__['display_debug'] == 1:
            print('vString')
            print(vString)
            ##  print 'vString2'
            ##  print vString2
            ##  print 'DvOrder'
            ##  print DvOrder
        return mapString, mapString_R, vString

    def __initmodel__(self):
        """
        __initmodel__()

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

    def __initvar__(self):
        """
        __initvar__()

        Compile and initialise the model variable species and derivatives and associated mapping arrays.

        Arguments:
        None

        """
        mvarString = ''
        sOrder = []
        self.__remaps = ''
        DvarUpString = ''

        # self.__inspec__ = numpy.zeros((len(self.__species__)),'d') moved to ParseModel
        for x in range(0, len(self.__species__)):
            key = self.__species__[x]
            ##  self.__inspec__[x] = eval(key)
            self.__inspec__[x] = getattr(self, key)
            mvarString += 'self.__inspec__[{}] = float(self.{})\n'.format(x, key)
            self.__remaps += 'self.{} = self.{}_ss\n'.format(key, key)
            sOrder.append('self.' + key)
            DvarUpString += 'self.{} = input[{}]\n'.format(self.__species__[x], x)

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
        __initfixed__()

        Compile and initialise the fixed species and associated mapping arrays

        Arguments:
        None

        """
        runmapString = ''
        if self.__settings__['display_debug'] == 1:
            print('InitParams2')
            print(self.__parameters__)
            print('InitStrings2')
            ##  print self.__InitStrings

        # Initialise parameter elasticities
        par_map2store = ''
        par_remap = ''
        for x in range(len(self.__parameters__)):
            par_map2store += 'parVal_hold[{}] = self.{}\n'.format(
                x, self.__parameters__[x]
            )
            par_remap += 'self.{} = parVal_hold[{}]\n'.format(self.__parameters__[x], x)

        if self.__settings__['display_debug'] == 1:
            print('par_map2store')
            print(par_map2store)
            print('par_remap')
            print(par_remap)
        return (par_map2store, par_remap)

    # pysces core - the steady-state solver and integration routines
    def _SpeciesAmountToConc(self, s):
        """takes and returns s"""
        for x in range(len(self.__CsizeAllIdx__)):
            self.__CsizeAll__[x] = getattr(self, self.__CsizeAllIdx__[x])
        ##  print '__CALL__', self.__CsizeAll__
        return s / self.__CsizeAll__

    def _SpeciesConcToAmount(self, s):
        """takes and returns s"""
        for x in range(len(self.__CsizeAllIdx__)):
            self.__CsizeAll__[x] = getattr(self, self.__CsizeAllIdx__[x])
        ##  print '__CALL__', self.__CsizeAll__
        return s * self.__CsizeAll__

    def _FixedSpeciesAmountToConc(self):
        for x in range(len(self.__fixed_species__)):
            am = getattr(self, self.__fixed_species__[x])
            self.__CsizeFS__[x] = getattr(self, self.__CsizeFSIdx__[x])
            setattr(self, self.__fixed_species__[x], am / self.__CsizeFS__[x])

    def _FixedSpeciesConcToAmount(self):
        for x in range(len(self.__fixed_species__)):
            cn = getattr(self, self.__fixed_species__[x])
            self.__CsizeFS__[x] = getattr(self, self.__CsizeFSIdx__[x])
            setattr(self, self.__fixed_species__[x], cn * self.__CsizeFS__[x])

    # extract SI and "correct" SD in s solve for v
    # Vtemp is passed through the function to avoid an irritating bug
    # that appears if it isn't (ie defined as a global) - brett
    def _EvalREq2_alt(self, s, Vtemp):
        """
        _EvalREq2_alt(s,Vtemp)

        Evaluate the rate equations correcting for mass conservation.
        Takes full [s] returns full [s] --> moiety aware

        Arguments:

        s: a vector of species cosnservations
        Vtemp: the rate vector

        """
        # form an SI vector
        for x in range(0, len(self.lzeromatrix_col)):
            self.__SI__[x] = s[self.lzeromatrix_col[x]]

        # replace SD with SI calculated values
        for x in range(0, len(self.lzeromatrix_row)):
            s[self.lzeromatrix_row[x]] = self.__tvec_a__[x] + numpy.add.reduce(
                self.__lzeromatrix__[x, :] * self.__SI__
            )
        return self._EvalREq(s, Vtemp)

    def _EvalREq(self, s, Vtemp):
        if self.__HAS_RATE_RULES__:
            exec(self.__CODE_raterule_map)
        if self.__HAS_COMPARTMENTS__ and self.__KeyWords__['Species_In_Conc']:
            s = self._SpeciesAmountToConc(s)
        exec(self.__mapFunc__)
        try:
            self.Forcing_Function()
        except Exception as de:
            print('INFO: forcing function failure', de)
        try:
            exec(self.__vFunc__)
        except (
            ArithmeticError,
            AttributeError,
            NameError,
            ZeroDivisionError,
            ValueError,
        ) as detail:
            print('INFO: REq evaluation failure:', detail)
            Vtemp[:] = self.__settings__['mach_floateps']
        if self.__HAS_RATE_RULES__:
            try:
                exec(self.__CODE_raterule)
            except (
                ArithmeticError,
                AttributeError,
                NameError,
                ZeroDivisionError,
                ValueError,
            ) as detail:
                print('INFO: RateRule evaluation failure:', detail)
        return Vtemp

    def _EvalREq2(self, s, Vtemp):
        """
        _EvalREq2(s,Vtemp)

        Evaluate the rate equations, as PySCeS uses the reduced set of ODE's for its core operations, this method
        takes mass conservation into account and regenerates the full species vector
        takes [si] returns full [s]

        Arguments:

        s: species vector
        Vtemp: rate vector

        """
        # stick SI into s using Nrrow order
        for x in range(len(s)):
            self.__SALL__[self.nrmatrix_row[x]] = s[x]
        # stick SD into s using lzerorow order
        for x in range(len(self.lzeromatrix_row)):
            self.__SALL__[self.lzeromatrix_row[x]] = self.__tvec_a__[
                x
            ] + numpy.add.reduce(self.__lzeromatrix__[x, :] * s)
        return self._EvalREq(self.__SALL__, Vtemp)

    def _EvalODE(self, s, Vtemp):
        """
        _EvalODE(s,Vtemp)

        Core ODE evaluation routine evaluating the reduced set of ODE's. Depending on whether mass conservation is present or not either N*v (_EvalREq) or Nr*v (_EvalREq2) is used to automatically generate the ODE's.
        Returns sdot.

        Arguments:

        s: species vector
        Vtemp: rate vector

        """
        if self.__HAS_RATE_RULES__:
            s, self.__rrule__ = numpy.split(s, [self.Nrmatrix.shape[0]])
        if self.__HAS_MOIETY_CONSERVATION__ == True:
            for x in range(len(s)):
                self.__SALL__[self.nrmatrix_row[x]] = s[x]
            for x in range(len(self.lzeromatrix_row)):
                self.__SALL__[self.lzeromatrix_row[x]] = self.__tvec_a__[
                    x
                ] + numpy.add.reduce(self.__lzeromatrix__[x, :] * s)
            self.__vvec__ = self._EvalREq(self.__SALL__, Vtemp)
            s = numpy.add.reduce(self.__nrmatrix__ * self.__vvec__, 1)
        else:
            self.__vvec__ = self._EvalREq(s, Vtemp)
            s = numpy.add.reduce(self.__nmatrix__ * self.__vvec__, 1)
        if self.__HAS_RATE_RULES__:
            s = numpy.concatenate([s, self.__rrule__]).copy()
        return s

    # the following are user functions defined in the input file or assigned after instantiated
    def _EvalODE_CVODE(self, s, Vtemp):
        """
        _EvalODE_CVODE(s,Vtemp)

        Evaluate the full set of ODE's. Depending on whether mass conservation is present or not
        either N*v (_EvalREq) or Nr*v (_EvalREq2_alt) is used.

        Arguments:

        s: species vector
        Vtemp: rate vector

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

        User defined forcing function either defined in the PSC input file as !F or by overwriting this method.
        This method is evaluated prior to every rate equation evaluation.

        Arguments:
        None

        """
        # initialized in __initFunction__() - brett 20050621
        exec(self.__CODE_forced)

    def User_Function(self):
        """
        **Deprecated**
        """
        exec(self.__CODE_user)

    # pysundials

    def _EvalExtraData(self, xdata):
        """
        Takes a list of model attributes and returns an array of values
        """
        return numpy.array([getattr(self, d) for d in xdata])

    __CVODE_initialise__ = True
    CVODE_continuous_result = None
    __CVODE_initial_num__ = None

    def CVODE_continue(self, tvec):
        """
        **Experimental:** continues a simulation over a new time vector, the
        CVODE memobj is reused and not reinitialised and model parameters can be
        changed between calls to this method. The mod.data_sim objects from
        the initial simulation and all calls to this method are stored in the list
        *mod.CVODE_continuous_result*.

         - *tvec* a numpy array of time points

        """
        assert (
            self.mode_integrator == 'CVODE'
        ), "\nFor what should be rather obvious reasons, this method requires CVODE to be used as the default integration algorithm.\n"
        if self.CVODE_continuous_result == None:
            self.CVODE_continuous_result = [self.data_sim]

        self.__CVODE_initialise__ = False
        self.sim_time = tvec

        Tsim0 = time.time()
        sim_res, rates, simOK = self.CVODE(None)
        Tsim1 = time.time()
        print(
            "{} time for {} points: {}".format(
                self.mode_integrator, len(self.sim_time), Tsim1 - Tsim0
            )
        )

        if self.__HAS_RATE_RULES__:
            sim_res, rrules = numpy.split(sim_res, [len(self.__species__)], axis=1)
            print('RateRules evaluated and added to mod.data_sim.')

        # TODO: split this off into a method shared by this and Simulate()
        self.data_sim = IntegrationDataObj()
        self.IS_VALID = simOK
        self.data_sim.setTime(self.sim_time)
        self.data_sim.setSpecies(sim_res, self.__species__)
        self.data_sim.setRates(rates, self.__reactions__)
        if self.__HAS_RATE_RULES__:
            self.data_sim.setRules(rrules, self.__rate_rules__)
        if len(self.CVODE_extra_output) > 0:
            self.data_sim.setXData(self.CVODE_xdata, lbls=self.CVODE_extra_output)
            self.CVODE_xdata = None
        if not simOK:
            print('Simulation failure')
        del sim_res

        self.CVODE_continuous_result.append(self.data_sim)
        self.__CVODE_initialise__ = True

    def CVODE(self, initial):
        """
        CVODE(initial)

        PySCeS interface to the CVODE integration algorithm.

        Arguments:
        initial: vector containing initial species concentrations

        """
        assert (
            _HAVE_PYSUNDIALS
        ), '\nPySundials is not installed or did not import correctly\n{}'.format(
            _PYSUNDIALS_LOAD_ERROR
        )
        Vtemp = numpy.zeros((self.__Nshape__[1]))

        def findi(t, y, ydot, f_data):
            self._TIME_ = t
            ydota = self._EvalODE(numpy.array(y), Vtemp)
            ydot[:] = ydota[:]
            return 0  # non-zero return indicates error state

        def ffull(t, y, ydot, f_data):
            self._TIME_ = t
            ##  ya = numpy.array(y)
            ydota = self._EvalODE_CVODE(numpy.array(y), Vtemp)  # unreduced ODE's
            ydot[:] = ydota[:]
            return 0

        func = None
        if self.mode_integrate_all_odes:
            func = ffull
        else:
            func = findi

        if self.__CVODE_initialise__:
            tZero = initial.copy()
            if self.__HAS_RATE_RULES__:
                initial, rrules = numpy.split(initial, [self.Nrmatrix.shape[0]])
                tZero = initial.copy()
            if self.__HAS_MOIETY_CONSERVATION__ and self.mode_integrate_all_odes:
                initial = self.Fix_S_indinput(initial, amounts=True)
                tZero = initial.copy()
            elif self.__HAS_MOIETY_CONSERVATION__:
                tZero = self.Fix_S_indinput(tZero, amounts=True)
            if self.__HAS_RATE_RULES__:
                initial = numpy.concatenate([initial, rrules])
                tZero = numpy.concatenate([tZero, rrules])
        else:
            tZero = None

        # the following block initialises the cvode integrator, and sets various options
        if self.__CVODE_initialise__:
            self.__CVODE_y__ = cvode.NVector(initial.tolist())
            self.__CVODE_mem__ = cvode.CVodeCreate(
                cvode.CV_BDF, cvode.CV_NEWTON
            )  # initialisation with basic options, newtonian solver etc...
            self.__CVODE_initial_num__ = len(initial)
        del initial
        t = cvode.realtype(0.0)
        rates = numpy.zeros((len(self.sim_time), len(self.__reactions__)))
        output = None
        if self.__HAS_RATE_RULES__:
            output = numpy.zeros(
                (len(self.sim_time), len(self.__rrule__) + len(self.__species__))
            )
        else:
            output = numpy.zeros((len(self.sim_time), len(self.__species__)))

        CVODE_XOUT = False
        if len(self.CVODE_extra_output) > 0:
            out = []
            for d in self.CVODE_extra_output:
                if (
                    hasattr(self, d)
                    and d
                    not in self.__species__ + self.__reactions__ + self.__rate_rules__
                ):
                    out.append(d)
                else:
                    print(
                        '\nWarning: CVODE is ignoring extra data ({}), it either doesn\'t exist or it\'s a species or rate.\n'.format(
                            d
                        )
                    )
            if len(out) > 0:
                self.CVODE_extra_output = out
                CVODE_XOUT = True
            del out

        if CVODE_XOUT:
            self.CVODE_xdata = numpy.zeros(
                (len(self.sim_time), len(self.CVODE_extra_output))
            )

        sim_st = 0
        if self.sim_time[0] == 0.0:

            if self.__CVODE_initialise__:
                out0 = tZero[:].copy()
            else:
                out0 = numpy.array(self.__CVODE_y__[:])
            output[0, :] = out0

            if not self.mode_integrate_all_odes:
                self._EvalODE(out0.copy(), self._CVODE_Vtemp)
            else:
                self._EvalODE_CVODE(out0.copy(), self._CVODE_Vtemp)
            rates[0] = self.__vvec__

            if CVODE_XOUT:
                self.CVODE_xdata[0, :] = self._EvalExtraData(self.CVODE_extra_output)
            sim_st = 1
        del tZero

        var_store = {}
        if self.__HAS_EVENTS__:
            for ev in self.__events__:
                ev.reset()
                for ass in ev.assignments:
                    var_store.update({ass.variable: getattr(self, ass.variable)})
        TOL_ADJUSTER = 0
        MAX_TOL_CNT = 5
        MAX_REL_TOLERANCE = 1.0e-3
        RELTOL_ADJUST_FACTOR = 1.0e3
        MIN_ABS_TOL = self.__settings__["cvode_abstol"]  # 1.0e-15
        ##  MAX_ABS_TOL = self.__settings__["cvode_abstol_max"] #1.0e-3 not used anymore
        ABSTOL_ADJUST_FACTOR = self.__settings__["cvode_abstol_factor"]  # 1.0e-6
        cvode_sim_range = list(range(sim_st, len(self.sim_time)))
        cvode_scale_range = list(
            range(
                sim_st,
                len(self.sim_time),
                len(self.sim_time) // 4 or len(self.sim_time),
            )
        )
        cvode_scale_range = cvode_scale_range[1:]
        reltol = cvode.realtype(
            self.__settings__["cvode_reltol"]
        )  # relative tolerance must be a realtype
        abstol = cvode.NVector(
            self.__CVODE_initial_num__ * [self.__settings__["cvode_abstol"]]
        )
        for s in range(len(self.__CVODE_y__)):
            newVal = abs(self.__CVODE_y__[s]) * ABSTOL_ADJUST_FACTOR
            if newVal < MIN_ABS_TOL:
                abstol[s] = MIN_ABS_TOL
            else:
                abstol[s] = newVal
        if self.__CVODE_initialise__:
            cvode.CVodeMalloc(
                self.__CVODE_mem__,
                func,
                0.0,
                self.__CVODE_y__,
                cvode.CV_SV,
                reltol,
                abstol,
            )  # set tolerances, specify vector of initial conditions, set integrtion function etc...
            cvode.CVDense(
                self.__CVODE_mem__, self.__CVODE_initial_num__
            )  # set dense option, specify dimension of problem (3)
        if self.__HAS_EVENTS__:
            cvode.CVodeRootInit(
                self.__CVODE_mem__, len(self.__events__), self.CVODE_EVENTS, None
            )  # specify to 'root' conditions, and function that calculates them

        for st in cvode_sim_range:
            tout = self.sim_time[st]
            errcount = 0
            ##  print TOL_ADJUSTER, MAX_TOL_CNT, abstol, MAX_REL_TOLERANCE
            if (
                self.__settings__["cvode_auto_tol_adjust"]
                and TOL_ADJUSTER >= MAX_TOL_CNT
                and reltol.value < MAX_REL_TOLERANCE
            ):
                reltol.value = reltol.value * RELTOL_ADJUST_FACTOR
                cvode.CVodeSetTolerances(
                    self.__CVODE_mem__, cvode.CV_SV, reltol, abstol
                )
                ##  print '\nCVODE: new tolerance set:\nreltol={}'.format(reltol.value)
                ##  print '\nAbs tolerance:\n{}'.format(abstol)
                TOL_ADJUSTER = 0
            if (st in cvode_scale_range) and self.__settings__["cvode_auto_tol_adjust"]:
                for s in range(len(self.__CVODE_y__)):
                    newVal = abs(self.__CVODE_y__[s]) * ABSTOL_ADJUST_FACTOR
                    if newVal < MIN_ABS_TOL:
                        abstol[s] = MIN_ABS_TOL
                    else:
                        abstol[s] = newVal
                ##  print '\nCVODE: new tolerance set, abstol:\n{}'.format(abstol)
                cvode.CVodeSetTolerances(
                    self.__CVODE_mem__, cvode.CV_SV, reltol, abstol
                )
                ##  cvode.CVodeReInit(self.__CVODE_mem__, func, self.sim_time[st-1], self.__CVODE_y__, cvode.CV_SV, reltol, abstol)
            while True:
                try:
                    flag = cvode.CVode(
                        self.__CVODE_mem__,
                        tout,
                        self.__CVODE_y__,
                        cvode.ctypes.byref(t),
                        cvode.CV_NORMAL,
                    )
                except AssertionError as ex:
                    print('cvode error1', ex)
                    flag = None
                self._TIME_ = tout
                if (
                    flag == cvode.CV_ROOT_RETURN
                ):  # if a root was found before desired time point, output it
                    ya = numpy.array(self.__CVODE_y__)
                    rootsfound = cvode.CVodeGetRootInfo(
                        self.__CVODE_mem__, len(self.__events__)
                    )
                    reInit = False
                    for ev in range(len(self.__events__)):
                        if rootsfound[ev] == 1:
                            for ass in self.__events__[ev].assignments:
                                # only can assign to independent species vector
                                if ass.variable in self.L0matrix.getLabels()[1] or (
                                    self.mode_integrate_all_odes
                                    and ass.variable in self.__species__
                                ):
                                    assVal = ass.getValue()
                                    assIdx = self.__species__.index(ass.variable)
                                    if self.__KeyWords__['Species_In_Conc']:
                                        ##  print self.__CVODE_y__
                                        self.__CVODE_y__[assIdx] = assVal * getattr(
                                            self, self.__CsizeAllIdx__[assIdx]
                                        )
                                        ##  raw_input(self.__CVODE_y__)
                                    else:
                                        self.__CVODE_y__[assIdx] = assVal
                                    reInit = True
                                elif (
                                    not self.mode_integrate_all_odes
                                    and ass.variable in self.L0matrix.getLabels()[0]
                                ):
                                    print(
                                        'Event assignment to dependent species consider setting \"mod.mode_integrate_all_odes = True\"'
                                    )
                                elif (
                                    self.__HAS_RATE_RULES__
                                    and ass.variable in self.__rate_rules__
                                ):
                                    ##  print 'Event is assigning to rate rule'
                                    assVal = ass.getValue()
                                    rrIdx = self.__rate_rules__.index(ass.variable)
                                    ##  print ass.variable, assVal
                                    self.__rrule__[rrIdx] = assVal
                                    ##  print self.L0matrix.shape[1], rrIdx, len(self.__CVODE_y__)
                                    self.__CVODE_y__[
                                        self.L0matrix.shape[1] + rrIdx
                                    ] = assVal
                                    setattr(self, ass.variable, assVal)
                                    reInit = True
                                else:
                                    try:
                                        setattr(self, ass.variable, ass.getValue())
                                        reInit = True
                                    except:
                                        print(
                                            'ERROR: Updating model attribute from event: ',
                                            ass.variable,
                                        )

                    if reInit:
                        cvode.CVodeReInit(
                            self.__CVODE_mem__,
                            func,
                            tout,
                            self.__CVODE_y__,
                            cvode.CV_SV,
                            reltol,
                            abstol,
                        )

                    # this gets everything into the current tout state
                    tmp = None
                    if not self.mode_integrate_all_odes:
                        tmp = self._EvalODE(ya.copy(), self._CVODE_Vtemp)
                    else:
                        tmp = self._EvalODE_CVODE(ya.copy(), self._CVODE_Vtemp)
                    del tmp

                    # here we regenerate Sd's and fix concentrations
                    rrules = None
                    if (
                        self.__HAS_MOIETY_CONSERVATION__
                        and not self.mode_integrate_all_odes
                    ):
                        if self.__HAS_RATE_RULES__:
                            ya, rrules = numpy.split(ya, [self.Nrmatrix.shape[0]])
                        ya = self.Fix_S_indinput(ya, amounts=True)
                        # convert to concentrations
                        if (
                            self.__HAS_COMPARTMENTS__
                            and self.__KeyWords__['Output_In_Conc']
                        ):
                            ya = self._SpeciesAmountToConc(ya)
                        if self.__HAS_RATE_RULES__:
                            ya = numpy.concatenate([ya, rrules])
                    else:
                        if self.__HAS_RATE_RULES__:
                            ya, rrules = numpy.split(ya, [self.Nmatrix.shape[0]])
                        if (
                            self.__HAS_COMPARTMENTS__
                            and self.__KeyWords__['Output_In_Conc']
                        ):
                            ya = self._SpeciesAmountToConc(ya)
                        if self.__HAS_RATE_RULES__:
                            ya = numpy.concatenate([ya, rrules])

                    output[st] = ya
                    # set with self._EvalODE above
                    rates[st] = self.__vvec__

                    if CVODE_XOUT:
                        self.CVODE_xdata[st, :] = self._EvalExtraData(
                            self.CVODE_extra_output
                        )
                    # this should adjust the expected time to the new output time time
                    self.sim_time[st] = float(tout)
                    # dont need anymore i think
                    if _HAVE_VPYTHON:
                        self.CVODE_VPYTHON(ya)
                    del ya, rrules
                    break
                if flag == cvode.CV_SUCCESS:
                    ya = numpy.array(self.__CVODE_y__)
                    # this gets everything into the current tout state
                    tmp = None
                    if not self.mode_integrate_all_odes:
                        tmp = self._EvalODE(ya.copy(), self._CVODE_Vtemp)
                    else:
                        tmp = self._EvalODE_CVODE(ya.copy(), self._CVODE_Vtemp)
                    del tmp
                    # here we regenerate Sd's and fix concentrations
                    rrules = None
                    if (
                        self.__HAS_MOIETY_CONSERVATION__
                        and not self.mode_integrate_all_odes
                    ):
                        if self.__HAS_RATE_RULES__:
                            ya, rrules = numpy.split(ya, [self.Nrmatrix.shape[0]])
                        ya = self.Fix_S_indinput(ya, amounts=True)
                        # convert to concentrations
                        if (
                            self.__HAS_COMPARTMENTS__
                            and self.__KeyWords__['Output_In_Conc']
                        ):
                            ya = self._SpeciesAmountToConc(ya)
                        if self.__HAS_RATE_RULES__:
                            ya = numpy.concatenate([ya, rrules])
                    else:
                        if self.__HAS_RATE_RULES__:
                            ya, rrules = numpy.split(ya, [self.Nmatrix.shape[0]])
                        if (
                            self.__HAS_COMPARTMENTS__
                            and self.__KeyWords__['Output_In_Conc']
                        ):
                            ya = self._SpeciesAmountToConc(ya)
                        if self.__HAS_RATE_RULES__:
                            ya = numpy.concatenate([ya, rrules])

                    output[st] = ya
                    # set with self._EvalODE above
                    rates[st] = self.__vvec__

                    if CVODE_XOUT:
                        self.CVODE_xdata[st, :] = self._EvalExtraData(
                            self.CVODE_extra_output
                        )
                    if _HAVE_VPYTHON:
                        self.CVODE_VPYTHON(ya)
                    del ya, rrules
                    break
                elif flag == -1:
                    if self.__settings__["cvode_mxstep"] == 1000:
                        self.__settings__["cvode_mxstep"] = 3000
                        TOL_ADJUSTER += 1
                    elif self.__settings__["cvode_mxstep"] == 3000:
                        self.__settings__["cvode_mxstep"] = 10000
                        TOL_ADJUSTER += 2
                    elif self.__settings__["cvode_mxstep"] == 10000:
                        ##  TOL_ADJUSTER += 1
                        output[st] = numpy.NaN
                        break
                    print(
                        'mxstep warning ({}) mxstep set to {}'.format(
                            flag, self.__settings__["cvode_mxstep"]
                        )
                    )
                    cvode.CVodeSetMaxNumSteps(
                        self.__CVODE_mem__, self.__settings__["cvode_mxstep"]
                    )
                elif flag < -3:
                    print('CVODE error:', flag)
                    print('At ', tout)
                    output[st] = numpy.NaN
                    rates[st] = numpy.NaN
                    if CVODE_XOUT:
                        self.CVODE_xdata[st] = numpy.NaN
                    break
            self.__settings__["cvode_mxstep"] = 1000
            cvode.CVodeSetMaxNumSteps(
                self.__CVODE_mem__, self.__settings__["cvode_mxstep"]
            )

        if self.__HAS_EVENTS__:
            for ass in list(var_store.keys()):
                # print 'old value', ass, getattr(self, ass)
                setattr(self, ass, var_store[ass])
                # print 'new value', getattr(self, ass)

        if self.__settings__["cvode_stats"]:
            # print some stats from the intgrator
            nst = cvode.CVodeGetNumSteps(self.__CVODE_mem__)
            nfe = cvode.CVodeGetNumRhsEvals(self.__CVODE_mem__)
            nsetups = cvode.CVodeGetNumLinSolvSetups(self.__CVODE_mem__)
            netf = cvode.CVodeGetNumErrTestFails(self.__CVODE_mem__)
            nni = cvode.CVodeGetNumNonlinSolvIters(self.__CVODE_mem__)
            ncfn = cvode.CVodeGetNumNonlinSolvConvFails(self.__CVODE_mem__)
            nje = cvode.CVDenseGetNumJacEvals(self.__CVODE_mem__)
            nfeLS = cvode.CVDenseGetNumRhsEvals(self.__CVODE_mem__)
            nge = cvode.CVodeGetNumGEvals(self.__CVODE_mem__)

            print("\nFinal Statistics:")
            print(
                "nst = {} nfe  = {} nsetups = {} nfeLS = {} nje = {}".format(
                    nst, nfe, nsetups, nfeLS, nje
                )
            )
            print(
                "nni = {} ncfn = {} netf = {} nge = {}\n ".format(nni, ncfn, netf, nge)
            )
            print('reltol = {}'.format(reltol))
            print('abstol:\n{}'.format(abstol))

        if cvode.CV_SUCCESS >= 0:
            return output, rates, True
        else:
            return output, rates, False

    def CVODE_EVENTS(self, t, svec, eout, f_data):
        svec = numpy.array(svec)
        self._TIME_ = t
        rrules = []
        tmp = None
        if not self.mode_integrate_all_odes:
            tmp = self._EvalODE(svec.copy(), self._CVODE_Vtemp)
        else:
            tmp = self._EvalODE_CVODE(svec.copy(), self._CVODE_Vtemp)
        del tmp
        eventFired = False
        for ev in range(len(self.__events__)):
            if self.__events__[ev](t):
                eout[ev] = 0
                eventFired = True
            else:
                eout[ev] = 1
        if self.__HAS_RATE_RULES__:
            exec(self.__CODE_raterule)
        return 0

    def CVODE_VPYTHON(self, s):
        """Future VPython hook for CVODE"""
        pass

    def LSODA(self, initial):
        """
        LSODA(initial)

        PySCeS interface to the LSODA integration algorithm. Given a set of initial
        conditions LSODA returns an array of species concentrations and a status flag.
        LSODA controls are accessible as mod.lsoda_<control>

        Arguments:

        initial: vector containing initial species concentrations

        """
        Vtemp = numpy.zeros((self.__Nshape__[1]), 'd')

        def function_sim(s, t):
            self._TIME_ = t
            return self._EvalODE(s, Vtemp)

        iter = 0
        go = True

        status = 0
        while go:
            sim_res, infodict = scipy.integrate.odeint(
                function_sim,
                initial,
                self.sim_time,
                atol=self.__settings__['lsoda_atol'],
                rtol=self.__settings__['lsoda_rtol'],
                mxstep=self.__settings__["lsoda_mxstep"],
                h0=self.__settings__["lsoda_h0"],
                hmax=self.__settings__["lsoda_hmax"],
                hmin=self.__settings__["lsoda_hmin"],
                mxordn=self.__settings__["lsoda_mxordn"],
                mxords=self.__settings__["lsoda_mxords"],
                printmessg=self.__settings__["lsoda_mesg"],
                full_output=1,
            )
            if infodict['message'] == 'Integration successful.':
                status = 0
            else:
                status = 1
            if status > 0 and iter < self.__settings__['mode_sim_max_iter']:
                if self.__settings__["lsoda_mxstep"] == 0:
                    print(
                        '\nIntegration error\n\nSetting self.__settings__["lsoda_mxstep"] = 1000 and reSimulating ...'
                    )
                    self.__settings__["lsoda_mxstep"] = 1000
                else:
                    print(
                        'Integration error\n\nSetting self.__settings__["lsoda_mxstep"] = '
                        + repr(self.__settings__["lsoda_mxstep"] * 3)
                        + ' and reSimulating ...'
                    )
                    self.__settings__["lsoda_mxstep"] = (
                        self.__settings__["lsoda_mxstep"] * 3
                    )
                iter += 1
            elif status > 0 and iter == self.__settings__['mode_sim_max_iter']:
                print(
                    '\nThis simulation is going nowhere fast\nConsider trying CVODE (mod.mode_integrator = \'CVODE\')\n'
                )
                print(
                    'self.__settings__["lsoda_mxstep"] = '
                    + repr(self.__settings__["lsoda_mxstep"])
                )
                print('__settings__[\'mode_sim_max_iter\'] = ' + repr(iter))
                go = False
            else:
                go = False
        self.__settings__["lsoda_mxstep"] = 0

        rates = numpy.zeros((sim_res.shape[0], len(self.__reactions__)))
        if status == 0:

            tmp = None
            for r in range(sim_res.shape[0]):
                tmp = self._EvalODE(sim_res[r].copy(), self._CVODE_Vtemp)
                # set with self._EvalODE above
                rates[r] = self.__vvec__
            del tmp

            # regenerate dependent variables
            if self.__HAS_RATE_RULES__:
                sim_res, rrules = numpy.split(sim_res, [self.Nrmatrix.shape[0]], axis=1)
            if self.__HAS_MOIETY_CONSERVATION__ == True:
                res = numpy.zeros((sim_res.shape[0], len(self.__species__)))
                for x in range(0, sim_res.shape[0]):
                    res[x, :] = self.Fix_S_indinput(sim_res[x, :], amounts=True)
                sim_res = res
                del res
            if self.__HAS_COMPARTMENTS__ and self.__KeyWords__['Output_In_Conc']:
                for x in range(0, sim_res.shape[0]):
                    sim_res[x] = self._SpeciesAmountToConc(sim_res[x])

            if self.__HAS_RATE_RULES__:
                sim_res = numpy.concatenate([sim_res, rrules], axis=1)
            return sim_res, rates, True
        else:
            if self.__HAS_MOIETY_CONSERVATION__ == True and self.__HAS_RATE_RULES__:
                sim_res = numpy.zeros(
                    (sim_res.shape[0], len(self.__species__) + len(self.__rrule__)), 'd'
                )
            elif self.__HAS_MOIETY_CONSERVATION__ == True:
                sim_res = numpy.zeros((sim_res.shape[0], len(self.__species__)), 'd')
            sim_res[:] = scipy.NaN
            return sim_res, rates, False

    def HYBRD(self, initial):
        """
        HYBRD(initial)

        PySCeS interface to the HYBRD solver. Returns a steady-state solution and
        error flag. Good general purpose solver.
        Algorithm controls are available as mod.hybrd_<control>

        Arguments:

        initial: vector of initial species concentrations

        """
        Vtemp = numpy.zeros((self.__Nshape__[1]), 'd')

        def function_state(s):
            return self._EvalODE(s, Vtemp)

        state_out = scipy.optimize.fsolve(
            function_state,
            initial,
            args=(),
            xtol=self.__settings__['hybrd_xtol'],
            maxfev=self.__settings__['hybrd_maxfev'],
            epsfcn=self.__settings__['hybrd_epsfcn'],
            factor=self.__settings__['hybrd_factor'],
            col_deriv=0,
            full_output=1,
        )
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

        Forward integration steady-state solver. Finds a steady state when the
        maximum change in species concentration falls within a specified tolerance.
        Returns the steady-state solution and a error flag.
        Algorithm controls are available as mod.fintslv_<control>

        Arguments:

        initial: vector of initial concentrations

        """
        sim_time = self.__fintslv_range__ * self.__settings__['fintslv_rmult']
        # print sim_time
        Vtemp = numpy.zeros((self.__Nshape__[1]), 'd')

        def function_sim(s, t):
            return self._EvalODE(s, Vtemp)

        res, infodict = scipy.integrate.odeint(
            function_sim,
            initial.copy(),
            sim_time,
            atol=self.__settings__['lsoda_atol'],
            rtol=self.__settings__[
                'lsoda_rtol'
            ],  ##  mxstep=self.__settings__["lsoda_mxstep"],\
            mxstep=10000,
            h0=self.__settings__["lsoda_h0"],
            hmax=self.__settings__["lsoda_hmax"],
            hmin=self.__settings__["lsoda_hmin"],
            mxordn=self.__settings__["lsoda_mxordn"],
            mxords=self.__settings__["lsoda_mxords"],
            printmessg=self.__settings__["lsoda_mesg"],
            full_output=1,
        )
        if infodict['message'] == 'Integration successful.':
            status = True
        else:
            status = False

        # run through results if max(abs([x]-[x-1])) < self.__settings__['fintslv_tol'] score +1
        # if you get 5 points by seq end ... happiness
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
        NLEQ2(initial)

        PySCeS interface to the (optional) NLEQ2 algorithm. This is a powerful steady-state
        solver that can usually find a solution for when HYBRD() fails. Algorithm
        controls are available as: mod.nleq2_<control>
        Returns as steady-state solution and error flag.

        Arguments:

        initial: vector of initial species concentrations

        """
        Vtemp = numpy.zeros((self.__Nshape__[1]), 'd')
        initial0 = initial.copy()
        s_scale = initial.copy()

        N = len(initial)
        iwk = numpy.zeros((N + 52), 'i')
        rwk = numpy.zeros(((N + max(N, 10) + 15) * N + 61), 'd')
        iopt = numpy.zeros((50), 'i')

        if self.__settings__['nleq2_jacgen'] == 1:
            print(
                '(nleq2)User supplied Jacobian not supported yet ... setting __settings__[\'nleq2_jacgen\'] = 0'
            )
            self.__settings__['nleq2_jacgen'] = 0

        rtol = mach_spec.eps * 10.0 * N
        iopt[2] = self.__settings__['nleq2_jacgen']  # 2
        iopt[8] = self.__settings__['nleq2_iscaln']  # 0
        iopt[10] = 0
        iopt[11] = 6
        iopt[12] = 0
        iopt[13] = 6
        iopt[14] = 0
        iopt[15] = 6
        iopt[30] = self.__settings__['nleq2_nonlin']  # 4
        iopt[31] = self.__settings__['nleq2_qrank1']  # 1
        iopt[34] = self.__settings__['nleq2_qnscal']  # 0
        iopt[37] = self.__settings__['nleq2_ibdamp']  # 0
        iopt[38] = self.__settings__['nleq2_iormon']  # 2

        iwk[30] = self.nleq2_nitmax

        def func(s, ifail):
            s = self._EvalODE(s, Vtemp)
            if numpy.isnan(s).any():
                ifail = -1
            elif (numpy.abs(s) > 1.0e150).any():
                ifail = -1
            return s, ifail

        def jacfunc(s):
            return s

        ierr = 0

        BRETT_DEBUG_MODE = False
        ADVANCED_MODE = self.__settings__['nleq2_advanced_mode']
        if ADVANCED_MODE:
            max_iter = self.__settings__['nleq2_iter']  # 3
            max_iter_ceiling = self.__settings__['nleq2_iter_max']  # 10
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

            res, s_scale, rtol, iopt, ierr = nleq2.nleq2(
                func, jacfunc, initial, s_scale, rtol, iopt, iwk, rwk
            )

            if ierr == 0:
                # success
                GO = False
            elif ierr == 21:
                # negative rtol
                GO = False
                ##  rtol = abs(rtol)
                ##  ierr = 0
            elif ierr == 2:
                # nitmax reached
                iwk[30] = iwk[30] * self.__settings__['nleq2_growth_factor']  # 10

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
        else:
            if self.__settings__['nleq2_mesg']:
                print('(nleq2) The solution converged.')
        if ierr > 0:
            return res, False
        else:
            return res, True

    def PITCON(self, scanpar, scanpar3d=None):
        """
        PITCON(scanpar,scanpar3d=None)

        PySCeS interface to the PITCON continuation algorithm. Single parameter
        continuation has been implemented as a "scan" with the continuation
        being initialised in mod.pitcon_par_space. The second argument does not
        affect the continuation but can be used to insert a third axis parameter into
        the results. Returns an array containing the results.
        Algorithm controls are available as mod.pitcon_<control>

        Arguments:

        scanpar: the model parameter to scan (x5)
        scanpar3d [default=None]: additional output parameter for 3D plots

        """
        if self.__HAS_RATE_RULES__:
            raise NotImplementedError(
                '\nBifurcation analysis not currently available for models containing RateRules'
            )

        assert (
            type(scanpar) == str
        ), '\nscanpar must be a <string> representing a model parameter'
        modpar = list(self.__parameters__)
        try:
            a = modpar.index(scanpar)
        except:
            raise NameError(repr(scanpar) + ' is not a parameter of this model')
        if scanpar3d != None:
            if type(scanpar3d) == str:
                scanpar3d = float(scanpar3d)
                assert type(scanpar3d) == float, 'scanpar3d must be a <float>'

        par_hold = getattr(self, scanpar)

        if self.__settings__["pitcon_jac_opt"] < 1:
            self.__settings__["pitcon_jac_opt"] = 1
            print(
                '\nINFO: .__settings__["pitcon_jac_opt"] set to 1 - user defined jacobian function not yet supported'
            )

        # DONE!
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
        Vtemp = numpy.zeros((self.__Nshape__[1]), 'd')
        xr2 = numpy.zeros((parm), 'd')
        sdot = numpy.zeros((parm), 'd')
        iwork = numpy.zeros((30 + parm), 'i')
        rwork = numpy.zeros((30 + (6 * (parm)) * (parm)), 'd')
        ipar = numpy.zeros((parm), 'i')
        fpar = numpy.zeros((parm), 'd')

        res = []
        for xscan in self.pitcon_par_space:
            Vtemp[:] = 0.0
            xr2[:] = 0.0
            sdot[:] = 0.0
            ipar[:] = 0
            fpar[:] = 0.0

            iwork[0] = 0  # This is a startup
            iwork[1] = self.__settings__[
                "pitcon_init_par"
            ]  # Use X(1) for initial parameter
            iwork[2] = self.__settings__[
                "pitcon_par_opt"
            ]  # Parameterization option 0:allows program
            iwork[3] = self.__settings__[
                "pitcon_jac_upd"
            ]  # Update jacobian every newton step
            iwork[4] = self.__settings__[
                "pitcon_targ_val_idx"
            ]  # Seek target values for X(n)
            iwork[5] = self.__settings__[
                "pitcon_limit_point_idx"
            ]  # Seek limit points in X(n)
            iwork[6] = self.__settings__[
                "pitcon_output_lvl"
            ]  # Control amount of output.
            iwork[7] = 6  # Output unit 6=PC
            iwork[8] = self.__settings__[
                "pitcon_jac_opt"
            ]  # Jacobian choice. 0:supply jacobian,1:use forward difference,2:central difference
            iwork[9] = 0
            iwork[10] = 0
            iwork[11] = 0
            iwork[12] = 30
            iwork[13] = len(iwork)
            iwork[14] = 30 + (4 * parm)
            iwork[15] = len(rwork)
            iwork[16] = self.__settings__["pitcon_max_step"]  # max corrector steps
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

            rwork[0] = self.__settings__["pitcon_abs_tol"]  # Absolute error tolerance
            rwork[1] = self.__settings__["pitcon_rel_tol"]  # Relative error tolerance
            rwork[2] = self.__settings__["pitcon_min_step"]  # Minimum stepsize
            rwork[3] = self.__settings__["pitcon_max_step"]  # Maximum stepsize
            rwork[4] = self.__settings__["pitcon_start_step"]  # Starting stepsize
            rwork[5] = self.__settings__[
                "pitcon_start_dir"
            ]  # Starting direction +1.0/-1.0
            rwork[
                6
            ] = self.pitcon_targ_val  # Target value (Seek solution with iwork[4]=)
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
            rwork[19] = self.__settings__["pitcon_max_grow"]  # maximum growth factor
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

            if self.__HAS_COMPARTMENTS__ and self.__KeyWords__['Output_In_Conc']:
                ##  if self.__KeyWords__['Species_In_Conc']:
                self.__inspec__ = copy.copy(
                    self._SpeciesConcToAmount(self.state_species)
                )
            else:
                self.__inspec__ = copy.copy(self.state_species)

            go = False
            if self.__StateOK__:
                go = True
            elif not self.__StateOK__ and self.__settings__["pitcon_allow_badstate"]:
                go = True
            else:
                go = False

            if go:
                if self.__HAS_MOIETY_CONSERVATION__ == True:
                    temp = numpy.zeros((len(self.__SI__)), 'd')
                    # sI0_sim_init
                    for x in range(0, len(self.lzeromatrix_col)):
                        xr2[x] = self.__inspec__[self.nrmatrix_row[x]]
                else:
                    xr2[:-1] = copy.copy(self.state_species[:])
                xr2[-1] = copy.copy(xscan)
                for x in range(int(self.pitcon_iter)):
                    ierror, iwork, rwork, xr2 = pitcon.pitcon1(
                        fx, fpar, fx, ipar, iwork, rwork, xr2
                    )
                    if iwork[0] == 2:
                        if self.__settings__["pitcon_flux_gen"]:
                            if (
                                min(xr2[:-1]) < 0.0
                                and self.__settings__["pitcon_filter_neg"]
                            ):
                                pass
                            else:
                                xout = xr2.tolist()
                                a = xout.pop(-1)
                                if self.__HAS_MOIETY_CONSERVATION__ == True:
                                    xout = self.Fix_S_indinput(xout, amounts=True)
                                else:
                                    xout = numpy.array(xout)
                                if (
                                    self.__HAS_COMPARTMENTS__
                                    and self.__KeyWords__['Output_In_Conc']
                                ):
                                    xout = self._SpeciesAmountToConc(xout)
                                xout2 = (self.__FluxGen__(xout)).tolist()
                                xout = xout.tolist()
                                xout.insert(0, a)
                                if scanpar3d != None:
                                    xout.insert(0, scanpar3d)
                                xout2 = xout + xout2
                                res.append(xout2)
                        else:
                            if (
                                min(xr2[:-1]) < 0.0
                                and self.__settings__["pitcon_filter_neg"]
                            ):
                                pass
                            else:
                                xout = xr2.tolist()
                                a = xout.pop(-1)
                                if self.__HAS_MOIETY_CONSERVATION__ == True:
                                    xout = self.Fix_S_indinput(xout, amounts=True)
                                else:
                                    xout = numpy.array(xout)
                                if (
                                    self.__HAS_COMPARTMENTS__
                                    and self.__KeyWords__['Output_In_Conc']
                                ):
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
                        if (
                            self.__HAS_COMPARTMENTS__
                            and self.__KeyWords__['Output_In_Conc']
                        ):
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
                        if (
                            self.__HAS_COMPARTMENTS__
                            and self.__KeyWords__['Output_In_Conc']
                        ):
                            xout = self._SpeciesAmountToConc(xout)
                        xout = xout.tolist()
                        xout.insert(0, a)

                        print(xout)
                        self.pitcon_limit_points.append(xout)
                    elif iwork[0] == 1:
                        pass
                    else:
                        print(iwork[0])
                        # raw_input()
            else:
                print('\nInvalid steady state, skipping ...')

        if self.__settings__["pitcon_filter_neg_res"]:
            for result in range(len(res) - 1, -1, -1):
                if min(res[result][: len(self.species) + 1]) < 0.0:
                    res.pop(result)
        setattr(self, scanpar, par_hold)
        return numpy.array(res)

    def Simulate(self, userinit=0):
        """
        PySCeS integration driver routine that evolves the system over the time.
        Resulting array of species concentrations is stored in the **mod.data_sim** object
        Initial concentrations can be selected using *mod.__settings__['mode_sim_init']*
        (default=0):

        - 0  initialise with intial concentrations
        - 1  initialise with a very small (close to zero) value
        - 2  initialise with results of previously calculated steady state
        - 3  initialise with final point of previous simulation

        *userinit* values can be (default=0):

        - 0: initial species concentrations intitialised from (mod.S_init), time array calculated from sim_start/sim_end/sim_points
        - 1: intial species concentrations intitialised from (mod.S_init) existing "mod.sim_time" used directly
        - 2: initial species concentrations read from "mod.__inspec__", "mod.sim_time" used directly
        """
        # check if the stoichiometry has been adjusted using Stoich_nmatrix_SetValue - brett 20050719
        if not self.__StoichOK:
            self.Stoichiometry_ReAnalyse()

        # check for zero first point in user-supplied mod.sim_time, add if needed
        self._sim_time_bak = None
        if userinit != 0:
            if self.sim_time[0] != 0:
                self._sim_time_bak = copy.copy(self.sim_time)
                self.sim_time = [0.0] + list(self._sim_time_bak)

        # initialises self.__inspec__[x] with self.sXi
        if userinit == 1:
            eval(self.__mapFunc_R__)
        elif userinit == 2:
            try:
                assert len(self.__inspec__) == len(self.__species__)
            except:
                print(
                    '\nINFO: sim_sinit is the incorrect length initialising with .sX_init'
                )
                self.__inspec__ = numpy.zeros(len(self.__species__))
                eval(self.__mapFunc_R__)
        else:
            self.sim_start = float(self.sim_start)
            self.sim_end = float(self.sim_end)
            self.sim_points = int(self.sim_points)
            if self.sim_points == 1.0:
                print(
                    '*****\nWARNING: simulations require a minimum of 2 points, setting sim_points = 2.0\n*****'
                )
                self.sim_points = 2.0
            self.sim_time = numpy.linspace(
                self.sim_start, self.sim_end, self.sim_points, endpoint=1, retstep=0
            )
            eval(self.__mapFunc_R__)

        # initialise __rrule__ to mod.<rule>_init
        if self.__HAS_RATE_RULES__:
            for r in range(len(self.__rate_rules__)):
                self.__rrule__[r] = getattr(
                    self, '{}_init'.format(self.__rate_rules__[r])
                )
        for c in range(len(self.__CsizeAllIdx__)):
            cval = getattr(self, '{}_init'.format(self.__CsizeAllIdx__[c]))
            setattr(self, self.__CsizeAllIdx__[c], cval)
            self.__CsizeAll__[c] = cval
        # set initialisation array to amounts
        if self.__HAS_COMPARTMENTS__ and self.__KeyWords__['Species_In_Conc']:
            self.__inspec__ = self._SpeciesConcToAmount(self.__inspec__)

        # set mod.<species> to corrected mod.<species>_init value
        for s in range(len(self.__species__)):
            setattr(self, self.__species__[s], self.__inspec__[s])

        # This should work ok __Build_Tvec__ uses .self.__inspec__ to create Tvec
        if self.__HAS_MOIETY_CONSERVATION__ == True:
            self.__Build_Tvec__(amounts=True)
            self.showConserved(screenwrite=0)
            if self.__settings__['display_debug'] == 1:
                print(self.conserved_sums)

        # Initialise the simulation ...
        if self.__settings__['mode_sim_init'] == 0:
            # Start with s[x] at their initial values
            s0_sim_init = copy.copy(self.__inspec__)
        elif self.__settings__['mode_sim_init'] == 1:
            # Start with s[x] at zero ... well close to it anyway
            s0_sim_init = copy.copy(self.__inspec__)
            s0_sim_init[:] = self.__settings__['small_concentration']
        elif self.__settings__['mode_sim_init'] == 2:
            # Start with s[x] at the previous steady state if it exists and check if s[x] < 0.0
            # if so set that value to 1.0e-10
            if self.state_species != None:
                s0_sim_init = copy.copy(self.state_species) * 1.0
                for x in range(len(s0_sim_init)):
                    if s0_sim_init[x] < 0.0:
                        s0_sim_init[x] = self.__settings__['small_concentration']
                        print(
                            'Negative concentration detected in SimInit: s['
                            + repr(x)
                            + '] set to '
                            + repr(self.__settings__['small_concentration'])
                        )
            else:
                s0_sim_init = copy.copy(self.__inspec__)
        elif self.__settings__['mode_sim_init'] == 3:
            # Start with s[x] at the final point of the previous simulation if exists -- johann 20050220
            try:
                s0_sim_init = copy.copy(self.data_sim.species[-1])
            except:
                s0_sim_init = copy.copy(self.__inspec__)
        else:
            s0_sim_init = copy.copy(self.__inspec__)

        # print 'Sim debug'
        # print s0_sim_init
        if self.__HAS_MOIETY_CONSERVATION__ == True:
            temp = numpy.zeros((len(self.__SI__)), 'd')
            for x in range(0, len(self.lzeromatrix_col)):
                temp[x] = s0_sim_init[self.nrmatrix_row[x]]
            s0_sim_init = temp
            del temp
        # print s0_sim_init

        # ok so now we add RateRules to the initialisation_vec and see what happens
        if self.__HAS_RATE_RULES__:
            s0_sim_init = numpy.concatenate([s0_sim_init, self.__rrule__])

        # re-set self._sim recarray (otherwise self.sim is not updated)
        self._sim = None

        # real pluggable integration routines - brett 2007
        # the array copy is important ... brett leave it alone!
        Tsim0 = time.time()
        if self.mode_integrator == 'LSODA':
            sim_res, rates, simOK = self.LSODA(copy.copy(s0_sim_init))
        elif self.mode_integrator == 'CVODE':
            sim_res, rates, simOK = self.CVODE(copy.copy(s0_sim_init))
        # remove zero point from reported simulation data if necessary
        if self._sim_time_bak is not None:
            self.sim_time = copy.copy(self._sim_time_bak)
            sim_res = sim_res[1:]
            rates = rates[1:]
        Tsim1 = time.time()
        if self.__settings__['lsoda_mesg']:
            print(
                "{} time for {} points: {}".format(
                    self.mode_integrator, len(self.sim_time), Tsim1 - Tsim0
                )
            )

        if self.__HAS_RATE_RULES__:
            sim_res, rrules = numpy.split(sim_res, [len(self.__species__)], axis=1)
            print('RateRules evaluated and added to mod.data_sim.')

        self.data_sim = IntegrationDataObj()
        self.IS_VALID = simOK
        self.data_sim.setTime(self.sim_time)
        self.data_sim.setSpecies(sim_res, self.__species__)
        self.data_sim.setRates(rates, self.__reactions__)
        if self.__HAS_RATE_RULES__:
            self.data_sim.setRules(rrules, self.__rate_rules__)
        if len(self.CVODE_extra_output) > 0:
            self.data_sim.setXData(self.CVODE_xdata, lbls=self.CVODE_extra_output)
            self.CVODE_xdata = None
        if not simOK:
            print('Simulation failure')
        del sim_res

    @property
    def sim(self):
        if self._sim is None and self.data_sim is not None:
            data = self.data_sim.getAllSimData(lbls=True)
            self._sim = numpy.rec.fromrecords(data[0], names=data[1])
        return self._sim

    def State(self):
        """
        State()

        PySCeS non-linear solver driver routine. Solve for a steady state using HYBRD/NLEQ2/FINTSLV
        algorithms. Results are stored in mod.state_species and mod.state_flux. The results
        of a steady-state analysis can be viewed with the mod.showState() method.

        The solver can be initialised in 3 ways using the mode_state_init switch.
        mod.mode_state_init = 0 initialize with species initial values
        mod.mode_state_init = 1 initialize with small values
        mod.mode_state_init = 2 initialize with the final value of a 10-logstep simulation numpy.logspace(0,5,18)

        Arguments:
        None

        """
        # check if the stoichiometry has been adjusted using Stoich_nmatrix_SetValue - brett 20050719
        if not self.__StoichOK:
            self.Stoichiometry_ReAnalyse()

        self.__StateOK__ = True

        # function to feed the simulation routine if used to initialise the solver
        Vtemp = numpy.zeros((self.__Nshape__[1]), 'd')

        def function_sim(s, t):
            return self._EvalODE(s, Vtemp)

        # set self.__inspec__ to current Xi values
        eval(self.__mapFunc_R__)

        # initialise __rrule__ to mod.<rule>_init
        if self.__HAS_RATE_RULES__:
            for r in range(len(self.__rate_rules__)):
                self.__rrule__[r] = getattr(
                    self, '{}_init'.format(self.__rate_rules__[r])
                )
        # set compartment values to initial values
        for c in range(len(self.__CsizeAllIdx__)):
            cval = getattr(self, '{}_init'.format(self.__CsizeAllIdx__[c]))
            setattr(self, self.__CsizeAllIdx__[c], cval)
            self.__CsizeAll__[c] = cval
        # set initialisation array to amounts
        if self.__HAS_COMPARTMENTS__ and self.__KeyWords__['Species_In_Conc']:
            self.__inspec__ = self._SpeciesConcToAmount(self.__inspec__)
        # set mod.<species> to corrected mod.<species>_init value
        for s in range(len(self.__species__)):
            setattr(self, self.__species__[s], self.__inspec__[s])
        # This should work ok __Build_Tvec__ uses .self.__inspec__ to create Tvec
        if self.__HAS_MOIETY_CONSERVATION__ == True:
            self.__Build_Tvec__(amounts=True)
            self.showConserved(screenwrite=0)
            if self.__settings__['display_debug'] == 1:
                print(self.conserved_sums)

        # clear the solver initialisation array
        s0_ss_init = None
        # save the __settings__['hybrd_factor']
        hybrd_factor_temp = copy.copy(self.__settings__['hybrd_factor'])

        # use StateInit to initialise the solver with either 0:initval 1:zeroval 2:sim 3:1%state
        # in all cases fallback to init
        # initialising to previous ss seems problematic so I use |10%| of previous - brett
        if self.mode_state_init == 0:
            # Use initial S values
            s0_ss_init = self.__inspec__.copy()
        elif self.mode_state_init == 1:
            # Use almost zero values 1.0e-8 any smaller seems to cause a problem
            # this is user defineable option --> model.ZeroVal - brett 20040121
            s0_ss_init = self.__inspec__.copy()
            s0_ss_init[:] = self.__settings__['small_concentration']
        elif self.mode_state_init == 2:
            print('This initialisation mode has been disabled, using initial values.')
            s0_ss_init = self.__inspec__.copy()
        ##  # Perform a 10-logstep simulation numpy.logspace(0,5,18) and initialise solver with final value
        ##  # This is guesstimate ... if anyone has a better idea please let me know
        ##  # it should and be a user defineable array (min/max/len)- brett 20031215
        ##  self.__settings__['hybrd_factor'] = 10
        ##  try:
        ##  s0_ss_init = scipy.integrate.odeint(function_sim,self.__inspec__.tolist(),self.__mode_state_init2_array__,10000,full_output=0)
        ##  if self.__HAS_MOIETY_CONSERVATION__ == True:
        ##  s0_ss_init = self.Fix_S_fullinput(s0_ss_init[-1], amounts=True)
        ##  else:
        ##  s0_ss_init = s0_ss_init[-1]
        ##  except:
        ##  s0_ss_init = copy.copy(self.__inspec__)
        elif self.mode_state_init == 3:
            # Use |10%| of previous steady-state - seems to be safe, needs more testing and probably professional help
            # My rationale for this is "a set of approximate values where all variables are relatively close"
            # without having the overhead of a simulation ...
            # will eventually be a user option so that it can be +/- x% previous steadystate - brett 20031215
            # Note: hybrid doesn't seem to like to start too close to its solution more than 10% is dodgy
            if (
                self.state_species != None
                and self.__HAS_COMPARTMENTS__
                and self.__KeyWords__['Output_In_Conc']
            ):
                s0_ss_init = self._SpeciesConcToAmount(
                    abs(copy.copy(self.state_species))
                    * float(self.__settings__['mode_state_init3_factor'])
                )
            else:
                s0_ss_init = copy.copy(self.__inspec__)
        else:
            s0_ss_init = copy.copy(self.__inspec__)
        # reset __settings__['hybrd_factor']
        self.__settings__['hybrd_factor'] = hybrd_factor_temp

        # print 'State debug'
        # print s0_ss_init
        # form an initialisation vector of independent species
        if self.__HAS_MOIETY_CONSERVATION__ == True:
            temp = numpy.zeros((len(self.__SI__)), 'd')
            for x in range(0, len(self.lzeromatrix_col)):
                temp[x] = s0_ss_init[self.nrmatrix_row[x]]
            # add RateRules to the initialisation_vec and see what happens
            if self.__HAS_RATE_RULES__:
                temp = numpy.concatenate([temp, self.__rrule__])
            s0_ss_init = temp
            del temp
        # print s0_ss_init

        available_solvers = ['HYBRD']

        # set mode_solver to new syntax for backwards compatibility only
        if self.mode_solver == 0:
            self.mode_solver = 'HYBRD'
        elif self.mode_solver == 1:
            self.mode_solver = 'FINTSLV'
        elif self.mode_solver == 2:
            self.mode_solver = 'NLEQ2'

        # check for nleq2 and add if available this should be first
        if nleq2_switch == 1:
            available_solvers.append('NLEQ2')
        else:
            if self.mode_solver == 'NLEQ2':
                self.mode_solver = 'HYBRD'
            print(
                'INFO: switching to HYBRD.\nNleq2 solver not available see /nleq/readme.txt for details'
            )

        # ******* OTHER SOLVERS GO IN HERE *******

        # ******* OTHER SOLVERS GO IN HERE *******

        # shall we add HYBRD to fallback? this should be last
        if self.__settings__['mode_solver_fallback_integration']:
            available_solvers.append('FINTSLV')

        if not self.mode_solver_fallback == 1:
            assert (
                self.mode_solver in available_solvers
            ), '\nERROR: {} is not a valid ({}) solver!'.format(
                solver, str(available_solvers)
            )
            available_solvers = [self.mode_solver]

        # if solver other than HYBRD is selected move it to front of list so it is run first
        if self.mode_solver != 'HYBRD':
            available_solvers.insert(
                0, available_solvers.pop(available_solvers.index(self.mode_solver))
            )

        STATE_XOUT = False
        STATE_xdata = None
        if len(self.STATE_extra_output) > 0:
            out = []
            for d in self.STATE_extra_output:
                if (
                    hasattr(self, d)
                    and d
                    not in self.__species__ + self.__reactions__ + self.__rate_rules__
                ):
                    out.append(d)
                else:
                    print(
                        '\nWARNING: STATE is ignoring extra data ({}), it either doesn\'t exist or it\'s a species, rate or rule.\n'.format(
                            d
                        )
                    )
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
                state_species, self.__StateOK__ = self.FINTSLV(s0_ss_init.copy())
            # In case of a number (scalar) output from the solver
            if numpy.isscalar(state_species):
                state_species = numpy.array([state_species], 'd')

            # this gets everything into the current tout state
            tmp = self._EvalODE(state_species.copy(), Vtemp)
            state_flux = self.__vvec__

            if self.__HAS_RATE_RULES__:
                state_species, rrules = numpy.split(
                    state_species, [self.Nrmatrix.shape[0]]
                )
            # test for negative concentrations
            if (state_species < 0.0).any():
                self.__StateOK__ = False
                if self.__settings__['mode_state_mesg']:
                    print('WARNING!! Negative concentrations detected.')
            if self.__StateOK__:
                break
            else:
                if (
                    self.mode_solver_fallback
                    and self.__settings__['solver_switch_warning']
                ):
                    slv_idx = available_solvers.index(solver)
                    if slv_idx != len(available_solvers) - 1:
                        print(
                            'INFO: STATE is switching to {} solver.'.format(
                                available_solvers[slv_idx + 1]
                            )
                        )
                    else:
                        print('INFO: STATE calculation failed!')
        if STATE_XOUT:
            STATE_xdata = self._EvalExtraData(self.STATE_extra_output)

        # the current status quo is all state algorithms will use and produce results
        # in amounts and autoconversion to and from species takes place in the calling
        # algorithms -- brett2008
        # THIS IS OPPOSITE TO THE SIMULATE METHOD brett - again

        if self.__HAS_MOIETY_CONSERVATION__ == True:
            state_species = self.Fix_S_indinput(state_species, amounts=True)
        if self.__HAS_COMPARTMENTS__ and self.__KeyWords__['Output_In_Conc']:
            self.state_species = self._SpeciesAmountToConc(state_species.copy())
        else:
            self.state_species = state_species

        self.state_flux = state_flux

        del state_species, state_flux

        # final check for a bad state set check if fluxes are == mach_eps
        # this is almost never going to be true
        if (self.state_flux == self.__settings__['mach_floateps']).any():
            print('\nWARNING: extremely small flux detected! proceed with caution:')
            print(self.state_flux)
            ##  self.__StateOK__ = False

        if not self.__StateOK__:
            print(
                '\n***\nWARNING: invalid steady state solution (species concentrations and fluxes)\n***\n'
            )
            if self.__settings__['mode_state_nan_on_fail']:
                self.state_species[:] = numpy.NaN
                self.state_flux[:] = numpy.NaN

        # set the instance steady state flux and species attributes
        self.SetStateSymb(self.state_flux, self.state_species)

        # coming soon to a terminal near you
        self.data_sstate = StateDataObj()
        self.data_sstate.setSpecies(self.state_species, self.__species__)
        self.data_sstate.setFluxes(self.state_flux, self.__reactions__)
        if self.__HAS_RATE_RULES__:
            self.data_sstate.setRules(rrules, self.__rate_rules__)
        if STATE_XOUT:
            self.data_sstate.setXData(STATE_xdata, lbls=self.STATE_extra_output)
            del STATE_xdata

        self.data_sstate.IS_VALID = self.__StateOK__

        if self.__settings__['display_debug'] == 1:
            print('self.state_species')
            print(self.state_species)
            print('self.state_flux')
            print(self.state_flux)

    # driver routines that support the core routines
    #   core support

    # takes the flux and species arrays and
    # assigns them as instances variable self.Jx and self.Xss
    # I'm not sure if anyone wants to use this in real life so it might be hidden at some point - brett 20040122
    # gone - brett 20040506 ... back brett 20040720
    def SetStateSymb(self, flux, metab):
        """
        SetStateSymb(flux,metab)

        Sets the individual steady-state flux and concentration attributes as
        mod.J_<reaction> and mod.<species>_ss

        Arguments:

        flux: the steady-state flux array
        metab: the steady-state concentration array

        """
        for x in range(0, len(self.state_species)):
            setattr(self, self.__species__[x] + '_ss', metab[x])

        for x in range(0, len(self.state_flux)):
            setattr(self, 'J_' + self.__reactions__[x], flux[x])

    # uses __inspec__ to build Tvec
    def __Build_Tvec__(self, amounts=True):
        """
        __Build_Tvec_(concs=False)

        Creates vectors of conserved moiety totals from __inspec__
        self.__tvec_a__ = amounts
        self.__tvec_c__ = concentrations

        Arguments:
        None

        """
        if self.__HAS_MOIETY_CONSERVATION__ == True:
            if amounts:
                self.__tvec_a__ = numpy.add.reduce(
                    copy.copy(self.__reordered_lcons) * self.__inspec__, 1
                )
                self.__tvec_c__ = numpy.add.reduce(
                    copy.copy(self.__reordered_lcons)
                    * self._SpeciesAmountToConc(self.__inspec__),
                    1,
                )
            else:
                self.__tvec_c__ = numpy.add.reduce(
                    copy.copy(self.__reordered_lcons) * self.__inspec__, 1
                )
                self.__tvec_a__ = numpy.add.reduce(
                    copy.copy(self.__reordered_lcons)
                    * self._SpeciesConcToAmount(self.__inspec__),
                    1,
                )

    def showConserved(self, File=None, screenwrite=1, fmt='%2.3f'):
        """
        showConserved(File=None,screenwrite=1,fmt='%2.3f')

        Print the moiety conserved cycles present in the system.

        Arguments:

        File [default=None]: an open writable Python file object
        screenwrite [default=1]: write results to console (0 means no reponse)
        fmt [default='%2.3f']: the output number format string

        """
        if self.__HAS_MOIETY_CONSERVATION__ == True:
            Tlist = list(range(0, len(self.__tvec_a__)))
            if Tlist != []:
                ConSumPstr = ''
                for x in range(0, len(Tlist)):
                    for y in range(0, len(self.conservation_matrix_col)):
                        if self.conservation_matrix[Tlist[x], y] > 0.0:
                            # coeff = self.__settings__['mode_number_format'] % (s_init[y]*abs(conservation_matrix[Tlist[x],y]))
                            coeff = fmt % abs(self.conservation_matrix[Tlist[x], y])
                            met = self._D_s_Order[self.conservation_matrix_col[y]]
                            ConSumPstr += (
                                ' + {' + coeff + '}' + met.replace('self.', '')
                            )
                        elif self.conservation_matrix[Tlist[x], y] < 0.0:
                            # coeff = self.__settings__['mode_number_format'] % (s_init[y]*abs(conservation_matrix[Tlist[x],y]))
                            coeff = fmt % abs(self.conservation_matrix[Tlist[x], y])
                            met = self._D_s_Order[self.conservation_matrix_col[y]]
                            ConSumPstr += (
                                ' - {' + coeff + '}' + met.replace('self.', '')
                            )
                    if (
                        self.__HAS_COMPARTMENTS__
                        and self.__KeyWords__['Output_In_Conc']
                    ):
                        ConSumPstr += (
                            ' = '
                            + self.__settings__['mode_number_format']
                            % self.__tvec_c__[Tlist[x]]
                            + '\n'
                        )
                    else:
                        ConSumPstr += (
                            ' = '
                            + self.__settings__['mode_number_format']
                            % self.__tvec_a__[Tlist[x]]
                            + '\n'
                        )
            self.conserved_sums = ConSumPstr
        else:
            self.conserved_sums = 'No moiety conservation'

        if File != None:
            print('\nConserved relationships')
            # assert type(File) == file, 'showConserved() needs an open file object'
            File.write('\n## Conserved relationships\n')
            File.write(self.conserved_sums)
        elif screenwrite:
            print('\nConserved relationships')
            print(self.conserved_sums)

    def showFluxRelationships(self, File=None):
        """
        showConserved(File=None)

        Print the flux relationships present in the system.

        Arguments:

        File [default=None]: an open writable Python file object

        """
        Ostr = ''
        for row in range(self.__kzeromatrix__.shape[0]):
            Ostr += "{} =".format(self.reactions[self.kzeromatrix_row[row]])
            for col in range(self.__kzeromatrix__.shape[1]):
                if self.__kzeromatrix__[row, col] != 0.0:
                    if self.__kzeromatrix__[row, col] > 0.0:
                        Ostr += " + {%2.2f}%s" % (
                            abs(self.__kzeromatrix__[row, col]),
                            self.reactions[self.kzeromatrix_col[col]],
                        )
                    else:
                        Ostr += " - {%2.2f}%s" % (
                            abs(self.__kzeromatrix__[row, col]),
                            self.reactions[self.kzeromatrix_col[col]],
                        )
            Ostr += '\n'

        if File != None:
            print('\nFlux relationships')
            # assert type(File) == file, 'showConserved() needs an open file object'
            File.write('\n## Flux relationships\n')
            File.write(Ostr)
        else:
            print('\nFlux relationships')
            print(Ostr)

    # Calculate dependant variables done directly in EvalREq2 this is a utility version
    def Fix_S_fullinput(self, s_vec, amounts=True):
        """
        Fix_S_fullinput(s_vec)

        Using the full concentration vector evaluate the dependent species

        Arguments:

        s_vec: a full length concentration vector

        """
        # s_vec = copy.copy(s)
        for x in range(0, len(self.lzeromatrix_col)):
            self.__SI__[x] = s_vec[self.lzeromatrix_col[x]]

        for x in range(0, len(self.lzeromatrix_row)):
            if amounts:
                s_vec[self.lzeromatrix_row[x]] = self.__tvec_a__[x] + numpy.add.reduce(
                    self.__lzeromatrix__[x, :] * self.__SI__
                )  # there might be a way to reduce the 2 for loops
            else:
                s_vec[self.lzeromatrix_row[x]] = self.__tvec_c__[x] + numpy.add.reduce(
                    self.__lzeromatrix__[x, :] * self.__SI__
                )  # there might be a way to reduce the 2 for loops

        return s_vec

    # Calculate dependant variables done directly in EvalREq2B this is a utility version
    def Fix_S_indinput(self, s_vec, amounts=True):
        """
        Fix_S_indinput(s_vec, amounts=True)
        whether to use self.__tvec_a__ (default)
        or self.__tvec_c__

        Given a vector of independent species evaluate and return a full concentration vector.

        Arguments:

        s_vec: vector of independent species

        """
        # stick SI into s using Nrrow order
        for x in range(len(s_vec)):
            self.__SALL__[self.nrmatrix_row[x]] = s_vec[x]
        # stick SD into s using lzerorow order
        for x in range(len(self.lzeromatrix_row)):
            if amounts:
                self.__SALL__[self.lzeromatrix_row[x]] = self.__tvec_a__[
                    x
                ] + numpy.add.reduce(self.__lzeromatrix__[x, :] * s_vec)
            else:
                self.__SALL__[self.lzeromatrix_row[x]] = self.__tvec_c__[
                    x
                ] + numpy.add.reduce(self.__lzeromatrix__[x, :] * s_vec)
        return copy.copy(self.__SALL__)

    #   output support routines

    # Quick and dirty flux regeneration uses the Rx form of the RE's to cater for potential conservation
    # caters for both vectors and arrays and is conservation aware
    def __FluxGen__(self, s):
        """
        **Deprecating** only used by **PITCON**

        s: species vector/array
        """
        Vtemp = numpy.zeros((self.__Nshape__[1]), 'd')
        try:
            s.shape[0]
            Vout = numpy.zeros((len(s), len(self.__vvec__)), 'd')
            for x in range(len(s)):
                if (
                    self.__HAS_MOIETY_CONSERVATION__ == True
                ):  # depending if there is conservation -- N.L_switch is: 0 for none or 1 for present
                    Vout[x, :] = self._EvalREq2_alt(s[x, :], Vtemp)
                elif self.__HAS_MOIETY_CONSERVATION__ == False:
                    Vout[x, :] = self._EvalREq(s[x, :], Vtemp)
        except:
            if (
                self.__HAS_MOIETY_CONSERVATION__ == True
            ):  # depending if there is conservation -- N.L_switch is: 0 for none or 1 for present
                Vout = self._EvalREq2_alt(s, Vtemp)
            elif self.__HAS_MOIETY_CONSERVATION__ == False:
                Vout = self._EvalREq(s, Vtemp)
        return Vout

    def FluxGenSim(self, s):
        """
        **Deprecated**
        """
        pass

    def ParGenSim(self):
        """
        **Deprecated**
        """
        pass

    def Fix_Sim(self, metab, flux=0, par=0):
        """
        **Deprecated**
        """
        pass

    # MCA routines
    #   elasticity routines

    def EvalEvar(self, input=None, input2=None):
        """
        EvalEvar(input=None,input2=None)

        Calculate reaction elasticities towards the variable species.

        Both inputs (input1=species,input2=rates) should be valid (steady state for MCA) solutions and given in the correct order for them to be used. If either or both are missing the last state values are used automatically.
        Elasticities are scaled using input 1 and 2.

        Arguments::

         - input [default=None]:  species concentration vector
         - input2 [default=None]: reaction rate vector

        Settings, set in mod.__settings__::

         - elas_evar_upsymb   [default = 1] attach individual elasticity symbols to model instance
         - elas_zero_conc_fix [default=False] if zero concentrations are detected in a steady-state solution make it a very small number
         - elas_zero_flux_fix [default=False] if zero fluxes are detected in a steady-state solution make it a very small number
         - elas_scaling_div0_fix [default=False] if INf's are detected after scaling set to zero

        """
        if input == None or input2 == None:
            input = self.state_species
            input2 = self.state_flux
            # print 'INFO: Using state_species and state_flux as input'
        else:
            assert len(input) == len(self.species), (
                'length error this array must have '
                + str(len(self.species))
                + ' elements'
            )
            assert len(input2) == len(self.reactions), (
                'length error this array must have '
                + str(len(self.reactions))
                + ' elements'
            )
            # not really necessary but in here just in case - brett 20040930
            # map input back to self.Sx
            exec(self.__CDvarUpString__)

        if self.__settings__['elas_zero_flux_fix']:
            val = self.__settings__['mach_floateps'] ** 2
            for x in range(len(input2)):
                if abs(input2[x]) <= val:
                    print(
                        'Info: zero flux detected: J_{} set to {}'.format(
                            self.reactions[x], val
                        )
                    )
                    input2[x] = val
        if self.__settings__['elas_zero_conc_fix']:
            val = self.__settings__['mach_floateps'] ** 2
            for x in range(len(input)):
                if abs(input[x]) <= val:
                    print(
                        'Info: zero concentration detected: {} set to {}'.format(
                            self.species[x], val
                        )
                    )
                    input[x] = val

        if self.__settings__['display_debug'] == 1:
            print('\nVarinput = ' + repr(input) + '\n')

        self.__evmatrix__ = numpy.zeros(
            (len(self.__reactions__), len(self.__species__)), 'd'
        )

        # scale KL (future refactoring allow user input)
        self.ScaleKL(input, input2)

        # create tuples of the rows and columns for the Ev matrix
        self.elas_var_row = tuple(copy.copy(self.__reactions__))
        col = copy.copy(self._D_s_Order)
        for x in range(len(self._D_s_Order)):
            col[x] = col[x].replace('self.', '')
        self.elas_var_col = tuple(col)
        del col

        # attempt to evaluate every variable against every flux: E's > mach_eps are assumed to exist
        for react in range(len(self.__reactions__)):
            if self.__settings__['display_debug'] == 1:
                print('\nReaction: ' + self.__reactions__[react])
            for met in range(len(self.__species__)):
                countV = 0
                countMet = 0
                # this modification should make this independant of a steady state - finally implimented july2004
                # eval('self.'+ self.__species__[met] + '_ss'),order=self.__settings__['mode_elas_deriv_order'],dx=hstep,n=1,\
                try:
                    """
                    I borrowed the idea of scaling the stepsize to So from Herbert Sauro's Jarnac TModel.uEEOp function
                    brett - 2004-08-18
                    """

                    hstep = (
                        input[met] * self.__settings__['mode_elas_deriv_factor']
                    )  # self.__settings__['mode_elas_deriv_factor'] = 0.0001
                    if (
                        abs(hstep) < self.__settings__['mode_elas_deriv_min']
                    ):  # self.__settings__['mode_elas_deriv_min'] = 1.0e-12
                        hstep = self.__settings__['mode_elas_deriv_min']

                    a = ScipyDerivative(
                        self.__num_deriv_function__,
                        input[met],
                        order=self.__settings__['mode_elas_deriv_order'],
                        dx=hstep,
                        n=1,
                        args=(self.__reactions__[react], self.__species__[met]),
                    )
                except Exception as ex:
                    print(ex)
                    print(
                        '\nINFO: Elasticity evaluation failure in ',
                        self.__reactions__[react],
                        self.__species__[met],
                    )
                    print('Elasticity has been set to zero')
                    print(
                        'A stepsize that is too large might cause this ... try decreasing the factor and or min stepsize'
                    )
                    print(
                        'Keep in mind machine precision, is',
                        self.__settings__['mach_floateps'],
                        ' and if min stepsize',
                    )
                    print('becomes too small numeric error can become significant')
                    a = 0.0
                if abs(a) >= self.__settings__['mach_floateps']:
                    if self.__settings__['display_debug'] == 1:
                        print('species: ' + self.__species__[met])
                        print(
                            '--> d('
                            + self.__reactions__[react]
                            + ')d('
                            + self.__species__[met]
                            + ') = '
                            + str(a)
                        )
                    self.__evmatrix__[react, met] = a

        # restore variables to ss values only if steady state used
        if self.__settings__['elas_evar_remap']:
            eval(compile(self.__remaps, '__remaps', 'exec'))
            # print "\nREMAPPING\n"

        self.elas_var_u = self.__evmatrix__

        # If scaled mca is requested
        if self.__settings__['mode_mca_scaled'] == 1:
            Ds = numpy.zeros((len(input), len(input)), 'd')
            Dj = numpy.zeros((len(input2), len(input2)), 'd')

            for x in range(0, len(input)):
                Ds[x, x] = input[x]
            # print Ds

            for x in range(0, len(input2)):
                Dj[x, x] = 1.0 / input2[x]
                if self.__settings__['elas_scaling_div0_fix'] and numpy.isinf(Dj[x, x]):
                    print(
                        'Infinite elasticity detected during scaling setting to zero ({})'.format(
                            self.reactions[x]
                        )
                    )
                    Dj[x, x] = 1.0e-16
            # print Dj

            Dj_e = numpy.dot(Dj, self.__evmatrix__)
            self.__evmatrix__ = numpy.dot(Dj_e, Ds)
            self.elas_var = self.__evmatrix__
        else:
            self.elas_var = None

        if self.__settings__["elas_evar_upsymb"] == 1:
            # Upload elasticity symbols into namespace
            r, c = self.__evmatrix__.shape
            output2 = ''
            for x in range(0, r):
                react = self.__reactions__[x]
                for y in range(0, c):
                    met = self._D_s_Order[y]
                    if self.__settings__['mode_mca_scaled']:
                        setattr(
                            self,
                            'ec' + react + '_' + met.replace('self.', ''),
                            self.__evmatrix__[x, y],
                        )
                    else:
                        setattr(
                            self,
                            'uec' + react + '_' + met.replace('self.', ''),
                            self.__evmatrix__[x, y],
                        )

            # the new way of doing things (test phase) brett - 2007
            self.ec = BagOfStuff(
                self.__evmatrix__, self.elas_var_row, self.elas_var_col
            )
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
            # print 'INFO: variable elasticity symbols not attached - .__settings__["elas_evar_upsymb"] = ' + `self.__settings__["elas_evar_upsymb"]`

        if self.__settings__['display_debug'] == 1:
            print('\ne_vmatrix')
            print(repr(self._D_s_Order).replace('self.', ''))
            print(self.__reactions__)
            print(self.__evmatrix__)

    def CleanNaNsFromArray(self, arr, replace_val=0.0):
        """
        Scan a matrix for NaN's and replace with zeros:

         - *arr* the array to be cleaned

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
        EvalEpar(input=None,input2=None)

        Calculate reaction elasticities towards the parameters.

        Both inputs (input1=species,input2=rates) should be valid (steady state for MCA) solutions and given in the correct order for them to be used. If either or both are missing the last state values are used automatically. Elasticities are scaled using input 1 and 2.

        Arguments:

         - input [default=None]: species concentration vector
         - input2 [default=None]: reaction rate vector

        Settings, set in mod.__settings__::

         - elas_epar_upsymb   [default = 1] attach individual elasticity symbols to model instance
         - elas_scaling_div0_fix [default=False] if NaN's are detected in the variable and parameter elasticity matrix replace with zero

        """
        if input == None or input2 == None:
            input = self.state_species
            input2 = self.state_flux
        else:
            assert len(input) == len(self.species), (
                'length error this array must have '
                + str(len(self.species))
                + ' elements'
            )
            assert len(input2) == len(self.reactions), (
                'length error this array must have '
                + str(len(self.reactions))
                + ' elements'
            )
            # put in to fix the juicy bug - brett 20040930
            # iow automatic derivatives use the current values of self.Sx to operate
            exec(self.__CDvarUpString__)

        # create parameter holding array
        parVal_hold = numpy.zeros((len(self.__parameters__)), 'd')
        if self.__settings__['display_debug'] == 1:
            print('\nParinput = ' + repr(input) + '\n')
            print('\nparVal_hold1')
            print(parVal_hold)

        # Store parameter values into the storage array and copy them to the working array (parVal2)
        exec(self.__par_map2storeC)
        parVal2 = copy.copy(parVal_hold)

        # create the matrix of parameter elasticities
        self.__epmatrix__ = numpy.zeros(
            (len(self.__reactions__), len(self.__parameters__)), 'd'
        )

        # create tuples of the rows and columns for the Ep matrix
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

                    hstep = (
                        getattr(self, self.__parameters__[par])
                        * self.__settings__['mode_elas_deriv_factor']
                    )  # self.__settings__['mode_elas_deriv_factor'] = 0.0001
                    if (
                        abs(hstep) < self.__settings__['mode_elas_deriv_min']
                    ):  # self.__settings__['mode_elas_deriv_min'] = 1.0e-12
                        hstep = self.__settings__['mode_elas_deriv_min']

                    a = ScipyDerivative(
                        self.__num_deriv_function__,
                        getattr(self, self.__parameters__[par]),
                        order=self.__settings__['mode_elas_deriv_order'],
                        dx=hstep,
                        n=1,
                        args=(self.__reactions__[react], self.__parameters__[par]),
                    )
                except Exception as ex:
                    print(
                        '\nNumeric derivative evaluation failure in ',
                        self.__reactions__[react],
                        self.__parameters__[par],
                    )
                    print('Elasticity has been set to NaN')
                    print(
                        'A stepsize that is too large might cause this ... try decreasing the factor and or min stepsize'
                    )
                    print(
                        'Keep in mind machine precision, is',
                        self.__settings__['mach_floateps'],
                        ' and if min stepsize',
                    )
                    print('becomes too small numeric error can become significant')
                    print(ex)
                    a = numpy.NaN

                if numpy.isnan(a) or abs(a) > self.__settings__['mach_floateps']:
                    if self.__settings__['display_debug'] == 1:
                        print('parameter: ' + self.__parameters__[par])
                        print(
                            '--> d('
                            + self.__reactions__[react]
                            + ')d('
                            + self.__parameters__[par]
                            + ') = '
                            + repr(a)
                        )
                    self.__epmatrix__[react, par] = a

        self.elas_par_u = self.__epmatrix__

        # Retrieve parameter values from the storage array
        exec(self.__par_remapC)

        if self.__settings__['display_debug'] == 1:
            print('\nparVal_hold2')
            print(parVal_hold)

        # Parameters are scaled by [1/J]*[Ep]*[P]
        # If scaled mca is requested scale self.__epmatrix__
        if self.__settings__['mode_mca_scaled'] == 1:
            Dp = numpy.zeros((len(self.__parameters__), len(self.__parameters__)), 'd')
            Dj = numpy.zeros((len(input2), len(input2)), 'd')

            for x in range(0, len(self.__parameters__)):
                Dp[x, x] = parVal_hold[x]
            # print Dp

            for x in range(0, len(input2)):
                Dj[x, x] = 1.0 / input2[x]
                if self.__settings__['elas_scaling_div0_fix'] and numpy.isinf(Dj[x, x]):
                    print(
                        'Infinite elasticity detected during scaling setting to zero ({})'.format(
                            self.reactions[x]
                        )
                    )
                    Dj[x, x] = 1.0e-16
            # print Dj

            Dj_e = numpy.dot(Dj, self.__epmatrix__)
            self.__epmatrix__ = numpy.dot(Dj_e, Dp)
            self.elas_par = self.__epmatrix__
        else:
            self.elas_par = None

        del parVal_hold

        if self.__settings__["elas_epar_upsymb"] == 1:
            # Upload elasticity symbols into namespace
            r, c = self.__epmatrix__.shape
            output2 = ''
            for x in range(0, r):
                react = self.__reactions__[x]
                for y in range(0, c):
                    met = self.__parameters__[y]
                    if self.__settings__['mode_mca_scaled']:
                        setattr(
                            self,
                            'ec' + react + '_' + met.replace('self.', ''),
                            self.__epmatrix__[x, y],
                        )
                    else:
                        setattr(
                            self,
                            'uec' + react + '_' + met.replace('self.', ''),
                            self.__epmatrix__[x, y],
                        )

            # the new way of doing things (test phase) brett - 2007
            self.ecp = BagOfStuff(
                self.__epmatrix__, self.elas_par_row, self.elas_par_col
            )
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
        showEvar(File=None)

        Write out all variable elasticities as \'LaTeX\' formatted strings, alternatively write results to a file.

        Arguments:

        File [default=None]: an open writable Python file object

        """
        # The Variable elasticities
        try:
            r, c = self.__evmatrix__.shape
            evar_output = ''
            for x in range(0, r):
                react = self.__reactions__[x]
                evar_output += '\n' + repr(self.__reactions__[x]) + '\n'
                for y in range(0, c):
                    rtemp = (
                        self.__settings__['mode_number_format']
                        % self.__evmatrix__[x, y]
                    )
                    met = self._D_s_Order[y]
                    if self.__evmatrix__[x, y] != 0.0:
                        if self.__settings__['mode_mca_scaled']:
                            elas = '\\ec{' + react + '}{' + met + '} = ' + rtemp
                        else:
                            elas = '\\uec{' + react + '}{' + met + '} = ' + rtemp
                        evar_output += elas + '\n'
            evar_output = evar_output.replace('self.', '')
        except:
            evar_output = 'No variable elasticities - run EvalEvar() to calculate'

        if File != None:
            print('\nspecies elasticities')
            File.write('\n## species elasticities\n')
            File.write(evar_output)
        else:
            print('\nspecies elasticities')
            print(evar_output)

    def showEpar(self, File=None):
        """
        showEpar(File=None)

        Write out all nonzero parameter elasticities as \'LaTeX\' formatted strings, alternatively write to file.

        Arguments:

        File [default=None]: an open writable Python file object

        """
        # The parameter elasticities
        try:
            r, c = self.__epmatrix__.shape
            epar_output = ''
            for x in range(0, r):
                react = self.__reactions__[x]
                epar_output += '\n' + repr(self.__reactions__[x]) + '\n'
                for y in range(0, c):
                    rtemp = (
                        self.__settings__['mode_number_format']
                        % self.__epmatrix__[x, y]
                    )
                    met = self.__parameters__[y]
                    # brett (paranoia removal) October 2004
                    # if abs(self.__epmatrix__[x,y]) >= 1.e-15:
                    if self.__settings__['mode_mca_scaled']:
                        elas = '\\ec{' + react + '}{' + met + '} = ' + rtemp
                    else:
                        elas = '\\uec{' + react + '}{' + met + '} = ' + rtemp
                        # print elas
                    epar_output += elas + '\n'
            epar_output = epar_output.replace('self.', '')
        except:
            epar_output = 'No parameter elasticities - run EvalEpar() to calculate'

        if File != None:
            print('\nParameter elasticities')
            File.write('\n## Parameter elasticities\n')
            File.write(epar_output)
        else:
            print('\nParameter elasticities')
            print(epar_output)

    def showElas(self, File=None):
        """
        showElas(File=None)

        Print all elasticities to screen or file as \'LaTeX\' compatible strings.
        Calls showEvar() and showEpar()

        Arguments:

        File [default=None]: an open writable Python file object

        """
        if File != None:
            self.showEvar(File)
            self.showEpar(File)
        else:
            self.showEvar()
            self.showEpar()

    def ScaleKL(self, input, input2):
        """
        ScaleKL(input,input2)

        Scale the K and L matrices with current steady state (if either input1 or 2 == None) or user input.

        Arguments:

        input: vector of species concentrations
        input2: vector of reaction rates

        """
        # called by EvalEvar
        # Scale L matrix
        lmat = copy.copy(self.__lmatrix__)

        if type(input) == type(None) or type(input2) == type(None):
            input = self.state_species
            input2 = self.state_flux
            # print 'INFO: Using state_species and state_flux as input'
        else:
            assert len(input) == len(self.species), (
                'length error this array must have '
                + str(len(self.species))
                + ' elements'
            )
            assert len(input2) == len(self.reactions), (
                'length error this array must have '
                + str(len(self.reactions))
                + ' elements'
            )

        s_1_scale = numpy.zeros((len(input), len(input)), 'd')

        for x in range(0, len(input)):
            try:
                s_1_scale[x, x] = (
                    1.0 / input[self.lmatrix_row[x]]
                )  # create 1/D Using the steady-state met's ordered to the rows
            except:
                print(
                    'Zero species detected: '
                    + self.__species__[self.lmatrix_row[x]]
                    + '_ss = '
                    + repr(input[self.lmatrix_row[x]])
                )
                s_1_scale[x, x] = 0.0

        Si_order = numpy.zeros(len(self.lmatrix_col))  # check
        Si_scale = numpy.zeros((len(self.lmatrix_col), len(self.lmatrix_col)), 'd')

        for x in range(0, len(self.lmatrix_col)):
            Si_order[x] = self.lmatrix_col[x]  # check
            Si_scale[x, x] = input[self.lmatrix_col[x]]

        L_Si = numpy.dot(lmat, Si_scale)
        L_scaled = numpy.dot(s_1_scale, L_Si)

        self.lmatrix_scaled = L_scaled

        # Scale K matrix
        kmat = copy.copy(self.__kmatrix__)

        j_1_scale = numpy.zeros((len(input2), len(input2)), 'd')

        for x in range(0, len(input2)):
            try:
                j_1_scale[x, x] = 1.0 / input2[self.kmatrix_row[x]]
            except:
                print(
                    '\nNull flux detected: '
                    + self.__reactions__[self.kmatrix_row[x]]
                    + '_ss = '
                    + repr(input2[self.kmatrix_row[x]])
                )
                j_1_scale[x, x] = 0.0

        Ji_order = numpy.zeros(len(self.kmatrix_col))  # check
        Ji_scale = numpy.zeros((len(self.kmatrix_col), len(self.kmatrix_col)), 'd')

        for x in range(0, len(self.kmatrix_col)):
            Ji_order[x] = self.kmatrix_col[x]  # check
            Ji_scale[x, x] = input2[self.kmatrix_col[x]]

        K_Ji = numpy.dot(kmat, Ji_scale)
        K_scaled = numpy.dot(j_1_scale, K_Ji)

        self.kmatrix_scaled = K_scaled

    def __num_deriv_function__(self, x, react, met):
        """
        __num_deriv_function__(x,react,met)

        System function that evaluates the rate equation, used by numeric perturbation methods to derive
        elasticities. It uses a specific format assigning x to met and evaluating react for v and is tailored for the ScipyDerivative() function using precompiled function strings.

        Arguments:

        x: value to assign to met
        react: reaction
        met: species

        """
        # exec(met)
        # self.Forcing_Function()
        # exec(react)
        # return v
        bk = getattr(self, met)  # backup current metabolite value
        setattr(self, met, x)
        self.Forcing_Function()
        v = getattr(self, react)
        vout = v()
        setattr(self, met, bk)  # reset metabolite value
        return vout

    #   Control coefficient routines
    def EvalCC(self):
        """
        EvalCC()

        Calculate the MCA control coefficients using the current steady-state solution.

        mod.__settings__["mca_ccj_upsymb"] = 1  attach the flux control coefficients to the model instance
        mod.__settings__["mca_ccs_upsymb"] = 1 attach the concentration control coefficients to the model instance

        Arguments:
        None

        """
        # sort E to match K and L
        e2 = numpy.zeros((len(self.kmatrix_row), len(self.lmatrix_row)), 'd')
        # print e2

        for x in range(0, len(self.kmatrix_row)):
            for y in range(0, len(self.lmatrix_row)):
                e2[x, y] = self.elas_var_u[self.kmatrix_row[x], self.lmatrix_row[y]]

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
            # fortran has a very fp idea of zero I use abs(val)<1.0e-15
            self.__structural__.MatrixFloatFix(Ci_u, val=1.0e-15)
            go = 1
        except Exception as ex:
            print(ex)
            print(
                '\nINFO: K-EL matrix inversion failed this is possibly due to NaN values in the Elasticity matrix'
            )
            print(
                'NaN elasticities can be caused by zero fluxes or concentrations look at the settings (in mod.__settings__)'
            )
            print('\'elas_scaling_div0_fix\' and \'elas_zero_flux_fix\'')
            go = 0

        if go == 1 and not self.__settings__['mode_mca_scaled']:
            self.mca_ci = Ci_u
            self.__mca_CCi_unscaled = Ci_u
            del Ci_u
        elif go == 1 and self.__settings__['mode_mca_scaled']:
            Vi_l = self.__kmatrix__.shape[1]  # Ind fluxes
            Vd_l = abs(
                self.__kmatrix__.shape[0] - self.__kmatrix__.shape[1]
            )  # Can never be the case but b.paranoid
            Si_l = self.__lmatrix__.shape[1]  # Ind species
            Sd_l = abs(
                self.__lmatrix__.shape[0] - self.__lmatrix__.shape[1]
            )  # Can never be the case but b.paranoid

            scale1_l = Vi_l + Si_l
            scale1 = numpy.zeros((scale1_l, scale1_l), 'd')

            for x in range(0, Vi_l):
                scale1[x, x] = 1.0 / self.state_flux[self.kmatrix_row[x]]
            for x in range(Vi_l, scale1_l):
                scale1[x, x] = 1.0 / self.state_species[self.lmatrix_row[x - Vi_l]]

            scale2 = numpy.zeros((Vi_l + Vd_l, Vi_l + Vd_l), 'd')

            for x in range(0, len(self.state_flux)):
                scale2[x, x] = self.state_flux[self.kmatrix_row[x]]

            left = numpy.dot(scale1, Ci_u)
            Ci = numpy.dot(left, scale2)

            self.mca_ci = Ci
            del scale1_l, scale1, scale2, left, Ci, Ci_u

        Cirow = []
        Cicol = []

        # the wierdness continues
        for x in range(0, len(self.kmatrix_row)):
            Cicol.append(self.__reactions__[self.kmatrix_row[x]])

        for x in range(0, len(self.kmatrix_col)):
            Cirow.append(self.__reactions__[self.kmatrix_col[x]])

        for x in range(0, len(self.lmatrix_col)):
            Cirow.append(self._D_s_Order[self.lmatrix_col[x]].replace('self.', ''))

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
        CSi = numpy.zeros(
            (self.mca_ci.shape[0] - len(self.kmatrix_col), self.mca_ci.shape[1]), 'd'
        )

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

        # we need the scaled k/l matrices to calculate the dependent cc's
        # self.ScaleKL()

        for x in range(self.__kzeromatrix__.shape[0] - 1, -1, -1):
            Ko_row.append(self.__reactions__[self.kzeromatrix_row[x]])
            for y in range(
                self.__kzeromatrix__.shape[1] - 1, -1, -1
            ):  # this can be replaced by a row slice operation
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

        # new
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
            Lo_row.append(
                self._D_s_Order[self.lzeromatrix_row[x]].replace('self.', '')
            )
            for y in range(
                self.__lzeromatrix__.shape[1] - 1, -1, -1
            ):  # this can be replaced by a row slice operation
                sLo[x, y] = self.lmatrix_scaled[x + xFactor, y]
                if x == 0:
                    Lo_col.append(self._D_s_Order[self.lzeromatrix_col[y]])

        Lo_row.reverse()
        Lo_col.reverse()

        if (
            self.__HAS_MOIETY_CONSERVATION__ == True
        ):  # Only do this if there is S dependency
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
            self.cc_flux = numpy.concatenate(
                (self.mca_ci[: self.__kmatrix__.shape[1], :], self.mca_cjd)
            )
            self.cc_flux_row = copy.copy(
                self.mca_ci_row[: self.__kmatrix__.shape[1]] + self.mca_cjd_row
            )
            self.cc_flux_col = copy.copy(self.mca_ci_col)
            self.cc_conc = numpy.concatenate(
                (self.mca_ci[self.__kmatrix__.shape[1] :, :], self.mca_csd)
            )
            self.cc_conc_row = copy.copy(
                self.mca_ci_row[self.__kmatrix__.shape[1] :] + self.mca_csd_row
            )
            self.cc_conc_col = copy.copy(self.mca_ci_col)
        elif self.__HAS_FLUX_CONSERVATION__:
            self.cc_flux = numpy.concatenate(
                (self.mca_ci[: self.__kmatrix__.shape[1], :], self.mca_cjd)
            )
            self.cc_flux_row = copy.copy(
                self.mca_ci_row[: self.__kmatrix__.shape[1]] + self.mca_cjd_row
            )
            self.cc_flux_col = copy.copy(self.mca_ci_col)
            self.cc_conc = copy.copy(self.mca_ci[self.__kmatrix__.shape[1] :, :])
            self.cc_conc_row = copy.copy(self.mca_ci_row[self.__kmatrix__.shape[1] :])
            self.cc_conc_col = copy.copy(self.mca_ci_col)
        elif self.__HAS_MOIETY_CONSERVATION__:
            print(
                'INFO: this is interesting no dependent flux cc\'s only dependent conc cc\'s!'
            )
            self.cc_flux = copy.copy(self.mca_ci[: self.__kmatrix__.shape[1], :])
            self.cc_flux_row = copy.copy(self.mca_ci_row[: self.__kmatrix__.shape[1]])
            self.cc_flux_col = copy.copy(self.mca_ci_col)
            self.cc_conc = numpy.concatenate(
                (self.mca_ci[self.__kmatrix__.shape[1] :, :], self.mca_csd)
            )
            self.cc_conc_row = copy.copy(
                self.mca_ci_row[self.__kmatrix__.shape[1] :] + self.mca_csd_row
            )
            self.cc_conc_col = copy.copy(self.mca_ci_col)
        else:
            print('INFO: this is interesting no dependent flux/conc coefficients!')
            self.cc_flux = copy.copy(self.mca_ci[: self.__kmatrix__.shape[1], :])
            self.cc_flux_row = copy.copy(self.mca_ci_row[: self.__kmatrix__.shape[1]])
            self.cc_flux_col = copy.copy(self.mca_ci_col)
            self.cc_conc = copy.copy(self.mca_ci[self.__kmatrix__.shape[1] :, :])
            self.cc_conc_row = copy.copy(self.mca_ci_row[self.__kmatrix__.shape[1] :])
            self.cc_conc_col = copy.copy(self.mca_ci_col)

        self.cc_all = numpy.concatenate((self.cc_flux, self.cc_conc))
        self.cc_all_row = self.cc_flux_row + self.cc_conc_row
        self.cc_all_col = copy.copy(self.mca_ci_col)

        # the new way of doing things (test phase) brett - 2007
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

        if self.__settings__["mca_ccj_upsymb"]:
            # CJ
            r, c = self.cc_all.shape
            CJoutput = ''
            for x in range(0, len(self.__reactions__)):
                for y in range(0, c):
                    if self.__settings__['mode_mca_scaled']:
                        setattr(
                            self,
                            'ccJ' + self.cc_all_row[x] + '_' + self.cc_all_col[y],
                            self.cc_all[x, y],
                        )
                    else:
                        setattr(
                            self,
                            'uccJ' + self.cc_all_row[x] + '_' + self.cc_all_col[y],
                            self.cc_all[x, y],
                        )

        if self.__settings__["mca_ccs_upsymb"]:
            # CS
            CSoutput = ''
            for x in range(len(self.__reactions__), r):
                for y in range(0, c):
                    if self.__settings__['mode_mca_scaled']:
                        setattr(
                            self,
                            'cc' + self.cc_all_row[x] + '_' + self.cc_all_col[y],
                            self.cc_all[x, y],
                        )
                    else:
                        setattr(
                            self,
                            'ucc' + self.cc_all_row[x] + '_' + self.cc_all_col[y],
                            self.cc_all[x, y],
                        )

        if self.__settings__['display_debug'] == 1:
            print('CJoutput')
            print(CJoutput)
            print('CSoutput')
            print(CSoutput)

    def showCC(self, File=None):
        """
        showCC(File=None)

        Print all control coefficients as \'LaTex\' formatted strings to the screen or file.

        Arguments:

        File [default=None]: an open, writable Python file object

        """
        r, c = self.cc_all.shape
        if self.__settings__["mca_ccall_altout"]:
            CAltoutput = ''
            for x in range(0, c):
                col = self.cc_all_col[x]
                CAltoutput += '\n`Reaction ' + col + '\'\n'
                for y in range(0, r):
                    if y == 0:
                        CAltoutput += '\"Flux control coefficients\"\n'
                    if y == len(self.__reactions__):
                        CAltoutput += '\"Concentration control coefficients\"\n'
                    row = self.cc_all_row[y]
                    rtemp = self.__settings__['mode_number_format'] % self.cc_all[y, x]
                    if y < len(self.__reactions__):
                        if self.__settings__['mode_mca_scaled']:
                            cc = '\\cc{J' + row + '}{' + col + '} = ' + rtemp
                        else:
                            cc = '\\ucc{J' + row + '}{' + col + '} = ' + rtemp
                    else:
                        if self.__settings__['mode_mca_scaled']:
                            cc = '\\cc{' + row + '}{' + col + '} = ' + rtemp
                        else:
                            cc = '\\ucc{' + row + '}{' + col + '} = ' + rtemp
                    CAltoutput += cc + '\n'
        else:
            if self.__settings__["mca_ccall_fluxout"]:
                # CJ
                CJoutput = ''
                for x in range(0, len(self.__reactions__)):
                    row = self.cc_all_row[x]
                    CJoutput += '\n`J' + row + '\'\n'
                    for y in range(0, c):
                        col = self.cc_all_col[y]
                        rtemp = (
                            self.__settings__['mode_number_format'] % self.cc_all[x, y]
                        )
                        if self.__settings__['mode_mca_scaled']:
                            cc = '\\cc{J' + row + '}{' + col + '} = ' + rtemp
                        else:
                            cc = '\\ucc{J' + row + '}{' + col + '} = ' + rtemp
                        CJoutput += cc + '\n'
            if self.__settings__["mca_ccall_concout"]:
                # CS
                CSoutput = ''
                for x in range(len(self.__reactions__), r):
                    row = self.cc_all_row[x]
                    CSoutput += '\n`' + row + '\'\n'
                    for y in range(0, c):
                        col = self.cc_all_col[y]
                        rtemp = (
                            self.__settings__['mode_number_format'] % self.cc_all[x, y]
                        )
                        if self.__settings__['mode_mca_scaled']:
                            cc = '\\cc{' + row + '}{' + col + '} = ' + rtemp
                        else:
                            cc = '\\ucc{' + row + '}{' + col + '} = ' + rtemp
                        CSoutput += cc + '\n'
        if File != None:
            # assert type(File) == file, 'showCC() needs an open file object'
            if self.__settings__["mca_ccall_altout"]:
                print('\nControl coefficients grouped by reaction')
                File.write('\n## Control coefficients grouped by reaction\n')
                File.write(CAltoutput)
            else:
                if self.__settings__["mca_ccall_fluxout"]:
                    print('\nFlux control coefficients')
                    File.write('\n## Flux control coefficients\n')
                    File.write(CJoutput)
                if self.__settings__["mca_ccall_concout"]:
                    print('\nConcentration control coefficients')
                    File.write('\n## Concentration control coefficients\n')
                    File.write(CSoutput)
        else:
            if self.__settings__["mca_ccall_altout"]:
                print('\nControl coefficients grouped by reaction')
                print(CAltoutput)
            else:
                if self.__settings__["mca_ccall_fluxout"]:
                    print('\nFlux control coefficients')
                    print(CJoutput)
                if self.__settings__["mca_ccall_concout"]:
                    print('\nConcentration control coefficients')
                    print(CSoutput)

    def EvalRC(self):
        """
        EvalRC()

        Calculate the MCA response coefficients using the current steady-state solution.

        Arguments:
        None

        """

        reordered_cc_all = numpy.zeros(self.cc_all.shape, 'd')
        for reac in range(len(self.elas_par_row)):
            reordered_cc_all[:, reac] = self.cc_all[
                :, list(self.cc_all_col).index(self.elas_par_row[reac])
            ]

        self.mca_rc_par = numpy.dot(reordered_cc_all, self.elas_par)
        self.mca_rc_par_row = self.cc_all_row
        self.mca_rc_par_col = self.elas_par_col
        del reordered_cc_all
        self.__structural__.MatrixFloatFix(self.mca_rc_par)
        self.rc = BagOfStuff(self.mca_rc_par, self.mca_rc_par_row, self.mca_rc_par_col)
        self.rc.load()
        if self.__settings__['mode_mca_scaled'] == 1:
            self.rc.scaled = True
        else:
            self.rc.scaled = False

    def EvalRCT(self):
        """
        EvalRCT()

        Calculate the MCA response coefficients using the current steady-state solution.

        Responses to changes in the sums of moiety conserved cycles are also calculated.

        Arguments:
        None

        """

        # We arbitrarily choose the order of reactions as that of elas_var_row
        # and the order of species as that of elas_var_col. The order of dependent
        # species (there is one for each moiety conserved cycle) is from the
        # conservation matrix (mod.Consmatrix.row). The final output is the
        # (m x (m-r)) R^S_T matrix and the (n x (m-r)) R^J_T matrix.
        #
        # See "Metabolic control analysis in a nutshell" by Hofmeyr and also
        # Kholodenko, Sauro, Westerhoff 1994 for the matrix formulations.
        #
        # Danie Palm - 2011-10-26

        # Reorder reactions in Cs and CJ to match elas_var_row
        reaction_reordered_cc_conc = numpy.zeros(self.cc_conc.shape, 'd')
        for reac in range(len(self.elas_var_row)):
            reaction_reordered_cc_conc[:, reac] = self.cc_conc[
                :, list(self.cc_conc_col).index(self.elas_var_row[reac])
            ]
        reaction_reordered_cc_flux = numpy.zeros(self.cc_flux.shape, 'd')
        for reac in range(len(self.elas_var_row)):
            reaction_reordered_cc_flux[:, reac] = self.cc_flux[
                :, list(self.cc_flux_col).index(self.elas_var_row[reac])
            ]

        # Reorder metabolites in Cs to match elas_var_col
        reordered_cc_conc = numpy.zeros(self.cc_conc.shape, 'd')
        for metab in range(len(self.elas_var_col)):
            reordered_cc_conc[metab, :] = reaction_reordered_cc_conc[
                list(self.cc_conc_row).index(self.elas_var_col[metab]), :
            ]

        # Reorder fluxes in CJ to match elas_var_row
        reordered_cc_flux = numpy.zeros(self.cc_flux.shape, 'd')
        for flux in range(len(self.elas_var_row)):
            reordered_cc_flux[flux, :] = reaction_reordered_cc_flux[
                list(self.cc_flux_row).index(self.elas_var_row[flux]), :
            ]

        # Some basic dimensions and matrices
        m = len(self.species)
        r = len(self.Lmatrix.col)
        Im = numpy.matrix(numpy.eye(m))

        # Construct the (right) pseudoinverse of the the conservation matrix in the order of elas_var_col
        reordered_zero_I = numpy.zeros((m, m - r), 'd')
        for species in self.elas_var_col:
            if species in self.Consmatrix.row:
                row = self.elas_var_col.index(species)
                col = self.Consmatrix.row.index(species)
                reordered_zero_I[row, col] = 1.0

        # Unlike normal RCs, separate calculations are required for scaled
        # and unscaled RCT's. We also need to pick the right elasticity
        # matrix: elas_var for scaled and elas_var_u for unscaled.
        if self.__settings__['mode_mca_scaled'] == 1:
            # Some scaling matrices
            diag_T = None
            if self.__KeyWords__['Species_In_Conc']:
                diag_T = numpy.matrix(numpy.diag(self.__tvec_c__))
            else:
                diag_T = numpy.matrix(numpy.diag(self.__tvec_a__))
            diag_S = numpy.matrix(
                numpy.diag(self.data_sstate.getStateData(*self.elas_var_col))
            )

            # Calculate concentration responses to changes in T
            # (Cs*es + Im) * diag_S.I * reordered_zero_I * diag_T
            self.mca_rct_conc = numpy.dot(
                numpy.dot(
                    numpy.dot(
                        numpy.dot(reordered_cc_conc, self.elas_var) + Im,
                        numpy.linalg.inv(diag_S),
                    ),
                    reordered_zero_I,
                ),
                diag_T,
            )

            # Calculate flux responses to changes in T
            # CJ*es * diag_S.I * reordered_zero_I * diag_T
            self.mca_rct_flux = numpy.dot(
                numpy.dot(
                    numpy.dot(
                        numpy.dot(reordered_cc_flux, self.elas_var),
                        numpy.linalg.inv(diag_S),
                    ),
                    reordered_zero_I,
                ),
                diag_T,
            )
        else:
            # Calculate concentration responses to changes in T
            # (Cs*es + Im) * reordered_zero_I
            self.mca_rct_conc = numpy.dot(
                numpy.dot(reordered_cc_conc, self.elas_var_u) + Im, reordered_zero_I
            )

            # Calculate flux responses to changes in T
            # CJ*es * reordered_zero_I
            self.mca_rct_flux = numpy.dot(
                numpy.dot(reordered_cc_flux, self.elas_var_u), reordered_zero_I
            )

        # Cleanup
        del reordered_cc_conc
        del reordered_cc_flux

        # We simply prepend 'T_' to the name of the dependent species to prevent
        # confusion with the real species. This is a parameter, not variable.
        moiety_sum_names = [
            'T_{0}'.format(dependent_species)
            for dependent_species in self.Consmatrix.row
        ]

        self.__structural__.MatrixFloatFix(self.mca_rct_conc)
        self.mca_rct_conc_row = self.elas_var_col
        self.mca_rct_conc_col = moiety_sum_names
        self.__structural__.MatrixFloatFix(self.mca_rct_flux)
        self.mca_rct_flux_row = self.elas_var_row
        self.mca_rct_flux_col = moiety_sum_names

        # Augment the existing m.rc object (this could be more efficient)

        # Reordering as a matter of principle: too many burnt fingers
        vstacked_col = list(self.elas_var_row) + list(self.elas_var_col)
        vstacked = numpy.vstack((self.mca_rct_flux, self.mca_rct_conc))
        reordered_vstack = numpy.zeros(vstacked.shape, 'd')
        for row_index in range(len(self.rc.row)):
            reordered_vstack[row_index, :] = vstacked[
                list(vstacked_col).index(self.rc.row[row_index]), :
            ]
        hstacked = numpy.hstack((self.rc.matrix, vstacked))

        self.rc = BagOfStuff(hstacked, self.rc.row, self.rc.col + moiety_sum_names)
        self.rc.load()
        if self.__settings__['mode_mca_scaled'] == 1:
            self.rc.scaled = True
        else:
            self.rc.scaled = False

    def EvalEigen(self):
        """
        EvalEigen()

        Calculate the eigenvalues or vectors of the unscaled Jacobian matrix and thereby
        analyse the stability of a system

        Arguments:
        None

        """
        Nr = copy.copy(self.__nrmatrix__)

        # sort E to match K and L
        Es = numpy.zeros((self.elas_var_u.shape), 'd')
        # print e2

        for x in range(0, len(self.nmatrix_col)):
            for y in range(0, len(self.lmatrix_row)):
                Es[x, y] = self.elas_var_u[self.nmatrix_col[x], self.lmatrix_row[y]]
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
                # returns eigenvalues as well as left and right eigenvectors
                eigenval, self.eigen_vecleft, self.eigen_vecright = scipy.linalg.eig(
                    jacobian, numpy.identity(jacobian.shape[0], 'd'), left=1, right=1
                )
                if self.__settings__['display_debug'] == 1:
                    print('\nEigenvalues')
                    print(eigenval)
                    print('\nLeft Eigenvector')
                    print(self.eigen_vecleft)
                    print('\nRight Eigenvector')
                    print(self.eigen_vecright)
            else:
                # returns eigenvalues
                eigenval = scipy.linalg.eigvals(jacobian)
                if self.__settings__['display_debug'] == 1:
                    print('\nEigenvalues')
                    print(eigenval)

        self.eigen_values = eigenval
        # self.eigen_order = tuple(eigenorder)

        for x in range(len(self.eigen_values)):
            setattr(self, 'lambda' + str(x + 1), self.eigen_values[x])

        if self.__settings__['display_debug'] == 1:
            print(
                '\nEigenvalues attached as lambda1 ... lambda'
                + repr(len(eigenorder))
                + '\n'
            )

    def showEigen(self, File=None):
        """
        showEigen(File=None)

        Print the eigenvalues and stability analysis of a system generated with EvalEigen()
        to the screen or file.

        Arguments:

        File [default=None]: an open, writable Python file object

        """
        eigenstats = ''
        eigenstats += (
            'Max real part: '
            + self.__settings__['mode_number_format'] % max(self.eigen_values.real)
            + '\n'
        )
        eigenstats += (
            'Min real part: '
            + self.__settings__['mode_number_format'] % min(self.eigen_values.real)
            + '\n'
        )
        eigenstats += (
            'Max absolute imaginary part: '
            + self.__settings__['mode_number_format'] % max(abs(self.eigen_values.imag))
            + '\n'
        )
        eigenstats += (
            'Min absolute imaginary part: '
            + self.__settings__['mode_number_format'] % min(abs(self.eigen_values.imag))
            + '\n'
        )
        eigenstats += (
            'Stiffness: '
            + self.__settings__['mode_number_format']
            % (max(abs(self.eigen_values.real)) / min(abs(self.eigen_values.real)))
            + '\n'
        )

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
            # assert type(File) == file, 'showEigen() needs an open file object'
            print('\nEigen values')
            File.write('\n## Eigen values\n')
            scipy.io.write_array(File, self.eigen_values, precision=2, keep_open=1)
            print('\nEigen statistics')
            File.write('\n## Eigen statistics\n')
            File.write(eigenstats)
            if self.__settings__['mode_eigen_output']:
                print('\nLeft eigen vector')
                File.write('\n## Left eigen vector\n')
                scipy.io.write_array(File, self.eigen_vecleft, precision=2, keep_open=1)
                print('\nRight eigen vector')
                File.write('\n## Right eigen vector\n')
                scipy.io.write_array(
                    File, self.eigen_vecright, precision=2, keep_open=1
                )
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

    # Utility functions
    # new generation metafunctions
    ##  def doLoad(self,stoich_load=0):
    ##  """
    ##  doLoad(stoich_load=0)

    ##  Load and instantiate a PySCeS model so that it can be used for further analyses.

    ##  Calls model loading subroutines:
    ##  Stoichiometry_Analyse() [override=0,load=stoich_load]
    ##  InitialiseModel()

    ##  Arguments:

    ##  stoich_load [default=0]: try to load a stoichiometry saved with Stoichiometry_Save_Serial()

    ##  """
    ##  self.InitialiseInputFile()
    ##  assert self.__parseOK, '\nError in input file, parsing could not complete'
    ##  self.Stoichiometry_Analyse(override=0,load=stoich_load)
    ##  # add new Style functions to model
    ##  self.InitialiseFunctions()
    ##  self.InitialiseCompartments()
    ##  self.InitialiseRules()
    ##  self.InitialiseEvents()
    ##  self.InitialiseOldFunctions() # TODO replace this with initialisation functions
    ##  self.InitialiseModel()
    ##  self.InitialiseRuleChecks()

    def doSim(self, end=10.0, points=21):
        """
        doSim(end=10.0,points=20.0)

        Run a time simulation from t=0 to t=sim_end with sim_points.

        Calls:
        Simulate()

        Arguments:

        end [default=10.0]: simulation end time
        points [default=20.0]: number of points in the simulation

        """
        self.sim_end = end
        self.sim_points = points
        self.Simulate()

    def doSimPlot(
        self, end=10.0, points=21, plot='species', fmt='lines', filename=None
    ):
        """
        Run a time simulation from t=0 to t=sim_end with sim_points and plot the results.
        The required output data and format can be set:

        - *end* the end time (default=10.0)
        - *points* the number of points in the simulation (default=20.0)
        - *plot* (default='species') select output data

         - 'species'
         - 'rates'
         - 'all' both species and rates

        - *fmt* plot format, UPI backend dependent (default='') or the *CommonStyle* 'lines' or 'points'.
        - *filename* if not None (default) then the plot is exported as *filename*.png

        Calls:
        - **Simulate()**
        - **SimPlot()**
        """
        self.sim_end = end
        self.sim_points = points
        self.Simulate()
        self.SimPlot(plot=plot, format=fmt, filename=filename)

    def exportSimAsSedML(
        self,
        output='files',
        return_sed=False,
        vc_given='PySCeS',
        vc_family='Software',
        vc_email='',
        vc_org='pysces.sourceforge.net',
    ):
        """
        Exports the current simulation as SED-ML in various ways it creates and stores the SED-ML files in a folder
        generated from the model name.

         - *output* [default='files'] the SED-ML export type can be one or more comma separated e.g. 'files,combine'
         - *files* export the plain SBML and SEDML XML files
         - *archive* export as a SED-ML archive *<file>.sedx* containing the SBML and SEDML xml files
         - *combine* export as a COMBINE archive *<file>.omex* containing the SBML, SEDML, manifest (XML) and metadata (RDF)
           - *vc_given* [default='PySCeS']
           - *vc_family* [default='Software']
           - *vc_email* [default='bgoli@users.sourceforge.net']
           - *vc_org* [default='<pysces.sourceforge.net>']

        """
        sedname = self.ModelFile.replace('.psc', '')

        sedout = os.path.join(OUTPUT_DIR, 'sedout')
        if not os.path.exists(sedout):
            os.makedirs(sedout)
        S = SED.SED(sedname + '_sed', sedout)
        S.addModel(sedname, self)
        S.addSimulation(
            'sim0',
            self.sim_start,
            self.sim_end,
            self.sim_points,
            self.__SIMPLOT_OUT__,
            initial=None,
            algorithm='KISAO:0000019',
        )
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
                S.writeCOMBINEArchive(
                    vc_given=vc_given,
                    vc_family=vc_family,
                    vc_email=vc_email,
                    vc_org=vc_org,
                )
        S._SED_CURRENT_ = False
        if return_sed:
            return S
        else:
            del S

    def doSimPerturb(self, pl, end):
        """
        **Deprecated**: use events instead
        """
        pass

    def doState(self):
        """
        doState()

        Calculate the steady-state solution of the system.

        Calls:
        State()

        Arguments:
        None

        """
        self.State()

    def doStateShow(self):
        """
        doStateShow()

        Calculate the steady-state solution of a system and show the results.

        Calls:
        State()
        showState()

        Arguments:
        None

        """
        self.State()
        assert (
            self.__StateOK__ == True
        ), '\n\nINFO: Invalid steady state: run mod.doState() and check for errors'
        self.showState()

    def doElas(self):
        """
        doElas()

        Calculate the model elasticities, this method automatically calculates a steady state.

        Calls:
        State()
        EvalEvar()
        EvalEpar()

        Arguments:
        None

        """
        self.State()
        assert (
            self.__StateOK__ == True
        ), '\n\nINFO: Invalid steady state: run mod.doState() and check for errors'
        self.EvalEvar()
        self.EvalEpar()

    def doEigen(self):
        """
        doEigen()

        Calculate the eigenvalues, automatically performs a steady state and elasticity analysis.

        Calls:
        State()
        EvalEvar()
        Evaleigen()

        Arguments:
        None

        """
        self.State()
        assert (
            self.__StateOK__ == True
        ), '\n\nINFO: Invalid steady state: run mod.doState() and check for errors'
        self.__settings__["elas_evar_upsymb"] = 1
        self.EvalEvar()
        self.__settings__["elas_evar_upsymb"] = 1
        self.EvalEigen()

    def doEigenShow(self):
        """
        doEigenShow()

        Calculate the eigenvalues, automatically performs a steady state and elasticity analysis
        and displays the results.

        Calls:
        doEigen()
        showEigen()

        Arguments:
        None

        """
        self.doEigen()
        self.showEigen()

    def doEigenMca(self):
        """
        doEigenMca()

        Calculate a full Control Analysis and eigenvalues, automatically performs a steady state, elasticity, control analysis.

        Calls:
        State()
        EvalEvar()
        EvalCC()
        Evaleigen()

        Arguments:
        None

        """
        self.State()
        assert (
            self.__StateOK__ == True
        ), '\n\nINFO: Invalid steady state: run mod.doState() and check for errors'
        self.EvalEvar()
        self.EvalCC()
        self.EvalEigen()

    def doMca(self):
        """
        doMca()

        Perform a complete Metabolic Control Analysis on the model, automatically calculates a steady state.

        Calls:
        State()
        EvalEvar()
        EvalEpar()
        EvalCC()

        Arguments:
        None

        """
        self.State()
        assert (
            self.__StateOK__ == True
        ), '\n\nINFO: Invalid steady state run: mod.doState() and check for errors'
        self.EvalEvar()
        self.EvalEpar()
        self.EvalCC()

    def doMcaRC(self):
        """
        doMca()

        Perform a complete Metabolic Control Analysis on the model, automatically calculates a steady state.

        Calls:
        State()
        EvalEvar()
        EvalEpar()
        EvalCC()
        EvalRC()

        Arguments:
        None

        """
        self.State()
        assert (
            self.__StateOK__ == True
        ), '\n\nINFO: Invalid steady state run: mod.doState() and check for errors'
        self.EvalEvar()
        self.EvalEpar()
        self.EvalCC()
        self.EvalRC()

    def doMcaRCT(self):
        """
        doMcaRCT()

        Perform a complete Metabolic Control Analysis on the model, automatically calculates a steady state.

        In additional, response coefficients to the sums of moiety-conserved cycles are calculated.

        Calls:
        State()
        EvalEvar()
        EvalEpar()
        EvalCC()
        EvalRC()
        EvalRCT()

        Arguments:
        None

        """
        self.State()
        assert (
            self.__StateOK__ == True
        ), '\n\nINFO: Invalid steady state run: mod.doState() and check for errors'
        self.EvalEvar()
        self.EvalEpar()
        self.EvalCC()
        self.EvalRC()
        if self.__HAS_MOIETY_CONSERVATION__:
            self.EvalRCT()
        else:
            print('No moiety conservation detected, reverting to simple doMcaRC().')

    # show/save function prototypes
    def showSpecies(self, File=None):
        """
        showSpecies(File=None)

        Prints the current value of the model's variable species (mod.X) to screen or file.

        Arguments:

        File [default=None]: an open, writable Python file object

        """
        out_list = []

        print('\nSpecies values')
        out_list.append('\n## species values\n')
        for x in range(len(self.__species__)):
            if File == None:
                ##  print self.__species__[x] + ' = ' +  self.__settings__['mode_number_format'] % eval('self.' + self.__species__[x])
                print(
                    self.__species__[x]
                    + ' = '
                    + self.__settings__['mode_number_format']
                    % getattr(self, self.__species__[x])
                )
            else:
                ##  out_list.append(self.__species__[x] + ' = ' +  self.__settings__['mode_number_format'] % eval('self.' + self.__species__[x]) + '\n')
                out_list.append(
                    self.__species__[x]
                    + ' = '
                    + self.__settings__['mode_number_format']
                    % getattr(self, self.__species__[x])
                    + '\n'
                )
        if File != None:
            for x in out_list:
                File.write(x)

    def showSpeciesI(self, File=None):
        """
        showSpeciesI(File=None)

        Prints the current value of the model's variable species initial values (mod.X_init) to screen or file.

        Arguments:

        File [default=None]: an open, writable Python file object

        """
        out_list = []

        print('\nSpecies initial values')
        out_list.append('\n## species initial values\n')
        for x in range(len(self.__species__)):
            if File == None:
                print(
                    self.__species__[x]
                    + '_init = '
                    + self.__settings__['mode_number_format']
                    % getattr(self, self.__species__[x] + '_init')
                )
            else:
                out_list.append(
                    self.__species__[x]
                    + ' = '
                    + self.__settings__['mode_number_format']
                    % getattr(self, self.__species__[x] + '_init')
                    + '\n'
                )
        if File != None:
            for x in out_list:
                File.write(x)

    def showSpeciesFixed(self, File=None):
        """
        showSpeciesFixed(File=None)

        Prints the current value of the model's fixed species values (mod.X) to screen or file.

        Arguments:

        File [default=None]: an open, writable Python file object

        """
        out_list = []

        print('\nFixed species')
        out_list.append('\n## fixed species\n')
        for x in range(len(self.__fixed_species__)):
            if File == None:
                print(
                    self.__fixed_species__[x]
                    + ' = '
                    + self.__settings__['mode_number_format']
                    % getattr(self, self.__fixed_species__[x])
                )
            else:
                out_list.append(
                    self.__fixed_species__[x]
                    + ' = '
                    + self.__settings__['mode_number_format']
                    % getattr(self, self.__fixed_species__[x])
                    + '\n'
                )
        if File != None:
            for x in out_list:
                File.write(x)

    def showPar(self, File=None):
        """
        showPar(File=None)

        Prints the current value of the model's parameter values (mod.P) to screen or file.

        Arguments:

        File [default=None]: an open, writable Python file object

        """
        out_list = []

        print('\nParameters')
        out_list.append('\n## parameters\n')
        for x in range(len(self.__parameters__)):
            if File == None:
                print(
                    self.__parameters__[x]
                    + ' = '
                    + self.__settings__['mode_number_format']
                    % getattr(self, self.__parameters__[x])
                )
            else:
                out_list.append(
                    self.__parameters__[x]
                    + ' = '
                    + self.__settings__['mode_number_format']
                    % getattr(self, self.__parameters__[x])
                    + '\n'
                )
        if File != None:
            for x in out_list:
                File.write(x)

    def showModifiers(self, File=None):
        """
        showModifiers(File=None)

        Prints the current value of the model's modifiers per reaction to screen or file.

        Arguments:

        File [default=None]: an open, writable Python file object

        """
        noMod = []
        out_list = []

        print('\nModifiers per reaction:')
        out_list.append('\n## modifiers per reaction\n')
        for reac in self.__modifiers__:
            rstr = ''
            if len(reac[1]) == 1:
                if File == None:
                    print(reac[0] + ' has modifier:', end=' ')
                rstr = rstr + reac[0] + ' has modifier: '
                for x in reac[1]:
                    if File == None:
                        print(x, end=' ')
                    rstr = rstr + x + ' '
                if File == None:
                    print(' ')
                rstr += '\n'
            elif len(reac[1]) > 1:
                if File == None:
                    print(reac[0] + ' has modifiers: ', end=' ')
                rstr = rstr + reac[0] + ' has modifiers: '
                for x in reac[1]:
                    if File == None:
                        print(x, end=' ')
                    rstr = rstr + x + ' '
                if File == None:
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
                    if File == None:
                        print(' ')
                    rstr += '\n'
                    cntr = 0
                if File == None:
                    print(noMod[n], end=' ')
                rstr = rstr + noMod[n] + ' '
            if File == None:
                print(' ')
            rstr += '\n'
            out_list.append(rstr)

            if File != None:
                for x in out_list:
                    File.write(x)

    def showState(self, File=None):
        """
        showState(File=None)

        Prints the result of the last steady-state analyses. Both steady-state flux's and species concentrations are shown.

        Arguments:

        File [default=None]: an open, writable Python file object

        """
        out_list = []

        print('\nSteady-state species concentrations')
        out_list.append('\n## Current steady-state species concentrations\n')
        if self.__StateOK__:
            for x in range(len(self.state_species)):
                if File == None:
                    print(
                        self.__species__[x]
                        + '_ss = '
                        + self.__settings__['mode_number_format']
                        % self.state_species[x]
                    )
                else:
                    out_list.append(
                        self.__species__[x]
                        + '_ss = '
                        + self.__settings__['mode_number_format']
                        % self.state_species[x]
                        + '\n'
                    )
        else:
            print('No valid steady state found')
            out_list.append('No valid steady state found.\n')

        print('\nSteady-state fluxes')
        out_list.append('\n## Steady-state fluxes\n')
        if self.__StateOK__:
            for x in range(len(self.state_flux)):
                if File == None:
                    print(
                        'J_'
                        + self.__reactions__[x]
                        + ' = '
                        + self.__settings__['mode_number_format'] % self.state_flux[x]
                    )
                else:
                    out_list.append(
                        'J_'
                        + self.__reactions__[x]
                        + ' = '
                        + self.__settings__['mode_number_format'] % self.state_flux[x]
                        + '\n'
                    )
        else:
            print('No valid steady state found')
            out_list.append('No valid steady state found.\n')

        if File != None:
            for x in out_list:
                File.write(x)

    def showRate(self, File=None):
        """
        Prints the current rates of all the reactions using the current parameter values and species concentrations

        - *File* an open, writable Python file object (default=None)

        """
        Vtemp = [r() for r in getattr(self, self.__reactions__[r])]
        outrate = ''
        for x in range(len(self.__reactions__)):
            outrate += (
                self.__reactions__[x]
                + ' = '
                + self.__settings__['mode_number_format'] % Vtemp[x]
                + '\n'
            )
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

        Prints the reaction stoichiometry and rate equations to screen or File.

        Arguments:

        File [default=None]: an open, writable Python file object
        """
        out_list = []

        tnform = self.__settings__['mode_number_format']

        self.__settings__['mode_number_format'] = '%2.5f'

        print('\nReaction stoichiometry and rate equations')
        out_list.append('\n## Reaction stoichiometry and rate equations\n')

        # writes these out in a better order
        for key in self.Kmatrix.row:
            if File == None:
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
                        reagR.append(
                            '{'
                            + self.__settings__['mode_number_format']
                            % abs(self.__nDict__[key]['Reagents'][reagent])
                            + '}'
                            + reagent.replace('self.', '')
                        )
                elif self.__nDict__[key]['Reagents'][reagent] < 0:
                    if self.__nDict__[key]['Reagents'][reagent] == -1.0:
                        reagL.append(reagent.replace('self.', ''))
                    else:
                        reagL.append(
                            '{'
                            + self.__settings__['mode_number_format']
                            % abs(self.__nDict__[key]['Reagents'][reagent])
                            + '}'
                            + reagent.replace('self.', '')
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

            if self.__nDict__[key]['Type'] == 'Rever':
                symbol = ' = '
            else:
                symbol = ' > '
            if File == None:
                print('\t' + substring + symbol + prodstring)
                print('\t' + self.__nDict__[key]['RateEq'].replace('self.', ''))
            else:
                out_list.append('\t' + substring + symbol + prodstring + '\n')
                out_list.append(
                    '\t' + self.__nDict__[key]['RateEq'].replace('self.', '') + '\n\n'
                )
        if len(list(self.__rules__.keys())) > 0:
            out_list.append('\n# Assignment rules\n')
            for ass in self.__rules__:
                out_list.append(
                    '!F {} = {}\n'.format(
                        self.__rules__[ass]['name'], self.__rules__[ass]['formula']
                    )
                )
            out_list.append('\n')
        if File != None:
            for x in out_list:
                File.write(x)
        self.__settings__['mode_number_format'] = tnform

    def showODE(self, File=None, fmt='%2.3f'):
        """
        showODE(File=None,fmt='%2.3f')

        Print a representation of the full set of ODE's generated by PySCeS to screen or file.

        Arguments:

        File [default=None]: an open, writable Python file object
        fmt [default='%2.3f']: output number format

        """
        maxmetlen = 0
        for x in self.__species__:
            if len(x) > maxmetlen:
                maxmetlen = len(x)

        maxreaclen = 0
        for x in self.__reactions__:
            if len(x) > maxreaclen:
                maxreaclen = len(x)
        odes = '\n## ODE\'s (unreduced)\n'
        for x in range(self.__nmatrix__.shape[0]):
            odes += (
                'd'
                + self.__species__[x]
                + (maxmetlen - len(self.__species__[x])) * ' '
                + '|'
            )
            beginline = 0
            for y in range(self.__nmatrix__.shape[1]):
                if abs(self.__nmatrix__[x, y]) > 0.0:
                    if self.__nmatrix__[x, y] > 0.0:
                        if beginline == 0:
                            odes += (
                                '  '
                                + fmt % abs(self.__nmatrix__[x, y])
                                + '*'
                                + self.__reactions__[y]
                                + (maxreaclen - len(self.__reactions__[y])) * ' '
                            )
                            beginline = 1
                        else:
                            odes += (
                                ' + '
                                + fmt % abs(self.__nmatrix__[x, y])
                                + '*'
                                + self.__reactions__[y]
                                + (maxreaclen - len(self.__reactions__[y])) * ' '
                            )
                    else:
                        if beginline == 0:
                            odes += (
                                ' -'
                                + fmt % abs(self.__nmatrix__[x, y])
                                + '*'
                                + self.__reactions__[y]
                                + (maxreaclen - len(self.__reactions__[y])) * ' '
                            )
                        else:
                            odes += (
                                ' - '
                                + fmt % abs(self.__nmatrix__[x, y])
                                + '*'
                                + self.__reactions__[y]
                                + (maxreaclen - len(self.__reactions__[y])) * ' '
                            )
                        beginline = 1
            odes += '\n'

        if self.__HAS_RATE_RULES__:
            odes += "\n## Rate Rules:\n"
            for rule in self.__rate_rules__:
                odes += "d{} | {}\n".format(
                    rule, self.__rules__[rule]['formula'].replace('()', ''),
                )
            odes += '\n'

        if File != None:
            print('\nODE\'s (unreduced)\n')
            File.write(odes)
        else:
            print(odes)

    def showODEr(self, File=None, fmt='%2.3f'):
        """
        showODEr(File=None,fmt='%2.3f')

        Print a representation of the reduced set of ODE's generated by PySCeS to screen or file.

        Arguments:

        File [default=None]: an open, writable Python file object
        fmt [default='%2.3f']: output number format

        """
        maxmetlen = 0
        for x in self.__species__:
            if len(x) > maxmetlen:
                maxmetlen = len(x)
        maxreaclen = 0
        for x in self.__reactions__:
            if len(x) > maxreaclen:
                maxreaclen = len(x)

        odes = '\n## ODE\'s (reduced)\n'
        for x in range(self.nrmatrix.shape[0]):
            odes += (
                'd'
                + self.__species__[self.nrmatrix_row[x]]
                + (maxmetlen - len(self.__species__[self.nrmatrix_row[x]])) * ' '
                + '|'
            )
            beginline = 0
            for y in range(self.__nrmatrix__.shape[1]):
                if abs(self.__nrmatrix__[x, y]) > 0.0:
                    if self.__nrmatrix__[x, y] > 0.0:
                        if beginline == 0:
                            odes += (
                                '  '
                                + fmt % abs(self.__nrmatrix__[x, y])
                                + '*'
                                + self.__reactions__[self.nrmatrix_col[y]]
                                + (
                                    maxreaclen
                                    - len(self.__reactions__[self.nrmatrix_col[y]])
                                )
                                * ' '
                            )
                            beginline = 1
                        else:
                            odes += (
                                ' + '
                                + fmt % abs(self.__nrmatrix__[x, y])
                                + '*'
                                + self.__reactions__[self.nrmatrix_col[y]]
                                + (
                                    maxreaclen
                                    - len(self.__reactions__[self.nrmatrix_col[y]])
                                )
                                * ' '
                            )
                    else:
                        if beginline == 0:
                            odes += (
                                ' -'
                                + fmt % abs(self.__nrmatrix__[x, y])
                                + '*'
                                + self.__reactions__[self.nrmatrix_col[y]]
                                + (
                                    maxreaclen
                                    - len(self.__reactions__[self.nrmatrix_col[y]])
                                )
                                * ' '
                            )
                        else:
                            odes += (
                                ' - '
                                + fmt % abs(self.__nrmatrix__[x, y])
                                + '*'
                                + self.__reactions__[self.nrmatrix_col[y]]
                                + (
                                    maxreaclen
                                    - len(self.__reactions__[self.nrmatrix_col[y]])
                                )
                                * ' '
                            )
                        beginline = 1
            odes += '\n'

        if self.__HAS_RATE_RULES__:
            odes += "\n## Rate Rules:\n"
            for rule in self.__rate_rules__:
                odes += "d{} | {}\n".format(
                    rule, self.__rules__[rule]['formula'].replace('()', ''),
                )
            odes += '\n'

        if File != None:
            print('\nODE\'s (reduced)\n')
            File.write(odes)
            self.showConserved(File)
        else:
            print(odes)
            self.showConserved()

    def showModel(self, filename=None, filepath=None, skipcheck=0):
        """
        showModel(filename=None,filepath=None,skipcheck=0)

        The PySCeS 'save' command, prints the entire model to screen or File in a PSC format.
        (Currently this only applies to basic model attributes, ! functions are not saved).

        Arguments:

        filename [default=None]: the output PSC file
        filepath [default=None]: the output directory
        skipcheck [default=0]: skip check to see if the file exists (1) auto-averwrite

        """
        if filepath == None:
            filepath = self.ModelDir

        if filename == None:
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
            # assert type(filename) == str, 'showModel() takes a string filename argument'
            FBuf = 0
            if type(filename) == str:
                # filename = str(filename)
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
                            inp = input(
                                '\nfile ' + filex + ' exists.\nOverwrite? ([y]/n) '
                            )
                        else:
                            raise IOError
                        if inp == 'y' or inp == '':
                            go = 1
                            loop = 1
                        elif inp == 'n':
                            filename = input(
                                '\nfile "'
                                + filename
                                + '" exists. Enter a new filename: '
                            )
                            go = 1
                            filex = os.path.join(filepath, filename)
                            filename = chkpsc(filename)
                        else:
                            print('\nInvalid input')
                    except IOError:
                        print('\nfile "' + filex + '" does not exist, proceeding')
                        loop = 1
                        go = 1
            else:
                print('\nI hope we have a filebuffer')
                if (
                    type(filename).__name__ == 'file'
                    or type(filename).__name__ == 'StringIO'
                    or type(filename).__name__ == 'TextIOWrapper'
                ):
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

                # header =   '# \n'
                header = '# PySCeS (' + __version__ + ') model input file \n'
                # header += '# Copyright Brett G. Olivier, 2004 (bgoli at sun dot ac dot za)  \n'
                header += '# http://pysces.sourceforge.net/ \n'
                header += '# \n'
                header += '# Original input file: ' + self.ModelFile + '\n'
                header += (
                    '# This file generated: '
                    + time.strftime("%a, %d %b %Y %H:%M:%S")
                    + '\n\n'
                )
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
                if (
                    type(filename) == str
                    or type(filename).__name__ == 'file'
                    or type(filename).__name__ == 'TextIOWrapper'
                ):
                    # don't close StringIO objects
                    outFile.close()

    def showN(self, File=None, fmt='%2.3f'):
        """
        showN(File=None,fmt='%2.3f')

        Print the stoichiometric matrix (N), including row and column labels to screen or File.

        Arguments:

        File [default=None]: an open, writable Python file object
        fmt [default='%2.3f']: output number format
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
                km += self.__reactions__[a] + (9 - len(self.__reactions__[a]) - 1) * ' '
        for x in range(len(self.nmatrix_row)):
            if len(self.__species__[x]) > 6:
                km += '\n' + self.__species__[x][:5] + '*' + 1 * ' '
                km_trunc2 = 1
                rtr.append(x)
            else:
                km += (
                    '\n'
                    + self.__species__[x]
                    + (8 - len(self.__species__[x]) - 1) * ' '
                )
                # km += '\n' + self.__reactions__[self.nmatrix_row[x]] + 4*' '
            for y in range(len(self.nmatrix_col)):
                if y != 0:
                    km += 3 * ' '
                if self.__nmatrix__[x, y] >= 10.0:
                    km += ' ' + fmt % self.__nmatrix__[x, y]
                elif self.__nmatrix__[x, y] >= 0.0:
                    km += '  ' + fmt % self.__nmatrix__[x, y]
                else:
                    if abs(self.__nmatrix__[x, y]) >= 10.0:
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
        showNr(File=None,fmt='%2.3f')

        Print the reduced stoichiometric matrix (Nr), including row and column labels to screen or File.

        Arguments:

        File [default=None]: an open, writable Python file object
        fmt [default='%2.3f']: output number format
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
                km += self.__reactions__[a] + (9 - len(self.__reactions__[a]) - 1) * ' '
        for x in range(len(self.nrmatrix_row)):
            if len(self.__species__[self.nrmatrix_row[x]]) > 6:
                km += '\n' + self.__species__[self.nrmatrix_row[x]][:5] + '*' + 1 * ' '
                km_trunc2 = 1
                rtr.append(x)
            else:
                km += (
                    '\n'
                    + self.__species__[self.nrmatrix_row[x]]
                    + (8 - len(self.__species__[self.nrmatrix_row[x]]) - 1) * ' '
                )
                # km += '\n' + self.__reactions__[self.nrmatrix_row[x]] + 4*' '
            for y in range(len(self.nrmatrix_col)):
                if y != 0:
                    km += 3 * ' '
                if self.__nrmatrix__[x, y] >= 10.0:
                    km += ' ' + fmt % self.__nrmatrix__[x, y]
                elif self.__nrmatrix__[x, y] >= 0.0:
                    km += '  ' + fmt % self.__nrmatrix__[x, y]
                else:
                    if abs(self.__nrmatrix__[x, y]) >= 10.0:
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
        showK(File=None,fmt='%2.3f')

        Print the Kernel matrix (K), including row and column labels to screen or File.

        Arguments:

        File [default=None]: an open, writable Python file object
        fmt [default='%2.3f']: output number format
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
                km += (
                    self.__reactions__[self.kmatrix_col[a]]
                    + (9 - len(self.__reactions__[self.kmatrix_col[a]]) - 1) * ' '
                )

        for x in range(len(self.kmatrix_row)):
            if len(self.__reactions__[self.kmatrix_row[x]]) > 6:
                km += '\n' + self.__reactions__[self.kmatrix_row[x]][:5] + '*' + 1 * ' '
                km_trunc2 = 1
                rtr.append(x)
            else:
                km += (
                    '\n'
                    + self.__reactions__[self.kmatrix_row[x]]
                    + (8 - len(self.__reactions__[self.kmatrix_row[x]]) - 1) * ' '
                )
                # km += '\n' + self.__reactions__[self.kmatrix_row[x]] + 4*' '
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
        showL(File=None,fmt='%2.3f')

        Print the Link matrix (L), including row and column labels to screen or File.

        Arguments:

        File [default=None]: an open, writable Python file object
        fmt [default='%2.3f']: output number format
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
                km += (
                    self.__species__[self.lmatrix_col[a]]
                    + (9 - len(self.__species__[self.lmatrix_col[a]]) - 1) * ' '
                )

        for x in range(len(self.lmatrix_row)):
            if len(self.__species__[self.lmatrix_row[x]]) > 6:
                km += '\n' + self.__species__[self.lmatrix_row[x]][:5] + '* '
                km_trunc2 = 1
                rtr.append(x)
            else:
                km += (
                    '\n'
                    + self.__species__[self.lmatrix_row[x]]
                    + (8 - len(self.__species__[self.lmatrix_row[x]]) - 1) * ' '
                )
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
                print('\"No moiety conservation\"\n')
                print(km)
            else:
                print(km)

    # Internal test functions
    def TestSimState(self, endTime=10000, points=101, diff=1.0e-5):
        """
        **Deprecated**
        """
        pass

    def ResetNumberFormat(self):
        """
        ResetNumberFormat()

        Reset PySCeS default number format stored as mod.mode_number format to %2.4e

        Arguments:
        None

        """
        self.__settings__['mode_number_format'] = '%2.4e'

    def Write_array(
        self, inputd, File=None, Row=None, Col=None, close_file=0, separator='  '
    ):
        """
        Write_array(input,File=None,Row=None,Col=None,close_file=0,separator='  ')

        Write an array to File with optional row/col labels. A ',' separator can be specified to create
        a CSV style file.

        mod.__settings__['write_array_header']: add <filename> as a header line (1 = yes, 0 = no)
        mod.__settings__['write_array_spacer']: add a space after the header line (1 = yes, 0 = no)
        mod.__settings__['write_arr_lflush']: set the flush rate for large file writes

        Arguments:

        input: the array to be written
        File [default=None]: an open, writable Python file object
        Row [default=None]: a list of row labels
        Col [default=None]: a list of column labels
        close_file [default=0]: close the file after write (1) or leave open (0)
        separator [default='  ']: the column separator to use

        """
        fname = (
            'Write_array_'
            + self.ModelFile.replace('.psc', '')
            + '_'
            + time.strftime("%H:%M:%S")
        )

        # seeing as row indexes are now tuples and not arrays - brett 20050713
        if type(inputd) == tuple or type(inputd) == list:
            inputd = numpy.array(inputd)

        if Row != None:
            assert (
                len(Row) == inputd.shape[0]
            ), 'len(Row) must be equal to len(input_row)'
        if Col != None:
            assert (
                len(Col) == inputd.shape[1]
            ), 'len(Col) must be equal to len(input_col)'

        if File != None:
            # assert File != file, 'WriteArray(input,File=None,Row=None,Col=None,close_file=0)'
            if self.__settings__['write_array_header']:
                File.write('\n## ' + fname + '\n')
            if self.__settings__['write_array_spacer']:
                File.write('\n')
            if Col != None:
                print('Writing column')
                File.write('#')
                try:
                    input_width = (
                        len(self.__settings__['mode_number_format'] % inputd[0, 0])
                        + 1
                        + len(str(separator))
                    )
                except:
                    input_width = (
                        len(self.__settings__['mode_number_format'] % inputd[0])
                        + 1
                        + len(str(separator))
                    )
                for x in Col:
                    if len(str(x)) <= len(str(inputd[0, 0])):
                        spacer = (input_width - len(str(x))) * ' '
                        File.write(str(x) + spacer)
                    else:
                        spacer = (len(str(separator))) * ' '
                        File.write(str(x[:input_width]) + spacer)
                File.write('\n')
            try:
                print('\nWriting array (normal) to file')
                # scipy.io.write_array(File,input,separator=' ',linesep='\n',precision=3,keep_open=1)
                flush_count = 0
                if self.__settings__['write_arr_lflush'] < 2:
                    print('INFO: LineFlush must be >= 2')
                    self.__settings__['write_arr_lflush'] = 2
                for x in range(inputd.shape[0]):
                    flush_count += 1
                    for y in range(inputd.shape[1]):
                        if inputd[x, y] < 0.0:
                            File.write(
                                self.__settings__['mode_number_format'] % inputd[x, y]
                            )
                        else:
                            File.write(
                                ' '
                                + self.__settings__['mode_number_format'] % inputd[x, y]
                            )
                        if y < inputd.shape[1] - 1:
                            File.write(separator)
                    File.write('\n')
                    if flush_count == self.__settings__['write_arr_lflush']:
                        File.flush()
                        flush_count = 0
                        # print 'flushing'
                # File.write('\n')
            except IndexError as e:
                print('\nWriting vector (normal) to file')
                for x in range(len(inputd)):
                    File.write(self.__settings__['mode_number_format'] % inputd[x])
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

    def Write_array_latex(self, inputa, File=None, Row=None, Col=None, close_file=0):
        """
        Write_array_latex(input,File=None,Row=None,Col=None,close_file=0)

        Write an array to an open file as a \'LaTeX\' {array}

        Arguments:

        input: the array to be written
        File [default=None]: an open, writable Python file object
        Row [default=None]: a list of row labels
        Col [default=None]: a list of column labels
        close_file [default=0]: close the file after write (1) or leave open (0)

        """
        # seeing as row indexes are now tuples and not arrays - brett 20050713
        if type(inputa) == tuple or type(inputa) == list:
            inputa = numpy.array(inputa)

        try:
            a = inputa.shape[1]
        except:
            inputa.shape = (1, inputa.shape[0])

        if Row != None:
            assert (
                len(Row) == inputa.shape[0]
            ), 'len(Row) must be equal to len(input_row)'
            print('Writing row')
        if Col != None:
            assert (
                len(Col) == inputa.shape[1]
            ), 'len(Col) must be equal to len(input_col)'
            print('Writing column')

        if self.__settings__['write_arr_lflush'] < 2:
            print('INFO: LineFlush must be >= 2')
            self.__settings__['write_arr_lflush'] = 2

        fname = (
            'Write_array_latex_'
            + self.ModelFile.replace('.psc', '')
            + '_'
            + time.strftime("%H:%M:%S")
        )
        if File != None:
            # assert File != file, 'WriteArray(input,File=None,Row=None,Col=None,close_file=0)'
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
                            elx = str(Col[el]).replace('_', '\_')
                            if Row != None:
                                File.write(' & $\\small{' + elx + '}$')
                            else:
                                if el == 0:
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
                        # print 'flushing'
                # File.write('\n')
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

    def Write_array_html(
        self, inputa, File=None, Row=None, Col=None, name=None, close_file=0
    ):
        """
        Write_array_html(input,File=None,Row=None,Col=None,name=None,close_file=0)

        Write an array as an HTML table (no header/footer) or complete document. Tables
        are formatted with coloured columns if they exceed a specified size.

        mod.__settings__['write_array_html_header']: write the HTML document header
        mod.__settings__['write_array_html_footer']: write the HTML document footer

        Arguments:

        input: the array to be written
        File [default=None]: an open, writable Python file object
        Row [default=None]: a list of row labels
        Col [default=None]: a list of column labels
        name [default=None]: an HTML table description line
        close_file [default=0]: close the file after write (1) or leave open (0)

        """
        # seeing as row indexes are now tuples and not arrays - brett 20050713
        if type(inputa) == tuple or type(inputa) == list:
            inputa = numpy.array(inputa)
        try:
            a = inputa.shape[1]
        except:
            inputa.shape = (1, inputa.shape[0])

        if Row != None:
            assert (
                len(Row) == inputa.shape[0]
            ), 'len(Row) must be equal to len(input_row)'
            print('Writing row')
        if Col != None:
            assert (
                len(Col) == inputa.shape[1]
            ), 'len(Col) must be equal to len(input_col)'
            print('Writing column')
        if self.__settings__['write_arr_lflush'] < 2:
            print('INFO: LineFlush must be >= 2')
            self.__settings__['write_arr_lflush'] = 2
        if name != None:
            fname = (
                'PySCeS data "'
                + name
                + '" generated from model file: '
                + self.ModelFile
                + ' '
                + time.strftime("%H:%M:%S")
            )
        else:
            fname = (
                'PySCeS data generated from model file: '
                + self.ModelFile
                + ' '
                + time.strftime("%H:%M:%S")
            )
        if File != None:
            # assert File != file, 'WriteArray(input,File=None,Row=None,Col=None,close_file=0)'
            header = '\n'
            header += (
                '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">\n'
            )
            header += '<html>\n'
            header += '<head>\n'
            header += (
                '<title>PySCeS data generated from: '
                + self.ModelFile
                + ' - '
                + time.strftime("%H:%M:%S (%Z)")
                + '</title>\n'
            )
            header += '<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">\n'
            header += '</head>\n'
            header += '<body>\n\n'
            header += (
                '<h4><a href="http://pysces.sourceforge.net">PySCeS</a> data generated from: '
                + self.ModelFile
                + '</h4>\n\n'
            )

            if self.__settings__['write_array_html_header']:
                File.write(header)
            File.write('<!-- ' + fname + '-->\n')
            if name != None:
                File.write(
                    '\n<p>\n<font face="Arial, Helvetica, sans-serif"><strong>'
                    + name
                    + '</strong></font>'
                )
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
                            '  <td bgcolor="#CCCCCC"><div align="center"><b>'
                            + str(col)
                            + '</b></div></td>\n'
                        )
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
                    if (
                        rowcntr == colour_count_row
                        and inputa.shape[0] > 3 * colour_count_row
                    ):
                        File.write('\n<tr bgcolor="' + html_colour + '">\n')
                        rowcntr = 0
                    else:
                        File.write('\n<tr>\n')
                    if Row != None:
                        File.write(
                            '  <td bgcolor="#CCCCCC"><div align="center"><b>'
                            + str(Row[x])
                            + '</b></div></td>\n'
                        )
                    colcntr = 0
                    for y in range(inputa.shape[1]):
                        colcntr += 1
                        if (
                            colcntr == colour_count_col
                            and inputa.shape[1] > 3 * colour_count_row
                        ):
                            File.write(
                                '  <td nowrap bgcolor="'
                                + html_colour
                                + '">'
                                + self.__settings__['write_array_html_format']
                                % inputa[x, y]
                                + '</td>\n'
                            )
                            colcntr = 0
                        else:
                            File.write(
                                '  <td nowrap>'
                                + self.__settings__['write_array_html_format']
                                % inputa[x, y]
                                + '</td>\n'
                            )
                    ##                        File.write('  <td nowrap>' + self.__settings__['write_array_html_format'] % input[x,y] + '</td>\n')
                    if Row != None and inputa.shape[1] + 1 >= double_index:
                        File.write(
                            '  <td bgcolor="#CCCCCC"><div align="center"><b>'
                            + Row[x]
                            + '</b></div></td>\n'
                        )
                    File.write('</tr>')

                    if flush_count == self.__settings__['write_arr_lflush']:
                        File.flush()
                        flush_count = 0
                        # print 'flushing'
                if Col != None and inputa.shape[0] + 1 >= double_index:
                    File.write('\n<tr>\n')
                    if Row != None:
                        File.write('  <td>&nbsp;</td>\n')
                    for col in Col:
                        File.write(
                            '  <td bgcolor="#CCCCCC"><div align="center"><b>'
                            + col
                            + '</b></div></td>\n'
                        )
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
                        + __version__
                        + '</font></a><font size="2"> HTML output (model <i>'
                        + self.ModelFile
                        + '</i> analysed at '
                        + time.strftime("%H:%M:%S")
                        + ' by <i>'
                        + getuser()
                        + '</i>)</font></p>\n'
                    )
                except:
                    File.write(
                        '<p><a href="http://pysces.sourceforge.net"><font size="3">PySCeS '
                        + __version__
                        + '</font></a><font size="2"> HTML output (model <i>'
                        + self.ModelFile
                        + '</i> analysed at '
                        + time.strftime("%H:%M:%S - %Z")
                        + ')</font></p>\n'
                    )
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

    def SimPlot(
        self, plot='species', filename=None, title=None, log=None, format='lines'
    ):
        """
        Plot the simulation results, uses the new UPI pysces.plt interface:

        - *plot*: output to plot (default='species')

         - 'all' rates and species
         - 'species' species
         - 'rates' reaction rates
         - `['S1', 'R1', ]` a list of model attributes (species, rates)

        - *filename* if not None file is exported to filename (default=None)
        - *title* the plot title (default=None)
        - *log* use log axis for 'x', 'y', 'xy' (default=None)
        - *format* line format, backend dependant (default='')
        """
        data = None
        labels = None
        allowedplots = ['all', 'species', 'rates']
        if type(plot) != list and plot not in allowedplots:
            raise RuntimeError(
                '\nPlot must be one of {} not \"{}\"'.format(str(allowedplots), plot)
            )
        if plot == 'all':
            data, labels = self.data_sim.getAllSimData(lbls=True)
        elif plot == 'species':
            data, labels = self.data_sim.getSpecies(lbls=True)
        elif plot == 'rates':
            data, labels = self.data_sim.getRates(lbls=True)
        else:
            plot = [
                at
                for at in plot
                if at
                in self.__species__ + self.__reactions__ + [self.data_sim.time_label]
            ]
            kwargs = {'lbls': True}
            if len(plot) > 0:
                data, labels = self.data_sim.getSimData(*plot, **kwargs)
        del allowedplots
        self.__SIMPLOT_OUT__ = labels
        plt.plotLines(
            data, 0, list(range(1, data.shape[1])), titles=labels, formats=[format]
        )
        # set the x-axis range so that it is original range + 0.2*sim_end
        # this is a sceintifcally dtermned amount of space that is needed for the title at the
        # end of the line :-) - brett 20040209
        RngTime = self.data_sim.getTime()
        end = RngTime[-1] + 0.2 * RngTime[-1]
        plt.setRange('x', RngTime[0], end)
        del RngTime

        if self.__KeyWords__['Output_In_Conc']:
            M = 'Concentration'
        else:
            M = (
                'Amount (%(multiplier)s x %(kind)s x 10**%(scale)s)**%(exponent)s'
                % self.__uDict__['substance']
            )
        xu = (
            'Time (%(multiplier)s x %(kind)s x 10**%(scale)s)**%(exponent)s'
            % self.__uDict__['time']
        )
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
        if title == None:
            plt.setGraphTitle(
                'PySCeS Simulation ('
                + self.ModelFile
                + ') '
                + time.strftime("%a, %d %b %Y %H:%M:%S")
            )
        else:
            plt.setGraphTitle(title)
        plt.replot()
        if filename != None:
            plt.export(filename, directory=self.ModelOutput, type='png')

    def Scan1(self, range1=[], runUF=0):
        """
        Scan1(range1=[],runUF=0)

        Perform a single dimension parameter scan using the steady-state solvers. The parameter to be scanned is defined (as a model attribute "P") in mod.scan_in while the required output is entered into the list mod.scan_out.
        Results of a parameter scan can be easilly viewed with Scan1Plot().

        mod.scan_in - a model attribute written as in the input file (eg. P, Vmax1 etc)
        mod.scan_out -  a list of required output ['A','T2', 'ecR1_s1', 'ccJR1_R1', 'rcJR1_s1', ...]
        mod.scan_res -  the results of a parameter scan
        mod.scan - numpy record array with the scan results (scan_in and scan_out), call as mod.scan.Vmax, mod.scan.A_ss, mod.scan.J_R1, etc.
        mod.__settings__["scan1_mca_mode"] - force the scan algorithm to evaluate the elasticities (1) and control coefficients (2)
        (this should also be auto-detected by the Scan1 method).

        Arguments:

        range1 [default=[]]: a predefined range over which to scan.
        runUF [default=0]: run (1) the user defined function mod.User_Function (!U) before evaluating the steady state.

        """
        if self.__settings__["scan1_dropbad"] != 0:
            print(
                '\n****\nINFO: Dropping invalid steady states can mask interesting behaviour\n****\n'
            )

        # check the legitimacy of the input and output lists
        # self.__settings__["scan1_mca_mode"] = scanMCA
        run = 1
        # we now accept parameters and intial species values as scanable values - brett 20060519
        parameters = list(self.__parameters__) + [a + '_init' for a in self.__species__]
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
            if x[:2] == 'ec' and self.__settings__["scan1_mca_mode"] < 2:
                self.__settings__["scan1_mca_mode"] = 1
            elif x[:2] == 'cc':
                self.__settings__["scan1_mca_mode"] = 2
            elif x[:2] == 'rc':
                self.__settings__["scan1_mca_mode"] = 3

        if self.__settings__["scan1_mca_mode"] == 1:
            self.doElas()
        elif self.__settings__["scan1_mca_mode"] == 2:
            self.doMca()
        elif self.__settings__["scan1_mca_mode"] == 3:
            self.doMcaRC()
        else:
            self.doState()  # we are going to run a whole bunch anyway

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
                        "INFO: using steady-state flux for reaction ({} --> J_%s)".format(
                            x, x
                        )
                    )
                    self.scan_out[self.scan_out.index(x)] = 'J_{}'.format(x)
                elif x.startswith('J_'):
                    getattr(self.data_sstate, x[2:])
                elif x in self.species:
                    getattr(self.data_sstate, x)
                    print(
                        "INFO: using steady-state concentration for species ({} --> {}_ss)".format(
                            x, x
                        )
                    )
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

        assert (
            len(self.scan_out) != 0
        ), "Output parameter list (mod.scan_out) empty - do the model attributes exist ... see manual for details."

        # print run, self.scan_in, self.scan_out

        self._scan = None
        result = []
        cntr = 0
        cntr2 = 1
        cntr3 = 0

        # reset scan error controls
        self.__scan_errors_par__ = None
        self.__scan_errors_idx__ = None

        if self.__settings__['scan1_mesg']:
            print('\nScanning ...')
        if len(self.scan_in) > 0 and run == 1:
            badList = []
            badList_idx = []
            if self.__settings__['scan1_mesg']:
                print(len(range1) - (cntr * cntr2), end=' ')
            for xi in range(len(range1)):
                x = range1[xi]
                setattr(self, self.scan_in, x)
                ##  exec('self.' + self.scan_in + ' = ' + `x`)
                if self.__settings__["scan1_mca_mode"] == 1:
                    self.doElas()
                elif self.__settings__["scan1_mca_mode"] == 2:
                    self.doMca()
                elif self.__settings__["scan1_mca_mode"] == 3:
                    self.doMcaRC()
                else:
                    self.State()
                rawres = [x]
                # these two lines are going to be terminated - brett
                # rawres.append(x)
                # rawres.insert(0,x)
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

                # The following is for user friendly reporting:
                # next we check if the state is ok : if bad report it and add it : if good add it
                if not self.__StateOK__ and self.__settings__["scan1_dropbad"] != 0:
                    pass
                elif not self.__StateOK__:
                    if self.__settings__["scan1_nan_on_bad"]:
                        result.append([numpy.NaN] * len(rawres))
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
                        print(len(range1) - (cntr * cntr2), end=' ')
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
            print(
                '\nINFO: '
                + str(len(badList))
                + ' invalid steady states detected at '
                + self.scan_in
                + ' values:'
            )
            print(badList)
        if len(result) == 0:
            self.scan_res = numpy.zeros((len(range1), len(self.scan_out) + 1))
        else:
            self.scan_res = numpy.array(result)

    @property
    def scan(self):
        if self._scan is None and self.scan_res is not None:
            self._scan = numpy.rec.fromrecords(
                self.scan_res, names=[self.scan_in] + self.scan_out
            )
        return self._scan

    def Scan1Plot(self, plot=[], title=None, log=None, format='lines', filename=None):
        """
        Plot the results of a parameter scan generated with **Scan1()**

        - *plot* if empty mod.scan_out is used, otherwise any subset of mod.scan_out (default=[])
        - *filename* the filename of the PNG file (default=None, no export)
        - *title* the plot title (default=None)
        - *log* if None a linear axis is assumed otherwise one of ['x','xy','xyz'] (default=None)
        - *format* the backend dependent line format (default='lines')  or the *CommonStyle* 'lines' or 'points'.
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

        plotidx = [self.scan_out.index(c) + 1 for c in plot]
        labels = [self.scan_in] + self.scan_out
        plt.plotLines(data, 0, plotidx, titles=labels, formats=[format])

        end = data[-1, 0] + 0.2 * data[-1, 0]
        plt.setRange('x', data[0, 0], end)
        plt.setAxisLabel('x', self.scan_in)

        yl = 'Steady-state variable'
        plt.setAxisLabel('y', yl)
        if log != None:
            plt.setLogScale(log)
        if title == None:
            plt.setGraphTitle(
                'PySCeS Scan1 ('
                + self.ModelFile
                + ') '
                + time.strftime("%a, %d %b %Y %H:%M:%S")
            )
        else:
            plt.setGraphTitle(title)
        plt.replot()
        if filename != None:
            plt.export(filename, directory=self.ModelOutput, type='png')

    def Scan2D(self, p1, p2, output, log=False):
        """
        Generate a 2 dimensional parameter scan using the steady-state solvers.

        - *p1* is a list of [parameter1, start, end, points]
        - *p2* is a list of [parameter2, start, end, points]
        - *output* steady-state variable/properties e.g. 'J_R1', 'A_ss', 'ecR1_s1'
        - *log* scan using log ranges for both axes
        """

        for p in [p1, p2]:
            if not hasattr(self, p[0]):
                raise RuntimeError('\"{}\" is not a valid model attribute'.format(p[0]))

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
        Plot the results of a 2D scan generated with Scan2D

        - *filename* the filename of the PNG file (default=None, no export)
        - *title* the plot title (default=None)
        - *log* if None a linear axis is assumed otherwise one of ['x','xy','xyz'] (default=None)
        - *format* the backend dependent line format (default='lines')  or the *CommonStyle* 'lines' or 'points'.
        """

        plt.splot(
            self.__scan2d_results__, 0, 1, 2, titles=self.__scan2d_pars__, format=format
        )
        plt.setRange('xyz')
        plt.setAxisLabel('x', self.__scan2d_pars__[0])
        plt.setAxisLabel('y', self.__scan2d_pars__[1])
        plt.setAxisLabel('z', 'Steady-state variable')
        if log != None:
            plt.setLogScale(log)
        if title == None:
            plt.setGraphTitle(
                'PySCeS Scan2D ('
                + self.ModelFile
                + ') '
                + time.strftime("4a, %d %b %Y %H:%M:%S")
            )
        else:
            plt.setGraphTitle(title)
        plt.replot()
        if filename != None:
            plt.export(filename, directory=self.ModelOutput, type='png')

    def SetQuiet(self):
        """
        SetQuiet()

        Turn off as much solver reporting noise as possible:
        mod.__settings__['hybrd_mesg'] = 0
        mod.__settings__['nleq2_mesg'] = 0
        mod.__settings__["lsoda_mesg"] = 0
        mod.__settings__['mode_state_mesg'] = 0
        mod.__settings__['scan1_mesg'] = 0
        mod.__settings__['solver_switch_warning'] = False

        Arguments:
        None

        """
        self.__settings__['hybrd_mesg'] = 0
        self.__settings__['nleq2_mesg'] = 0
        self.__settings__["lsoda_mesg"] = 0
        self.__settings__['mode_state_mesg'] = 0
        self.__settings__['scan1_mesg'] = 0
        self.__settings__['solver_switch_warning'] = False

    def SetLoud(self):
        """
        SetLoud()

        Turn on as much solver reporting noise as possible:
        mod.__settings__['hybrd_mesg'] = 1
        mod.__settings__['nleq2_mesg'] = 1
        mod.__settings__["lsoda_mesg"] = 1
        mod.__settings__['mode_state_mesg'] = 1
        mod.__settings__['scan1_mesg'] = 1
        mod.__settings__['solver_switch_warning'] = True

        Arguments:
        None

        """
        self.__settings__['hybrd_mesg'] = 1
        self.__settings__['nleq2_mesg'] = 1
        self.__settings__["lsoda_mesg"] = 1
        self.__settings__['mode_state_mesg'] = 1
        self.__settings__['scan1_mesg'] = 1
        self.__settings__['solver_switch_warning'] = True

    def clone(self):
        """
        Returns a deep copy of this model object (experimental!)

        """

        print('\nINFO: Cloning function is currently experimental ... use with care.')
        return copy.deepcopy(self)


class WasteManagement(object):
    _gc_delay_hack = (
        gplt  # dont f$%^&*)ing even ask ... only 5 hrs to hack this workaround
    )


__waste_manager = WasteManagement()

# if __psyco_active__:
# psyco.bind(PysMod)

if __name__ == '__main__':
    print('\nTo use PySCeS import it from a Python Shell: \n\timport pysces\n')
