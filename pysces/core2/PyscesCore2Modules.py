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

# Structural Analysis module
from pysces.PyscesStoich import Stoich
from .PyscesCore2 import StructMatrix


class PyscesEnhancedStoich(Stoich):
    """PySCeS stoichiometry class for use with core2"""

    N = None
    Nr = None
    K = None
    K0 = None
    L = None
    L0 = None
    Gamma = None

    def __init__(self, core):
        Stoich.__init__(self, core.stoichiometric_matrix.array)
        self.species = core.stoichiometric_matrix.row
        self.reactions = core.stoichiometric_matrix.col

    def getNullSpaces(self):
        self.AnalyseK()
        self.AnalyseL()

    def testNullSpaces(self):
        # TODO: build in nullspace validity checks from PyscesModel
        pass

    def setStructMatrices(self):
        self.N = StructMatrix(self.nmatrix, self.nmatrix_row, self.nmatrix_col)
        self.N.setRow(self.species)
        self.N.setCol(self.reactions)

        self.Nr = StructMatrix(self.nrmatrix, self.nrmatrix_row, self.nrmatrix_col)
        self.Nr.setRow(self.species)
        self.Nr.setCol(self.reactions)

        self.K = StructMatrix(self.kmatrix, self.kmatrix_row, self.kmatrix_col)
        self.K.setRow(self.reactions)
        self.K.setCol(self.reactions)

        self.K0 = StructMatrix(
            self.kzeromatrix, self.kzeromatrix_row, self.kzeromatrix_col
        )
        self.K0.setRow(self.reactions)
        self.K0.setCol(self.reactions)

        self.L = StructMatrix(self.lmatrix, self.lmatrix_row, self.lmatrix_col)
        self.L.setRow(self.species)
        self.L.setCol(self.species)

        self.L0 = StructMatrix(
            self.lzeromatrix, self.lzeromatrix_row, self.lzeromatrix_col
        )
        self.L0.setRow(self.species)
        self.L0.setCol(self.species)

        if self.info_moiety_conserve:
            self.Gamma = StructMatrix(
                self.conservation_matrix,
                self.conservation_matrix_row,
                self.conservation_matrix_col,
            )
            self.Gamma.setRow(self.species)
            self.Gamma.setCol(self.species)


class StructuralModule(object):
    core = None
    struct = None

    def setCore(self, core):
        self.core = core
        self.struct = None
        if self.core.stoichiometric_matrix == None:
            print("StructuralModule building stoichiometric matrix ...")
            self.core.setStoichiometricMatrix()

    def getCore(self):
        self.core.struct = self.struct
        return self.core

    def analyseStoichiometry(self):
        self.struct = PyscesEnhancedStoich(self.core)
        self.struct.getNullSpaces()
        self.struct.setStructMatrices()


# Integration Module
import numpy

"""
class StateDataObj(object):
    flux = None
    flux_labels = None
    species = None
    species_labels = None
    valid = True
    _suffix = None
    _prefix = None

    def __init__(self):
        self.species_labels = []
        self.flux_labels = []

    def setSpecies(self, name, value, suffix=None):
        if suffix != None:
            name = name + suffix
        if name not in self.species_labels:
            self.species_labels.append(name)
        self._suffix = suffix
        setattr(self, name, value)

    def setFlux(self, name, value, prefix=None):
        if prefix != None:
            name = prefix + name
        if name not in self.flux_labels:
            self.flux_labels.append(name)
        self._prefix = prefix
        setattr(self, name, value)

    def setAllSpecies(self, species_labels, species, suffix=None):
        assert len(species_labels) == len(species), '\nThis aint gonna work1'
        self.species_labels = []
        #  self.species_labels = tuple(species_labels)
        self.species = species.copy()
        for S in range(len(species_labels)):
            self.setSpecies(species_labels[S], species[S], suffix)

    def setAllFluxes(self, flux_labels, flux, prefix=None):
        assert len(flux_labels) == len(flux), '\nThis aint gonna work2'
        self.flux_labels = []
        #  self.flux_labels = tuple(flux_labels)
        self.flux = flux.copy()
        for J in range(len(flux_labels)):
            self.setFlux(flux_labels[J], flux[J], prefix)

    def getFlux(self, name):
        if prefix != None:
            name = prefix + name
        return getattr(self, name)

    def getSpecies(self, name):
        if suffix != None:
            name = name + suffix
        return getattr(self, name)
"""


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
                print('I don\'t have an attribute %s ... ignoring.' % roc)
        if not lbls:
            return output
        else:
            return numpy.array(output), lout


'''
class IntegrationDataObj(object):
    """
    This class is specifically designed to store the results of a time simulation
    It has methods for setting the Time, Labels, Species and Rate data and
    getting Time, Species and Rate (including time) arrays. However, of more use:
    - getOutput(*arg) feed this method species/rate labels and it will return
      an array of [time, sp1, r1, ....]
    - getDataAtTime(time) the data generated at time point "time".
    - getDataInTimeInterval(time, bounds=None) more intelligent version of the above
      returns an array of all data points where: time-bounds <= time <= time+bounds
      where bounds defaults to stepsize.

    time = None
    rates = None
    species = None
    rate_labels = None
    species_labels = None

    def setLabels(self, species, rates):
        """set the species and rate label lists"""
        self.species_labels = species
        self.rate_labels = rates

    def setTime(self, time):
        """Set the time vector"""
        self.time = time.reshape(len(time), 1)

    def setSpecies(self, species):
        """Set the species array"""
        self.species = species

    def setRates(self, rates):
        """set the rate array"""
        self.rates = rates

    def getTime(self):
        """return the time vector"""
        assert self.time != None, "\nNo time"
        return self.time.reshape(len(self.time),)

    def getSpecies(self):
        """return time+species array"""
        assert self.species != None, "\nNo species"
        return numpy.hstack((self.time, self.species))

    def getRates(self):
        """return time+rate array"""
        assert self.rates != None, "\nNo rates"
        return numpy.hstack((self.time, self.rates))

    def getDataAtTime(self, time):
        """Return all data generated at "time" """
        t = None
        sp = None
        ra = None
        temp_t = self.time.reshape(len(self.time),)
        for tt in range(len(temp_t)):
            if temp_t[tt] == time:
                t = tt
                if self.species is not None:
                    sp = self.species.take([tt], axis=0)
                if self.rates is not None:
                    ra = self.rates.take([tt], axis=0)
                break
        output = None
        if t is not None:
            output = numpy.array([[temp_t[t]]])
            if sp is not None:
                output = numpy.hstack((output,sp))
            if ra is not None:
                output = numpy.hstack((output,ra))
        return output

    def getDataInTimeInterval(self, time, bounds=None):
        """
         getDataInTimeInterval(time, bounds=None) returns an array of all
         data points where: time-bounds <= time <= time+bounds
         where bound defaults to stepsize
        """
        temp_t = self.time.reshape(len(self.time),)
        if bounds == None:
            bounds = temp_t[1] - temp_t[0]
        c1 = (temp_t >= time-bounds)
        c2 = (temp_t <= time+bounds)
        print 'Searching (%s:%s:%s)' % (time-bounds, time, time+bounds)

        t = []
        sp = None
        ra = None

        for tt in range(len(c1)):
            if c1[tt] and c2[tt]:
                t.append(tt)
        output = None
        if len(t) > 0:
            output = self.time.take(t)
            output = output.reshape(len(output),1)
            if self.species is not None:
                output = numpy.hstack((output, self.species.take(t, axis=0)))
            if self.rates is not None:
                output = numpy.hstack((output, self.rates.take(t, axis=0)))
        return output

    def getOutput(self, *args):
        """getOutput(*arg) feed this method species/rate labels and it
        will return an array of [time, sp1, r1, ....]
        """
        output = self.time
        for roc in args:
            if roc in self.species_labels:
                assert self.species != None, "\nNo species"
                output = numpy.hstack((output, self.species.take([self.species_labels.index(roc)], axis=-1)))
            if roc in self.rate_labels:
                assert self.rates != None, "\nNo rates"
                output = numpy.hstack((output, self.rates.take([self.rate_labels.index(roc)], axis=-1)))
        return output
'''


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
        print('Searching (%s:%s:%s)' % (time - bounds, time, time + bounds))

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

    def getOutput(self, *args):
        """
        Old alias for getSimData()
        getOutput(\*args) feed this method species/rate labels and it
        will return an array of [time, sp1, r1, ....]
        """
        return self.getSimData(*args)

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


# TODO:
class IntegrationBase(object):
    name = None
    core = None
    data = None

    sim_start = None
    sim_end = None
    sim_point = None
    initial_value_vector = None

    def setName(self, name):
        self.name = name

    def getName(self):
        return self.name

    def setCore(self, core):
        self.core = core
        self.data = IntegrationData()
        self.data.setLabels(self.core.hasVariableSpecies(), self.core.hasReactions())

    def getCore(self):
        return self.core

    def getData(self):
        return self.data
