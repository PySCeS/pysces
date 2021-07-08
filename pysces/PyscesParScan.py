"""
PySCeS - Python Simulator for Cellular Systems (http://pysces.sourceforge.net)

Copyright (C) 2004-2020 B.G. Olivier, J.M. Rohwer, J.-H.S Hofmeyr all rights reserved,

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
import sys, os, pickle

flush = sys.stdout.flush
from time import sleep, time
import subprocess

try:
    import ipyparallel
except ImportError as ex:
    print()
    print(ex)
    raise ImportError(
        'PARSCANNER: Requires ipyparallel version >=4.0 (http://ipython.org).'
    )

MULTISCANFILE = __file__.replace(
    'PyscesParScan', '_multicorescan'
)
# print 'MULTISCANFILE is', MULTISCANFILE
__psyco_active__ = 0


class ParScanner(Scanner):
    """
    Arbitrary dimension generic distributed scanner.
    Subclassed from pysces.PyscesScan.Scanner.
    This class is initiated with a loaded PySCeS model and then allows
    the user to define scan parameters, see self.addScanParameter()
    and user output, see self.addUserOutput().
    Steady-state results are always stored in self.SteadyStateResults while
    user output can be found in self.UserOutputResults.
    Distributed (parallel) execution is achieved with the clustering capability
    of IPython. See "ipcluster --help".
    """

    # --johann 20101206
    genOn = True
    _MODE_ = 'state'
    HAS_USER_OUTPUT = False
    nan_on_bad_state = True
    MSG_PRINT_INTERVAL = 500
    __AnalysisModes__ = ('state', 'elasticity', 'mca', 'stability')
    invalid_state_list = None
    invalid_state_list_idx = None
    scans_per_run = 100

    def __init__(self, mod, engine='multiproc'):
        """
        Instantiate the parallel scanner class with a PySCeS model instance
        and an optional 'engine' argument specifying the parallel engine:
        'multiproc' -- multiprocessing (default)
        'ipcluster' -- IPython cluster
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
                raise OSError(
                    str(ex)
                    + '\nPARSCANNER: Requires a running IPython cluster. See "ipcluster --help".\n'
                )
            dv = rc[:]  # direct view
            lv = rc.load_balanced_view()
            self.dv = dv
            self.lv = lv
            dv.execute('from pysces.PyscesParScan import Analyze, setModValue')
        else:
            raise UserWarning(engine + " is not a valid parallel engine!")

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
        Generates the parameter scan space, partitioned according to self.scans_per_run
        """
        spr = self.scans_per_run
        Tsteps = 1
        for gen in self.GenOrder:
            if self.GenDict[gen][4] == False:  # don't increase Tsteps for followers
                Tsteps *= self.GenDict[gen][2]
        for step in range(Tsteps):
            pars = self.__nextValue__()
            # if pars not in self.ScanSpace:
            self.ScanSpace.append(pars)
        self.ScanSpace = np.array(self.ScanSpace)
        self.SeqArray = np.arange(1, Tsteps + 1)
        if Tsteps % spr == 0:
            numparts = Tsteps // spr
        else:
            numparts = Tsteps // spr + 1
        self.ScanPartition = [
            self.ScanSpace[n * spr : (n + 1) * spr] for n in range(numparts)
        ]
        self.SeqPartition = [
            self.SeqArray[n * spr : (n + 1) * spr] for n in range(numparts)
        ]
        self.Tsteps = Tsteps

    def Prepare(self, ReRun=False):
        """
        Internal method to prepare the parameters and generate ScanSpace.
        """
        print("\nPREPARATION\n-----------")
        self.scanT.normal_timer('PREP')
        self._MODE_ = self._MODE_.lower()
        assert self._MODE_ in self.__AnalysisModes__, (
            '\nSCANNER: \"%s\" is not a valid analysis mode!' % self._MODE_
        )
        if ReRun:
            self.ScanSpace = []
            self.UserOutputResults = []
            self.SteadyStateResults = []
        self.invalid_state_list = []
        self.invalid_state_list_idx = []
        # if self.nan_on_bad_state:
        # self.mod.__settings__["mode_state_nan_on_fail"] = True
        # self.dv.execute('mod.__settings__["mode_state_nan_on_fail"] = True')
        # generate the scan space
        self.genScanSpace()
        print("Generated ScanSpace:", next(self.scanT.PREP))
        print('PARSCANNER: Tsteps', self.Tsteps)
        flush()

    def Run(self, ReRun=False):
        """
        Run the parameter scan with a load balancing task client.
        """
        self.Prepare(ReRun)
        # this is where the parallel magic fun starts....
        arl = []  # asynchronous results list
        if self.engine == 'multiproc':
            fN = str(time()).split('.')[0]
            with open(fN, 'wb') as F:
                pickle.dump(
                    (
                        self.mod,
                        self.ScanPartition,
                        self.SeqPartition,
                        self.GenOrder,
                        self.UserOutputList,
                    ),
                    F,
                    protocol=-1,
                )
            fN = os.path.abspath(fN)
            print("Preparation completed:", next(self.scanT.PREP))
            self.scanT.normal_timer('RUN')
            subprocess.call([sys.executable, MULTISCANFILE, self._MODE_, fN])
            with open(fN, 'rb') as F:
                res_list = pickle.load(F)
            os.remove(fN)
            for result in res_list:
                self.StoreData(result)
        elif self.engine == 'ipcluster':
            for i in range(len(self.ScanPartition)):
                arl.append(
                    self.lv.apply(
                        Analyze,
                        self.ScanPartition[i],
                        self.SeqPartition[i],
                        self.GenOrder,
                        self._MODE_,
                        self.UserOutputList,
                        self.mod,
                    )
                )
            print("Submitted tasks:", len(arl))
            print("Preparation completed:", next(self.scanT.PREP))
            print("\nPARALLEL COMPUTATION\n--------------------")
            flush()
            self.scanT.normal_timer('RUN')
            while self.lv.queue_status()['unassigned'] > 0:
                sleep(5)
                print('Tasks to go... ', self.lv.queue_status()['unassigned'])
            # wait until all tasks are completed
            self.lv.wait()
            flush()
            print("\nGATHER RESULT\n-------------")
            flush()
            for ar in arl:
                result = ar.get()
                # tuple: 0 - state_species
                #   1 - state_flux
                #   2 - user_output_results
                #   3 - invalid_state_list
                #   4 - invalid_state_list_idx
                self.StoreData(result)
        print("Parallel calculation completed:", next(self.scanT.RUN))
        self.GatherScanResult()

    def RunScatter(self, ReRun=False):
        """
        Run the parameter scan by using scatter and gather for the ScanSpace.
        Not load balanced, equal number of scan runs per node.
        """
        if self.engine != 'ipcluster':
            print("RunScatter() only supported with ipcluster!")
            return
        self.Prepare(ReRun)
        # this is where the parallel magic fun starts....
        # push details into client namespace
        self.dv.push(
            {
                'GenOrder': self.GenOrder,
                'mode': self._MODE_,
                'mod': self.mod,
                'UserOutputList': self.UserOutputList,
            }
        )
        # scatter ScanSpace and SeqArray
        self.dv.scatter('partition', self.ScanSpace)
        self.dv.scatter('seqarray', self.SeqArray)
        print("Scattered ScanSpace on number of engines:", len(self.dv))
        print("Preparation completed:", next(self.scanT.PREP))
        print("\nPARALLEL COMPUTATION\n--------------------")
        flush()
        self.scanT.normal_timer('RUN')
        # executes scan on partitioned ScanSpace on every node
        self.dv.execute(
            'y=Analyze(partition,seqarray,GenOrder,mode,UserOutputList,mod)', block=True
        )
        print("Parallel calculation completed:", next(self.scanT.RUN))
        flush()

        ## this is analysis stuff
        print("\nGATHER RESULT\n-------------")
        flush()
        ar = self.dv.gather('y')
        results = ar.get()
        results = [results[i : i + 5] for i in range(0, len(results), 5)]
        for result in results:
            # tuple: 0 - state_species
            #   1 - state_flux
            #   2 - user_output_results
            #   3 - invalid_state_list
            #   4 - invalid_state_list_idx
            self.StoreData(result)
        print("Parallel calculation completed:", next(self.scanT.RUN))
        self.GatherScanResult()

    def StoreData(self, result):
        """
        Internal function which concatenates and stores single result generated by Analyze.

        - *result* IPython client result object
        """
        self.SteadyStateResults.append(
            np.hstack((np.array(result[0]), np.array(result[1])))
        )
        if self.HAS_USER_OUTPUT:
            self.UserOutputResults.append(np.array(result[2]))
        self.invalid_state_list += result[3]
        self.invalid_state_list_idx += result[4]

    def GatherScanResult(self):
        """
        Concatenates and combines output result fragments from the parallel scan.
        """
        # from here on we have the complete results
        self.SteadyStateResults = np.vstack([i for i in self.SteadyStateResults])
        if self.HAS_USER_OUTPUT:
            self.UserOutputResults = np.vstack([i for i in self.UserOutputResults])
        self.resetInputParameters()
        # print "Gather completed:", self.scanT.GATHER.next()
        print("\nPARSCANNER: %s states analysed" % len(self.SteadyStateResults))
        print("Total time taken: ", next(self.scanT.PREP))
        self.scanT.PREP.close()  # close timer
        self.scanT.RUN.close()
        if len(self.invalid_state_list) > 0:
            print('\nBad steady states encountered at:\n')
            print("Sequence: ", self.invalid_state_list_idx)
            print("Parameters: ", self.invalid_state_list)


# Utility methods for scatter/gather
# was in engineCode.py but refactored into main script
# --johann 20120813

# set model attribute values
def setModValue(mod, name, value):
    assert hasattr(mod, name), 'Model does not have an attribute: %s ' % name
    setattr(mod, name, float(value))


# analyze steady states
def Analyze(partition, seqarray, GenOrder, mode, UserOutputList, mod):
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
        user_output_results.append([getattr(mod, res) for res in UserOutputList])
    return (
        state_species,
        state_flux,
        user_output_results,
        invalid_state_list,
        invalid_state_list_idx,
    )
