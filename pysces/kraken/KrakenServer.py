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

# begin header
import pysces
import numpy, scipy

from .KrakenNET import os, time, socket
from .KrakenNET import StatusServer, BasicServer
from .KrakenNET import HOSTNAME, BLOCK_SIZE, STATUS_PORT, PYSCES_PORT, CFSERVE_PORT

print('I AM:', HOSTNAME, 'using BLOCK_SIZE: %s' % BLOCK_SIZE)

# end header

__psyco_active__ = 0

from pysces import TimerBox


class PyscesServer(BasicServer):
    COMMAND_LIST = (
        'P_INIT',
        'P_LOAD',
        'P_SET_STR',
        'P_SET_FLOAT',
        'P_STATE',
        'P_SIM',
        'P_ELAS',
        'P_MCA',
        'P_SCAN',
        'P_CONTINUATION',
        'P_CONTINUATION_WITH_EIGEN',
        'P_SET_SETTINGS_STRING',
        'P_SET_SETTINGS_FLOAT',
        'P_CONTINUATION_WITH_EIGEN_USER',
        'P_CONTINUATION2D_WITH_EIGEN_USER',
    )
    model_name = None
    model = None

    def __init__(self, port, block_size, status_server=None, myname=None):
        BasicServer.__init__(self, port, block_size, status_server, myname)

    def P_INIT(self, *args):
        args = args[0]
        print('Args', args)
        self.RESULT = None
        self.setStatus('INITIALISING')
        client = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        client.connect((self.client_address[0], CFSERVE_PORT))
        client.send('GET')
        GO = True
        model = ''
        while GO:
            data = client.recv(self.block_size)
            model += data
            if data == '':
                GO = False
        self.model_name = args[0]
        self.model = pysces.model(self.model_name, loader="string", fString=model)
        if self.debug:
            print('P_INIT ... ok.')
        self.setStatus('DONE_INITIALISE')
        return True

    def P_LOAD(self, *args):
        print('Loading:', self.model_name)
        self.RESULT = None
        self.model.doLoad()
        self.model.SetQuiet()
        if self.debug:
            print('P_LOAD ... ok.')
        self.setStatus('DONE_LOAD')
        return True

    def P_SET_SETTINGS_FLOAT(self, *args):
        args = args[0]
        print('Setting:', args)
        self.RESULT = None
        self.model.__settings__.update({args[0]: float(args[1])})
        ##  setattr(self.model, args[0], args[1])
        print(self.model.__settings__[args[0]])
        return True

    def P_SET_SETTINGS_STRING(self, *args):
        args = args[0]
        print('Setting:', args)
        self.RESULT = None
        self.model.__settings__.update({args[0]: args[1]})
        ##  setattr(self.model, args[0], args[1])
        print(self.model.__settings__[args[0]])
        return True

    def P_SET_STR(self, *args):
        args = args[0]
        print('Setting:', args)
        self.RESULT = None
        setattr(self.model, args[0], args[1])
        print(getattr(self.model, args[0]))
        return True

    def P_SET_FLOAT(self, *args):
        # note to self in future investigate numpy.float64 for this
        args = args[0]
        print('Setting:', args)
        self.RESULT = None
        setattr(self.model, args[0], float(args[1].strip()))
        print(getattr(self.model, args[0]))
        return True

    def P_STATE(self, *args):
        self.model.doState()
        self.RESULT = (self.model.state_species, self.model.state_flux)
        if self.debug:
            print('P_STATE ... ok.')
        self.setStatus('DONE_STATE')
        return True

    def P_ELAS(self, *args):
        self.model.doElas()
        self.RESULT = (self.model.elas_var, self.model.elas_par)
        if self.debug:
            print('P_ELAS ... ok.')
        self.setStatus('DONE_ELAS')
        return True

    def P_MCA(self, *args):
        self.model.doMca()
        self.RESULT = (self.model.mca_ci, self.model.mca_cjd, self.model.mca_csd)
        if self.debug:
            print('P_MCA ... ok.')
        self.setStatus('DONE_MCA')
        return True

    def P_SIM(self, *args):
        """Args: (float(time_end), int(points))"""
        self.setStatus('SIMULATING')
        args = args[0]
        print('Args', args)
        self.model.doSim(float(args[0]), float(args[1]))
        self.RESULT = self.model.sim_res
        if self.debug:
            print('P_STATE ... ok.')
        self.setStatus('DONE_SIM')
        return True

    def P_SCAN(self, *args):
        """Args: (par, out, par*[str(p),float(start),float(end),int(points),
                 log=(int(0/1)),follower=(int(0/1))], [out*str(o)])
        """
        self.setStatus('SCANNING')
        args = args[0]
        print('Args', args)
        commands = []
        cpoint = 1
        shift = 0
        for set in range(int(args[0])):
            ctmp = []
            for t in range(2, 8):
                cpoint += 1
                ctmp.append(args[t + shift].strip())
            commands.append(ctmp)
            shift += 6
        output = []
        for out in range(int(args[1])):
            cpoint += 1
            output.append(args[cpoint].strip())

        print('\n', commands)
        print(output, '\n')

        scan = pysces.Scanner(self.model)
        scan.quietRun = True
        for gen in commands:
            print(gen)
            print(
                gen[0],
                float(gen[1]),
                float(gen[2]),
                int(gen[3]),
                bool(int(gen[4])),
                bool(int(gen[5])),
            )

            scan.addGenerator(
                gen[0],
                float(gen[1]),
                float(gen[2]),
                int(gen[3]),
                log=bool(int(gen[4])),
                follower=bool(int(gen[5])),
            )
        scan.addUserOutput(*output)
        scan.Run()
        self.RESULT = scan.UserOutputResults
        self.setStatus('DONE_SCAN')
        del scan
        return True

    def P_CONTINUATION(self, *args):
        """
        [str(p), float(min), float(max), int(points), *float(y-value)]
        """
        self.setStatus('CONTINUING')
        args = args[0]
        print('Args', args)

        self.model.pitcon_par_space = scipy.logspace(
            scipy.log10(float(args[1].strip())),
            scipy.log10(float(args[2].strip())),
            int(args[3].strip()),
        )

        if len(args) == 4:
            self.RESULT = self.model.PITCON(args[0].strip())
        elif len(args) == 5:
            self.RESULT = self.model.PITCON(args[0].strip(), float(args[4].strip()))
        return True

    def P_CONTINUATION_WITH_EIGEN(self, *args):
        """
        [str(p), float(min), float(max), int(points), *float(y-value)]
        """
        self.setStatus('CONTINUING')
        args = args[0]
        print('Args', args)

        pscan = pysces.PITCONScanUtils(self.model)
        if len(args) == 4:
            pscan.runContinuation(
                args[0].strip(),
                float(args[1].strip()),
                float(args[2].strip()),
                int(args[3].strip()),
            )
        elif len(args) == 5:
            pscan.runContinuation(
                args[0].strip(),
                float(args[1].strip()),
                float(args[2].strip()),
                int(args[3].strip()),
                float(args[4].strip()),
            )
        pscan.analyseData(analysis='eigen')
        self.RESULT = (
            pscan.res_idx,
            pscan.res_metab,
            pscan.res_flux,
            pscan.getArrayListAsArray(pscan.res_eigen),
        )
        del pscan
        return True

    def P_CONTINUATION_WITH_EIGEN_USER(self, *args):
        """
        [str(p), float(min), float(max), int(points), float(y-value)]
        """
        self.setStatus('CONTINUING')
        args = args[0]
        print('Args', args)

        pscan = pysces.PITCONScanUtils(self.model)
        # TODO HARD CODED FOR NOW BUT NEEDS GENERALISATION
        pscan.setUserOuput('ecR1_A', 'ecR2_A')
        if len(args) == 4:
            pscan.runContinuation(
                args[0].strip(),
                float(args[1].strip()),
                float(args[2].strip()),
                int(args[3].strip()),
            )
        elif len(args) == 5:
            pscan.runContinuation(
                args[0].strip(),
                float(args[1].strip()),
                float(args[2].strip()),
                int(args[3].strip()),
                float(args[4].strip()),
            )
        pscan.analyseData(analysis='eigen')
        self.RESULT = (
            pscan.res_idx,
            pscan.res_metab,
            pscan.res_flux,
            pscan.getArrayListAsArray(pscan.res_eigen),
            pscan.res_user,
        )
        del pscan
        return True

    def P_CONTINUATION2D_WITH_EIGEN_USER(self, *args):
        """
        [str(p), float(min), float(max), int(points), float(y-value), str(userNx)]
        """
        self.setStatus('CONTINUING')
        args = args[0]
        print('Args', args)
        userout = []
        for a in range(len(args)):
            args[a] = args[a].strip()
            if a > 4:
                userout.append(args[a])

        pscan = pysces.PITCONScanUtils(self.model)
        if len(userout) > 0:
            pscan.setUserOuput(*userout)
        pscan.runContinuation(
            args[0], float(args[1]), float(args[2]), int(args[3]), float(args[4])
        )
        pscan.analyseData(analysis='all')
        self.RESULT = (
            pscan.res_idx,
            pscan.res_metab,
            pscan.res_flux,
            pscan.getArrayListAsArray(pscan.res_eigen),
            pscan.res_user,
        )
        del pscan
        return True


def StartPysCeSKrakenServer(
    STATUS_PORT=STATUS_PORT, PYSCES_PORT=PYSCES_PORT, BLOCK_SIZE=BLOCK_SIZE
):
    STATUSserve = StatusServer(STATUS_PORT, BLOCK_SIZE, 'KrakenStatus')
    STATUSserve.start()
    PYSCESserve = PyscesServer(PYSCES_PORT, BLOCK_SIZE, STATUSserve, 'KrakenPySCeS')
    PYSCESserve.start()
    PYSCESserve.join()


if __name__ == '__main__':
    ##  mod = pysces.model('isola_thesis')
    ##  mod.doLoad()
    ##  pscan = pysces.PITCONScanUtils(mod)
    ##  pscan.do2Dcontinuation('P',0.01,3000,5,par3d=0.1)
    ##  pscan.doEigenAnalysis()
    ##  print 'PITCON', pscan.pitcon_res[0]
    ##  print 'FLUX',   pscan.flux_res[0]
    ##  print 'EIGEN',  pscan.eigen_res[0]
    StartPysCeSKrakenServer()
    print('\nAll threads terminated ... exiting now.\n')
