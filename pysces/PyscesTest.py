"""
PySCeS - Python Simulator for Cellular Systems (http://sourceforge.net)

Copyright (C) 2004-2025 B.G. Olivier, J.M. Rohwer, J.-H.S Hofmeyr all rights reserved,

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

__doc__ = """PySCeS unittest module"""

import os
import sys
import shutil
import unittest
from time import sleep

from . import PyscesStoich, pitcon_switch

if pitcon_switch:
    from . import pitcon
from . import install_dir as INSTALL_DIR
from . import model_dir as MODEL_DIR
from . import model as PSCMODEL
from .version import __version__
from . import interface

# evaluation comparison stuff
from numpy.testing import assert_array_equal, assert_array_almost_equal, assert_allclose
from numpy import array, zeros, linspace, copy, allclose
from scipy.optimize import fsolve
from scipy.integrate import odeint


def CopyTestModels(
    dirIn=os.path.join(INSTALL_DIR, 'pscmodels'),
    dirOut=MODEL_DIR,
    overwrite=0,
):
    """
    CopyTestModels(dirIn=os.path.join(pysces_install,'pscmodels'),dirOut=pysces_model,overwrite=0)

    'Copy all PySCeS test model files in dirIn to dirOut, defaults to:
    in: pysces/pscmodels
    out: model_dir

    Arguments:
    =========
    dirIn [default=os.path.join(pysces_install,'pscmodels')]: target directory
    dirOut [default=pysces_model]: destination directory
    overwrite [default=0]: automatically (1) overwrite target files

    """
    if not os.path.exists(dirOut):
        os.makedirs(dirOut)

    if os.path.exists(dirIn):
        if os.path.exists(dirOut):
            print('src : ' + dirIn)
            print('dest: ' + dirOut)
            flist = os.listdir(dirIn)
            if overwrite == 0:
                flist2 = os.listdir(dirOut)
            maxlen = 0
            for x in flist:
                if len(x) > maxlen:
                    maxlen = len(x)
            if len(flist) != 0:
                for File in flist:
                    if File[:12] == 'pysces_test_' and (File[-4:] == '.psc' or File[-4:] == '.xml'):
                        if overwrite == 0:
                            try:
                                a = flist2.index(File)
                                # print File +  (maxlen-len(File))*'.' + ' skipped'
                            except:
                                shutil.copy(
                                    os.path.join(dirIn, File),
                                    os.path.join(dirOut, File),
                                )
                                print(File + (maxlen - len(File)) * '.' + ' ok')
                        else:
                            shutil.copy(
                                os.path.join(dirIn, File), os.path.join(dirOut, File)
                            )
                            print(File + (maxlen - len(File)) * '.' + ' ok')
            else:
                print('Empty directory?')
        else:
            print(dirOut + ' does not exist')
    else:
        print(dirIn + ' does not exist')


class PyscesTest:
    """PySCeS test suite: takes a test level as an argument"""

    __version__ = __version__

    def __init__(self, lvl=2, std2file=0):
        # copy models from server to local model store
        print('\nCopying pysces_test models if necessary ...')
        CopyTestModels()
        print('done.')
        self.basic_runner = unittest.TextTestRunner()
        self.BasicTest = unittest.defaultTestLoader.loadTestsFromTestCase(
            PyscesBasicTest
        )
        self.ExtendedTest = unittest.defaultTestLoader.loadTestsFromTestCase(
            PyscesExtendedTest
        )
        self.ExternalTest = unittest.defaultTestLoader.loadTestsFromTestCase(
            PyscesExternalTest
        )
        self.OptionalDependencyTest = unittest.defaultTestLoader.loadTestsFromTestCase(
            PyscesOptionalDependencyTest
        )
        if lvl > 4:
            lvl = 4

        class NullWriter:
            def __init__(self, dres=0):
                if dres:
                    self.Fout = open(
                        os.path.join(os.getcwd(), 'pysces_test_results.txt'), 'w'
                    )
                    self.fof = 1
                else:
                    self.fof = 0

            def write(self, s):
                if self.fof:
                    self.Fout.write(s)

            def close(self):
                if self.fof:
                    self.Fout.flush()
                    self.Fout.close()

        self.__dfi__ = NullWriter(std2file)
        tmpSTDOUT = sys.stdout
        tmpSTDERR = sys.stderr
        for x in range(lvl):
            print('\nLevel ' + str(x + 1) + ' tests')
            sys.stdout = self.__dfi__
            sys.stderr = self.__dfi__
            getattr(self, 'lvl_' + str(x + 1))()
            sys.stdout = tmpSTDOUT
            sys.stderr = tmpSTDERR
            sleep(1)
        self.__dfi__.close()

    def lvl_1(self):
        self.basic_runner.run(self.BasicTest)

    def lvl_2(self):
        self.basic_runner.run(self.ExtendedTest)

    def lvl_3(self):
        self.basic_runner.run(self.ExternalTest)

    def lvl_4(self):
        self.basic_runner.run(self.OptionalDependencyTest)


class PyscesBasicTest(unittest.TestCase):
    """Basic test class, tests low-level numerical algorithms"""

    __version__ = __version__
    model_dir = MODEL_DIR

    MathArrayTest = PyscesStoich.MathArrayFunc()

    def test_swaprow_d(self):
        arr = array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], 'd')
        karr = array([[7, 8, 9], [4, 5, 6], [1, 2, 3]], 'd')
        arr = self.MathArrayTest.SwapRowd(arr, 0, 2)
        assert_array_equal(arr, karr)

    def test_swapcol_d(self):
        arr = array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], 'd')
        karr = array([[3, 2, 1], [6, 5, 4], [9, 8, 7]], 'd')
        arr = self.MathArrayTest.SwapCold(arr, 0, 2)
        assert_array_equal(arr, karr)

    def test_swaprow_z(self):
        arr = array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], 'D')
        karr = array([[7, 8, 9], [4, 5, 6], [1, 2, 3]], 'D')
        arr = self.MathArrayTest.SwapRowz(arr, 0, 2)
        assert_array_equal(arr, karr)

    def test_swapcol_z(self):
        arr = array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], 'D')
        karr = array([[3, 2, 1], [6, 5, 4], [9, 8, 7]], 'D')
        arr = self.MathArrayTest.SwapColz(arr, 0, 2)
        assert_array_equal(arr, karr)

    def test_matrixfloatfix(self):
        arr = array([[1.0e-15, 2, 3], [4, -0.0, 6], [7, 8, 1.0e-15]])
        karr = array([[0.0, 2, 3], [4, 0.0, 6], [7, 8, 0.0]])
        self.MathArrayTest.MatrixFloatFix(arr, val=1.0e-13)
        assert_array_equal(arr, karr)

    def test_FSOLVE_linear(self):
        ds = zeros((3), 'd')
        s = zeros((3), 'd')

        def linear1_ode(s):
            ds[0] = 100.0 - 6.0 * s[0] + s[1]
            ds[1] = 5.0 * s[0] - 4.0 * s[1] + s[2]
            ds[2] = 3.0 * s[1] - 3.0 * s[2] + 1.0
            return ds

        res = fsolve(linear1_ode, [1.0, 1.0, 1.0])
        known = [23.1025641, 38.61538462, 38.94871795]
        assert_array_almost_equal(res, known)

    def test_ODEINT_linear(self):
        ds = zeros((3), 'd')
        s = zeros((3), 'd')

        def linear1_ode(s, t):
            ds[0] = 100.0 - 6.0 * s[0] + s[1]
            ds[1] = 5.0 * s[0] - 4.0 * s[1] + s[2]
            ds[2] = 3.0 * s[1] - 3.0 * s[2] + 1.0
            return ds

        sim_res = odeint(
            linear1_ode, [1.0, 1.0, 1.0], array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        )
        known = array(
            [
                [1.0, 1.0, 1.0],
                [20.70698909, 27.55450322, 20.9819434],
                [22.4418045, 35.47695146, 33.60751644],
                [22.91227529, 37.71013683, 37.40388317],
                [23.04761842, 38.35410015, 38.50268374],
                [23.08668173, 38.53997688, 38.81993704],
                [23.09809306, 38.59337738, 38.91160251],
                [23.10124276, 38.60909619, 38.93798616],
                [23.10217397, 38.61358525, 38.94561234],
                [23.10245088, 38.61486396, 38.94781836],
                [23.10253812, 38.61526501, 38.94851136],
            ]
        )

        assert_array_almost_equal(sim_res, known, 3)

    def test_ODEINT_moiety(self):
        """This set of ODE's is numerically unstable and might fail without error"""
        ds = zeros((3), 'd')
        s = zeros((3), 'd')

        def moibranch_ode(s, t):
            s[2] = 1.0 - s[1]
            ds[0] = s[1] / (2.0 + 2.0 * s[1]) - 2.0 * s[0] / (1.0 + s[0])
            ds[1] = 10.0 * s[2] / (2.0 + 2.0 * s[2]) - s[1] / (2.0 + 2.0 * s[1])
            ds[2] = s[1] / (2.0 + 2.0 * s[1]) - 10.0 * s[2] / (2.0 + 2.0 * s[2])
            return ds

        sim_res = odeint(
            moibranch_ode, [0.7, 1.0, 0.3], array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        )
        known = array(
            [
                [0.7, 1.0, 0.3],
                [0.43521812, 0.94587168, 0.05412832],
                [0.21547096, 0.94879365, 0.05120635],
                [0.15590709, 0.94882109, 0.05117891],
                [0.14233313, 0.9488214, 0.0511786],
                [0.13938884, 0.94882111, 0.05117889],
                [0.13875741, 0.94882111, 0.05117889],
                [0.1386223, 0.94882133, 0.05117867],
                [0.1385936, 0.94882134, 0.05117866],
                [0.13858742, 0.94882133, 0.05117867],
                [0.13858601, 0.94882133, 0.05117867],
            ]
        )

        assert_array_almost_equal(sim_res[-5:, :], known[-5:, :], 3)


class PyscesExtendedTest(unittest.TestCase):
    """Extended test class, tests modelling related methods"""

    __version__ = __version__
    model_dir = MODEL_DIR

    def test_statemetab_linear1(self):
        lin = PSCMODEL('pysces_test_linear1.psc', self.model_dir)
        lin.State()
        linmet = array([23.1025641, 38.61538462, 38.94871795], 'd')
        assert_array_almost_equal(lin.state_species, linmet)

    def test_statemetab_branch1(self):
        bra = PSCMODEL('pysces_test_branch1.psc', self.model_dir)
        bra.State()
        bramet = array([4.8583996, 1.88547254, 1.49124431, 1.49124431], 'd')
        assert_array_almost_equal(bra.state_species, bramet)

    def test_statemetab_moiety1(self):
        moi = PSCMODEL('pysces_test_moiety1.psc', self.model_dir)
        moi.State()
        moimet = array(
            [
                3.6886875,
                16.25569882,
                7.3113125,
                4.39229787,
                41.02504596,
                2.60770213,
                0.42718994,
                2.57281006,
                2.44791155,
                17.0012171,
            ],
            'd',
        )
        assert_array_almost_equal(moi.state_species, moimet)

    def test_stateflux_linear1(self):
        lin = PSCMODEL('pysces_test_linear1.psc', self.model_dir)
        lin.State()
        linflux = array([76.8974359, 76.8974359, 76.8974359, 76.8974359], 'd')
        assert_array_almost_equal(lin.state_flux, linflux)

    def test_stateflux_branch1(self):
        bra = PSCMODEL('pysces_test_branch1.psc', self.model_dir)
        bra.State()
        braflux = array(
            [2.42139889, 2.42139889, 1.21069945, 1.21069945, 1.21069945, 1.21069945],
            'd',
        )
        assert_array_almost_equal(bra.state_flux, braflux)

    def test_stateflux_moiety1(self):
        moi = PSCMODEL('pysces_test_moiety1.psc', self.model_dir)
        moi.State()
        moiflux = array(
            [
                250.01825652,
                250.01825652,
                250.01825652,
                250.01825652,
                250.01825652,
                250.01825652,
                250.01825652,
            ],
            'd',
        )
        assert_array_almost_equal(moi.state_flux, moiflux)

    def test_elas_linear1(self):
        lin = PSCMODEL('pysces_test_linear1.psc', self.model_dir)
        lin.State()
        lin.EvalEvar()
        line = [
            -0.30043348,
            0.0,
            0.0,
            1.50216739,
            -0.50216739,
            0.0,
            0.0,
            1.50650217,
            -0.50650217,
            0.0,
            0.0,
            1.01300433,
        ]
        line_new = []
        for x in range(lin.elas_var.shape[0]):
            for y in range(lin.elas_var.shape[1]):
                line_new.append(round(lin.elas_var[x, y], 2))
        for x in range(len(line)):
            line[x] = round(line[x], 2)

        assert_array_almost_equal(line, line_new)

    def test_elas_branch1(self):
        bra = PSCMODEL('pysces_test_branch1.psc', self.model_dir)
        bra.State()
        bra.EvalEvar()
        brae = [
            -0.66930781448548349,
            0.0,
            0.0,
            0.0,
            0.78845903183743249,
            -0.52920041809870422,
            0.0,
            0.0,
            0.0,
            0.95441603882382564,
            -0.60578219088834051,
            0.0,
            0.0,
            0.95441603882382564,
            0.0,
            -0.60578219088834051,
            0.0,
            0.0,
            0.9421058792599355,
            0.0,
            0.0,
            0.0,
            0.0,
            0.9421058792599355,
        ]
        brae_new = []
        for x in range(bra.elas_var.shape[0]):
            for y in range(bra.elas_var.shape[1]):
                brae_new.append(round(bra.elas_var[x, y], 2))
        for x in range(len(brae)):
            brae[x] = round(brae[x], 2)
        assert_array_almost_equal(brae, brae_new)

    def test_elas_moiety1(self):
        moi = PSCMODEL('pysces_test_moiety1.psc', self.model_dir)
        moi.State()
        moi.EvalEvar()
        moie = [
            1.4753672613878632,
            -0.47536726138786306,
            -0.47536726138786312,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.4278931518740323,
            0.0,
            1.4278931518740325,
            -0.42789315187403271,
            -0.42789315187403271,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0514524437327983,
            0.0,
            1.0514524437327983,
            -0.051452443732798468,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            -0.04300468618304179,
            0.0,
            1.0430046861830418,
            0.0,
            0.0,
            -0.043004686183041797,
            0.0,
            -0.073768363069393092,
            0.0,
            1.0737683630693931,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0737683630693931,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            -0.029048874655981143,
            1.0290488746559812,
            0.0,
            -0.029048874655981143,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0199985395848277,
        ]
        moie_new = []
        for x in range(moi.elas_var.shape[0]):
            for y in range(moi.elas_var.shape[1]):
                moie_new.append(round(moi.elas_var[x, y], 2))
        for x in range(len(moie)):
            moie[x] = round(moie[x], 2)
        assert_array_almost_equal(moie, moie_new)

    def test_cc_linear1(self):
        lin = PSCMODEL('pysces_test_linear1.psc', self.model_dir)
        lin.State()
        lin.EvalEvar()
        lin.EvalCC()
        lincc = [
            0.02564102564102570300,
            0.05128205128205126600,
            0.15384615384615383000,
            0.76923076923076916000,
            0.02564102564102571000,
            0.05128205128205128000,
            0.15384615384615385000,
            0.76923076923076938000,
            0.02564102564102571000,
            0.05128205128205128000,
            0.15384615384615385000,
            0.76923076923076938000,
            0.02564102564102570300,
            0.05128205128205126600,
            0.15384615384615383000,
            0.76923076923076916000,
            -0.30636428644396779000,
            -0.61272857288793536000,
            0.15318214322198387000,
            0.76591071610991934000,
            -0.08534676570192656400,
            -0.17069353140385327000,
            -0.51208059421155983000,
            0.76812089131733974000,
            -0.96185074526088343000,
            0.05062372343478333000,
            0.15187117030435002000,
            0.75935585152175000000,
        ]
        lincc_new = []
        for x in range(lin.cc_all.shape[0]):
            for y in range(lin.cc_all.shape[1]):
                lincc_new.append(round(lin.cc_all[x, y], 2))
        for x in range(len(lincc)):
            lincc[x] = round(lincc[x], 2)
        assert_array_almost_equal(lincc, lincc_new)

    def test_cc_branch1(self):
        bra = PSCMODEL('pysces_test_branch1.psc', self.model_dir)
        bra.State()
        bra.EvalEvar()
        bra.EvalCC()
        bracc = [
            0.25338970963029361000,
            -0.13797075277724413000,
            -0.21457061516904127000,
            0.32372624345386269000,
            0.39406892242342095000,
            0.38135649243870806000,
            -0.13797075277724413000,
            0.25338970963029356000,
            0.39406892242342095000,
            0.32372624345386258000,
            -0.21457061516904130000,
            0.38135649243870812000,
            -0.13797075277724413000,
            0.25338970963029356000,
            0.39406892242342095000,
            0.32372624345386258000,
            -0.21457061516904130000,
            0.38135649243870812000,
            0.05770947842652475500,
            0.05770947842652472700,
            0.08974915362718985400,
            0.32372624345386269000,
            0.08974915362718984000,
            0.38135649243870817000,
            0.25338970963029361000,
            -0.13797075277724413000,
            -0.21457061516904127000,
            0.32372624345386269000,
            0.39406892242342095000,
            0.38135649243870806000,
            0.05770947842652474100,
            0.05770947842652471300,
            0.08974915362718984000,
            0.32372624345386264000,
            0.08974915362718982600,
            0.38135649243870806000,
            -0.23751396180748979000,
            -0.23751396180748979000,
            -0.36937913195668803000,
            0.55728841856739109000,
            -0.36937913195668803000,
            0.65649776896096446000,
            -0.08622262758262820600,
            -0.08622262758262820600,
            -0.13409249329650016000,
            -0.48367318660804903000,
            -0.13409249329650016000,
            0.92430342836630586000,
            -0.79249085140642661000,
            -0.14644930661681688000,
            -0.22775636995026044000,
            0.34361981023636490000,
            0.41828517483934929000,
            0.40479154289778968000,
            -0.14644930661681685000,
            -0.79249085140642661000,
            0.41828517483934929000,
            0.34361981023636484000,
            -0.22775636995026050000,
            0.40479154289778979000,
        ]
        bracc_new = []
        for x in range(bra.cc_all.shape[0]):
            for y in range(bra.cc_all.shape[1]):
                bracc_new.append(round(bra.cc_all[x, y], 2))
        for x in range(len(bracc)):
            bracc[x] = round(bracc[x], 2)
            # os.sys.stderr.write(str(bracc[x]) + ' --> ' + str(bracc_new[x])+'\n')
        assert_array_almost_equal(bracc, bracc_new)

    def test_cc_moiety1(self):
        moi = PSCMODEL('pysces_test_moiety1.psc', self.model_dir)
        moi.State()
        moi.EvalEvar()
        moi.EvalCC()
        moicc = [
            0.01114694087227875000,
            0.07381787576415140000,
            0.18139104475836643000,
            0.04685617079784056700,
            0.25208717198140801000,
            0.04329618737852313600,
            0.39140460844743269000,
            0.01114694087227862500,
            0.07381787576415056700,
            0.18139104475836440000,
            0.04685617079784004000,
            0.25208717198140518000,
            0.04329618737852265100,
            0.39140460844742830000,
            0.01114694087227874000,
            0.07381787576415133100,
            0.18139104475836623000,
            0.04685617079784051800,
            0.25208717198140773000,
            0.04329618737852309500,
            0.39140460844743230000,
            0.01114694087227875000,
            0.07381787576415140000,
            0.18139104475836643000,
            0.04685617079784056700,
            0.25208717198140801000,
            0.04329618737852313600,
            0.39140460844743269000,
            0.01114694087227875000,
            0.07381787576415140000,
            0.18139104475836643000,
            0.04685617079784056700,
            0.25208717198140801000,
            0.04329618737852313600,
            0.39140460844743269000,
            0.01114694087227876100,
            0.07381787576415146900,
            0.18139104475836659000,
            0.04685617079784060900,
            0.25208717198140823000,
            0.04329618737852317800,
            0.39140460844743302000,
            0.01114694087227875200,
            0.07381787576415141400,
            0.18139104475836645000,
            0.04685617079784057400,
            0.25208717198140806000,
            0.04329618737852314300,
            0.39140460844743280000,
            -0.96946517151898726000,
            0.07237057005414701500,
            0.17783461222620814000,
            0.04593748812318156800,
            0.24714464011293155000,
            0.04244730330314599300,
            0.38373055769937375000,
            -0.00493043463469603270,
            -0.03265059135931840800,
            -0.08023158100034956400,
            0.14913606725290945000,
            0.02650472034900546900,
            0.11529508789995771000,
            -0.17312326850750914000,
            -0.00258941341724495350,
            -0.01714775382110165000,
            -0.04213679882643606200,
            0.25950964810000904000,
            -0.07785637143685617000,
            -0.02885675558547331700,
            -0.09092255501289707400,
            -0.00651179467430067560,
            -0.04312275948862360300,
            -0.10596460973080281000,
            -0.02003169777867189200,
            0.40783108111031047000,
            -0.00355036168876400080,
            -0.22864985774914801000,
            -0.07520202175595115700,
            -0.49800690277268156000,
            1.11329060303376930000,
            0.28758054022385904000,
            1.54718927875470390000,
            0.26573108181778565000,
            -2.64058257930148920000,
            0.08511194652599511600,
            -0.37976718223865702000,
            -0.93319355560031703000,
            -0.24105862936563618000,
            -1.29690043219020780000,
            -0.22274375836758617000,
            2.98855161123640300000,
            -0.01413200623441426100,
            0.06305662607990611400,
            0.15494766227242857000,
            0.04002542759392844300,
            0.21533763168639128000,
            0.03698442240380646300,
            -0.49621976380204558000,
            0.01332315621836589900,
            0.08822932693215754200,
            0.21680398717628285000,
            -0.25121007495788511000,
            0.32322678367369212000,
            -0.85819164176559659000,
            0.46781846272298522000,
            0.01096817827331699600,
            0.07263406439629339900,
            0.17848209108569948000,
            0.03374050370795465800,
            -0.68693259335574453000,
            0.00598007183654026410,
            0.38512768405594056000,
            0.00513245176028354240,
            0.03398840011328187900,
            0.08351894906745090100,
            -0.51437161070192172000,
            0.15431837495286327000,
            0.05719670138973560700,
            0.18021673341830685000,
        ]
        moicc_new = []
        for x in range(moi.cc_all.shape[0]):
            for y in range(moi.cc_all.shape[1]):
                moicc_new.append(round(moi.cc_all[x, y], 2))
        for x in range(len(moicc)):
            moicc[x] = round(moicc[x], 2)
            # os.sys.stderr.write(str(moicc[x]) + ' --> ' + str(moicc_new[x])+'\n')
        assert_array_almost_equal(moicc, moicc_new)

    def test_scan_linear1(self):
        lin = PSCMODEL('pysces_test_linear1.psc', self.model_dir)
        lin.scan_in = 'x0'
        lin.scan_out = ['J_R1', 's0_ss', 's1_ss', 's2_ss']
        scan_range = linspace(2, 20, 10)
        lin.Scan1(scan_range)
        linscan = array(
            [
                [2.0, 15.35897436, 4.64102564, 7.84615385, 8.17948718],
                [4.0, 30.74358974, 9.25641026, 15.53846154, 15.87179487],
                [6.0, 46.12820513, 13.87179487, 23.23076923, 23.56410256],
                [8.0, 61.51282051, 18.48717949, 30.92307692, 31.25641026],
                [10.0, 76.8974359, 23.1025641, 38.61538462, 38.94871795],
                [12.0, 92.28205128, 27.71794872, 46.30769231, 46.64102564],
                [14.0, 107.66666667, 32.33333333, 54.0, 54.33333333],
                [16.0, 123.05128205, 36.94871795, 61.69230769, 62.02564103],
                [18.0, 138.43589744, 41.56410256, 69.38461538, 69.71794872],
                [20.0, 153.82051282, 46.17948718, 77.07692308, 77.41025641],
            ]
        )
        assert_allclose(linscan, lin.scan_res)

    def test_scan_branch1(self):
        bra = PSCMODEL('pysces_test_branch1.psc', self.model_dir)
        bra.scan_in = 'Vf5'
        bra.scan_out = ['J_R1', 'J_R3', 's1_ss', 's2_ss', 's4_ss']
        scan_range = linspace(5, 50, 10)
        bra.Scan1(scan_range)
        brascan = array(
            [
                [5.0, 2.31904671, 0.98414893, 5.17769673, 2.23691074, 2.57516754],
                [10.0, 2.42139889, 1.21069945, 4.8583996, 1.88547254, 1.49124431],
                [15.0, 2.47522288, 1.32840199, 4.70029299, 1.71892539, 1.08136736],
                [20.0, 2.50940431, 1.4026989, 4.60314727, 1.61904914, 0.86179112],
                [25.0, 2.53330794, 1.45446464, 4.53665539, 1.55176582, 0.72390142],
                [30.0, 2.55106248, 1.49281822, 4.48801604, 1.50310207, 0.62890054],
                [35.0, 2.56481301, 1.52246889, 4.45077512, 1.4661592, 0.55932108],
                [40.0, 2.57579746, 1.5461228, 4.42129073, 1.43710561, 0.50609274],
                [45.0, 2.58478513, 1.56545628, 4.39733911, 1.41363069, 0.46402152],
                [50.0, 2.5922813, 1.58156758, 4.37748021, 1.39425319, 0.42991222],
            ]
        )
        assert_allclose(brascan, bra.scan_res)

    def test_scan_moiety1(self):
        moi = PSCMODEL('pysces_test_moiety1.psc', self.model_dir)
        moi.scan_in = 's0_init'
        moi.scan_out = ['J_R1', 's0_ss', 's1_ss', 's4_ss', 's5_ss']
        scan_range = linspace(1, 10, 10)
        moi.Scan1(scan_range)
        moiscan = array(
            [
                [1.0, 250.01825652, 0.42718994, 2.57281006, 3.6886875, 7.3113125],
                [2.0, 273.03336666, 1.07112697, 2.92887303, 3.84554096, 7.15445904],
                [3.0, 281.12962035, 1.83808367, 3.16191633, 3.9009623, 7.0990377],
                [4.0, 284.88428449, 2.6408052, 3.3591948, 3.9267074, 7.0732926],
                [5.0, 287.00414515, 3.45698856, 3.54301144, 3.94125508, 7.05874492],
                [6.0, 288.35651302, 4.27949058, 3.72050942, 3.95054037, 7.04945963],
                [7.0, 289.29162429, 5.10542741, 3.89457259, 3.95696289, 7.04303711],
                [8.0, 289.97583488, 5.93342945, 4.06657055, 3.96166326, 7.03833674],
                [9.0, 290.49779112, 6.76276696, 4.23723304, 3.9652496, 7.0347504],
                [10.0, 290.90891151, 7.59301668, 4.40698332, 3.96807476, 7.03192524],
            ]
        )
        assert_allclose(moiscan, moi.scan_res)


class PyscesExternalTest(unittest.TestCase):
    """Extended test class, tests external/add-in numerical algorithms"""

    __version__ = __version__
    model_dir = MODEL_DIR

    def test_PITCON1(self):
        assert pitcon_switch, (
            'PITCON is not installed! Continuation analysis not available.'
        )

        print(
            """
        C  PCPRB1.FOR  The Freudenstein-Roth function.
        C  Limit points in the first variable occur at:
        C    (14.28309, -1.741377,  0.2585779)
        C    (61.66936,  1.983801, -0.6638797)
        """
        )

        def fx(X):
            FX[0] = (
                X[0] - ((X[1] - 5.0) * X[1] + 2.0) * X[1] - 13.0 + 34.0 * (X[2] - 1.0)
            )
            FX[1] = (
                X[0] + ((X[1] + 1.0) * X[1] - 14.0) * X[1] - 29.0 + 10.0 * (X[2] - 1.0)
            )
            return FX

        def fjac(s):
            return

        parm = 3
        iwork = zeros((30 + parm), 'i')
        rwork = zeros((30 + (6 * parm) * parm), 'd')
        ipar = zeros((parm), 'i')
        fpar = zeros((parm), 'd')
        xr = zeros((parm), 'd')  # output array
        FX = zeros((3), 'd')  # function array 'v'

        iwork[0] = 0
        iwork[1] = 2
        iwork[2] = 0
        iwork[3] = 0
        iwork[4] = 0
        iwork[5] = 1
        iwork[6] = 0
        iwork[7] = 6
        iwork[8] = 2
        iwork[9] = 0
        iwork[10] = 0
        iwork[11] = 0
        iwork[12] = 30
        iwork[13] = len(iwork)
        iwork[14] = 30 + 4 * parm
        iwork[15] = len(rwork)
        iwork[16] = 20

        rwork[0] = 0.00001
        rwork[1] = 0.00001
        rwork[2] = 0.01
        rwork[3] = 20.0
        rwork[4] = 0.3
        rwork[5] = 1.0
        rwork[6] = 1.0
        rwork[7] = 0.0
        rwork[19] = 3.0

        xr[0] = 15.0
        xr[1] = -2.0
        xr[2] = 0.0

        output = []
        limits = []
        iterations = 10

        for run in range(iterations):
            ierror, iwork, rwork, xrout = pitcon.pitcon1(
                fjac, fpar, fx, ipar, iwork, rwork, xr
            )
            if iwork[0] == 4:
                print('\nLimit point in run = ' + repr(run + 1))
                print(xrout)
                limits.append(xrout)
            output.append(xrout)

        output = array(output)
        limit0 = [14.28309125, -1.74137688, 0.25857788]
        limit1 = [61.66936258, 1.98380112, -0.66387974]
        out6 = [4.48781323e001, 4.87659925e-001, 5.95322940e-002]

        assert_array_almost_equal(limits[0], limit0, 4)
        assert_array_almost_equal(limits[1], limit1, 4)
        assert_array_almost_equal(output[6], out6, 4)

    def test_NLEQ2_statemetab_linear1(self):
        lin = PSCMODEL('pysces_test_linear1.psc', self.model_dir)
        lin.mode_solver = 'NLEQ2'
        lin.State()
        linmet = array([23.1025641, 38.61538462, 38.94871795], 'd')
        assert_array_almost_equal(lin.state_species, linmet)

    def test_NLEQ2_statemetab_branch1(self):
        bra = PSCMODEL('pysces_test_branch1.psc', self.model_dir)
        bra.mode_solver = 'NLEQ2'
        bra.State()
        bramet = array([4.8583996, 1.88547254, 1.49124431, 1.49124431], 'd')
        assert_array_almost_equal(bra.state_species, bramet)

    def test_NLEQ2_statemetab_moiety1(self):
        moi = PSCMODEL('pysces_test_moiety1.psc', self.model_dir)
        moi.mode_solver = 'NLEQ2'
        moi.State()
        moimet = array(
            [
                3.6886875,
                16.25569882,
                7.3113125,
                4.39229787,
                41.02504596,
                2.60770213,
                0.42718994,
                2.57281006,
                2.44791155,
                17.0012171,
            ],
            'd',
        )
        assert_array_almost_equal(moi.state_species, moimet)

    def test_NLEQ2_stateflux_linear1(self):
        lin = PSCMODEL('pysces_test_linear1.psc', self.model_dir)
        lin.mode_solver = 'NLEQ2'
        lin.State()
        linflux = array([76.8974359, 76.8974359, 76.8974359, 76.8974359], 'd')
        assert_array_almost_equal(lin.state_flux, linflux)

    def test_NLEQ2_stateflux_branch1(self):
        bra = PSCMODEL('pysces_test_branch1.psc', self.model_dir)
        bra.mode_solver = 'NLEQ2'
        bra.State()
        braflux = array(
            [2.42139889, 2.42139889, 1.21069945, 1.21069945, 1.21069945, 1.21069945],
            'd',
        )
        assert_array_almost_equal(bra.state_flux, braflux)

    def test_NLEQ2_stateflux_moiety1(self):
        moi = PSCMODEL('pysces_test_moiety1.psc', self.model_dir)
        moi.mode_solver = 'NLEQ2'
        moi.State()
        moiflux = array(
            [
                250.01825652,
                250.01825652,
                250.01825652,
                250.01825652,
                250.01825652,
                250.01825652,
                250.01825652,
            ],
            'd',
        )
        assert_array_almost_equal(moi.state_flux, moiflux)

    # def test_PITCON2(self): # unreliable
    # mod = model('pysces_test_pitcon.psc', self.model_dir)
    # mod.doLoad()
    # mod.V2 = mod.V3 = 100.0
    # mod.A05 = 10.0
    # mod.pitcon_iter = 5
    # mod.pitcon_par_space = logspace(0,1,10)
    # res = mod.PITCON('P')
    ## pysces.plt.plot2D(res, 0, style='.'); pysces.plt.logxy()

    # os.sys.stderr.write('\n\"Iteration terminates after NITMAX and Newton method fails to converge\" warnings can be ignored\n')
    # benchmark = array([
    # [8.39949877e-01,1.61211501e+01,4.84068424e+00,7.12011640e+01,7.12011640e+01,7.12011640e+01],
    # [5.40984421e-01,1.57720807e+01,4.31351073e+00,7.27545617e+01,7.27545617e+01,7.27545617e+01],
    # [-5.14362451e-02,1.45210385e+01,3.06361620e+00,7.64860963e+01,7.64860963e+01,7.64860963e+01],
    # [-4.66310864e-01,1.31962735e+01,1.95980785e+00,8.04668684e+01,8.04668684e+01,8.04668683e+01],
    # [1.13190191e+00,1.61011691e+01,5.26326909e+00,6.96411063e+01,6.96411063e+01,6.96411063e+01],
    # [8.30873009e-01,1.61149804e+01,4.82588194e+00,7.12478479e+01,7.12478479e+01,7.12478479e+01],
    # [4.39606946e-01,1.56001789e+01,4.11842961e+00,7.33077027e+01,7.33077027e+01,7.33077027e+01],
    # [-3.07339668e-01,1.37808508e+01,2.41799659e+00,7.87206889e+01,7.87206889e+01,7.87206889e+01],
    # [1.64608357e+00,1.37635551e+01,5.46388538e+00,6.53427535e+01,6.53427535e+01,6.53427535e+01],
    # [1.57039736e+00,1.45248830e+01,5.53382469e+00,6.63454792e+01,6.63454792e+01,6.63454792e+01],
    # [1.36033858e+00,1.56423070e+01,5.48938893e+00,6.81979734e+01,6.81979734e+01,6.81979734e+01],
    # [3.90033764e-01,1.55086998e+01,4.02030887e+00,7.35869382e+01,7.35869382e+01,7.35869382e+01],
    # [1.90300355e+00,4.61218050e-01,8.11483528e-01,1.67232557e+01,1.67232557e+01,1.67232557e+01],
    # [1.65604156e+00,1.14527548e+00,1.32711044e+00,2.91604808e+01,2.91604810e+01,2.91604810e+01],
    # [1.60792245e+00,3.36846985e+00,2.47903648e+00,4.55723011e+01,4.55723011e+01,4.55723011e+01],
    # [1.76260920e+00,1.02529273e+01,4.77285351e+00,6.09994737e+01,6.09994737e+01,6.09994737e+01],
    # [2.48410121e+00,1.76577346e-01,5.36816396e-01,7.17264824e+00,7.17264824e+00,7.17264824e+00],
    # [2.04216203e+00,3.36094901e-01,6.94802141e-01,1.31279273e+01,1.31279273e+01,1.31279273e+01],
    # [1.68318189e+00,9.87364532e-01,1.22120154e+00,2.69666997e+01,2.69666997e+01,2.69666997e+01],
    # [1.60195518e+00,3.09604637e+00,2.35831199e+00,4.43144773e+01,4.43144773e+01,4.43144773e+01],
    # [3.29451248e+00,1.10989590e-01,5.07996700e-01,3.71775348e+00,3.71775351e+00,3.71775348e+00],
    # [2.53835638e+00,1.67389490e-01,5.28506741e-01,6.75388116e+00,6.75388117e+00,6.75388117e+00],
    # [2.13066473e+00,2.84747678e-01,6.44229065e-01,1.14218470e+01,1.14218470e+01,1.14218470e+01],
    # [1.85930898e+00,5.18944892e-01,8.62231652e-01,1.81726016e+01,1.81726016e+01,1.81726016e+01],
    # [4.93995542e+00,1.06865480e-01,6.50716299e-01,2.37791780e+00,2.37791781e+00,2.37791780e+00],
    # [5.51998987e+00,1.12930767e-01,7.15671948e-01,2.26203167e+00,2.26203167e+00,2.26203167e+00],
    # [7.20197168e+00,1.35245848e-01,9.14818799e-01,2.13476038e+00,2.13476038e+00,2.13476038e+00],
    # [1.21464994e+01,2.10315787e-01,1.52485159e+00,2.11433598e+00,2.11433598e+00,2.11433598e+00],
    # [6.29277796e+00,1.22633783e-01,8.05857247e-01,2.18036060e+00,2.18036060e+00,2.18036060e+00],
    # [6.71463987e+00,1.28373921e-01,8.56137906e-01,2.15469266e+00,2.15469266e+00,2.15469266e+00],
    # [7.96638113e+00,1.46372121e-01,1.00775931e+00,2.11668562e+00,2.11668562e+00,2.11668562e+00],
    # [1.16768618e+01,2.02994936e-01,1.46632082e+00,2.11150946e+00,2.11150946e+00,2.11150946e+00],
    # [8.04040339e+00,1.47466151e-01,1.01680363e+00,2.11553054e+00,2.11553054e+00,2.11553055e+00],
    # [8.46158253e+00,1.53734006e-01,1.06838323e+00,2.11040547e+00,2.11040547e+00,2.11040547e+00],
    # [9.21463978e+00,1.65086727e-01,1.16102019e+00,2.10586653e+00,2.10586653e+00,2.10586653e+00],
    # [1.07956247e+01,1.89324600e-01,1.35672050e+00,2.10728985e+00,2.10728985e+00,2.10728985e+00],
    # [1.02976931e+01,1.81645169e-01,1.29494069e+00,2.10576586e+00,2.10576586e+00,2.10576586e+00],
    # [1.07186771e+01,1.88135560e-01,1.34716589e+00,2.10700665e+00,2.10700665e+00,2.10700665e+00],
    # [1.13140312e+01,1.97355266e-01,1.42115741e+00,2.10957637e+00,2.10957637e+00,2.10957637e+00],
    # [1.21560089e+01,2.10464256e-01,1.52603757e+00,2.11439651e+00,2.11439651e+00,2.11439651e+00]])
    # assert_array_almost_equal(benchmark[:5,:],res[:5,:],3)
    # assert_array_almost_equal(benchmark[-5:,:],res[-5:,:],1)


class PyscesOptionalDependencyTest(unittest.TestCase):
    """Test functions that require optional dependencies: assimulo and libsbml"""

    __version__ = __version__
    model_dir = MODEL_DIR

    def test_sbml_log_ln(self):
        refdata = array(
            [
                [0.0, 0.0, 0.0, 0.0],
                [40.04004004, 549.52160273, 103.66716917, 405.29619564],
                [80.08008008, 73.89946487, 950.21148336, 394.62757153],
                [120.12012012, 1497.90258854, 308.38181296, 71.12217444],
                [160.16016016, 272.97335258, 73.98220441, 1853.80464473],
                [200.2002002, 73.89557669, 2059.78641999, 265.25805896],
                [240.24024024, 2171.23593707, 273.7411225, 70.26188816],
                [280.28028028, 292.95477084, 64.99312756, 2211.92916261],
                [320.32032032, 59.76067809, 2192.09743141, 320.24799747],
                [360.36036036, 2120.23371348, 354.37388765, 55.45540257],
                [400.4004004, 394.85037388, 52.40634153, 2006.27397873],
                [440.44044044, 50.6559264, 1861.19005023, 441.64759609],
                [480.48048048, 1695.84305088, 495.02140025, 50.1333647],
                [520.52052052, 555.41632185, 50.74137638, 1520.10307512],
                [560.56056056, 52.39360167, 1342.37127729, 623.401937],
                [600.6006006, 1169.38689863, 699.62285956, 55.02816614],
                [640.64064064, 784.74938261, 58.61092737, 1006.20448896],
                [680.68068068, 63.13480612, 856.2800681, 879.41826124],
                [720.72072072, 721.63761082, 984.15368027, 68.61801675],
                [760.76076076, 1099.25830557, 75.10234636, 603.09440286],
                [800.8008008, 82.65190152, 500.51982087, 1224.66467103],
                [840.84084084, 413.09911447, 1359.73964232, 91.35243743],
                [880.88088088, 1503.04174969, 101.31126888, 339.57665537],
                [920.92092092, 112.65772564, 278.46096926, 1652.04599086],
                [960.96096096, 228.18341849, 1802.87632682, 125.54410944],
            ]
        )
        # input
        model_file = 'pysces_test_BIOMD012.xml'
        # convert to psc
        interface.convertSBML2PSC(model_file, sbmldir=self.model_dir, pscdir=self.model_dir)
        mod = PSCMODEL(model_file + '.psc', dir=self.model_dir)
        # explicitly call LSODA as CVODE gets invoked as soon as there are assignment rules
        mod.mode_integrator = 'LSODA'
        mod.doSim(1000, 1000)
        lsoda_data = copy(mod.data_sim.getSimData('PX', 'PY', 'PZ'))[::40]
        mod.reLoad()
        mod.mode_integrator = 'CVODE'
        mod.__settings__['cvode_track_assignment_rules'] = False
        mod.doSim(1000, 1000)
        cvode_data = copy(mod.data_sim.getSimData('PX', 'PY', 'PZ'))[::40]
        assert allclose(lsoda_data, refdata), (
            "LSODA data doesn't match reference!"
        )
        assert allclose(cvode_data, refdata), (
            "CVODE data doesn't match reference!"
        )

    def test_sbml_rules_functions(self):
        refdata = array(
            [
                [0.0, 0.0, 0.0],
                [400.0, 0.0, 0.0],
                [800.0, 0.0, 0.0],
                [1200.0, 0.0, 0.0],
                [1600.0, 0.0, 0.0],
                [2000.0, 0.0, 0.0],
                [2400.0, 0.0, 0.0],
                [2800.0, 0.0, 0.0],
                [3200.0, 0.0, 0.0],
                [3600.0, 0.0, 0.0],
                [3996.0, 36.28124582, 518.81003737],
                [4396.0, 67.88089104, 523.23386418],
                [4796.0, 94.80833245, 396.36492755],
                [5196.0, 117.7543844, 312.57310603],
                [5596.0, 137.30772005, 265.67078904],
                [5996.0, 153.96997358, 236.6711118],
                [6396.0, 168.16860943, 216.96627888],
                [6796.0, 180.26788878, 202.77162667],
                [7196.0, 190.57821453, 192.14973106],
                [7596.0, 199.36409458, 183.98460585],
                [7996.0, 206.85092769, 177.58121669],
                [8396.0, 213.23078603, 172.48178027],
                [8796.0, 218.66734269, 168.37158625],
                [9196.0, 223.30007067, 165.02678928],
                [9596.0, 227.24782103, 162.28369752],
                [9996.0, 230.61187199, 160.01982735],
                [10396.0, 233.47852712, 158.14175147],
                [10796.0, 235.92132949, 156.57703989],
                [11196.0, 238.00294835, 155.26876267],
                [11596.0, 239.77678693, 154.1716452],
                [11996.0, 241.28835239, 153.24931975],
                [12396.0, 242.57642354, 152.47232202],
                [12796.0, 243.67404537, 151.81660551],
                [13196.0, 244.60937701, 151.26242178],
                [13596.0, 245.40641406, 150.79346492],
                [13996.0, 246.08560424, 150.39620951],
                [14396.0, 246.66437194, 150.05939229],
            ]
        )
        # input
        model_file = 'pysces_test_BIOMD860.xml'
        # convert to psc
        interface.convertSBML2PSC(model_file, sbmldir=self.model_dir, pscdir=self.model_dir)
        mod = PSCMODEL(model_file + '.psc', dir=self.model_dir)
        mod.__settings__['cvode_track_assignment_rules'] = False
        mod.doSim(4 * 3600, 3601)
        sim_data = copy(mod.data_sim.getSpecies())[::100]
        assert (allclose(sim_data, refdata)), "Data doesn't match reference!"

    def test_sbml_time_events(self):
        refdata = array(
            [
                [0.00000000e00, 0.00000000e00, 0.00000000e00, 0.00000000e00],
                [5.00000000e-01, 3.27535396e-03, 1.07158876e-04, 1.44643354e-06],
                [1.00000000e00, 6.51815647e-03, 2.14760296e-04, 2.90857654e-06],
                [1.50000000e00, 9.73085184e-03, 3.20599229e-04, 4.31859107e-06],
                [2.00000000e00, 1.29129602e-02, 4.25418590e-04, 5.71454558e-06],
                [2.50000000e00, 1.60647608e-02, 5.29239440e-04, 7.09719935e-06],
                [3.00000000e00, 1.91865422e-02, 6.32071471e-04, 8.46668899e-06],
                [3.50000000e00, 2.22785903e-02, 7.33924082e-04, 9.82313881e-06],
                [4.00000000e00, 2.53411882e-02, 8.34806603e-04, 1.11666732e-05],
                [4.50000000e00, 2.83746162e-02, 9.34728272e-04, 1.24974152e-05],
                [5.00000000e00, 3.13791525e-02, 1.03369825e-03, 1.38154865e-05],
                [5.50000000e00, 6.31897720e-03, 1.92799645e-02, 8.14828305e-02],
                [6.00000000e00, 7.66221177e-03, 5.53368763e-02, 2.76566872e-01],
                [6.50000000e00, 1.56751633e-02, 1.04563567e-01, 4.37783740e-01],
                [7.00000000e00, 2.37432451e-02, 1.44027928e-01, 5.39753117e-01],
                [7.50000000e00, 2.96549914e-02, 1.70232842e-01, 5.99674986e-01],
                [8.00000000e00, 3.34433158e-02, 1.86169817e-01, 6.33384265e-01],
                [8.50000000e00, 3.56877258e-02, 1.95234348e-01, 6.51034455e-01],
                [9.00000000e00, 3.69174144e-02, 1.99926808e-01, 6.58832833e-01],
                [9.50000000e00, 3.75056082e-02, 2.01891928e-01, 6.60580293e-01],
                [1.00000000e01, 3.76948097e-02, 2.02171128e-01, 6.58618153e-01],
                [1.05000000e01, 3.76379702e-02, 2.01414903e-01, 6.54393772e-01],
                [1.10000000e01, 3.74304993e-02, 2.00027549e-01, 6.48803209e-01],
                [1.15000000e01, 3.71316668e-02, 1.98259867e-01, 6.42401447e-01],
                [1.20000000e01, 3.67782297e-02, 1.96267297e-01, 6.35532146e-01],
                [1.25000000e01, 3.63929622e-02, 1.94146104e-01, 6.28407973e-01],
                [1.30000000e01, 3.59899641e-02, 1.91955855e-01, 6.21160377e-01],
                [1.35000000e01, 3.55779546e-02, 1.89733361e-01, 6.13870457e-01],
                [1.40000000e01, 3.51623167e-02, 1.87501318e-01, 6.06588084e-01],
                [1.45000000e01, 3.47463647e-02, 1.85273676e-01, 5.99343777e-01],
                [1.50000000e01, 3.43321303e-02, 1.83058960e-01, 5.92156057e-01],
                [1.55000000e01, 3.39208503e-02, 1.80862332e-01, 5.85036009e-01],
                [1.60000000e01, 3.35132690e-02, 1.78686874e-01, 5.77990117e-01],
                [1.65000000e01, 3.31098254e-02, 1.76534375e-01, 5.71022010e-01],
                [1.70000000e01, 3.27107692e-02, 1.74405825e-01, 5.64133554e-01],
                [1.75000000e01, 3.23162334e-02, 1.72301721e-01, 5.57325527e-01],
                [1.80000000e01, 3.19262786e-02, 1.70222256e-01, 5.50598035e-01],
                [1.85000000e01, 3.15409209e-02, 1.68167435e-01, 5.43950774e-01],
                [1.90000000e01, 3.11601491e-02, 1.66137146e-01, 5.37383189e-01],
                [1.95000000e01, 3.07839354e-02, 1.64131211e-01, 5.30894576e-01],
                [2.00000000e01, 3.04122416e-02, 1.62149407e-01, 5.24484139e-01],
                [2.05000000e01, 6.08004379e-01, 9.16295968e-02, 1.20706135e-02],
                [2.10000000e01, 6.84842857e-01, 2.36206540e-02, 6.10769675e-04],
                [2.15000000e01, 6.83064835e-01, 2.28426770e-02, 5.64968369e-04],
                [2.20000000e01, 6.80575906e-01, 2.27529115e-02, 5.61957724e-04],
                [2.25000000e01, 6.78100446e-01, 2.26692684e-02, 5.59251028e-04],
                [2.30000000e01, 6.75644017e-01, 2.25863097e-02, 5.56564861e-04],
                [2.35000000e01, 6.73206508e-01, 2.25039882e-02, 5.53896958e-04],
                [2.40000000e01, 6.70787771e-01, 2.24222944e-02, 5.51246953e-04],
                [2.45000000e01, 6.68387649e-01, 2.23412283e-02, 5.48614990e-04],
                [2.50000000e01, 6.66005995e-01, 2.22607829e-02, 5.46000840e-04],
                [2.55000000e01, 6.63642660e-01, 2.21809530e-02, 5.43404359e-04],
                [2.60000000e01, 6.61297496e-01, 2.21017337e-02, 5.40825435e-04],
                [2.65000000e01, 6.58970356e-01, 2.20231203e-02, 5.38263952e-04],
                [2.70000000e01, 6.56661096e-01, 2.19451080e-02, 5.35719785e-04],
                [2.75000000e01, 6.54369570e-01, 2.18676918e-02, 5.33192811e-04],
                [2.80000000e01, 6.52095638e-01, 2.17908670e-02, 5.30682909e-04],
                [2.85000000e01, 6.49839157e-01, 2.17146289e-02, 5.28189959e-04],
                [2.90000000e01, 6.47599987e-01, 2.16389727e-02, 5.25713843e-04],
                [2.95000000e01, 6.45377988e-01, 2.15638939e-02, 5.23254444e-04],
                [3.00000000e01, 6.43173023e-01, 2.14893876e-02, 5.20811643e-04],
            ]
        )
        # input
        model_file = 'pysces_test_BIOMD126.xml'
        # convert to psc
        interface.convertSBML2PSC(model_file, sbmldir=self.model_dir, pscdir=self.model_dir)
        mod = PSCMODEL(model_file+'.psc', dir=self.model_dir)
        mod.__settings__['cvode_return_event_timepoints'] = False
        mod.__settings__['cvode_track_assignment_rules'] = False
        mod.doSim(30, 61)
        sim_data = mod.data_sim.getSimData('IC3', 'IC2', 'IF')
        assert allclose(sim_data, refdata),  "Data doesn't match reference!"


if __name__ == '__main__':
    unittest.main()
