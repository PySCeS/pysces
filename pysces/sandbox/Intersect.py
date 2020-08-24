from __future__ import division, print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import os, sys
import pysces, scipy
import scipy.io
import gc
from time import sleep


class Intersect:
    dtype = scipy.double

    def __init__(self, mod, data_dir):
        self.mod = mod
        self.data_dir = data_dir
        self.STDOUT = sys.stdout

    def D2_setup_args(
        self,
        min1,
        max1,
        step1,
        min2,
        max2,
        step2,
        input,
        output,
        raw_data_file,
        surfaces,
        isect_data_file,
        isect_tol=1.0e-2,
        filename=None,
    ):
        self._min1 = min1
        self._max1 = max1
        self._step1 = step1
        self._min2 = min2
        self._max2 = max2
        self._step2 = step2
        self._input = input
        self._output = output
        self._raw_data_file = raw_data_file
        self._surfaces = surfaces
        self._isect_data_file = isect_data_file
        self._isect_tol = isect_tol
        self._filename = filename

    def D2_run(self):
        # try:
        self.D2_setup_scan(
            self._min1, self._max1, self._step1, self._min2, self._max2, self._step2
        )
        # self.load_raw_data('p4_raw.dat')
        self.D2_generate_raw_data(self._input, self._output, self._raw_data_file)
        self.D2_calculate_intersection_data(
            self._surfaces, self._isect_data_file, self._isect_tol
        )
        self.output_intersection_metrics(self._filename)
        # except Exception, EXC:
        #    print EXC
        #    sleep(5)
        #    sys.exit(1)

    def D2_setup_scan(self, min1, max1, step1, min2, max2, step2):
        self.P1min = min1
        self.P1max = max1
        self.P1steps = step1
        self.P2min = min2
        self.P2max = max2
        self.P2steps = step2
        self.P1range = scipy.logspace(scipy.log10(min1), scipy.log10(max1), step1)
        self.P2range = scipy.logspace(scipy.log10(min2), scipy.log10(max2), step2)

    def load_raw_data(self, raw_data_file):
        raw_file = open(os.path.join(self.data_dir, raw_data_file), 'r')
        self.data_raw = scipy.io.read_array(raw_file)
        raw_file.close()
        print('\nData file load complete.\n')

    def D2_generate_raw_data(self, input, output, raw_data_file):
        assert len(input) == 2
        assert len(output) > 0
        self.inputP = input
        self.outputP = output
        self.data_raw = scipy.zeros(
            (self.P1steps * self.P2steps, len(output) + 2), self.dtype
        )
        raw_file = open(os.path.join(self.data_dir, raw_data_file), 'w')
        iDx = 0
        flush = 0
        for P1 in range(len(self.P1range)):
            # not agaaaaaaaaaaaaaaaain
            gc.collect()
            del gc.garbage[:]
            data_dump = scipy.zeros((len(self.P2range), len(output) + 2), self.dtype)
            setattr(self.mod, input[0], self.P1range[P1])
            flush += 1
            # print ' '
            print(self.P1range[P1])
            # print input[1],
            for P2 in range(len(self.P2range)):
                # print "%2.2f" % self.P2range[P2],
                setattr(self.mod, input[1], self.P2range[P2])
                self.mod.doState()
                outArr = [self.P1range[P1], self.P2range[P2]] + [
                    getattr(self.mod, oP) for oP in output
                ]
                self.data_raw[iDx, :] = outArr
                data_dump[P2, :] = outArr
                iDx += 1
            ##  raw_file = scipy.io.write_array(raw_file, data_dump, keep_open=1)
            raw_file = scipy.io.write_array(
                os.path.join(self.data_dir, raw_data_file), data_dump, keep_open=1
            )
            del data_dump
            if flush >= 5:
                raw_file.flush()
                flush = 0
        raw_file.close()

    def D2_calculate_intersection_data(
        self, surfaces, isect_data_file, isect_tol=1.0e-2
    ):
        """Data is P1 P2 Surf1 Surf2"""
        self.isect_tol = isect_tol
        assert len(surfaces) == 2
        self.surfP = surfaces
        surfidx = [self._output.index(el) + 2 for el in surfaces]
        print(surfidx)

        isect = []
        isectCount = 0
        for row in range(self.data_raw.shape[0]):
            if (
                abs(self.data_raw[row, surfidx[0]] - self.data_raw[row, surfidx[1]])
                < isect_tol
            ):
                isect.append(self.data_raw[row])
                isectCount += 1
                self.STDOUT.write(str(isectCount) + ' intersections found!\n')

        self.isect_data = scipy.array(isect)
        ##  BB = open(os.path.join(self.data_dir,isect_data_file),'w')
        scipy.io.write_array(
            os.path.join(self.data_dir, isect_data_file), self.isect_data, keep_open=0
        )

    def output_intersection_metrics(self, filename=None):
        print('\n**********\nIntersection metrics')
        print('Raw data space : ', self.data_raw.shape)
        print('Intersect space: ', self.isect_data.shape)
        print('Intersect tol  : ', self.isect_tol)
        print(
            'Sparsity %     : ',
            '%2.1e'
            % (
                float(self.isect_data.shape[0])
                / float(self.data_raw.shape[0] * self.data_raw.shape[1])
                * 100.0
            ),
        )
        print(
            'Efficiency %   : ',
            '%2.1e'
            % (float(self.isect_data.shape[0]) / float(self.data_raw.shape[0]) * 100.0),
        )
        print(pysces.session_time())

        if filename != None:
            CC = open(os.path.join(self.data_dir, filename), 'w')
            CC.write('\n**********\nIntersection metrics\n**********\n')
            CC.write('Raw data space : ' + repr(self.data_raw.shape) + '\n')
            CC.write('Intersect space: ' + repr(self.isect_data.shape) + '\n')
            CC.write('Intersect tol  : ' + repr(self.isect_tol) + '\n')
            CC.write(
                'Sparsity %     : '
                + '%2.1e'
                % (
                    float(self.isect_data.shape[0])
                    / float(self.data_raw.shape[0] * self.data_raw.shape[1])
                    * 100.0
                )
                + '\n'
            )
            CC.write(
                'Efficiency %   : '
                + '%2.1e'
                % (
                    float(self.isect_data.shape[0])
                    / float(self.data_raw.shape[0])
                    * 100.0
                )
                + '\n'
            )
            CC.write('\n' + pysces.session_time() + '\n')
            CC.close()


class IntersectD3(Intersect):
    def __init__(self, mod, data_dir):
        Intersect.__init__(self, mod, data_dir)

    def D3_setup_args(
        self,
        min1,
        max1,
        step1,
        min2,
        max2,
        step2,
        min3,
        max3,
        step3,
        input,
        output,
        raw_data_file,
        surfaces,
        isect_data_file,
        isect_tol=1.0e-2,
        filename=None,
    ):
        self._min1 = min1
        self._max1 = max1
        self._step1 = step1
        self._min2 = min2
        self._max2 = max2
        self._step2 = step2
        self._min3 = min3
        self._max3 = max3
        self._step3 = step3
        self._input = input
        self._output = output
        self._raw_data_file = raw_data_file
        self._surfaces = surfaces
        self._isect_data_file = isect_data_file
        self._isect_tol = isect_tol
        self._filename = filename

    def D3_run(self):
        # try:
        self.D3_setup_scan(
            self._min1,
            self._max1,
            self._step1,
            self._min2,
            self._max2,
            self._step2,
            self._min3,
            self._max3,
            self._step3,
        )
        # self.load_raw_data('p4_raw.dat')
        self.D3_generate_raw_data(self._input, self._output, self._raw_data_file)
        self.D3_calculate_intersection_data(
            self._surfaces, self._isect_data_file, self._isect_tol
        )
        self.output_intersection_metrics(self._filename)
        # except Exception, EXC:
        #    print EXC
        #    sleep(5)
        #    sys.exit(1)

    def D3_setup_scan(self, min1, max1, step1, min2, max2, step2, min3, max3, step3):
        self.P1min = min1
        self.P1max = max1
        self.P1steps = step1
        self.P2min = min2
        self.P2max = max2
        self.P2steps = step2
        self.P3min = min3
        self.P3max = max3
        self.P3steps = step3
        self.P1range = scipy.logspace(scipy.log10(min1), scipy.log10(max1), step1)
        self.P2range = scipy.logspace(scipy.log10(min2), scipy.log10(max2), step2)
        self.P3range = scipy.logspace(scipy.log10(min3), scipy.log10(max3), step3)

    def D3_generate_raw_data(self, input, output, raw_data_file):
        assert len(input) == 3
        assert len(output) > 0
        # needs cleaning
        self.inputP = input
        self.outputP = output
        self.data_raw = scipy.zeros(
            (self.P1steps * self.P2steps * self.P3steps, len(output) + 3), self.dtype
        )

        iDx = 0
        flush = 0
        SflushCount = 0
        raw_file = open(os.path.join(self.data_dir, raw_data_file), 'w')
        raw_file.close()
        for P1 in range(len(self.P1range)):
            # not agaaaaaaaaaaaaaaaain
            gc.collect()
            del gc.garbage[:]
            raw_file = open(os.path.join(self.data_dir, raw_data_file), 'a')
            # raw_file.seek(-1)
            setattr(self.mod, input[0], self.P1range[P1])
            for P2 in range(len(self.P2range)):
                setattr(self.mod, input[1], self.P2range[P2])
                data_dump = scipy.zeros(
                    (len(self.P3range), len(output) + 3), self.dtype
                )
                flush += 1
                for P3 in range(len(self.P3range)):
                    setattr(self.mod, input[2], self.P3range[P3])
                    self.mod.doState()
                    outArr = [self.P1range[P1], self.P2range[P2], self.P3range[P3]] + [
                        getattr(self.mod, oP) for oP in output
                    ]
                    self.data_raw[iDx, :] = outArr
                    data_dump[P3, :] = outArr
                    iDx += 1
                    SflushCount += 1
                    if SflushCount > self.data_raw.shape[0] / 50:
                        self.STDOUT.write('\n\n**************\n')
                        self.STDOUT.write(
                            '%3i' % (float(iDx) / float(self.data_raw.shape[0]) * 100)
                            + '% complete.\n'
                        )
                        self.STDOUT.write('**************\n\n')
                        self.STDOUT.flush()
                        SflushCount = 0
                scipy.io.write_array(raw_file, data_dump, keep_open=1)
                del data_dump
                if flush > 5:
                    raw_file.flush()
                    flush = 0
            raw_file.close()
        raw_file.close()

    def D3_calculate_intersection_data(
        self, surfaces, isect_data_file, isect_tol=1.0e-2
    ):
        self.isect_tol = isect_tol
        assert len(surfaces) == 2
        self.surfP = surfaces
        surfidx = [self._output.index(el) + 3 for el in surfaces]
        print(surfidx)

        isect = []
        isectCount = 0
        for row in range(self.data_raw.shape[0]):
            if (
                abs(self.data_raw[row, surfidx[0]] - self.data_raw[row, surfidx[1]])
                < isect_tol
            ):
                isect.append(self.data_raw[row])
                isectCount += 1
                self.STDOUT.write(str(isectCount) + ' intersections found!\n')

        self.isect_data = scipy.array(isect)
        BB = open(os.path.join(self.data_dir, isect_data_file), 'w')
        scipy.io.write_array(BB, self.isect_data, keep_open=0)
