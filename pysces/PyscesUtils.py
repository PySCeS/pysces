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

__doc__ = '''
            PyscesUtils
            -----------

            The PyscesUtils module holds a collection of methods of general use such as timers and
            array export fucntionality that can be accessed as **pysces.write.**

            '''

import fileinput
import os, sys
import time


def VersionCheck(ver='0.1.5'):
    print("dead")


def str2bool(s):
    """
    Tries to convert a string to a Python boolean

    - *s* True if 'True', 'true' or'1' else False

    """
    if s in ['True', 'true', '1']:
        return True
    else:
        return False


class TimerBox:
    """
    A timer "container" class that can be used to hold user defined timers

    """

    __go__ = 1

    def __init__(self):
        if os.sys.version_info[0] < 2 or os.sys.version_info[1] < 3:
            print(
                'Your version of Python ('
                + str(os.sys.version_info[:3])
                + ') might be too old\
            to use this class. Python 2.2 or newer is required'
            )

    def normal_timer(self, name):
        """
        normal_timer(name)

        Creates a normal timer method with <name> in the TimerBox instance. Normal timers
        print the elapsed time since creation when called.

         - *name* the timer name

        """
        assert type(name) == str or type(name) == unicode

        def timer(startTime=time.time()):
            while self.__go__:
                ms = divmod(time.time() - startTime, 60)
                nowTime = str(int(ms[0])) + ' min ' + str(int(ms[1])) + ' sec'
                yield nowTime

        setattr(self, name, timer())

    def step_timer(self, name, maxsteps):
        """
        step_timer(name,maxsteps)

        Creates a step timer method with <name> in the TimerBox instance. Step timers
        print the elapsed time as well as the next step out of maxsteps when called.

        - *name* the timer name
        - *maxsteps* the maximum number of steps associated with this timer

        """
        assert type(name) == str, 'Name is a string representing the timer name'
        assert type(maxsteps) == int or type(maxsteps) == float, 'int or float needed'

        setattr(self, name + '_cstep', 0)

        def timer(startTime=time.time(), maxS=maxsteps):
            while self.__go__:
                ms = divmod(time.time() - startTime, 60)
                nowStep = (
                    'step '
                    + str(getattr(self, name + '_cstep'))
                    + ' of '
                    + str(maxS)
                    + ' ('
                    + str(int(ms[0]))
                    + ' min '
                    + str(int(ms[1]))
                    + ' sec)'
                )
                yield nowStep

        setattr(self, name, timer())

    def stop(self, name):
        """
        Delete the timer <name> from the TimerBox instance

        - *name* the timer name to delete

        """
        delattr(self, name)
        try:
            getattr(self, name + '_cstep')
            delattr(self, name + '_cstep')
        except:
            pass

    def reset_step(self, name):
        """
        Reset the number of steps of timer <name> in the TimerBox to zero

        - *name* the step timer whose steps should be reset

        """
        try:
            getattr(self, name + '_cstep')
            setattr(self, name + '_cstep', 0)
        except:
            pass


def CopyTestModels(*args, **kwargs):
    print("moved to PyscesTest")


def CopyModels(*args, **kwargs):
    print("dead")


def ConvertFileD2U(Filelist):
    """
    Converts a [Filename] from rn to n inplace no effect if the line termination is correct

    - *Filelist* a file or list of files to convert

    """
    for line in fileinput.input(Filelist, inplace=1):
        try:
            if line[-2] == '\r':
                print(line[:-2] + '\n', end=' ')
            else:
                print(line, end=' ')
        except:
            print(line, end=' ')


def ConvertFileU2D(Filelist):
    """
    Converts a [Filename] from n to rn inplace no effect if the line termination is correct:

    - *Filelist* a file or list of files to convert

    """
    for line in fileinput.input(Filelist, inplace=1):
        try:
            if line[-2] != '\r':
                print(line[:-1] + '\n', end=' ')
            else:
                print(line, end=' ')
        except:
            print('\n', end=' ')


class WriteOutput(object):
    """
    This code is adapted from:

    CBMPy: CBTools module

    Constraint Based Modelling in Python (http://cbmpy.sourceforge.net)
    Copyright (C) 2009-2017 Brett G. Olivier, VU University Amsterdam, Amsterdam, The Netherlands

    """

    def exportLabelledArray(self, arr, names, fname, sep=',', sformat='%f'):
        """
        Write a 2D array type object to file

         - *arr* the an array like object
         - *names* the list of row names
         - *fname* the output filename
         - *sep* [default=','] the column separator
         - *format* [default='%s'] the output number format

        """
        if names != None:
            assert arr.shape[0] == len(
                names
            ), '\n ...  rows must equal number of names!'
        F = open(fname, 'w')
        cntr = 0
        for r in range(arr.shape[0]):
            if names != None:
                F.write(('{}' + sep).format(names[r]))
            for c in range(arr.shape[1]):
                if c < arr.shape[1] - 1:
                    F.write((sformat + sep) % arr[r, c])
                else:
                    F.write((sformat + '\n') % arr[r, c])
            cntr += 1
            if cntr >= 250:
                F.flush()
                cntr = 1
        F.write('\n')
        F.flush()
        F.close()
        print('exported to %s' % fname)

    def exportLabelledArrayWithHeader(
        self, arr, names, header, fname, sep=',', format='%f'
    ):
        """
        Export an array with row names and header

         - *arr* the an array like object
         - *names* the list of row names
         - *header* the list of column names
         - *fname* the output filename
         - *sep* [default=','] the column separator
         - *format* [default='%s'] the output number format
         - *appendlist* [default=False] if True append the array to *fname* otherwise create a new file

        """
        if names != None:
            assert arr.shape[0] == len(
                names
            ), '\n ...  rows must equal number of names!'
        if header != None:
            assert arr.shape[1] == len(
                header
            ), '\n ...  cols must equal number of header names!'
        F = open(fname, 'w')
        cntr = 0
        if header != None:
            if names != None:
                hstr = ' ' + sep
            else:
                hstr = ''
            for h in header:
                hstr += str(h) + sep
            hstr = hstr[:-1] + '\n'
            F.write(hstr)
            del hstr
        for r in range(arr.shape[0]):
            if names != None:
                F.write(('%s' + sep) % names[r])
            for c in range(arr.shape[1]):
                if c < arr.shape[1] - 1:
                    F.write((format + sep) % arr[r, c])
                else:
                    F.write((format + '\n') % arr[r, c])
            cntr += 1
            if cntr >= 250:
                F.flush()
                cntr = 1
        F.write('\n')
        F.flush()
        F.close()
        print('exported to %s' % fname)

    def exportLabelledLinkedList(
        self, arr, fname, names=None, sep=',', format='%s', appendlist=False
    ):
        """
        Write a 2D linked list [[...],[...],[...],[...]] and optionally a list of row labels to file:

         - *arr* the linked list
         - *fname* the output filename
         - *names* the list of row names
         - *sep* [default=','] the column separator
         - *format* [default='%s'] the output number format
         - *appendlist* [default=False] if True append the array to *fname* otherwise create a new file

        """
        if names != None:
            assert len(arr) == len(names), '\n ...  rows must equal number of names!'
        if not appendlist:
            F = open(fname, 'w')
        else:
            F = open(fname, 'a')
        cntr = 0
        for r in range(len(arr)):
            if names != None:
                F.write(('%s' + sep) % names[r])
            col_l = len(arr[0])
            for c in range(col_l):
                if c < col_l - 1:
                    if arr[r][c] == 0.0:
                        F.write('0.0' + sep)
                    else:
                        try:
                            F.write((format + sep) % arr[r][c])
                        except UnicodeEncodeError:
                            F.write((format + sep) % 'uError')
                else:
                    if arr[r][c] == 0.0:
                        F.write('0.0\n')
                    else:
                        try:
                            F.write((format + '\n') % arr[r][c])
                        except UnicodeEncodeError:
                            F.write((format + '\n') % 'uError')
            cntr += 1
            if cntr >= 250:
                F.flush()
                cntr = 1
        ##  F.write('\n')
        F.flush()
        F.close()
        del arr
        if not appendlist:
            print('exported to %s' % fname)

    def exportLabelledArrayWithHeader2CSV(self, arr, names, header, fname):
        """
        Export an array with row names and header to fname.csv

         - *arr* the an array like object
         - *names* the list of row names
         - *header* the list of column names
         - *fname* the output filename

        """
        fname += '.csv'
        self.exportLabelledArrayWithHeader(
            arr, names, header, fname, sep=',', format='%f'
        )

    def exportLabelledArray2CSV(self, arr, names, fname):
        """
        Export an array with row names to fname.csv

         - *arr* the an array like object
         - *names* the list of row names
         - *fname* the output filename

        """
        fname += '.csv'
        self.exportLabelledArray(arr, names, fname, sep=',', sformat='%f')

    def exportArray2CSV(self, arr, fname):
        """
        Export an array to fname.csv

         - *arr* the an array like object
         - *fname* the output filename
         - *sep* [default=','] the column separator

        """
        fname += '.csv'
        self.exportLabelledArray(arr, None, fname, sep=',', sformat='%f')

    def exportLabelledArrayWithHeader2TXT(self, arr, names, header, fname):
        """
        Export an array with row names and header to fname.txt

         - *arr* the an array like object
         - *names* the list of row names
         - *header* the list of column names
         - *fname* the output filename

        """
        fname += '.txt'
        self.exportLabelledArrayWithHeader(
            arr, names, header, fname, sep='\t', format='%f'
        )

    def exportLabelledArray2TXT(self, arr, names, fname):
        """
        Export an array with row names to fname.txt

         - *arr* the an array like object
         - *names* the list of row names
         - *fname* the output filename

        """
        fname += '.txt'
        self.exportLabelledArray(arr, names, fname, sep='\t', sformat='%f')

    def exportArray2TXT(self, arr, fname):
        """
        Export an array to fname.txt

         - *arr* the an array like object
         - *fname* the output filename
         - *sep* [default=','] the column separator

        """
        fname += '.txt'
        self.exportLabelledArray(arr, None, fname, sep='\t', sformat='%f')
