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

import sys, pickle
import multiprocessing
from pysces.PyscesParScan import Analyze
from time import sleep


if __name__ == '__main__':
    print("\nmultiprocessing script command line:")
    print(sys.argv)
    # called with
    # subprocess.call(['python', MULTISCANFILE, self._MODE_, fN])
    MODE = sys.argv[1]
    pool = multiprocessing.Pool()

    # load stuff from the pickle
    with open(sys.argv[2], 'rb') as F:
        mod, scanpartition, seqpartition, genorder, useroutputlist = pickle.load(F)

    mod.SetQuiet()  # kill verbose output during scan
    # append tasks to asynchronous results list
    arl = []
    for i in range(len(scanpartition)):
        arl.append(
            pool.apply_async(
                Analyze,
                (
                    scanpartition[i],
                    seqpartition[i],
                    genorder,
                    MODE,
                    useroutputlist,
                    mod,
                ),
            )
        )

    print("Submitted tasks:", len(arl))
    print("\nPARALLEL COMPUTATION\n--------------------")

    while arl[-1].ready() != True:
        sleep(5)
        print('Tasks to go... ', [r.ready() for r in arl].count(False))
        # wait until all tasks are completed
    arl[-1].wait()

    # get results
    print("\nGATHER RESULT\n-------------")
    res_list = []
    for ar in arl:
        res_list.append(ar.get())
    pool.close()
    # pickle results_list
    with open(sys.argv[2], 'wb') as F:
        pickle.dump(res_list, F, protocol=-1)
        F.flush()
