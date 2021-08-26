#!/usr/bin/env python
"""
execute all python files beginning with 'test...' in CWD
"""

import glob
from time import sleep

testfiles = glob.glob('test*.py')

for testfile in testfiles:
    exec(open(testfile, 'r').read())
    sleep(2)
