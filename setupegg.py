#!/usr/bin/env python

"""
PySCeS - Python Simulator for Cellular Systems (http://pysces.sourceforge.net)

Copyright (C) 2004-2020 B.G. Olivier, J.M. Rohwer, J.-H.S Hofmeyr all rights reserved,

Brett G. Olivier (bgoli@users.sourceforge.net)
Triple-J Group for Molecular Cell Physiology
Stellenbosch University, South Africa

Permission to use, modify, and distribute this software is given under the
terms of the PySceS (BSD style) license.  See LICENSE.txt that came with
this distribution for specifics.

NO WARRANTY IS EXPRESSED OR IMPLIED.  USE AT YOUR OWN RISK.
Brett G. Olivier
"""

"""
A setup.py script to use setuptools, which gives egg goodness, etc.

Adapted from the original NumPy src (numpy.scipy.org).
"""
FRYING_EGGS = True
from setuptools import setup
exec(compile(open('setup.py').read(), 'setup.py', 'exec'))
