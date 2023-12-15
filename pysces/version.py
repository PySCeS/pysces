#!/usr/bin/env python3

import re
import os

cwd = os.path.dirname(__file__)

with open(os.path.join(cwd, 'version.txt'), 'r') as f:
    v = f.read().strip().split('.')

MAJOR = int(v[0])
MINOR = int(v[1])
MICRO = int(re.split(r'(\d+)', v[2])[1])
STATUS = re.split(r'(\d+)', v[2])[2]

def current_version():
    return '%s.%s.%s%s' % (MAJOR, MINOR, MICRO, STATUS)

def current_version_tuple():
    return (MAJOR, MINOR, MICRO)

__version__ = current_version()

if __name__ == '__main__':
    print(__version__)
