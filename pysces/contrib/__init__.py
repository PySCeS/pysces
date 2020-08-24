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


__version__ = '0.4.0'
__doc__ = 'PySCeS contrib module loader'


import os
import pysces

tempdir = os.getcwd()
mod_path = os.path.join(pysces.install_dir, 'contrib')

os.chdir(mod_path)
mod_dir_list = [dir for dir in os.listdir(mod_path) if os.path.isdir(dir)]
mod_load_list = []
mod_dict = {}

mod_os = os.name
# mod_os = 'posix' # transform a windows machine into a posix one

print('\nAdding contrib modules ...')
for mod in mod_dir_list:
    status = 0
    author = ''
    email = ''
    base = ''
    web = ''
    try:
        exec('import ' + mod + ' as CurMod')
        author = getattr(CurMod, 'pysx_author')
        email = getattr(CurMod, 'pysx_email')
        base = getattr(CurMod, 'pysx_base_class')
        web = getattr(CurMod, 'pysx_web')
        if mod_os in getattr(CurMod, 'pysx_oscompat'):
            print('Including module \"' + mod + '\"')
            if email != '':
                print('\tAuthor(s): ' + author + ', (' + email + ')')
            else:
                print('\tAuthor(s): ' + author)
            if web != '':
                print('\t' + web)
            status = 1
            mod_load_list.append(mod)
        else:
            print(
                '\t'
                + getattr(CurMod, 'pysx_name')
                + ' only available on '
                + str(getattr(CurMod, 'pysx_oscompat'))
            )
            print('\tMaintainer ' + author + email)
            status = 0
        del CurMod
    except Exception as e:
        print('\nModule ' + mod + ' *not* included')
        print('\t', e, '\n')
    mod_dict.update({mod: {}})
    mod_dict[mod].update(
        {
            'status': status,
            'path': os.path.abspath(os.path.join(mod_path, mod, '__init__.py')),
            'author': author,
            'email': email,
            'base': base,
        }
    )
print(' ')

os.chdir(tempdir)

del pysces, os, tempdir
