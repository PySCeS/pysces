'''
This __init__ file is for a PySCeS "eXtension module" (pysx), the following variables can be set for each module
'''
pysx_name = 'PySCeS SBW interface'
pysx_oscompat = ['nt']
pysx_base_class = 'CONTRIB_sbw'

pysx_author = 'Herbert Sauro et al.'
pysx_email = ''
pysx_affiliation = 'Kegg Graduate Institute'
pysx_web = 'http://www.sys-bio.org'
pysx_notes = 'This is an interface to the Systems Biology Workbench'

from __future__ import absolute_import

from .sbw_func import *
