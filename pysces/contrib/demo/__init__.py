'''
This __init__ file is for a PySCeS "eXtension module" (pysx), the following variables can be set for each module
'''
pysx_name = 'Pysces module template'
pysx_oscompat = ['posix','nt']
pysx_base_class = 'CONTRIB_demo'

pysx_author = 'Brett G. Olivier'
pysx_email = 'bgoli@users.sourceforge.net'
pysx_affiliation = 'Triple-J Group for Molecular Cell Physiology'
pysx_web = 'http://pysces.sourceforge.net'
pysx_notes = ''

from .demotest import *
