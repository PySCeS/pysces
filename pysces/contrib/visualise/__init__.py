'''
This __init__ file is for a PySCeS "eXtension module" (pysx), the following variables should be set for each module
'''
pysx_name = 'VisualiseSBML'
pysx_oscompat = ['posix','nt']
pysx_base_class = 'CONTRIB_visualise'

pysx_author = 'Collin Gillespie'
pysx_email = ''
pysx_affiliation = ''
pysx_web = 'http://www.basis.ncl.ac.uk/team.html'
pysx_notes = ''

from __future__ import absolute_import

from .VisualiseSBML import *
