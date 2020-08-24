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

__doc__ = """
          PyscesInterfaces
          ----------------

          Interfaces converting to and from PySCeS models - makes use of Brett's Core2
          """

import os
from pysces import model_dir, output_dir
from pysces.core2.PyscesCore2 import NewCore
from pysces.core2.PyscesCore2Interfaces import (
    PscToCore,
    SbmlToCore,
    CoreToPsc,
    CoreToSBML,
)


class Core2interfaces(object):
    """
    Defines interfaces for translating PySCeS model objects into and from other formats.
    """

    core = None
    sbml = None
    core2sbml = None
    core2psc = None
    sbml2core = None
    sbml_level = 2
    sbml_version = 1

    def __buildPscComponents__(self):
        """
        Builds the PSC file
        """
        self.core2psc.setHeader()
        self.core2psc.setFixed()
        self.core2psc.setCompartments()
        self.core2psc.setFunctions()
        self.core2psc.setReactions()
        self.core2psc.setAssignmentRules()
        self.core2psc.setRateRules()
        self.core2psc.setEvents()
        self.core2psc.setSpecies()
        self.core2psc.setParameters()
        self.core2psc.setNotes()

    def __buildSbmlComponents__(self):
        """
        The SBML file build commands
        """
        self.core2sbml.createModel()
        self.core2sbml.setDescription()
        self.core2sbml.setUnits()
        self.core2sbml.setCompartments()
        self.core2sbml.setSpecies()
        self.core2sbml.setParameters()
        self.core2sbml.setReactions()
        self.core2sbml.setEvents()
        self.core2sbml.setRules()

    def writeMod2PSC(
        self, mod, filename=None, directory=None, iValues=True, getstrbuf=False
    ):
        """
        Writes a PySCeS model object to a PSC file.

        - *filename*: writes <filename>.psc or <model_name>.psc if None
        - *directory*: (optional) an output directory
        - *iValues*: if True then the models initial values are used (or the current values if False).
        - *getstrbuf*: if True a StringIO buffer is returned instead of writing to disk
        """
        self.core = NewCore(mod, iValues=iValues)
        assert self.core != None, "\nPlease set a PySCeS model or Core2 object"
        self.core2psc = CoreToPsc(self.core)
        if not iValues:
            self.core2psc.SPECIES_CURRENT_VALUE = True
            self.core2psc.FIXED_SPECIES_CURRENT_VALUE = True
            self.core2psc.PARAMETER_CURRENT_VALUE = True

        self.__buildPscComponents__()
        if filename == None:
            filename = self.core.name
        if filename[-4:] != '.psc':
            filename += '.psc'
        if not getstrbuf:
            self.core2psc.write(filename, directory, getstrbuf)
        else:
            return self.core2psc.write(filename, directory, getstrbuf)

    def writeMod2SBML(
        self,
        mod,
        filename=None,
        directory=None,
        iValues=True,
        getdocument=False,
        getstrbuf=False,
    ):
        """
        Writes a PySCeS model object to an SBML file.

        - *filename*: writes <filename>.xml or <model_name>.xml if None
        - *directory*: (optional) an output directory
        - *iValues*: if True then the models initial values are used (or the current values if False).
        - *getdocument*: if True an SBML document object is returned instead of writing to disk or
        - *getstrbuf*: if True a StringIO buffer is returned instead of writing to disk
        """
        self.core = NewCore(mod, iValues=iValues)
        ##  assert os.sys.platform != 'win32', '\nSBML translation currently only supported on Linux (Mac?)'
        assert self.core != None, "\nPlease set a PySCeS model or Core2 object"
        self.core2sbml = CoreToSBML(self.core)
        self.core2sbml.level = self.sbml_level
        self.core2sbml.version = self.sbml_version
        assert self.core2sbml.SBML != None, "\nProblems loading SBML/Python interface"
        self.__buildSbmlComponents__()
        if filename == None:
            filename = self.core.name
        if filename[-4:] != '.xml':
            filename += '.xml'
        if getdocument:
            return self.core2sbml.getSBMLdocument()
        elif getstrbuf:
            return self.core2sbml.getSBMLFileAsStrBuf()
        else:
            self.core2sbml.writeSBML2file(filename, directory)

    def readSBMLToCore(self, filename, directory=None):
        """
        Reads the SBML file specified with filename and converts it into
        a core2 object pysces.interface.core

        - *filename*: the SBML file
        - *directory*: (optional) the SBML file directory None means try the current working directory
        """
        if directory != None:
            fpath = os.path.join(directory, filename)
            assert os.path.exists(fpath), "\nFile \"%s\" does not exist!" % fpath

        self.sbml2core = SbmlToCore()
        self.sbml2core.getSbmlStringFromDisk(filename, Dir=directory)
        self.sbml2core.getSbmlModel()
        self.sbml2core.getParsedModel()
        self.sbml2core.getUnits()
        self.core = NewCore(self.sbml2core)

    def readMod2Core(self, mod, iValues=True):
        """
        Convert a PySCeS model object to core2

        - *iValues*: if True then the models initial values are used (or the current values if False).
        """
        self.core = NewCore(mod, iValues=iValues)

    def writeCore2PSC(self, filename=None, directory=None, getstrbuf=False):
        """
        Writes a Core2 object to a PSC file.

        - *filename*: writes <filename>.xml or <model_name>.xml if None
        - *directory*: (optional) an output directory
        - *getstrbuf*: if True a StringIO buffer is returned instead of writing to disk
        """
        assert self.core != None, "\nPlease set a PySCeS model or Core2 object"
        ##  print 'Using existing self.core'
        self.core2psc = CoreToPsc(self.core)
        self.__buildPscComponents__()
        if filename == None:
            filename = self.core.name
        if filename[-4:] != '.psc':
            filename += '.psc'
        if not getstrbuf:
            self.core2psc.write(filename, directory, getstrbuf)
        else:
            return self.core2psc.write(filename, directory, getstrbuf)

    def writeCore2SBML(self, filename=None, directory=None, getdocument=False):
        """
        Writes Core2 object to an SBML file.

        - *filename*: writes <filename>.xml or <model_name>.xml if None
        - *directory*: (optional) an output directory
        - *getdocument*: if True an SBML document object is returned instead of writing to disk or
        """
        assert self.core != None, "\nPlease set a PySCeS model or Core2 object"
        ##  print 'Using existing self.core'
        self.core2sbml = CoreToSBML(self.core)
        self.core2sbml.level = self.sbml_level
        self.core2sbml.version = self.sbml_version
        assert self.core2sbml.SBML != None, "\nProblems loading SBML/Python interface"
        self.__buildSbmlComponents__()
        if filename == None:
            filename = self.core.name
        if filename[-4:] != '.xml':
            filename += '.xml'
        if not getdocument:
            self.core2sbml.writeSBML2file(filename, directory)
        else:
            return self.core2sbml.getSBMLdocument()

    def convertSBML2PSC(self, sbmlfile, sbmldir=None, pscfile=None, pscdir=None):
        """
        Convert an SBML file to a PySCeS MDL input file.

        - *sbmlfile*: the SBML file name
        - *sbmldir*: the directory of SBML files (if None current working directory is assumed)
        - *pscfile*: the output PSC file name (if None *sbmlfile*.psc is used)
        - *pscdir*: the PSC output directory (if None the pysces.model_dir is used)
        """
        if sbmldir == None or not os.path.exists(sbmldir):
            sbmldir = os.getcwd()
        if not os.path.exists(os.path.join(sbmldir, sbmlfile)):
            raise RuntimeError(
                '\nSBML file \"%s\" does not exist!' % os.path.join(sbmldir, sbmlfile)
            )
        if pscdir == None or not os.path.exists(pscdir):
            pscdir = model_dir
        if pscfile == None:
            pscfile = '%s.psc' % sbmlfile
        self.readSBMLToCore(filename=sbmlfile, directory=sbmldir)
        self.writeCore2PSC(filename=pscfile, directory=pscdir, getstrbuf=False)
        print(
            '\nSBML2PSC\nin : %s\nout: %s'
            % (os.path.join(sbmldir, sbmlfile), os.path.join(pscdir, pscfile))
        )
