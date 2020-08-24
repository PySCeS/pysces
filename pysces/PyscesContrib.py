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

__doc__ = '''pysx contributed modules interface'''

# load any unofficial user contrib stuff
from .PyscesContribUser import *

# load official contrib stuff
import os, sys, fileinput
from . import contrib as CONTRIB_module
from . import PyscesWeb

try:
    from . import PyscesSBML

    sbmlOK = 1
    libsbml = PyscesSBML.libsbml
except Exception as e:
    print('contrib: sbml not available', e)
    sbmlOK = 0


class CONTRIB_demo:
    def run(self):
        demo = CONTRIB_module.demo.test()
        demo.run()


class CONTRIB_visualise(PyscesSBML.PyscesSBML, PyscesWeb.PyscesHTML):
    """Visualise PySCeS models with PySBML"""

    def __init__(self):
        assert sbmlOK, "visualisation requires access to (lib)sbml"
        if os.sys.platform == 'win32':
            self.image_format = 'png'
        else:
            self.image_format = 'jpg'

    def static(self, mod, filename, slvl=2, dir=None):
        """
        static(mod,filename,slvl=2,dir=None)

        Generate a png/jpg image of a PySCeS model object

        Arguments:
        =========
        mod: a loaded PySCeS model
        filename: image file name
        slvl [default=2]: sbml level to use
        dir [default=None]: output directory

        """
        assert self.image_format in [
            'png',
            'jpg',
            'gif',
        ], 'Currently supported formats are:\n jpg, png, gif'
        if dir != None:
            filename = os.path.join(dir, filename)
        self.SBML_createModel(mod, filename, slvl, dir)
        self.SBML_setCompartment(mod.ModelFile[:-4].lower())
        self.SBML_setSpecies()
        self.SBML_setReactions()
        self.SBML_setModel()
        visMod = CONTRIB_module.visualise.VisualiseSBML(self.sbml_document.getModel())
        jpg = visMod.visualiseModel().static(self.image_format)
        F = open(filename + '.' + self.image_format, 'wb')
        F.write(jpg)
        F.close()

    def static_image_return(self, mod, slvl=2):
        """
        returns a binary image
        """
        assert self.image_format in [
            'png',
            'jpg',
            'gif',
        ], 'Currently supported formats are:\n jpg, png, gif'
        filename = 'dbasetemp'
        dir = None
        self.SBML_createModel(mod, filename, slvl, dir)
        self.SBML_setCompartment(mod.ModelFile[:-4].lower())
        self.SBML_setSpecies()
        self.SBML_setReactions()
        self.SBML_setModel()
        visMod = CONTRIB_module.visualise.VisualiseSBML(self.sbml_document.getModel())
        jpg = visMod.visualiseModel().static(self.image_format)

        return jpg

    def web(self, mod, filename, slvl=2, dir=None):
        """
        web(mod,filename,slvl=2,dir=None)

        Generate a png/jpg image of a PySCeS model object embedded in an HTML file

        Arguments:
        =========
        mod: a loaded PySCeS model
        filename: image/html file name
        slvl [default=2]: sbml level to use
        dir [default=None]: output directory

        """
        self.SBML_createModel(mod, filename, slvl, dir)
        if dir != None:
            filename = os.path.join(dir, filename)
        self.SBML_setCompartment(mod.ModelFile[:-4].lower())
        self.SBML_setSpecies()
        self.SBML_setReactions()
        self.SBML_setModel()
        visMod = CONTRIB_module.visualise.VisualiseSBML(self.sbml_document.getModel())
        jpg = visMod.visualiseModel().web(self.image_format)
        F = open(filename + '.' + self.image_format, 'wb')
        F.write(jpg[0])
        F.close()

        F2 = open(filename + '.html', 'w')
        F2 = self.HTML_header(F2)

        F2.write(
            '<img src="'
            + filename
            + '.'
            + self.image_format
            + '" usemap="#'
            + filename
            + '" border="0" alt="PySces image map">\n'
        )
        F2.write('<map name="' + filename + '">\n')
        F2.write(jpg[1])
        F2.write('</map>\n')
        F2 = self.HTML_footer(F2)
        F2.close()


class CONTRIB_sbw:
    """PySCeS interface to the Systems Biology Workbench"""

    def __init__(self):
        assert sbmlOK, "SBW requires access to (lib)sbml"
        self.StructAnalysis = CONTRIB_module.sbw.StructAnalysis()
        # self.MultipleViewsOfTimeSeriesData = CONTRIB_module.sbw.MultipleViewsOfTimeSeriesData()

    def LoadDrawNetwork(self):
        self.DrawNetworkGUI = CONTRIB_module.sbw.DrawNetworkGUI()


class contrib:
    """Dynamic on-availability module loading framework"""

    def __init__(self):
        for module in CONTRIB_module.mod_dict:
            if CONTRIB_module.mod_dict[module]['status']:
                loadOK = 2
                # print mod
                print('Loading module \"' + module + '\"')
                # exec('code = '+ CONTRIB_module.mod_dict[mod]['base']+'()')
                try:
                    setattr(
                        self,
                        module,
                        eval(CONTRIB_module.mod_dict[module]['base'] + '()'),
                    )
                    loadOK -= 1
                except Exception as f:
                    print(' Unable to instantiate module:', f)
                try:
                    setattr(self, '_' + module, eval('CONTRIB_module.' + module))
                    loadOK -= 1
                except Exception as f:
                    print(' Unable to load module:', f)
                if not loadOK:
                    print('ok.')
                else:
                    print(' Module \"' + module + '\" loading incomplete')
        print(' ')


del PyscesSBML, PyscesWeb
