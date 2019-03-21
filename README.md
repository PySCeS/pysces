# PySCeS - Python Simulator for Cellular Systems

The official code repository (http://github.com/PySCeS) of The Python Simulator for Cellular Systems: PySCeS project: http://pysces.sourceforge.net release code will be mirrored on SF SVN and all official source archives and binary releases made available from SourceForge.

PySCeS is available as an Anaconda package for 64 bit Windows, Linux and OSX:

```
conda install -c sbmlteam -c pysces pysces

```

Copyright (c) 2004 - 2019, Brett G. Olivier, Johann M. Rohwer and Jan-Hendrik S. Hofmeyr
All rights reserved.

PySCeS is distributed under a BSD style licence, please see LICENCE.txt for details

# Referring to the latest release

PySCeS is archived on Zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2600906.svg)](https://doi.org/10.5281/zenodo.2600906):
 
 Brett G Olivier & Johann Rohwer. (2019, March 21). PySCeS/pysces: PySCeS release 0.9.7 (Version 0.9.7). Zenodo. http://doi.org/10.5281/zenodo.2600905


## Author information

Brett G. Olivier

Vrije Universiteit Amsterdam
Faculty of Earth and Life Sciences
Systems Bioinformatics
De Boelelaan 1085
NL-1081 HV AMSTERDAM
The Netherlands

previous address

Triple-J Group for Molecular Cell Physiology
Department of Biochemistry
University of Stellenbosch
Private Bag X1
Matieland
7602

email: bgoli@users.sourceforge.net
web:   http://pysces.sourceforge.net

# Installation

Please see INSTALL.txt for information on building PySCeS an its modules.

PySCeS is available as an Anaconda package for 64 bit Windows, Linux and OSX:

```
conda install -c sbmlteam -c pysces pysces
```

[![Research software impact](http://depsy.org/api/package/pypi/PySCeS/badge.svg)](http://depsy.org/package/python/PySCeS)

# Additional licence information

PySCeS has interfaces to and is distributed with additional packages that
are supplied under their own licence conditions.

## NLEQ2

The NLEQ2 algorithm from ZIB which has it's own licence and warranty
and does not fall under the PySCeS GPL licence. ZIB has kindly given us
permission to distribute NLEQ2 with PySCeS provided it's use falls
within their terms of usage (/nleq2/nleq2_readme):

You may use or modify this code for your own non commercial
purposes for an unlimited time.
In any case you should not deliver this code without a special
permission of ZIB.
In case you intend to use the code commercially, we oblige you
to sign an according licence agreement with ZIB.

If you are not using PySCeS for research or personal use you **must** disable
NLEQ2 by setting nleq2 = 0 in setup.py before installation. This is a plug-in
solver and its removal does not otherwise influence PySCeS.

# Other software (or parts there of) distributed with PySCeS

## PLY

David Beazley's PLY 3.3 (http://www.dabeaz.com/ply/). This package is included
unchanged and consists of two files lex.py and yacc.py in the lib/ and core2/
modules.

Brett G. Olivier, 21 March 2019
