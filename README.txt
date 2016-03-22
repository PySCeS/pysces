PySCeS - Python Simulator for Cellular Systems
----------------------------------------------

Copyright (c) 2004 - 2016, Brett G. Olivier, Johann M. Rohwer and Jan-Hendrik S. Hofmeyr
All rights reserved.

PySCeS is distributed under a BSD style licence, please see LICENCE.txt for details

Author information
------------------

Brett G. Olivier

Vrije Universiteit Amsterdam
Faculty of Earth and Life Sciences
Department of Molecular Cell Physiology
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

Installation
------------

Please see INSTALL.txt for information on building PySCeS an its modules.


Additional licence information
------------------------------

PySCeS has interfaces to and is distributed with additional packages that
are supplied under their own licence conditions.

NLEQ2
~~~~~

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

Other software (or parts there of) distributed with PySCeS
----------------------------------------------------------

PLY
~~~

David Beazley's PLY 3.3 (http://www.dabeaz.com/ply/). This package is included
unchanged and consists of two files lex.py and yacc.py in the lib/ and core2/
modules.


Brett G. Olivier, 8 June 2010
