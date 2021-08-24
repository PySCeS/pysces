# Additional licences

Licences for other software (or parts thereof) distributed with, or used by, PySCeS.

## NLEQ2

The NLEQ2 algorithm from ZIB which has its own licence and warranty and does not fall
under the PySCeS GPL licence. ZIB has kindly given us permission to distribute NLEQ2
with PySCeS provided its use falls within their terms of usage (/nleq2/nleq2_readme):

You may use or modify this code for your own non-commercial purposes for an unlimited
time. In any case you should not deliver this code without a special permission of ZIB.
In case you intend to use the code commercially, we oblige you to sign an according
licence agreement with ZIB.

If you are not using PySCeS for research or personal use you **must** disable NLEQ2 by
setting `nleq2 = 0` in setup.py before installation. This is a plug-in solver and its
removal does not otherwise influence PySCeS.

## PLY

David Beazley's PLY 3.3 (http://www.dabeaz.com/ply/). This package is included unchanged
and consists of two files *lex.py* and *yacc.py* in the `pysces/lib` and `pysces/core2`
modules.

© Brett G. Olivier & Johann M. Rohwer, August 2021
