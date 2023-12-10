This is an additional package that is distributed unchanged, with permission
from ZIB, under the terms of their own licence and warranty as described in nleq2.f
and *not* the PySCeS BSD style licence.

Please take careful notice of the usage conditions, licence and warranty as
described in nleq2_readme.txt and nleq2.f. If you do not agree with these usage
conditions you must either disable the installation of nleq2 or obtain a licence from ZIB.

The co-installation of the nleq2 non-linear solver can be disabled in one of two ways:

1. Use the environment variable PYSCES_SKIP, i.e.
  $ PYSCES_SKIP=nleq2 python -m build .    (to build)
or
  $ PYSCES_SKIP=nleq2 pip install .        (to install)

2. In the user configuration section at the top of config.py, set nleq to False:
-- config.py --
##### USER CONFIGURATION SECTION #####
  nleq2 = False

