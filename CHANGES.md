# PySCeS changes (per GitHub release)

## PySCeS release 0.9.8 (May 2020)

We are pleased to announce the release of the Python Simulator for Cellular Systems:
PySCeS (http://pysces.sourceforge.net) version 0.9.8

**What's new in this release:**

* The main change from a user perspective is that the default model directory and output
  directory on Windows have moved to `%USERPROFILE%\Pysces` and subfolders 
  (e.g. `C:\Users\<username>\Pysces`). These folders are created by default on a fresh
  install. The previous default location was `C:\Pysces`. This change brings the Windows
  version in line with the Linux and macOS versions, and moreover allows multiple users
  to maintain their individual PySCeS configurations and model files on the same Windows
  machine.

* If you are **upgrading** from a previous version: For backward compatibility, the
  original location is searched first, and used if found upon startup, so that existing
  installations should continue to function as previously.

* Users wishing to **migrate** the PySCeS directory to the new location without losing
  any work should move the folder `C:\Pysces` to `C:\Users\<username>\Pysces`, and then
  edit the `model_dir` and `output_dir` keys in the configuration
  file `C:\Users\<username>\Pysces\.pys_usercfg.ini` to reflect the new paths.

* Numerous bug fixes relating to Python 3 compatibility.

* Other bug fixes.

PySCeS has its own project on the Anaconda cloud: https://anaconda.org/pysces

Install the latest version from the new site (Python 3.6 - 3.8) using:

```bash    
conda install -c sbmlteam -c pysces pysces
```

If you are not using Anaconda, binary wheels are provided on PyPI for Windows, Linux and
macOS (Python 3.6 - 3.8), which can be installed using:

```bash
pip install pysces
```

With Python 2 having reached its end of life, binaries for Python 2.7 are no longer
provided. However, the PySCeS codebase continues to run under Python 2.7 and can be
compiled from source if needed.

## PySCeS release 0.9.7 (March 2019)

We are pleased to announce the release of the Python Simulator for Cellular Systems:
PySCeS (http://pysces.sourceforge.net) version 0.9.7

PySCeS is now available for Python 3, many thanks to @jmrohwer and now has it's own
project on the Anaconda cloud: https://anaconda.org/pysces

Install (0.9.7+) from the new site (Python 2.7, 3.6 and 3.7) using:

```bash
conda install -c sbmlteam -c pysces pysces 
```

## PySCeS release 0.9.6 (April 2018)

I am pleased to announce the release of the Python Simulator for Cellular Systems:
PySCeS (http://pysces.sourceforge.net) version 0.9.6.

This is mostly a bugfix release that enables SciPy 1.0 compatibility.

PySCeS has been updated to support the PyscesToolkit and is now available via the
Anaconda cloud.

## PySCeS release 0.9.5 (Aug 2017)

I am pleased to announce the release of the Python Simulator for Cellular Systems:
PySCeS (http://pysces.sourceforge.net) version 0.9.5.

PySCeS has been updated to support the PyscesToolkit and is now available via the
Anaconda cloud.

## PySCeS version 0.9.3 (March 2016)

I am pleased to announce the release of the Python Simulator for Cellular Systems:
PySCeS (http://pysces.sourceforge.net) version 0.9.3.

This release contains new functionality (response coefficients of moiety conserved
cycles), bugfixes (assignment rules) and updates (SED-ML and COMBINE archive export).

## PySCeS releases 2004 - 2016

These releases are available for download
from [SourceForge] (https://sourceforge.net/projects/pysces/files/pysces)

© Brett G. Olivier & Johann M. Rohwer, August 2021

