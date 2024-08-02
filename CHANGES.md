# PySCeS changes (per GitHub release)

## PySCeS release 1.2.2 (August 2024)

We are pleased to announce the release of the Python Simulator for Cellular Systems:
PySCeS (https://pysces.github.io/) version 1.2.2.

This is the first release to support `numpy 2.0`. All extension modules are now  
compiled against `numpy 2.0` to work with `numpy>=2.0`. The minimum supported `numpy` 
version is now 1.23.5 in line with `scipy` and other downstream projects.

### Fixes

- Fixed `numpy 2.0` compatibility-related bugs and updated build and runtime requirements.
- Added convenience functions for loading a model directly from and exporting a model 
  directly to SBML, as well as saving as PSC.
- Fixed parallel scanning (requires `dill` to pickle model object).
- Released Anaconda packages for Python 3.12 (now that `assimulo` has been ported to 
  Python 3.12).

## PySCeS release 1.2.1 (May 2024)

We are pleased to announce the release of the Python Simulator for Cellular Systems:
PySCeS (https://pysces.github.io/) version 1.2.1. This is a bugfix release.

### Fixes

- Fixed a bug in model evaluation where a model has events and assignment rules, and the
  assignment rule contains a variable name and a parameter with the variable a
  substring match of the parameter (e.g. `E1` and `E1tot`). Previously, this led to
  garbage code during string replacements.
- Cleaned up various error messages.

## PySCeS release 1.2.0 (Feb 2024)

We are pleased to announce the release of the Python Simulator for Cellular Systems:
PySCeS (https://pysces.github.io/) version 1.2.0. This is the first release in the 1.2
series and contains some new features and various bug fixes.

### New features

There have been major changes in the build system for PySCeS. This release
implements [meson-python](https://meson-python.readthedocs.io/),
a [PEP 517](https://peps.python.org/pep-0517/) compliant Python build backend that
is suitable for building C and Fortran extension modules by implementing
the [meson build system](https://mesonbuild.com/). There is no more `setup.py` and
project build options and metadata are now in `pyproject.toml`.

The main **motivation** for this has been support for Python 3.12, which deprecated
the use of `distutils`, and consequently `numpy` dropped `numpy.distutils`. The
previous build system had relied on `numpy.distutils` for building Fortran extension
modules, which would no longer work on Python 3.12. In the mean time, `scipy`, `numpy`,
as well as a number of other packages in the scientific Python ecosystem have moved
their build systems to `meson-python`. This version thus brings PySCeS in line with
the other packages.

In addition, **binary wheels for macOS on Apple silicon** (*arm64*) are provided
for the first time (Python 3.11 and 3.12).

### What has changed?

- From a regular user perspective, not much. Installation is still via `pip` or `conda`.

- From a developer perspective: there is no longer a `setup.py`. Metadata has been
  migrated to `pyproject.toml`. Build settings are spread between `pyproject.toml` and
  various `meson.build` files. The build is straightforward using one of:
  ```bash
  $ pip wheel -w dist .
  $ python -m build .
  ```
  as long as the relevant compiler toolchain is installed (gcc and gfortran on Linux and
  macOS, RTools on Windows). Installation can be done with
  ```bash
  $ pip install .
           or
  $ pip install --no-build-isolation -e .          (for an editable install)
  ```
  Refer to [INSTALL.md](https://github.com/PySCeS/pysces/blob/main/INSTALL.md).

- `Numpy` version >=1.26 is required for the build, as older versions of `f2py` (which
  is distributed with `numpy`) are not compatible. At runtime, any `numpy` version >
  =1.23 is supported, older versions don't work due to ABI incompatibility.
- Binaries (wheels) are distributed for Python versions 3.9-3.12. Python 3.8 support is
  dropped as the latest `numpy` and `scipy` versions also no longer support it.
- Documentation has been updated to reflect these changes.

### Bug fixes

- Fixed a bug with assignment rule evaluation when one species was an exact substring of
  another species.
- Use vectorised functions `numpy.log`, `numpy.log10`, `numpy.exp` instead of
  their `math.*` counterparts to support their application to arrays in e.g. assignment
  rules.
- Removed a bunch of deprecated `scipy.*` functions that have moved to the `numpy.*`
  namespace.
- Cleaned up unused imports.

## PySCeS release 1.1.1 (Jul 2023)

We are pleased to announce the release of the Python Simulator for Cellular Systems:
PySCeS (https://pysces.github.io/) version 1.1.1. This is the first release in the 1.1
series and contains some new features and various bug fixes.

### New features

- When used in a notebook environment, the `ipympl` backend for matplotlib is now
  enabled if installed. This allows use in JupyterLab (as opposed to classic notebook).
  If `ipympl` is not installed, fallback is to the standard `nbAgg`, which is part of
  matplotlib.
- Simulation results (`mod.sim` object) can now be returned as a pandas DataFrame if
  pandas is installed, otherwise a numpy recarray is returned as before. This option is
  configurable with `custom_datatype = pandas` in the user and system configuration
  files (
  see https://pyscesdocs.readthedocs.io/en/latest/userguide_doc.html#configuration), and
  can be enabled or disabled per session or per model:
  ```python
  pysces.enablePandas()
  pysces.enablePandas(False)
  mod.enableDataPandas()
  mod.enableDataPandas(False)
  ```

### Bug fixes

- Fixed a bug in simulations with `RateRules` using Assimulo, where a wrong solver
  variable was being assigned internally.
- Fixed SBML export when assignment rules were evaluating reaction rates.
- Enabled assignment rules (forcing functions) to track parameter changes when using
  CVODE. This is needed in case events change parameters during the course of the
  simulation.

## PySCeS release 1.1.0 (Apr 2023)

We are pleased to announce a new minor release (version 1.1.0) of the Python
Simulator for Cellular Systems: PySCeS (https://pysces.github.io/). This is the
first release in the 1.1 series.

### What's new?

The most significant new feature in Version 1.1 is a major upgrade in the way
**PySCeS handles events in simulations**. The definition of events follows the
framework described in the SBML Level 3 Version 2 specification, thus making the
event handling SBML-compliant. Specifically, event **persistence** (for events with
a delay) is now handled correctly, and simultaneous events can be executed
according to their assigned **priorities**.

The new **event specification** in the PySCeS input file reads:

```
Event: <name>, <trigger>, <optional_kwargs> { <assignments> }
```

To achieve this, three new optional keyword arguments have been added as a
comma-separated list to the event specification. The general syntax for these
arguments is `<attribute>=<value>`. The keywords are:

- *delay* (float): specifies a delay between when the trigger is fired (and the
  assignments are evaluated) and the eventual assignment to the model. If this
  keyword is not specified, a value of `0.0` is assumed.
- *priority* (integer or None): specifies a priority for events that trigger at the
  same simulation time. Events with a higher priority are executed before those
  with a lower priority, while events without a priority (`None`) are executed in
  random positions in the sequence. If this keyword is not specified, a value of
  `None` is assumed.
- *persistent* (boolean): is only relevant to events with a delay, where the
  situation may occur that the trigger condition no longer holds by the time the
  delay in the simulation has passed. The persistent attribute specifies how to
  deal with this situation: if `True`, the event executes nevertheless; if `False`,
  the event does not execute if the trigger condition is no longer valid. If the
  keyword is not specified, a default of `True` is assumed.

The following event illustrates the use of a delay of ten time units with a
non-persistent trigger and a priority of 3:

```
Event: event2, geq(_TIME_, 15.0), delay=10.0, persistent=False, priority=3 {
V3 = V3*vfact2
}
```

The legacy event syntax is still supported.

### Other changes

- A new setting has been added to the settings dictionary of the `PysMod` class
  with the following default:    
  `mod.__settings__["cvode_access_solver"] = True`    
  This specifies if the Assimulo solver object is available from within the
  `PysMod` instance to make low-level changes to the integration algorithm. The
  current default emulates previous behaviour, but the setting can be changed to
  `False`, which facilitates serialization of the `PysMod` class in e.g. parallel
  computations. Previously, the attached Assimulo solver object would prevent
  serialization. Thanks @c-barry
- Documentation has been updated to reflect the changes to the event syntax.
- Various bug fixes, including dealing with deprecations for Numpy 1.24.x
  compatibility.

## PySCeS release 1.0.3 (Sep 2022)

We are pleased to announce the release of the Python Simulator for Cellular
Systems: PySCeS (https://pysces.github.io/) version 1.0.3. This is the third
release in the 1.0 series.

### What's new?

- The build-system has been adapted to make use of `scikit-build`. This gets rid of
  the `distutils` and `numpy.distutils` dependencies, which are deprecated and will
  be removed with the release of Python 3.12.

### Bug Fixes:

- Fixed CVODE defaults and set default tolerances to more sane levels
- Fixed string replacement in parsing and construction of `PieceWise` functions

## PySCeS release 1.0.2 (May 2022)

We are pleased to announce the release of the Python Simulator for Cellular Systems:
PySCeS (https://pysces.github.io/) version 1.0.2.
This is the second bug-fix release in the 1.0 series.

### Fixes:

- Fixed a number of bugs with `RateRule` execution with CVODE
- Reintroduced the functionality to track additional items such as Assignment Rules
  during a simulation with CVODE
- The Assimulo `CVODE` implementation has been updated, and legacy PySUNDIALS CVODE
  code removed
- The `RateRule` and `AssignmentRule` implementations have been checked against the
  SBML Test Suite

## PySCeS release 1.0.1 (Feb 2022)

We are pleased to announce the release of the Python Simulator for Cellular Systems:
PySCeS (https://pysces.github.io/) version 1.0.1. This is the first bug-fix release
in the 1.0 series.

### Fixes:

- Fixed references to numpy/scipy
- Fixed a bug in `mod.Simulate(userinit=1)` with CVODE
- Fixed bug where the maximal number of steps in LSODA would not be honoured from the
  `mod.__settings__["lsoda_mxstep"]` dictionary entry
- General cleanup of license files, version files, and the way requirements are handled

## PySCeS release 1.0.0 (Sep 2021)

We are pleased to announce the release of the Python Simulator for Cellular
Systems: PySCeS (https://pysces.github.io/) version 1.0.0.

**What's new in this release:**

- Re-introduced support for CVODE as integrator under Python 3 (via Assimulo),
  which brings back support of events in models
- Improved import and export of SBML
- Automatic installation of optional dependencies with `pip install "pysces
  [optional_dep]"`
- Extensive update of documentation
- Numerous bug fixes

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
from [SourceForge](https://sourceforge.net/projects/pysces/files/pysces)

© Brett G. Olivier & Johann M. Rohwer, 2004-2024

