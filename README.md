[![CI build wheels](https://github.com/PySCeS/pysces/actions/workflows/cibuildwheel.yml/badge.svg)](https://github.com/PySCeS/pysces/actions/workflows/cibuildwheel.yml)
[![CI build Anaconda](https://github.com/PySCeS/pysces/actions/workflows/build-conda.yml/badge.svg)](https://github.com/PySCeS/pysces/actions/workflows/build-conda.yml)
[![Documentation Status](https://readthedocs.org/projects/pyscesdocs/badge/?version=latest)](https://pyscesdocs.readthedocs.io/en/latest/?badge=latest)

# PySCeS - Python Simulator for Cellular Systems

PySCeS is a flexible, user-friendly tool for the analysis of cellular systems. For more
information please see https://pysces.github.io/

## Table of contents

* [Introduction](#introduction)
* [Installation](#installation)
* [Getting help](#getting-help)
* [Citing PySCeS](#citing-pysces)
* [Contributing](#contributing)
* [Licence](#licence)
* [Authors and acknowledgments](#authors-and-acknowledgments)

## Introduction

Computer modelling has become an integral tool in the analysis and understanding of the
reaction networks that underlie cellular processes. PySCeS, first released in 2003, is
extremely flexible, user-extensible, open source, software actively used and developed
by a community of researchers and developers.

### Features

* **Simulate** your time courses with solvers like LSODA and CVODE.
* Efficiently determine **steady states** using a selection of non-linear, root-finding
  algorithms (e.g. HYBRD, NLEQ2).
* Use **Metabolic Control Analysis** (MCA) to investigate control and regulation of
  cellular systems (elasticities, flux- and concentration-control and response
  coefficients).
* Auto-generate the stoichiometric matrix and perform **structural analysis**, determine
  nullspaces and reduced stoichiometric and MCA matrices.
* Investigate systems which exhibit multiple (stable and unstable) steady-state
  solutions using a PITCON-based bifurcation analysis module.
* Define models with a human-readable **Model Description Language**.
* Use the full power of Python available to test, build and share your modelling
  experiments.
* User-friendly methods for generating and 1-, 2- and n-dimensional **parameter scans**.
* Distributed parameter scanning with the **ParScanner** module within the *ipyparallel*
  framework.
* **Visualise** results of simulations with flexible Matplotlib and Gnuplot interfaces.
* PySCeS supports open community standards for **model
  exchange** ([SBML](http://sbml.org/)) and archiving (SEDML, OMEX).

## Installation

PySCeS releases are available on [PyPI](https://pypi.org/project/pysces/#files), 
[Anaconda.org](https://anaconda.org/pysces/pysces) and 
[Github](https://github.com/PySCeS/pysces)

PySCeS binaries are available for Python 3.6+ on Windows, Linux and OSX and can be
installed using

**PyPI**

```bash
pip install pysces
```

**Anaconda.org**

```bash
conda install -c conda-forge -c pysces pysces
```

**Building from source**

PySCeS requires a Python 3 development environment as well as C++ and FORTRAN compilers
to build. For more information please see the [INSTALL](./INSTALL.md) instructions.

## Getting help

We currently encourage users to ask questions and request support using a
GitHub [issue](https://github.com/PySCeS/pysces/issues).

General information on PySCeS can be found on the [website](https://pysces.github.io/)
and extensive documentation on installing, configuring and using PySCeS is available
on [pyscesdocs.readthedocs.io](https://pyscesdocs.readthedocs.io/en/latest/).

## Citing PySCeS

**Publication**

Brett G. Olivier, Johann M. Rohwer and Jan-Hendrik S. Hofmeyr, **Modelling cellular
systems with PySCeS**,
*Bioinformatics*, Volume 21, Issue 4, 15 February 2005, Pages
560–561, https://doi.org/10.1093/bioinformatics/bti046

**Published code releases**

PySCeS releases are archived on [Zenodo](http://zenodo.org). 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2600905.svg)](https://doi.org/10.5281/zenodo.2600905)

## Contributing

Primary development is on GitHub, feel free to fork and submit pull requests! Want to
add new code to PySCeS, go 
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

## Licence

PySCeS is distributed as Open Source software under a BSD style licence, please
see [LICENCE.txt](./LICENCE.txt) for details. The licences of software that is bundled
with, or used by, PySCeS are available
in [ADDITIONAL_LICENCES.md](./ADDITIONAL_LICENCES.md).

## Authors and acknowledgments

PySCeS main authors are [Brett Olivier](https://research.vu.nl/en/persons/bg-olivier)
and [Johann Rohwer](https://github.com/jmrohwer) and in particular, PySCeS founding
father [Jannie Hofmeyr]().

We would also like to acknowledge everyone who has contributed to PySCeS and would like
to thank the following (in alphabetical order):

* Evert Bosdriesz
* James Dominy
* Jannie Hofmeyr
* Morgan Hough
* Jonathan Karr
* Ruchir Khandelwal
* Charl Moller
* Danie Palm
* Herbert Sauro
* Lafras Uys

This information can also be found in [CONTRIBUTORS.md](./CONTRIBUTORS.md)
In addition we would like to acknowledge financial support from the following funding
organisations:

* The South African National Bioinformatics Network (2003-2008)
* The South African National Research Foundation (2000-2004)

© Brett G. Olivier & Johann M. Rohwer, August 2021
