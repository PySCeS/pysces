# Installing PySCeS

Binary install packages for all three OSs and Python versions 3.6-3.9 are provided.
Anaconda users can conveniently install PySCeS with:

```bash
$ conda install -c conda-forge -c pysces pysces
```  

Any dependencies will be installed automatically. Alternatively, you can use *pip* to
install PySCeS from PyPI. Again, dependencies will be installed automatically.

```bash
$ pip install pysces
```

For more information on installing and configuring PySCeS please see the
[PySCeS User Guide](https://github.com/PySCeS/pysces-documentation/blob/main/source/userguide_doc.rst#installing-and-configuring)

## Compilation from source

As an alternative to a binary installation, you can also build your own PySCeS
installation from source. This requires Fortran and C compilers.

### Windows build

The fastest way to build your own copy of PySCeS is to use Anaconda Python.

* Download and install 
  [Anaconda for Python 3](https://www.anaconda.com/products/individual#Downloads>)
* Obtain [Git for Windows](https://git-scm.com/download/win)
* Create a PySCeS environment using conda and activate it:

```bash
$ conda create -n pyscesdev -c conda-forge python=3.8 numpy scipy \ 
  matplotlib sympy packaging pip wheel nose ipython python-libsbml \
  fortran-compiler assimulo 
$ conda activate pyscesdev
```

* Clone and enter the PySCeS code repository using git

```bash
(pyscesdev)$ git clone https://github.com/PySCeS/pysces.git pysces-src
(pyscesdev)$ cd pysces-src
```

* Now you can build and install PySCeS into the pyscesdev environment

```bash
(pyscesdev)$ python setup.py build
(pyscesdev)$ python setup.py install
```

### Linux build

All modern Linux distributions ship with gcc and gfortran. In addition, the Python
development headers (*python-dev* or *python-devel*, depending on your distro) need to
be installed.

Clone the source from Github as described above, change into the source directory and
run:

```bash
$ python setup.py install
```

### macOS build

The Anaconda build method, described above for Windows, should also work on macOS.

Alternatively, Python 3 may be obtained via [Homebrew](https://brew.sh) and the
compilers may be installed via [Xcode](https://developer.apple.com/xcode). Clone the
source from Github as described above, change into the source directory and run:

```bash
$ python setup.py install
```

© Brett G. Olivier & Johann M. Rohwer, August 2021
