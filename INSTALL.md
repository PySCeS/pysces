# Installing PySCeS

Binary install packages for all three OSs and recent Python versions are 
provided. Anaconda packages are available for 64-bit Windows and Linux, as well as 
macOS *x86_64* and Apple Silicon (*arm64*) architectures. Anaconda users can 
conveniently install PySCeS with:

```bash
$ conda install -c conda-forge -c pysces pysces
```  

Any dependencies will be installed automatically, including the optional dependencies 
*Assimulo*, *ipyparallel* and *libSBML*.
> **NOTE:** Anaconda packages are only provided for Python versions 3.9-3.11. The 
> reason is that Assimulo has not been ported to Python 3.12 and still depends on 
> `numpy.distutils`. As soon as this has happened, PySCeS Anaconda packages for 3.12 
> will be built.

Alternatively, you can use *pip* to
install PySCeS from PyPI. Core dependencies will be installed automatically. Wheels
are available for 64-bit Windows and Linux, as well as macOS architectures *x86_64* and
*arm64* (starting from PySCeS version 1.2.0, supporting Python 3.11 and 3.12).

```bash
$ pip install pysces
```

To install the optional dependences:

- `pip install "pysces[parscan]"` - for *ipyparallel*
- `pip install "pysces[sbml]"` - for *libSBML*
- `pip install "pysces[cvode]"` - for *Assimulo*
- `pip install "pysces[all]"` - for all of the above

>  **_NOTE:_**  Installation of *Assimulo* via `pip` may well require C and Fortran 
>  compilers to be
>  properly set up on your system, as binary packages are only provided for a
>  very limited number of Python versions and operating systems on PyPI. 
>  **This is not guaranteed to work!** If you require Assimulo, the conda
>  install is by far the easier option as up-to-date binaries are supplied
>  for all OSs and recent Python versions.

For more information on installing and configuring PySCeS please see the
[PySCeS User Guide](https://github.com/PySCeS/pysces-documentation/blob/main/source/userguide_doc.rst#installing-and-configuring)

## Compilation from source

As an alternative to a binary installation, you can also build your own PySCeS
installation from source. This requires Fortran and C compilers.

### Windows build

The fastest way to build your own copy of PySCeS is to use Anaconda for 
installing Python and the required libraries. In addition, the *RTools* compiler 
toolchain is required.

* Download and install 
  [Anaconda for Python 3](https://www.anaconda.com/download#Downloads)
* Obtain [Git for Windows](https://git-scm.com/download/win)
* Obtain the *RTools* compiler toolchain (**version 4.0.0.20220206**), either using
  [Chocolatey](https://chocolatey.org/install) 
  (`choco install rtools -y --version=4.0.0.20220206`) or by
  [direct download](https://github.com/r-windows/rtools-installer/releases/download/2022-02-06/rtools40-x86_64.exe)
  - Install in `C:\rtools40` (Chocolatey automatically installs to this path)
  - Add `c:\rtools40\ucrt64\bin` and `c:\rtools40\usr\bin` to the system `PATH`
* Create a PySCeS environment using conda and activate it:

```bash
$ conda create -n pyscesdev -c conda-forge python=3.11 numpy=1.26 scipy \ 
  matplotlib sympy packaging pip wheel ipython python-libsbml \
  assimulo meson meson-python ninja
$ conda activate pyscesdev
```

* Clone and enter the PySCeS code repository using git

```bash
(pyscesdev)$ git clone https://github.com/PySCeS/pysces.git pysces-src
(pyscesdev)$ cd pysces-src
```

* Now you can build and install PySCeS into the *pyscesdev* environment

```bash
(pyscesdev)$ pip install --no-deps --no-build-isolation .
```

### Linux build

All modern Linux distributions ship with gcc and gfortran. In addition, the Python
development headers (*python-dev* or *python-devel*, depending on your distro) need to
be installed.

Clone the source from Github as described above, change into the source directory and
run:

```bash
$ pip install .
```

### macOS build

The Anaconda build method, described above for Windows, should also work on macOS.

Alternatively, Python 3 may be obtained via [Homebrew](https://brew.sh) and the
compilers may be installed via [Xcode](https://developer.apple.com/xcode). Clone the
source from Github as described above, change into the source directory and run:

```bash
$ pip install .
```

© Johann M. Rohwer & Brett G. Olivier, January 2024
