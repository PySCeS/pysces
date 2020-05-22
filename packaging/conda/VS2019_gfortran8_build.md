Building with VS 2019-C++  and MINGW64-gfortran 8.1

* It appears that for newer Pythons one should use Visual C++ as a compiler, rather than MINGW64-GCC.
* For newer versions of numpy, you apparently need a recent version of MINGW64-gfortran. The one that comes with Anaconda (as part of the mingw package) is 4.7.0 which does not work!

Here is how I have managed to build a wheel for Python 3.7 that seems to work:

* I have Visual Studio 2017 installed for C++ stuff, however, the VS build tools are available for free. Try here for download details: https://wiki.python.org/moin/WindowsCompilers>
** Anaconda might also have the VC build tools (package vc) but I have never used it.
* Create a conda environment: conda create -n "pscbuild" python=3.7 ipython pip wheel numpy scipy matplotlib
** Alternatively, if using an existing environment, make sure the Anaconda mingw package is not installed (remove it now)
* I unpacked MINGW64 that I obtained here (https://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Win64/Personal%20Builds/mingw-builds/8.1.0/threads-posix/seh/x86_64-8.1.0-release-posix-seh-rt_v6-rev0.7z/download) based on the information I found here (https://stackoverflow.com/questions/48826283/compile-fortran-module-with-f2py-and-python-3-6-on-windows-10)
* In a terminal I then activated the conda environment and added the new GCC to the path (making sure there was no existing gcc)
** conda activate pscbuild
** path=c:\tools\mingw64\bin;%path%
* Build/install PySCeS with python setup.py build --compiler=msvc --fcompiler=gfortran install bdist_wheel

For conda build we need the following:

* Open a VS2015 terminal (v10) and change to your build directory
** path=c:\tools\mingw64\bin;%path%
** path="c:\Program Files (x86)\Windows Kits\10\bin\10.0.18362.0\x64";%path%
* 3.6 <= Python < 3.8
** conda build --python 3.6 --numpy 1.15 .
** conda build --python 3.7 --numpy 1.15 .
* Python >= 3.8
** conda build --python 3.8 --numpy 1.18 .

* cleanup and upload
** conda build purge
** conda clean --all
** anaconda upload -u pysces pysces-builds\*
 


Brett 2020