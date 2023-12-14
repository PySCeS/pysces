# IMPORTANT INFORMATION
## OUTDATED!

The information in this file is **outdated**! For updated building instructions 
refer to [INSTALL.md](../../INSTALL.md).

Note that *updated conda packages* are distributed via 
[anaconda.org](https://anaconda.org/pysces/pysces).

# Building Anaconda packages and Python wheels

This guide describes building **portable** Anaconda packages and binary Python wheels for 
distribution of PySCeS on all supported OSs (Linux, Windows, macOS). With the 
exception of Linux (see below), Python wheels are built in the Anaconda environment as 
well and a single build script provides both Anaconda packages and Python wheels.

Build variants (builds for different Python versions built with the correct associated 
NumPy version) are handled automatically by the build scripts making use of 
`conda_build_config.yaml` and Jinja 2 templates in the various `meta.yaml` files. 
Normally, it should not be necessary to mess with this other than to customise your build 
variants in `conda_build_config.yaml`.

## Contents

- [Prerequisites](#prerequisites)
- [Linux](#linux)
- [Windows](#windows)
- [macOS](#macos)

## Prerequisites

1. An updated Anaconda 64-bit environment for your OS. I recommend using 
[Miniconda3](https://docs.conda.io/en/latest/miniconda.html). Your base environment 
should have `conda-build>=3.0` installed.

2. Set conda-forge as the default channel:

    ```shell
    $ conda config --add channels conda-forge
    ```

3. (Optional) Install `mamba` for faster builds and dependency resolution of environments:

    ```shell
    $ conda install mamba
    ```

OS-specific instructions are provided below.

## Linux

1. Install `gcc` and `gfortran` via your Linux distribution.

2. Activate your Anaconda base environment and navigate to this directory.

3. Run the build script:

    ```shell
    $ bash ./build_conda_linux.sh
    ```

4. Depending on where you installed Anaconda/Miniconda, the binary packages will be in
`$HOME/miniconda3/conda-bld/linux-64/` or similar. You have to upload the packages to the 
Anaconda repo yourself, the script does not do this.

5. If desired, clean up disk space:

    ```shell
    $ conda build purge
    $ conda clean -a
    ```

6. Portable wheels on Linux have to be built with the `manylinux` specification; this is 
achieved with a docker image. See instructions under `packaging/many_linux/`.
    
## Windows

1. Get [Microsoft Visual C++ Build Tools for Visual Studio 
2019](https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2019). 
You can also get the complete VS Community 2019 edition, but it is not needed for this 
compile.

2. Download and install [MinGW-W64](https://sourceforge.net/projects/mingw-w64/) 
(x86_64-win32-seh build). Link to 
[online installer](https://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Win32/Personal%20Builds/mingw-builds/installer/mingw-w64-install.exe).
*Modify the batch file* that launches the MinGW-W64 Command Prompt to set the paths and also 
activate the conda base environment (modify the location of your Anaconda/Miniconda 
installation as needed):

    ```bat
    @echo off
    set PATH=C:\Program Files\mingw-w64\x86_64-8.1.0-win32-seh-rt_v6-rev0\mingw64\bin;%PATH%
    cd %USERPROFILE%
    C:\Windows\system32\cmd.exe /K %USERPROFILE%\miniconda3\condabin\activate.bat
    ```
3. Activate the environment with this modified shortcut, navigate to this directory, and run the build script:

    ```shell
    $ build_conda_windows.bat
    ```

4. Depending on where you installed Anaconda/Miniconda, the binary packages will be in 
`%USERPROFILE%\miniconda3\conda-bld\win-64\` or similar. You have to upload the packages 
to the Anaconda repo yourself, the script does not do this. Binary wheels will be in the 
`dist` subfolder of your PySCeS source directory.

5. If desired, clean up disk space:

    ```shell
    $ conda build purge
    $ conda clean -a
    ```

## macOS

1. Install the Command Line Tools for Xcode by running the following from a terminal and 
following the prompts:

    ```shell
    $ xcode-select --install
    ```

    You can of course install the complete Xcode, but it is a huge multi-Gb install, and 
    only the Command Line Tools are needed.

2. Create a file `conda_build_config.yaml` in your `$HOME` directory with the following 
content:

    ```yaml
    CONDA_BUILD_SYSROOT:
      - /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
    ```

3. Download and install 
[gfortran](https://github.com/fxcoudert/gfortran-for-macOS/releases/tag/8.2). This will 
set up the compiler executable automatically to be in the path
(`/usr/local/bin/`).

4. Open a terminal, activate the conda base environment,

    ```shell
    $ conda activate
    ```
    
    Install `delocate` (needed to create portable Python wheel files).
    
    ```shell
    $ conda install delocate
    ```
    
5. Change to this directory and run the build script:

    ```shell
    $ sh ./build_conda_macos.sh
    ```

6. Depending on where you installed Anaconda/Miniconda, the binary packages will be in
`$HOME/miniconda3/conda-bld/osx-64/` or similar. You have to upload the packages to the 
Anaconda repo yourself, the script does not do this. Binary wheels will be in the 
`dist/wheelhouse` subfolder of your PySCeS source directory.

5. If desired, clean up disk space:

    ```shell
    $ conda build purge
    $ conda clean -a
    ```

(C) Johann Rohwer, August 2021
