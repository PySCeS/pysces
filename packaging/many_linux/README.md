# Building `manylinux` compatible Python wheels for PySCeS

The various *manylinux* specifications were developed to set a minimum standard in terms 
of compiler compatibility for distribution of binary Python wheels. 
See https://github.com/pypa/manylinux for background.

This directory provides a set of scripts to build PySCeS *manylinux* wheels via a docker 
image.

## Prerequisites and preparation

1. Install `docker` via your Linux distribution's repositories.

2. Install the required docker image:

    ```shell
    $ docker pull quay.io/pypa/manylinux2014_x86_64
    ```

    ---
    **NOTE:**
    
    If compatibility with earlier compiler versions is desired, the `manylinux1` or 
    `manylinux2010` images can also be used (see https://github.com/pypa/manylinux).
    
    ---
    
## Build

1. Open a terminal and navigate to this folder.

2. If desired, update `build_pysces.sh` with your Unix UID to ensure write access to the 
files during the build and after the build has completed. 

3. If desired/needed, update the `DNS_SERVER` and `DOCKER_IMAGE` variables in 
`run_pysces.sh`.

4. Run the build script:

    ```shell
    $ sh run_pysces.sh
    ```

5. Portable patched wheels are in the `dist/wheelhouse` subdirectory of the PySCeS source.

(C) Johann Rohwer, August 2021
