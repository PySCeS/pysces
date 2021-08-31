#!/bin/bash

cp meta.posix.yaml meta.yaml

### building Anaconda packages for pysces under Python 3.6 - 3.9

mamba build .

### building wheels for Python 3.6 - 3.9
source $CONDA_PREFIX/etc/profile.d/conda.sh
cd ../..

conda create -y -n py36-build python=3.6 numpy=1.15
conda activate py36-build
pip wheel --no-deps -w dist .
conda deactivate
conda env remove -y -n py36-build

conda create -y -n py37-build python=3.7 numpy=1.15
conda activate py37-build
pip wheel --no-deps -w dist .
conda deactivate
conda env remove -y -n py37-build

conda create -y -n py38-build python=3.8 numpy=1.17
conda activate py38-build
pip wheel --no-deps -w dist .
conda deactivate
conda env remove -y -n py38-build

conda create -y -n py39-build python=3.9 numpy=1.20
conda activate py39-build
pip wheel --no-deps -w dist .
conda deactivate
conda env remove -y -n py39-build

### move linked shared libraries into the wheel with delocate
# install delocate in base environment with:
# (base) $ conda install -c conda-forge delocate
cd dist
for I in *macosx_10_9_x86_64.whl; do
    delocate-wheel -w wheelhouse $I
done
