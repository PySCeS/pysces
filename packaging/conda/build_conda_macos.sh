#!/bin/bash

cp meta.linux.yaml meta.yaml

# building Anaconda packages for pysces under Python 3.6 - 3.9

conda build . --python 3.6 --numpy 1.15
conda build . --python 3.7 --numpy 1.15
conda build . --python 3.8 --numpy 1.17
conda build . --python 3.9 --numpy 1.20

# building wheels for Python 3.6 - 3.9
# assumes appropriate conda environments are available (py3x-build)
# with appropriate numpy installed, create with:
# (base) $ conda create -n py36-build python=3.6 numpy=1.15
# (base) $ conda create -n py37-build python=3.7 numpy=1.15
# (base) $ conda create -n py38-build python=3.8 numpy=1.17
# (base) $ conda create -n py39-build python=3.9 numpy=1.20

source $CONDA_PREFIX/etc/profile.d/conda.sh
cd ../..
conda activate py36-build
python setup.py bdist_wheel
conda deactivate

conda activate py37-build
python setup.py bdist_wheel
conda deactivate

conda activate py38-build
python setup.py bdist_wheel
conda deactivate

conda activate py39-build
python setup.py bdist_wheel
conda deactivate

# move linked shared libraries into the wheel with delocate
# install delocate in base environment with:
# (base) $ conda install -c conda-forge delocate
cd dist
for I in *macosx_10_9_x86_64.whl; do
    delocate-wheel -w wheelhouse $I
done
