:: building Anaconda packages for pysces under Python 3.6, 3.7 and 3.8 (Windows)
@echo off

copy meta.windows.yaml meta.yaml

call conda build --python 3.6 --numpy 1.15 .
call conda build --python 3.7 --numpy 1.15 .
call conda build --python 3.8 --numpy 1.17 .
call conda build --python 3.9 --numpy 1.19 .

:: building wheels for Python 3.6, 3.7 and 3.8
:: assumes appropriate conda environments are available (py3x-build)
:: with appropriate numpy installed (1.15 for Py36 and Py37, 1.17 for Py38)
:: create with:
:: (base) C:\Users\user> conda create -n py36-build python=3.6 numpy=1.15
:: (base) C:\Users\user> conda create -n py37-build python=3.7 numpy=1.15
:: (base) C:\Users\user> conda create -n py38-build python=3.8 numpy=1.17
:: (base) C:\Users\user> conda create -n py39-build python=3.9 numpy=1.19
cd ..\..

call conda activate py36-build
pip wheel --no-deps -w dist .
call conda deactivate

call conda activate py37-build
pip wheel --no-deps -w dist .
call conda deactivate

call conda activate py38-build
pip wheel --no-deps -w dist .
call conda deactivate

call conda activate py39-build
pip wheel --no-deps -w dist .
call conda deactivate

cd packaging\conda
