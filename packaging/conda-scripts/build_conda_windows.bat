:: building Anaconda packages for pysces under Python 3.6 - 3.9 (Windows)
@echo off

echo;
echo Building Anaconda packages...
echo =============================
echo;

call mamba build .

:: building wheels for Python 3.6 - 3.9
echo;
echo Building Python wheels...
echo =========================
echo;

cd ..\..

call conda create -y -n py36-build python=3.6 numpy=1.15
call conda activate py36-build
pip wheel --no-deps -w dist .
call conda deactivate
call conda env remove -y -n py36-build
rmdir /s /q %USERPROFILE%\miniconda3\envs\py36-build

call conda create -y -n py37-build python=3.7 numpy=1.15
call conda activate py37-build
pip wheel --no-deps -w dist .
call conda deactivate
call conda env remove -y -n py37-build
rmdir /s /q %USERPROFILE%\miniconda3\envs\py37-build

call conda create -y -n py38-build python=3.8 numpy=1.17
call conda activate py38-build
pip wheel --no-deps -w dist .
call conda deactivate
call conda env remove -y -n py38-build
rmdir /s /q %USERPROFILE%\miniconda3\envs\py38-build

call conda create -y -n py39-build python=3.9 numpy=1.20
call conda activate py39-build
pip wheel --no-deps -w dist .
call conda deactivate
call conda env remove -y -n py39-build
rmdir /s /q %USERPROFILE%\miniconda3\envs\py39-build

cd packaging\conda
