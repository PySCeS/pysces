@echo off
rem This creates an Anaconda environment suitable
rem for the building of PySCeS binary wheels under Windows

if "%1"=="" (
	echo Provide name of new environment as first argument
	goto:end
)

conda create -y -n %1 certifi icc_rt intel-openmp libpython mkl numpy=1.11.3 pip python=2.7 setuptools vc vs2008_runtime wheel wincertstore

:end
endlocal
	