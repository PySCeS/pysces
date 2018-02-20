@echo off
rem This fixes a freshly created Anaconda environment
rem for the building of PySCeS binary wheels under Windows

echo Environment: %CONDA_PREFIX%
copy %CONDA_PREFIX%\libs\libmsvcr90.dll.a %CONDA_PREFIX%\libs\libmsvcr90.a
copy %CONDA_PREFIX%\libs\libpython27.dll.a %CONDA_PREFIX%\libs\libpython27.a
