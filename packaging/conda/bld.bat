::copy ..\_h_env\libs\libpython27.dll.a ..\_h_env\libs\libpython27.a
::copy ..\_h_env\libs\libmsvcr90.dll.a ..\_h_env\libs\libmsvcr90.a
"%PYTHON%" setup.py build --compiler=msvc --fcompiler=gfortran install --record=record.txt
if errorlevel 1 exit 1
