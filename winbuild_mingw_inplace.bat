::python setup.py build bdist_wheel sdist

::conda activate pscdev
set path=%path%;c:\tools\mingw64\bin

python setup.py develop

pause


