#!/bin/bash

for PYBIN in /opt/python/cp[23][67]*/bin; do
    ${PYBIN}/pip install numpy
    ${PYBIN}/python /io/setup.py bdist_wheel
done


# Bundle external shared libraries into the wheels
for whl in /io/dist/*.whl; do
    auditwheel repair $whl -w /io/dist/wheelhouse/
    mv $whl $whl.bak
done
