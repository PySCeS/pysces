#!/bin/bash

for PYBIN in /opt/python/cp27*/bin; do
    ${PYBIN}/pip install numpy
    ${PYBIN}/pip wheel pysces
done


# Bundle external shared libraries into the wheels
for whl in *.whl; do
    auditwheel repair $whl -w /io/wheelhouse/
done
