#!/bin/bash

for PYBIN in /opt/python/cp3[678]*/bin; do
    ${PYBIN}/pip --cache-dir /io/packaging/many_linux/pip_cache install numpy
    ${PYBIN}/python /io/setup.py bdist_wheel
done


# Bundle external shared libraries into the wheels
for whl in /io/dist/*.whl; do
    auditwheel repair $whl -w /io/dist/wheelhouse/
    mv $whl $whl.bak
done
