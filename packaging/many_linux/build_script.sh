#!/bin/bash
export PATH=$HOME/.local/bin:$PATH

for PYBIN in /opt/python/cp3[67]*/bin; do
    ${PYBIN}/pip --cache-dir /io/packaging/many_linux/pip_cache install numpy~=1.15.0 --user
    ${PYBIN}/python /io/setup.py bdist_wheel
done

for PYBIN in /opt/python/cp38*/bin; do
    ${PYBIN}/pip --cache-dir /io/packaging/many_linux/pip_cache install numpy~=1.17.0 --user
    ${PYBIN}/python /io/setup.py bdist_wheel
done

# Bundle external shared libraries into the wheels
for whl in /io/dist/pysces*linux*.whl; do
    auditwheel repair $whl -w /io/dist/wheelhouse/
    mv $whl $whl.bak
done
