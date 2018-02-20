#!/bin/bash
DNS_SERVER=8.8.4.4

export DOCKER_IMAGE=quay.io/pypa/manylinux1_x86_64
docker run --dns $DNS_SERVER --rm -v `pwd`:/io $DOCKER_IMAGE /io/build_pysces.sh
export DOCKER_IMAGE=quay.io/pypa/manylinux1_i686
docker run --dns $DNS_SERVER --rm -v `pwd`:/io $DOCKER_IMAGE /io/build_pysces.sh
