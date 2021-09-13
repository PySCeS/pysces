#!/bin/bash
DNS_SERVER=8.8.4.4
# DNS_SERVER=146.232.20.10

export DOCKER_IMAGE=quay.io/pypa/manylinux2014_x86_64
docker run --dns $DNS_SERVER --rm -v `pwd`/../..:/io $DOCKER_IMAGE /io/packaging/many_linux/build_pysces.sh
