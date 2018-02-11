export DOCKER_IMAGE=quay.io/pypa/manylinux1_x86_64
sudo docker run --rm -v `pwd`:/io $DOCKER_IMAGE /io/build_pysces.sh
export DOCKER_IMAGE=quay.io/pypa/manylinux1_i686
sudo docker run --rm -v `pwd`:/io $DOCKER_IMAGE /io/build_pysces.sh
