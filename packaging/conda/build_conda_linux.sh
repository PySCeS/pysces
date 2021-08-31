#!/bin/bash

cp meta.posix.yaml meta.yaml

# building Anaconda packages for pysces under Python 3.6 - 3.9
# new automated variants build using conda_build_config.yaml
# also use faster mamba instead of conda
mamba build .
