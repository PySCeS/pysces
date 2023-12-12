#!/bin/bash

# building Anaconda packages for pysces under Python 3.6 - 3.9
# new automated variants build using conda_build_config.yaml
# also use faster mamba instead of conda

echo "\nBuilding Anaconda packages..."
echo "=============================\n"

mamba build .
