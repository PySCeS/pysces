#!/bin/bash

cp meta.linux.yaml meta.yaml

# building Anaconda packages for pysces under Python 3.6, 3.7 and 3.8

conda build --python 3.6 --numpy 1.15 .
conda build --python 3.7 --numpy 1.15 .
conda build --python 3.8 --numpy 1.17 .
conda build --python 3.9 --numpy 1.19 .
