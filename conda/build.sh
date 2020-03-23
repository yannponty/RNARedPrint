#!/bin/sh

make && \
make install PREFIX=${CONDA_PREFIX}
