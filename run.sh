#!/bin/bash

set -e

SRC_DIR=src
OUT=main

FLAGS="-std=c++20 -rdc=true"

if [ "$1" = "debug" ]; then
    FLAGS="$FLAGS -g -G -O0"
else
    FLAGS="$FLAGS -O3"
fi

rm -f $OUT
mkdir -p $SRC_DIR/io/out

rm -f plot/vtk/*.vti

nvcc $FLAGS \
    -lineinfo \
    -Xptxas=-v \
    -o $OUT \
    $(find $SRC_DIR -name "*.cu")

./$OUT