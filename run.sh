#!/bin/bash
set -e

SRC="
main.cu
functions/*.cu
functions/boundary/*.cu
presets/stencil.cu
vtk.cu
"

rm -f ./main
rm -rf ./plot/vtk

nvcc -Xptxas -v -rdc=true $SRC -o main && ./main