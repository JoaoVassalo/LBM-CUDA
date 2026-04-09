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

nvcc -Xptxas -v -rdc=true -I. $SRC -o main && ./main

source animation/.venv/bin/activate
cd animation
python3 plot_graph.py