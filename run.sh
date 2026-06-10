#!/bin/bash
set -e

# ===== Diretórios =====
SRC_DIR="./src"
BUILD="./main"

# ===== Arquivos CUDA =====
SRC=$(find $SRC_DIR -name "*.cu")

# ===== Limpeza =====
rm -f $BUILD
rm -rf ./plot/vtk

# ===== Compilação =====
nvcc \
    -std=c++20 \
    -Xptxas -v \
    -rdc=true \
    -I$SRC_DIR \
    $SRC \
    -o $BUILD

# ===== Execução =====
./main

# ===== Pós-processamento =====
source animation/.venv/bin/activate

cd src/io
python3 plot_graph.py