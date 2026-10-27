#!/bin/bash

set -e

src_dir=src
out=main
flags=(-std=c++20 -rdc=true -O3)

echo "Removendo executável e arquivos VTI anteriores..."
rm -f "$out" plot/vtk/*.vti

echo "Compilando a simulação..."
mapfile -t sources < <(find "$src_dir" -name "*.cu" -print)

compile_start=$(date +%s%N)
nvcc "${flags[@]}" \
    -o "$out" \
    "${sources[@]}"
compile_end=$(date +%s%N)
compile_time_ms=$(((compile_end - compile_start) / 1000000))
printf 'Compilação concluída em %d.%03d s.\n' "$((compile_time_ms / 1000))" "$((compile_time_ms % 1000))"

echo "Executando a simulação..."
./"$out"
