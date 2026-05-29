#pragma once

#include "../../core/constexpr_for.cuh"
#include "../../core/indexing.cuh"
#include "../../config/layer_config.cuh"
#include "../build/build_mom.cuh"

__global__ void init_layers(D2Q9 sim, layer layer);

__global__ void swap_layers(D2Q9 sim, layer layer);

__device__ int grid_id_from_layer(int yref)
{
    int x = threadIdx.x + blockDim.x * blockIdx.x;
    int y = yref + layer::LNy;

    return x + Geometry::Nx * y;
}
