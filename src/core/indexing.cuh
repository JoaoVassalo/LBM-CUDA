#pragma once

#include "../config/geometry.h"
#include "../config/stencil.cuh"
#include "../config/layer_config.cuh"

#include "../lbm/build/build_mom.cuh"

#include "rest.cuh"

#include <cmath>

__host__ __device__ __forceinline__ int grid_id(int x, int y)
{
    return y * Geometry::Nx + x;
}

template <int I>
__host__ __device__ __forceinline__ int momIdx(int x, int y)
{
    return (y * Geometry::Nx + x) * D2Q9::num_var + I;
}

template <int I>
__host__ __device__ __forceinline__ int momIdx(int index)
{
    return index * D2Q9::num_var + I;
}

template <int I>
__host__ __device__ __forceinline__ int layerIdx(int x)
{
    return x * D2Q9::num_var + static_cast<int>(I);
}