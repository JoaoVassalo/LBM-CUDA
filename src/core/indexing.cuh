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
__device__ inline int from_id(int x, int y, int i)
{
    int x_from = x - (int)c_ix[i];
    int y_from = y - (int)c_iy[i];

    if (x_from < 0)
        x_from += Geometry::Nx;
    if (x_from >= Geometry::Nx)
        x_from -= Geometry::Nx;

    if (y_from < 0)
        y_from += Geometry::Ny;
    if (y_from >= Geometry::Ny)
        y_from -= Geometry::Ny;

    return grid_id(x_from, y_from);
}

template <int I>
__host__ __device__ __forceinline__ int momIdx(int x, int y)
{
    return (y * Geometry::Nx + x) * D2Q9::num_var + static_cast<int>(I);
}
template <int I>
__device__ __forceinline__ int layerIdx(int x)
{
    return x * D2Q9::num_var + static_cast<int>(I);
}