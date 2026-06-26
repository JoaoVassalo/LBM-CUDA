#pragma once

#include "../config/stencilConfig.cuh"

__host__ __device__ __forceinline__ int gridId(int x, int y)
{
    return y * Geometry::NX + x;
}

template <int I>
__host__ __device__ __forceinline__ int momIdx(int x, int y)
{
    return gridId(x, y) * D2Q9::momNum + I;
}

template <int I>
__host__ __device__ __forceinline__ int momIdx(int index)
{
    return index * D2Q9::momNum + I;
}

template <auto I>
__host__ __device__ __forceinline__ int layerIdx(auto x)
{
    return x * D2Q9::momNum + I;
}