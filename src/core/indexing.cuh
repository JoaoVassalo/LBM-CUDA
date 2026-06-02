#pragma once

#include "../config/geometry.h"
#include "../config/stencil.cuh"
#include "../config/layer_config.cuh"

#include "../lbm/build/build_mom.cuh"

#include "rest.cuh"

#include <cmath>

__host__ inline int grid_plot(int x, int y)
{

    // Precisa arrumar isso aq pq tem q ver a posição dos momentos agr q mudamos a forma de acessá-los
    int blockX = x / BX;
    int blockY = y / BY;

    int localX = rest(x, BX);
    int localY = rest(y, BY);

    int blockId = blockX + GX * blockY;
    int localId = localX + BX * localY;

    return blockId * (BX * BY) + localId;
}

__device__ inline int grid_id()
{
    return (
        // Índice do bloco atual dentro da grid 2D
        (blockIdx.x + gridDim.x * blockIdx.y)

            // Quantidade total de threads por bloco
            * blockDim.x * blockDim.y

        // Índice da thread dentro do bloco 2D
        + (threadIdx.x + blockDim.x * threadIdx.y));
}

__host__ __device__ __forceinline__ int wrap_y(int y)
{
    if (y < 0)
        return y + Geometry::Ny;
    if (y >= Geometry::Ny)
        return y - Geometry::Ny;
    return y;
}

__host__ __device__ __forceinline__ int grid_index_from_xy(int x, int y)
{
    constexpr int blockSize = BX * BY;

    const int blockX = x / BX;
    const int blockY = y / BY;
    const int localX = rest(x, BX);
    const int localY = rest(y, BY);

    return (blockX + GX * blockY) * blockSize + localX + BX * localY;
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

    int Tx = blockDim.x;
    int Ty = blockDim.y;

    int tx = rest(x_from, Tx);
    int ty = rest(y_from, Ty);

    int bx = x_from / Tx;
    int by = y_from / Ty;

    int Bx = gridDim.x;
    // int By = gridDim.y;

    return ((bx + Bx * by) * Tx * Ty + tx + Tx * ty);
}

__device__ inline int from_layer_row_major(int x, int y, int i)
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

    return x_from + Geometry::Nx * y_from;
}

__device__ inline int
pop_id(int g_id, int i)
{
    return (g_id * Q + i);
}

template <int I>
__device__ __forceinline__ int layerIndex(int index)
{
    constexpr int layerSize = layer::LNx * layer::LNy;

    int localId = rest(index, layerSize);

    return I * layerSize + localId;
}

template <int I>
__host__ __device__ __forceinline__ int layer_moment_index(int x)
{
    return I * layer::LNx + x;
}

__device__ __forceinline__ int layer_slot_for_y(layer current_layer, int y)
{
    const int yw = wrap_y(y);

    if (yw == wrap_y(current_layer.yref - 1))
        return 0;
    if (yw == wrap_y(current_layer.yref))
        return 1;

    return 2;
}

template <int I>
__device__ __forceinline__ float layer_moment_read(layer current_layer, int source_row_major)
{
    const int x = source_row_major % Geometry::Nx;
    const int y = source_row_major / Geometry::Nx;
    const int slot = layer_slot_for_y(current_layer, y);

    return current_layer.buffer[slot][layer_moment_index<I>(x)];
}

template <int I>
__device__ __forceinline__ void layer_moment_write(layer current_layer, int slot, int x, float value)
{
    current_layer.buffer[slot][layer_moment_index<I>(x)] = value;
}

template <int I>
__host__ __device__ __forceinline__ int momIdx(int index)
{
    constexpr int blockSize = BX * BY;
    constexpr int blockStride = blockSize * D2Q9::num_var;

    int blockId = index / blockSize;
    int localId = rest(index, blockSize);

    return blockId * blockStride + I * blockSize + localId;
}
