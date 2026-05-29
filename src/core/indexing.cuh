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
__host__ __device__ __forceinline__ int momIdx(int index)
{
    constexpr int blockSize = BX * BY;
    constexpr int blockStride = blockSize * D2Q9::num_var;

    int blockId = index / blockSize;
    int localId = rest(index, blockSize);

    return blockId * blockStride + I * blockSize + localId;
}
