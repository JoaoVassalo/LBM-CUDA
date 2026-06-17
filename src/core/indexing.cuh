#pragma once

#include "../config/geometry.h"
#include "../config/stencil.cuh"

#include "rest.cuh"

#include <cmath>

__device__ inline int from_id(int x, int y, int i)
{
    int x_from = x - (int)c_ix[i];
    int y_from = y - (int)c_iy[i];

    if (x_from < 0)
        x_from += Nx;
    if (x_from >= Nx)
        x_from -= Nx;

    if (y_from < 0)
        y_from += Ny;
    if (y_from >= Ny)
        y_from -= Ny;

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

__device__ inline int grid_id()
{
    return ((blockIdx.x + gridDim.x * blockIdx.y) * blockDim.x * blockDim.y +
            (threadIdx.x + blockDim.x * threadIdx.y));
}

__host__ inline int grid_plot(int x, int y)
{
    int blockX = x / BX;
    int blockY = y / BY;

    int localX = rest(x, BX);
    int localY = rest(y, BY);

    int blockId = blockX + GX * blockY;
    int localId = localX + BX * localY;

    return blockId * (BX * BY) + localId;
}

__device__ inline int
pop_id(int g_id, int i)
{
    return (g_id * Q + i);
}