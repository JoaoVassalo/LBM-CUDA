#pragma once

#include "../presets/geometry.h"
#include "../presets/stencil.cuh"

#include <cmath>

__device__ __forceinline__ int rest(int a, int b)
{
    return a - b * (a / b);
}

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

    return ((bx + Bx * by) * Tx * Ty + tx + Tx * ty) * 6;
}