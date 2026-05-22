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