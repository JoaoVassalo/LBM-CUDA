#pragma once

#include "../presets/geometry.h"
#include "../presets/stencil.cuh"

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

    int block_x = x_from / blockDim.x;
    int block_y = y_from / blockDim.y;

    int thread_x = x_from % blockDim.x;
    int thread_y = y_from % blockDim.y;

    int block_id = block_x + gridDim.x * block_y;
    int thread_id = thread_x + blockDim.x * thread_y;

    return block_id * blockDim.x * blockDim.y + thread_id;
}