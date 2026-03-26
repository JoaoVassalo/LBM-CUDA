#pragma once

#include "../presets/geometry.h"

__device__ inline int grid_id(int x, int y)
{
    //    return (blockIdx.x + gridDim.x * blockIdx.y) * blockDim.x * blockDim.y + (threadIdx.x + blockDim.x * threadIdx.y);
    return (x + Nx * y);
}