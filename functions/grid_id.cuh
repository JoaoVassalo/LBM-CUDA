#pragma once

#include "../presets/geometry.h"

__device__ inline int grid_id()
{
    return ((blockIdx.x + gridDim.x * blockIdx.y) * blockDim.x * blockDim.y + (threadIdx.x + blockDim.x * threadIdx.y));
}