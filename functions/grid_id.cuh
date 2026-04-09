#pragma once

#include "../presets/geometry.h"
#include "../cuda_config/var.cuh"

__device__ inline int grid_id()
{
    return ((blockIdx.x + gridDim.x * blockIdx.y) * blockDim.x * blockDim.y +
            (threadIdx.x + blockDim.x * threadIdx.y)) *
           6;
}