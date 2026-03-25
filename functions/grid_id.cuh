#pragma once

#include "../presets/geometry.h"

__host__ __device__ inline int grid_id(int x, int y)
{
    return (x + Nx * y);
}