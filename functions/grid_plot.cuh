#pragma once

#include "../presets/geometry.h"

__device__ int grid_plot(int x, int y)
{
    return x + Nx * y;
}