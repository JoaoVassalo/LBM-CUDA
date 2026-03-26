#pragma once

#include "../presets/geometry.h"

__host__ int grid_plot(int x, int y)
{
    return x + Nx * y;
}