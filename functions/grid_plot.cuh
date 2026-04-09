#pragma once

#include "../presets/geometry.h"
#include "../presets/config.h"

__host__ inline int grid_plot(int x, int y)
{
    return (((x / BX) + GX * (y / BY)) * (BX * BY) + ((x - BX * (x / BX)) + BX * (y - BY * (y / BY)))) * 6;
}