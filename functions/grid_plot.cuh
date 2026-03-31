#pragma once

#include "../presets/geometry.h"
#include "../presets/config.h"

__host__ inline int grid_plot(int x, int y)
{
    int block_x = x / BX;
    int block_y = y / BY;

    int thread_x = x % BX;
    int thread_y = y % BY;

    int block_id = block_x + GX * block_y;
    int thread_id = thread_x + BX * thread_y;

    return block_id * (BX * BY) + thread_id;
}