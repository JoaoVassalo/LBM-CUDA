#pragma once

#include "../config/geometry.h"
#include "rest.cuh"

__host__ inline int grid_plot(int x, int y)
{
    int blockX = x / BX;
    int blockY = y / BY;

    int localX = rest(x, BX);
    int localY = rest(y, BY);

    int blockId = blockX + GX * blockY;
    int localId = localX + BX * localY;

    return blockId * (BX * BY) + localId;
}