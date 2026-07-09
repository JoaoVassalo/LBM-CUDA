#pragma once

#include "../../config/gridConfig.cuh"
#include "../../config/momConfig.cuh"

#include "layerManager.cuh"
#include "streamCollide.cuh"

__host__ void initLayer(moments &sim, Grid2D &grid);

__host__ void midLayer(moments &sim, Grid2D &grid, int y);

__host__ void lastLayer(moments &sim, Grid2D &grid);