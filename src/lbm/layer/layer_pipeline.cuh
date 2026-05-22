#pragma once

#include "build_layer.cuh"
#include "../../config/geometry.h"
#include "../build/build_grid.cuh"
#include "../propagation.cuh"

__device__ void seed_layer(D2Q9 sim, Grid2D grid);
__device__ void advance_layer(D2Q9 sim, Grid2D grid, int y);
__device__ void end_layer(D2Q9 sim, Grid2D grid);