#pragma once

#include "../fi.cuh"

#include "../../../config/momConfig.cuh"
#include "../../../config/simulationConfig.cuh"
#include "../../../config/gridConfig.cuh"

#include "../../../core/to_u8.cuh"
#include "../../../core/constexprFor.cuh"

__device__ void applyBoundary(moments sim, Grid2D grid, int x, int y);
