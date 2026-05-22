#pragma once

#include "collision.cuh"
#include "propagation.cuh"
#include "grid_id.cuh"

#include "../config/geometry.h"

#include "build/build_grid.cuh"
#include "layer/build_layer.cuh"
#include "layer/layer_pipeline.cuh"
#include "build/build_mom.cuh"

__host__ void step(D2Q9 sim, Grid2D grid);