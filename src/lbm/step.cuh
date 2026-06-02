#pragma once

#include "collision.cuh"
#include "propagation.cuh"

#include "../config/geometry.h"

#include "build/build_grid.cuh"
#include "layer/layer_manager.cuh"
#include "layer/layer_pipeline.cuh"
#include "build/build_mom.cuh"

__host__ void step(D2Q9 sim, layer layer, Grid2D grid);
