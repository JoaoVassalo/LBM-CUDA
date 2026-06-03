#pragma once

#include "layer_manager.cuh"
#include "../../core/constexpr_for.cuh"
#include "../../config/geometry.h"
#include "../../config/layer_config.cuh"
#include "../build/build_grid.cuh"
#include "../propagation.cuh"

__host__ void seed_layer(D2Q9 sim, layer layer, Grid2D grid);
__host__ void advance_layer(D2Q9 sim, layer layer, Grid2D grid, int y);
__host__ void final_layers(D2Q9 sim, layer layer, Grid2D grid);
