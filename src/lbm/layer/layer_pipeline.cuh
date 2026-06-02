#pragma once

#include "layer_manager.cuh"
#include "../../core/constexpr_for.cuh"
#include "../../config/geometry.h"
#include "../../config/layer_config.cuh"
#include "../build/build_grid.cuh"
#include "../propagation.cuh"

__global__ void seed_layer(D2Q9 sim, layer layer, Grid2D grid);
__global__ void advance_layer(D2Q9 sim, layer layer, Grid2D grid, int y);
__global__ void end_layer(D2Q9 sim, layer layer, Grid2D grid);
__global__ void final_layers(D2Q9 sim, layer layer, Grid2D grid);
