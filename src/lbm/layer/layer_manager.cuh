#pragma once

#include "../../core/constexpr_for.cuh"
#include "../../core/indexing.cuh"
#include "../../config/layer_config.cuh"
#include "../build/build_mom.cuh"

__global__ void init_layers(D2Q9 sim, layer current_layer);

__global__ void swap_layers(D2Q9 sim, layer current_layer, int y);
