#pragma once

#include <cstdint>

#include "../../config/mom_config.cuh"
#include "../../config/layer_config.cuh"
#include "../build/build_grid.cuh"

#include "../../core/indexing.cuh"
#include "../../core/to_u8.cuh"
#include "../../core/def_bc.cuh"
#include "../../core/constexpr_for.cuh"

#include "../f_i.cuh"

__device__ void center(int x, int y,
                       D2Q9 sim, layer layer);

__device__ void boundary(uint8_t mask_in, uint8_t node_in, int x, int y,
                         D2Q9 sim, layer layer);
