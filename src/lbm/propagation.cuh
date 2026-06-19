#pragma once

#include <iostream>

#include "../config/geometry.h"
#include "../config/stencil.cuh"
#include "../config/physics.h"
#include "../config/mom_config.cuh"

#include "build/build_grid.cuh"

#include "../core/to_u8.cuh"
#include "../core/indexing.cuh"
#include "../core/def_bc.cuh"

#include "f_i.cuh"

__device__ void propagation(int x, int y, float **mom, float *layer, uint8_t *mask, uint8_t *node, int step_i);