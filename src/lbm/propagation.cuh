#pragma once

#include <iostream>

#include "../config/geometry.h"
#include "../config/stencil.cuh"
#include "../config/physics.h"
#include "../config/layer_config.cuh"

#include "build/build_grid.cuh"

#include "../core/to_u8.cuh"

#include "boundary/boundary.cuh"

#include "pop_id.cuh"
#include "grid_id.cuh"
#include "from_id.cuh"

__global__ void propagation(D2Q9 sim, layer layer, Grid2D grid);