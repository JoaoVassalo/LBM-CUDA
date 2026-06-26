#pragma once

#include "../../config/gridConfig.cuh"
#include "../../config/momConfig.cuh"

#include "streaming/boundary.cuh"
#include "collision/collision.cuh"

__global__ void streamCollide(moments sim, Grid2D grid, int y);