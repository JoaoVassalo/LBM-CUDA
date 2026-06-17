#pragma once

#include "collision.cuh"
#include "propagation.cuh"
#include "../core/indexing.cuh"

#include "../config/geometry.h"

__global__ void step(float **mom, float *layer, uint8_t *mask, uint8_t *node, int step_i);