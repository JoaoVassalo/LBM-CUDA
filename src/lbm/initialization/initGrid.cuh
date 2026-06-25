#pragma once

#include "../../config/stencilConfig.cuh"
#include "../../config/gridConfig.cuh"

#include "../../core/indexing.cuh"
#include "../../core/to_u8.cuh"

__global__ void initGrid(uint8_t *mask, uint8_t *node);