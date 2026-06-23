#pragma once

#include "../../config/momConfig.cuh"
#include "../../config/stencilConfig.cuh"

#include "../../core/indexing.cuh"
#include "../../core/constexprFor.cuh"

__global__ void defLayerMem(moments sim);

__global__ void swapLayerMem(moments sim, int y);