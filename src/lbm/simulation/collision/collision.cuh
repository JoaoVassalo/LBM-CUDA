#pragma once

#include "../../../config/momConfig.cuh"
#include "../../../config/simulationConfig.cuh"
#include "../../../config/stencilConfig.cuh"

#include "../../../core/indexing.cuh"

__device__ void collide(moments sim, int y);