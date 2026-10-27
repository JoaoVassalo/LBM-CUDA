#pragma once

#include "../../../config/momConfig.cuh"
#include "../../../config/simulationConfig.cuh"

#include "../../../core/indexing.cuh"

__device__ void collide(varType *mom, int index);
