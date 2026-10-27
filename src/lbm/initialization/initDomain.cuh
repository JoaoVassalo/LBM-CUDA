#pragma once

#include "../../config/simulationConfig.cuh"
#include "../../config/momConfig.cuh"
#include "../../core/indexing.cuh"

__global__ void initDomain(varType *mom);
