#pragma once

#include "../../config/stencilConfig.cuh"
#include "../../config/momConfig.cuh"
#include "../core/indexing.cuh"

__global__ void initDomain(float *mom);