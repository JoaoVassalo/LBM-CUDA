#pragma once

#include "../../config/stencil.cuh"
#include "../../config/geometry.h"
#include "../../config/mom_config.cuh"

#include "../build/build_mom.cuh"

#include "../../core/indexing.cuh"

__global__ void initDomain(float *mom);