#pragma once

#include "../config/geometry.h"
#include "../config/physics.h"

#include "../core/indexing.cuh"

#include "build/build_mom.cuh"

__device__ void collision(float *mom);