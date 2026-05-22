#pragma once

#include "../../core/grid_id.cuh"
#include "../../core/rest.cuh"
#include "../../lbm/build/build_mom.cuh"
#include "../../config/geometry.h"

enum layer_type
{
    bot = 0,
    mid,
    top
};

__global__ void first_layer(D2Q9 sim);

__device__ void other_layers(D2Q9 sim, int y, int i);

__global__ void last_layer(D2Q9 sim);