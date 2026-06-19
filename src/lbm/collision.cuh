#pragma once

#include "build/build_mom.cuh"

__global__ void collide_layer_at(D2Q9 sim, int y);