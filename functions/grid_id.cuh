#pragma once

#include "../presets/geometry.h"

__device__ inline int grid_id(int x, int y){
    return (x + Nx*y);
}