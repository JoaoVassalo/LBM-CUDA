#pragma once

#include "geometry.h"
#include "lbm/build/build_mom.cuh"

// Layer block size
#define LBX 64
#define LBY 1
#define LGX (Geometry::Nx / LBX)
#define LGY (Geometry::Ny / LBY)

struct layer
{
    static constexpr int LNx = Geometry::Nx;
    static const int LNy = 3;
    static const int layer_num = Geometry::Ny / LNy;

    dim3 layer_block = dim3(LBX, LBY);
    dim3 layer_Nblock = dim3(LGX, LGY);

    float *layer[3];
    static constexpr size_t layer_size = LNx * LNy * sizeof(float) * D2Q9::num_var;
    int yref = 0;
};
