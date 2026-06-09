#pragma once

#include "geometry.h"
#include "lbm/build/build_mom.cuh"

// Layer block size
#define LBX 64
#define LBY layer::LNy
#define LGX (Geometry::Nx / LBX)
#define LGY (layer::LNy - 2)

struct layer
{
    static constexpr int LNx = Geometry::Nx;
    static constexpr int LNy = 3;
    static constexpr int layer_num = Geometry::Ny;

    dim3 layer_block = dim3(LBX, LBY);
    dim3 layer_Nblock = dim3(LGX, LGY);

    float *buffer[3];
    static constexpr size_t buffer_size = LNx * D2Q9::num_var;
    static constexpr size_t buffer_bytesize = buffer_size * sizeof(float);
    static constexpr size_t layer_size = buffer_size * LNy;
    static constexpr size_t layer_bytesize = buffer_bytesize * LNy;
    int yref = 0;
};
