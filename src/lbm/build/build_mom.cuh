#pragma once

#include "../../config/geometry.h"
#include "../core/rest.cuh"

enum MomentId
{
    rho = 0,
    ux,
    uy,
    mxx,
    myy,
    mxy
};

struct D2Q9
{
    float *mom;
    float *layer[3];

    static constexpr int num_var = 6; // Number of moments in the stencil

    dim3 block = dim3(BX, BY);
    dim3 N_block = dim3(GX, GY);

    dim3 layer_block = dim3(LBX, LBY);
    dim3 layer_Nblock = dim3(LGX, LGY);

    size_t size = Geometry::Nx * Geometry::Ny * sizeof(float) * num_var;
    size_t layer_size = layer::LNx * layer::LNy * sizeof(float) * num_var;
};

template <int I>
__host__ __device__ __forceinline__ int momIdx(int index)
{
    constexpr int blockSize = BX * BY;
    constexpr int blockStride = blockSize * D2Q9::num_var;

    int blockId = index / blockSize;
    int localId = rest(index, blockSize);

    return blockId * blockStride + I * blockSize + localId;
}
