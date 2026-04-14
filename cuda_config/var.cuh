#pragma once

#include "../presets/config.h"
#include "../presets/geometry.h"

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
    float *momA;
    float *momB;

    static constexpr int num_var = 6; // Number of moments in the stencil

    dim3 block = dim3(BX, BY);
    dim3 N_block = dim3(GX, GY);

    size_t size = Nx * Ny * sizeof(float) * num_var;
};

template <int I>
__host__ __device__ __forceinline__ int momIdx(int index)
{
    int blockId = index - BX * BY * (index / (BX * BY));
    return index + BX * BY * (blockId * D2Q9::num_var + I);
}
