#pragma once

#include "../presets/config.h"
#include "../presets/geometry.h"

template <int I>
__host__ __device__ __forceinline__ int momIdx(int index)
{
    return index + I;
}

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

    size_t num_var = 6; // Number of moments in the stencil

    dim3 block = dim3(BX, BY);
    dim3 N_block = dim3(GX, GY);

    size_t size = Nx * Ny * sizeof(float) * num_var;
};