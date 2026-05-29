#pragma once

#include "../../config/geometry.h"
#include "../../core/rest.cuh"

enum MomentId
{
    rho = 0,
    ux,
    uy,
    mxx,
    myy,
    mxy,
    Count
};

// Block size
#define BX 32
#define BY 16
#define GX (Geometry::Nx / BX)
#define GY (Geometry::Ny / BY)

struct D2Q9
{
    float *mom;

    static constexpr int num_var = 6; // Number of moments in the stencil

    dim3 block = dim3(BX, BY);
    dim3 N_block = dim3(GX, GY);

    static constexpr size_t size = Geometry::Nx * Geometry::Ny * sizeof(float) * num_var;
};
