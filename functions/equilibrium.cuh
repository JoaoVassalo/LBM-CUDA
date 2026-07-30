#pragma once

#include "../presets/stencil.cuh"

__device__ inline varUnit equilibrium(varUnit rho, varUnit ux, varUnit uy, int i)
{
    return rho * w[i] *
           (1 +
            a_s2 * (ux * c_ix[i] + uy * c_iy[i]) +
            a_s4 * 0.5f *
                ((ux * c_ix[i] + uy * c_iy[i]) * (ux * c_ix[i] + uy * c_iy[i]) -
                 (ux * ux + uy * uy) / (a_s2)));
}