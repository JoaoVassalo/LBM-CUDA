#pragma once

#include "../../presets/stencil.cuh"

__device__ inline float f_i(int index, int i, float *rho, float *ux, float *uy, float *mxx, float *mxy, float *myy)
{
    return rho[index] * w[i] *
           (1 +
            a_s2 * ux[index] * c_ix[i] +
            a_s2 * uy[index] * c_iy[i] +
            a_s4 * 0.5f * mxx[index] * (c_ix[i] * c_ix[i] - inv_as2) +
            a_s4 * mxy[index] * (c_ix[i] * c_iy[i]) +
            a_s4 * 0.5f * myy[index] * (c_iy[i] * c_iy[i] - inv_as2));
}