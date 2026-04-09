#pragma once

#include "../../presets/stencil.cuh"

__device__ inline float f_i(int index, int i, float *mom)
{
    return mom[index] * w[i] * // rho
           (1.f +
            a_s2 * mom[index + 1] * c_ix[i] +                              // ux
            a_s2 * mom[index + 2] * c_iy[i] +                              // uy
            a_s4 * 0.5f * mom[index + 3] * (c_ix[i] * c_ix[i] - inv_as2) + // mxx
            a_s4 * mom[index + 5] * (c_ix[i] * c_iy[i]) +                  // mxy
            a_s4 * 0.5f * mom[index + 4] * (c_iy[i] * c_iy[i] - inv_as2)); // myy
}