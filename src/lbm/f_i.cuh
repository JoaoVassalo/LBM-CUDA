#pragma once

#include "build/build_mom.cuh"
#include "../config/stencil.cuh"

__device__ inline float f_i(int index, int i, float **mom)
{

    return mom[i][momIdx<MomentId::rho>(index)] * w[i] * // rho
           (1.f +
            a_s2 * mom[i][momIdx<MomentId::ux>(index)] * c_ix[i] +                               // ux
            a_s2 * mom[i][momIdx<MomentId::uy>(index)] * c_iy[i] +                               // uy
            a_s4 * 0.5f * mom[i][momIdx<MomentId::mxx>(index)] * (c_ix[i] * c_ix[i] - inv_as2) + // mxx
            a_s4 * mom[i][momIdx<MomentId::mxy>(index)] * (c_ix[i] * c_iy[i]) +                  // mxy
            a_s4 * 0.5f * mom[i][momIdx<MomentId::myy>(index)] * (c_iy[i] * c_iy[i] - inv_as2)); // myy
}