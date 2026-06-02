#pragma once

#include "build/build_mom.cuh"
#include "../core/indexing.cuh"
#include "../config/stencil.cuh"

template <int I>
__device__ __forceinline__ float f_i(int index, layer current_layer)
{
    return layer_moment_read<MomentId::rho>(current_layer, index) * w[I] *
           (1.f +
            a_s2 * layer_moment_read<MomentId::ux>(current_layer, index) * c_ix[I] +
            a_s2 * layer_moment_read<MomentId::uy>(current_layer, index) * c_iy[I] +
            a_s4 * 0.5f * layer_moment_read<MomentId::mxx>(current_layer, index) * (c_ix[I] * c_ix[I] - inv_as2) +
            a_s4 * layer_moment_read<MomentId::mxy>(current_layer, index) * (c_ix[I] * c_iy[I]) +
            a_s4 * 0.5f * layer_moment_read<MomentId::myy>(current_layer, index) * (c_iy[I] * c_iy[I] - inv_as2));
}
