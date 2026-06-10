#pragma once

#include "build/build_mom.cuh"
#include "../core/indexing.cuh"
#include "../config/stencil.cuh"

template <int I>
__device__ __forceinline__ float f_i(int index, layer current_layer)
{
    int x = rest(index, Geometry::Nx);
    int y = index / Geometry::Nx;

    return current_layer.buffer[y][layerIdx<MomentId::rho>(x)] * w[I] *
           (1.f +
            a_s2 * current_layer.buffer[y][layerIdx<MomentId::ux>(x)] * c_ix[I] +
            a_s2 * current_layer.buffer[y][layerIdx<MomentId::uy>(x)] * c_iy[I] +
            a_s4 * 0.5f * current_layer.buffer[y][layerIdx<MomentId::mxx>(x)] * (c_ix[I] * c_ix[I] - inv_as2) +
            a_s4 * current_layer.buffer[y][layerIdx<MomentId::mxy>(x)] * (c_ix[I] * c_iy[I]) +
            a_s4 * 0.5f * current_layer.buffer[y][layerIdx<MomentId::myy>(x)] * (c_iy[I] * c_iy[I] - inv_as2));
}
