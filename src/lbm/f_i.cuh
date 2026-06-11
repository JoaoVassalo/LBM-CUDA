#pragma once

#include "build/build_mom.cuh"
#include "../core/indexing.cuh"
#include "../config/stencil.cuh"

template <int I>
__device__ __forceinline__ float f_i(int x, layer current_layer)
{
    int x_from = x - (int)c_ix[I];

    if (x_from < 0)
        x_from += Geometry::Nx;
    if (x_from >= Geometry::Nx)
        x_from -= Geometry::Nx;

    int y = 1 - (int)c_iy[I];

    return current_layer.buffer[y][layerIdx<MomentId::rho>(x_from)] * w[I] *
           (1.f +
            a_s2 * current_layer.buffer[y][layerIdx<MomentId::ux>(x_from)] * c_ix[I] +
            a_s2 * current_layer.buffer[y][layerIdx<MomentId::uy>(x_from)] * c_iy[I] +
            a_s4 * 0.5f * current_layer.buffer[y][layerIdx<MomentId::mxx>(x_from)] * (c_ix[I] * c_ix[I] - inv_as2) +
            a_s4 * current_layer.buffer[y][layerIdx<MomentId::mxy>(x_from)] * (c_ix[I] * c_iy[I]) +
            a_s4 * 0.5f * current_layer.buffer[y][layerIdx<MomentId::myy>(x_from)] * (c_iy[I] * c_iy[I] - inv_as2));
}