#pragma once

#include "../../config/stencilConfig.cuh"
#include "../../config/momConfig.cuh"

#include "../../core/indexing.cuh"

__device__ __forceinline__ float fi(int x, int i, float **layer)
{
    int x_from = x - (int)D2Q9::cx(i);

    if (x_from < 0)
        x_from += Geometry::NX;
    if (x_from >= Geometry::NX)
        x_from -= Geometry::NX;

    int y = 1 - (int)D2Q9::cy(i);

    float rho = layer[y][layerIdx<momId::rho>(x_from)];
    float ux = layer[y][layerIdx<momId::ux>(x_from)];
    float uy = layer[y][layerIdx<momId::uy>(x_from)];
    float mxx = layer[y][layerIdx<momId::mxx>(x_from)];
    float mxy = layer[y][layerIdx<momId::mxy>(x_from)];
    float myy = layer[y][layerIdx<momId::myy>(x_from)];

    return rho * D2Q9::w(i) *
           (1.f +
            D2Q9::a_s2 * ux * D2Q9::cx(i) +
            D2Q9::a_s2 * uy * D2Q9::cy(i) +
            D2Q9::a_s4 * 0.5f * mxx * (D2Q9::cx(i) * D2Q9::cx(i) - D2Q9::inv_as2) +
            D2Q9::a_s4 * mxy * (D2Q9::cx(i) * D2Q9::cy(i)) +
            D2Q9::a_s4 * 0.5f * myy * (D2Q9::cy(i) * D2Q9::cy(i) - D2Q9::inv_as2));
}