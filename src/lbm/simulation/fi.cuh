#pragma once

#include "../../config/simulationConfig.cuh"
#include "../../config/momConfig.cuh"

#include "../../core/indexing.cuh"

__device__ __forceinline__ varType fi(int x, int i, varType **layer)
{
    int x_from = x - (int)D2Q9::cx(i);

    if (x_from < 0)
        x_from += Geometry::NX;
    if (x_from >= Geometry::NX)
        x_from -= Geometry::NX;

    int y = 1 - (int)D2Q9::cy(i);

    varType rho = layer[y][layerIdx<momId::rho>(x_from)];
    varType ux = layer[y][layerIdx<momId::ux>(x_from)];
    varType uy = layer[y][layerIdx<momId::uy>(x_from)];
    varType mxx = layer[y][layerIdx<momId::mxx>(x_from)];
    varType mxy = layer[y][layerIdx<momId::mxy>(x_from)];
    varType myy = layer[y][layerIdx<momId::myy>(x_from)];

    return rho * D2Q9::w(i) *
           (1. +
            D2Q9::a_s2 * ux * D2Q9::cx(i) +
            D2Q9::a_s2 * uy * D2Q9::cy(i) +
            D2Q9::a_s4 * 0.5 * mxx * (D2Q9::cx(i) * D2Q9::cx(i) - D2Q9::inv_as2) +
            D2Q9::a_s4 * mxy * (D2Q9::cx(i) * D2Q9::cy(i)) +
            D2Q9::a_s4 * 0.5 * myy * (D2Q9::cy(i) * D2Q9::cy(i) - D2Q9::inv_as2));
}
