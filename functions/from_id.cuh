#pragma once

#include "../presets/geometry.h"
#include "../presets/stencil.h"

__device__ inline int from_id(int x, int y, int i){

    int x_from = x - c_ix[i];
    int y_from = y - c_iy[i];

    if (x_from < 0) x_from += Nx;
    if (x_from >= Nx) x_from -= Nx;

    if (y_from < 0) y_from += Ny;
    if (y_from >= Ny) y_from -= Ny;

    return (x_from + Nx*y_from);
}