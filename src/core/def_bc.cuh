#pragma once

#include "../config/geometry.h"
#include "../config/mom_config.cuh"

__host__ __device__ __forceinline__ Boundary def_bc(int x, int y)
{
    if (x == 0 && y == 0)
    { // Sudoeste
        return Boundary::Southwest;
    }
    else if (x == Nx - 1 && y == 0)
    { // Sudeste
        return Boundary::Southeast;
    }
    else if (x == 0 && y == Ny - 1)
    { // Noroeste
        return Boundary::Northwest;
    }
    else if (x == Nx - 1 && y == Ny - 1)
    { // Nordeste
        return Boundary::Northeast;
    }
    else if (y == 0)
    { // Sul
        return Boundary::South;
    }
    else if (y == Ny - 1)
    { // Norte
        return Boundary::North;
    }
    else if (x == 0)
    { // Oeste
        return Boundary::West;
    }
    else if (x == Nx - 1)
    { // Leste
        return Boundary::East;
    }
    else
    { // Centro
        return Boundary::Center;
    }
}