#include "../presets/stencil.cuh"
#include "../presets/geometry.h"
#include "../presets/config.h"

#include "equilibrium.cuh"
#include "grid_id.cuh"
#include "pop_id.cuh"

__global__ void initDomain(float *mom)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    if (x >= Nx || y >= Ny)
        return;

    int index = grid_id();

    mom[index] = 1.f;     // rho
    mom[index + 1] = 0.f; // ux
    mom[index + 2] = 0.f; // uy
    mom[index + 3] = 0.f; // mxx
    mom[index + 4] = 0.f; // mxy
    mom[index + 5] = 0.f; // myy
}