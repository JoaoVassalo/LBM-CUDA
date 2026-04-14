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

    mom[momIdx<MomentId::rho>(index)] = 1.f; // rho
    mom[momIdx<MomentId::ux>(index)] = 0.f;  // ux
    mom[momIdx<MomentId::uy>(index)] = 0.f;  // uy
    mom[momIdx<MomentId::mxx>(index)] = 0.f; // mxx
    mom[momIdx<MomentId::mxy>(index)] = 0.f; // mxy
    mom[momIdx<MomentId::myy>(index)] = 0.f; // myy
}