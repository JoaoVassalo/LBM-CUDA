#include "../presets/stencil.cuh"
#include "../presets/geometry.h"
#include "../presets/config.h"

#include "equilibrium.cuh"
#include "grid_id.cuh"
#include "pop_id.cuh"

__global__ void initDomain(float *rho, float *ux, float *uy, float *mxx, float *mxy, float *myy)
{

    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    if (x >= Nx && y >= Ny)
        return;

    int index = grid_id();

    rho[index] = 1.f;
    ux[index] = 0.f;
    uy[index] = 0.f;
    mxx[index] = 0.f;
    mxy[index] = 0.f;
    myy[index] = 0.f;

    //     float f_i[Nx * Ny * Q] = {0};

    // #pragma unroll
    //     for (int i = 0; i < Q; i++)
    //     {
    //         f_i[Nx * Ny * Q] = equilibrium(rho[index], ux[index], uy[index], i);

    //         mxx[index] += f_i[pop_id(index, i)] * (c_ix[i] * c_ix[i] - inv_as2);
    //         mxy[index] += f_i[pop_id(index, i)] * (c_ix[i] * c_iy[i]);
    //         mxx[index] += f_i[pop_id(index, i)] * (c_iy[i] * c_iy[i] - inv_as2);
    //     }

    //     mxx[index] *= 1 / rho[index];
    //     mxy[index] *= 1 / rho[index];
    //     myy[index] *= 1 / rho[index];
}