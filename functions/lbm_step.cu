#include "collision.cuh"
#include "propagation.cuh"
#include "grid_id.cuh"

#include "../presets/geometry.h"

__global__ void lbm_step(varUnit *mom_in, varUnit *mom_out)
{

    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    if (x >= Nx || y >= Ny)
        return;

    int index = grid_id();

    propagation(mom_in, mom_out);

    collision(index, mom_out);
}