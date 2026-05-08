#include "collision.cuh"
#include "propagation.cuh"
#include "grid_id.cuh"

#include "../config/geometry.h"

__global__ void step(float *mom)
{

    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    if (x >= Nx || y >= Ny)
        return;

    int index = grid_id();

    propagation(mom);

    collision(index, mom);
}