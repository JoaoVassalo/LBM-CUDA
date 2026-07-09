#include "streamCollide.cuh"

__global__ void streamCollide(moments sim, Grid2D grid, int y)
{
    int x = threadIdx.x + blockDim.x * blockIdx.x;

    if (x >= Geometry::NX)
        return;

    applyBoundary(sim, grid, x, y);
    collide(sim, x, y);
    saveMom(sim, x, y);
}