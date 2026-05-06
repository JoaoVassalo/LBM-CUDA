#include "build_grid.cuh"

#include "../../config/geometry.h"

struct Grid2D
{
    uint32_t *mask;
    uint8_t *node;
};

__host__ void build_grid()
{

    Grid2D grid;

    cudaMalloc((void **)grid.mask, sizeof(uint32_t) * Nx * Ny);
    cudaMalloc((void **)grid.node, sizeof(uint8_t) * Nx * Ny);
}

__device__ void init_grid()
{

    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
}