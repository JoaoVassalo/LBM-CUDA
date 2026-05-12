#include "build_grid.cuh"

#include "../../config/geometry.h"
#include "../../config/stencil.cuh"

#include "../../core/grid_id.cuh"
#include "../../core/to_u8.cuh"

#include "build_mom.cuh"

// Got this from https://github.com/brunoyanjos/lbm/blob/ldc2D/

struct Grid2D
{
    uint8_t *mask;
    uint8_t *node;
};

__device__ void build_directions(uint8_t node, uint8_t mask)
{

    int *out;

    int out_cont = 0;

    out[out_cont++] = 0;

    for (int i = 1; i < Q; i++)
    {
        if (mask & (1u << i - 1))
            out[out_cont++] = i;
    }
}

__device__ void init_grid(uint8_t *node, uint8_t *mask)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    const int index = grid_id();

    const bool top = (y == Ny - 1);
    const bool left = (x == Nx - 1);
    const bool down = (y == 0);
    const bool right = (x == 0);

    const int bc_count = (int)top + (int)left + (int)down + (int)right;

    uint8_t wid = to_u8(domainTags::Fluid);

    if (bc_count > 0)
    {
        wid = to_u8(domainTags::Boundary);
    }

    node[index] = wid;

    uint8_t m = 0u;

#pragma unroll
    for (int i = 1; i < Q; i++)
    {
        unsigned const int xn = x + c_ix[i];
        unsigned const int yn = y + c_iy[i];
        if (xn > Nx || yn > Ny)
            continue;
        m |= (1u << (i - 1));
    }

    mask[index] = m;
}

__host__ void build_grid(D2Q9 sim)
{
    Grid2D grid;

    cudaMalloc((void **)grid.mask, sizeof(uint8_t) * Nx * Ny);
    cudaMalloc((void **)grid.node, sizeof(uint8_t) * Nx * Ny);

    init_grid<<<sim.N_block, sim.block>>>(grid.node, grid.mask);
}