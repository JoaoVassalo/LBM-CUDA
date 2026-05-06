#include "build_grid.cuh"

#include "../../config/geometry.h"

#include "../../core/grid_id.cuh"

#include "build_mom.cuh"

// Got this from https://github.com/brunoyanjos/lbm/blob/ldc2D/

struct Grid2D
{
    uint32_t *mask;
    uint8_t *node;
};

__host__ void build_grid(D2Q9 sim)
{

    Grid2D grid;

    cudaMalloc((void **)grid.mask, sizeof(uint32_t) * Nx * Ny);
    cudaMalloc((void **)grid.node, sizeof(uint8_t) * Nx * Ny);

    init_grid<<<sim.N_block, sim.block>>>();
}

__device__ void init_grid()
{

    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    const int index = grid_id();

    const bool top = (y == Ny - 1);
    const bool left = (x == Nx - 1);
    const bool down = (y == 0);
    const bool right = (x == 0);

    const int bc_count = (int)top + (int)left + (int)down + (int)right;

    // uint8_t wid = FUNCAO QUE TRANSFORMA INT EM u8 (FLUIDO)

    if (bc_count > 0)
    {
        return;

        // wid = FUNCAO QUE TRANSFORMA INT EM u8 (CONTORNO)
    }

    node[index] = wid;
}