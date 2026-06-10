#include "build_grid.cuh"

#include "../../config/geometry.h"
#include "../../config/stencil.cuh"

#include "../../core/indexing.cuh"
#include "../../core/to_u8.cuh"
#include "../../core/def_bc.cuh"

#include "build_mom.cuh"

// Got this from https://github.com/brunoyanjos/lbm/blob/ldc2D/

__global__ void init_grid(uint8_t *node, uint8_t *mask)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    const int index = grid_id(x, y);

    node[index] = to_u8(def_bc(x, y));

    uint8_t m = 0u;

#pragma unroll
    for (int i = 1; i < Q; i++)
    {

        int xn = x + c_ix[i];
        int yn = y + c_iy[i];

        if (xn < 0 || xn >= Geometry::Nx ||
            yn < 0 || yn >= Geometry::Ny)
        {
            continue;
        }
        // unsigned const int xn = x + c_ix[i];
        // unsigned const int yn = y + c_iy[i];
        // if (xn > Geometry::Nx || yn > Geometry::Ny)
        //     continue;
        m |= (1u << (i - 1));
    }

    mask[index] = m;
}

__host__ void build_grid(D2Q9 sim, Grid2D &grid)
{
    cudaMalloc((void **)&grid.mask, grid.size);
    cudaMalloc((void **)&grid.node, grid.size);

    init_grid<<<sim.N_block, sim.block>>>(grid.node, grid.mask);
}