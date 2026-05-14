#include "build_grid.cuh"

#include "../../config/geometry.h"
#include "../../config/stencil.cuh"

#include "../../core/grid_id.cuh"
#include "../../core/to_u8.cuh"
#include "../../core/def_bc.cuh"

#include "build_mom.cuh"

// Got this from https://github.com/brunoyanjos/lbm/blob/ldc2D/

__host__ __device_builtin__ __forceinline__ int income(int i)
{
    switch (i)
    {
    case 0:
        return 0;
        break;
    case 1:
        return 3;
        break;
    case 2:
        return 4;
        break;
    case 3:
        return 1;
        break;
    case 4:
        return 2;
        break;
    case 5:
        return 7;
        break;
    case 6:
        return 8;
        break;
    case 7:
        return 5;
        break;
    default: // case 8
        return 6;
        break;
    }
}

__global__ void init_grid(uint8_t *node, uint8_t *mask)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    const int index = grid_id();

    node[index] = to_u8(def_bc(x, y));

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

__host__ void build_grid(D2Q9 sim, Grid2D grid)
{
    cudaMalloc((void **)grid.mask, grid.size);
    cudaMalloc((void **)grid.node, grid.size);

    init_grid<<<sim.N_block, sim.block>>>(grid.node, grid.mask);
}