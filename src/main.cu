#include "config/stencilConfig.cuh"
#include "config/gridConfig.cuh"
#include "lbm/initialization/initDomain.cuh"
#include "lbm/initialization/initGrid.cuh"

int main()
{
    D2Q9 sim;

    cudaMalloc((void **)&sim.mom, sim.momByteSize);

    initDomain<<<block, blockNumber>>>(sim.mom);

    Grid2D grid;

    initGrid<<<block, blockNumber>>>(grid.mask, grid.node);
}