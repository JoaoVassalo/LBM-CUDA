#include "initialization.cuh"

__host__ __forceinline__ void initialization(moments sim, Grid2D grid)
{
    constexpr size_t momByteSize = D2Q9::momByteSize;
    constexpr size_t layerByteSize = D2Q9::layerByteSize;

    cudaMalloc((void **)&sim.mom, momByteSize);

    initDomain<<<block, blockNumber>>>(sim.mom);

    cudaError_t err = cudaGetLastError();
    printf("Kernel launch error: %s\n", cudaGetErrorString(err));

    cudaDeviceSynchronize();

    for (size_t i = 0; i < Geometry::LNY; i++)
    {
        cudaMalloc((void **)&sim.layer[i], layerByteSize);
    }

    cudaMalloc((void **)&grid.mask, grid.gridByteSize);
    cudaMalloc((void **)&grid.node, grid.gridByteSize);

    initGrid<<<block, blockNumber>>>(grid.mask, grid.node);

    cudaError_t err = cudaGetLastError();
    printf("Kernel launch error: %s\n", cudaGetErrorString(err));

    sim.mom_host = (float *)malloc(momByteSize);

    writeOutput(sim, grid, (size_t)0, momByteSize, output::vtkPath);
}