#include "initialization.cuh"

__host__ void initialization(moments &sim, Grid2D &grid)
{
    constexpr size_t momByteSize = D2Q9::momByteSize;
    constexpr size_t layerByteSize = D2Q9::layerByteSize;

    cudaMalloc((void **)&sim.mom, momByteSize);

    initDomain<<<blockNumber, block>>>(sim.mom);
    cudaDeviceSynchronize();

    for (size_t i = 0; i < Geometry::LNY; i++)
    {
        cudaMalloc((void **)&sim.layer[i], layerByteSize);
    }

    cudaMalloc((void **)&grid.mask, grid.gridByteSize);
    cudaMalloc((void **)&grid.node, grid.gridByteSize);

    initGrid<<<blockNumber, block>>>(grid.mask, grid.node);

    cudaMalloc((void **)&sim.vtk, D2Q9::outputByteSize);
    cudaHostAlloc((void **)&sim.mom_host, D2Q9::outputByteSize, cudaHostAllocDefault);

    writeOutput(sim, grid, (size_t)0, momByteSize, output::vtkPath);
}
