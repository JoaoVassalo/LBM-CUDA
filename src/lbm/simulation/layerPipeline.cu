#include "layerPipeline.cuh"

#include <cstdio>

__host__ void initLayer(moments &sim, Grid2D &grid)
{
    defLayerMem<<<blockNumberLayer, blockLayer>>>(sim);

    cudaError_t err = cudaDeviceSynchronize();
    printf("defLayerMem: %s\n", cudaGetErrorString(err));

    streamCollide<<<blockNumberLayer, blockLayer>>>(sim, grid, 0);

    err = cudaDeviceSynchronize();
    printf("streamCollide: %s\n", cudaGetErrorString(err));
}

__host__ void midLayer(moments &sim, Grid2D &grid, int y)
{
    swapLayerMem<<<blockNumberLayer, blockLayer>>>(sim, y);

    cudaError_t err = cudaDeviceSynchronize();
    printf("swapLayerMem: %s\n", cudaGetErrorString(err));

    streamCollide<<<blockNumberLayer, blockLayer>>>(sim, grid, y);

    err = cudaDeviceSynchronize();
    printf("streamCollide: %s\n", cudaGetErrorString(err));
}

__host__ void lastLayer(moments &sim, Grid2D &grid)
{
    swapLayerMem<<<blockNumberLayer, blockLayer>>>(sim, Geometry::NY - 1);

    cudaError_t err = cudaDeviceSynchronize();
    printf("swapLayerMem: %s\n", cudaGetErrorString(err));

    streamCollide<<<blockNumberLayer, blockLayer>>>(sim, grid, Geometry::NY - 1);

    err = cudaDeviceSynchronize();
    printf("streamCollide: %s\n", cudaGetErrorString(err));
}