#include "layerPipeline.cuh"

#include <cstdio>

__global__ void streamCollideStep(moments sim, Grid2D grid)
{
    const int x = threadIdx.x;
    const bool validX = x < Geometry::NX;

    if (validX)
    {
        for (int slot = 0; slot < Geometry::LNY; ++slot)
        {
            const int sourceY = slot == 0 ? 0 : slot - 1;
            for (int moment = 0; moment < D2Q9::momNum; ++moment)
            {
                const int layerIndex = x * D2Q9::momNum + moment;
                const int momIndex = (sourceY * Geometry::NX + x) * D2Q9::momNum + moment;
                sim.layer[slot][layerIndex] = sim.mom[momIndex];
            }
        }
    }

    __syncthreads();

    if (validX)
    {
        applyBoundary(sim, grid, x, 0);
        collide(sim.mom, gridId(x, 0));
    }

    __syncthreads();

    for (int y = 1; y < Geometry::NY; ++y)
    {
        if (validX)
        {
            const int sourceY = y + 1 < Geometry::NY ? y + 1 : Geometry::NY - 1;

            for (int moment = 0; moment < D2Q9::momNum; ++moment)
            {
                const int layerIndex = x * D2Q9::momNum + moment;
                const varType mid = sim.layer[1][layerIndex];
                const varType top = sim.layer[2][layerIndex];

                sim.layer[0][layerIndex] = mid;
                sim.layer[1][layerIndex] = top;
                sim.layer[2][layerIndex] = sim.mom[(sourceY * Geometry::NX + x) * D2Q9::momNum + moment];
            }
        }

        __syncthreads();

        if (validX)
        {
            applyBoundary(sim, grid, x, y);
            collide(sim.mom, gridId(x, y));
        }

        __syncthreads();
    }
}

__host__ void initLayer(moments &sim, Grid2D &grid)
{
    defLayerMem<<<blockNumberLayer, blockLayer>>>(sim);

    streamCollide<<<blockNumberLayer, blockLayer>>>(sim, grid, 0);
}

__host__ void midLayer(moments &sim, Grid2D &grid, int y)
{
    swapLayerMem<<<blockNumberLayer, blockLayer>>>(sim, y + 1);

    streamCollide<<<blockNumberLayer, blockLayer>>>(sim, grid, y);
}

__host__ void lastLayer(moments &sim, Grid2D &grid)
{
    swapLayerMem<<<blockNumberLayer, blockLayer>>>(sim, Geometry::NY - 1);

    streamCollide<<<blockNumberLayer, blockLayer>>>(sim, grid, Geometry::NY - 1);
}
