#include "layerPipeline.cuh"

__host__ void initLayer(moments sim, Grid2D grid)
{
    defLayerMem<<<blockNumberLayer, blockLayer>>>(sim);
    streamCollide<<<blockNumberLayer, blockLayer>>>(sim, grid, 0);
}

__host__ void midLayer(moments sim, Grid2D grid, int y)
{
    swapLayerMem<<<blockNumberLayer, blockLayer>>>(sim, y);
    streamCollide<<<blockNumberLayer, blockLayer>>>(sim, grid, y);
}

__host__ void lastLayer(moments sim, Grid2D grid)
{
    swapLayerMem<<<blockNumberLayer, blockLayer>>>(sim, Geometry::NY - 1);
    streamCollide<<<blockNumberLayer, blockLayer>>>(sim, grid, Geometry::NY - 1);
}