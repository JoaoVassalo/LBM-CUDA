#include "layerPipeline.cuh"

__host__ void initLayer(moments sim, Grid2D grid)
{
    defLayerMem<<<blockNumberLayer, blockLayer>>>(sim);
    streaming_at_y<<<blockNumberLayer, blockLayer>>>(sim, grid, 0);
}

__host__ void midLayer(moments sim, Grid2D grid, int y)
{
    swapLayerMem<<<blockNumberLayer, blockLayer>>>(sim, y);
    streaming_at_y<<<blockNumberLayer, blockLayer>>>(sim, grid, 0);
}

__host__ void lastLayer(moments sim, Grid2D grid, int y)
{
    swapLayerMem<<<blockNumberLayer, blockLayer>>>(sim, y);
    streaming_at_y<<<blockNumberLayer, blockLayer>>>(sim, grid, 0);
}