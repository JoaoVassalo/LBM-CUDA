#include "layerPipeline.cuh"

__host__ void initLayer(moments sim, Grid2D grid)
{
    defLayerMem<<<blockNumberLayer, blockLayer>>>(sim);
}