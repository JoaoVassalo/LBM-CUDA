#include "step.cuh"

__host__ void step(D2Q9 sim, layer layer, Grid2D grid)
{
    seed_layer<<<layer.layer_Nblock, layer.layer_block>>>(sim, layer, grid);

    for (int y = 1; y < Geometry::Ny - 1; y++)
    {
        advance_layer<<<layer.layer_Nblock, layer.layer_block>>>(sim, layer, grid, y);
    }

    end_layer<<<layer.layer_Nblock, layer.layer_block>>>(sim, layer, grid);
}