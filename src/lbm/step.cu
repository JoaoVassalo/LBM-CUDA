#include "step.cuh"

__host__ void step(D2Q9 sim, layer layer, Grid2D grid)
{
    layer.yref = 1;
    init_layers<<<layer.layer_Nblock, layer.layer_block>>>(sim, layer);
    seed_layer<<<layer.layer_Nblock, layer.layer_block>>>(sim, layer, grid);

    for (int y = 1; y < Geometry::Ny - 1; y++)
    {
        layer.yref = y;
        advance_layer<<<layer.layer_Nblock, layer.layer_block>>>(sim, layer, grid, y);
        swap_layers<<<layer.layer_Nblock, layer.layer_block>>>(sim, layer);
    }

    layer.yref = Geometry::Ny - 1;
    final_layers<<<layer.layer_Nblock, layer.layer_block>>>(sim, layer, grid);
}
