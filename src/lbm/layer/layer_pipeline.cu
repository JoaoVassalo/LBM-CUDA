#include "layer_pipeline.cuh"

__host__ void seed_layer(D2Q9 sim, layer layer, Grid2D grid)
{
    init_layers<<<layer.layer_Nblock, layer.layer_block>>>(sim, layer);
    propagate_layer_at<<<layer.layer_Nblock, layer.layer_block>>>(sim, layer, grid, 0);
    collide_layer_at<<<layer.layer_Nblock, layer.layer_block>>>(sim, 0);
}

__host__ void advance_layer(D2Q9 sim, layer layer, Grid2D grid, int y)
{
    propagate_layer_at<<<layer.layer_Nblock, layer.layer_block>>>(sim, layer, grid, y);
    swap_layers<<<layer.layer_Nblock, layer.layer_block>>>(sim, layer, y);
    collide_layer_at<<<layer.layer_Nblock, layer.layer_block>>>(sim, y);
}

__host__ void final_layers(D2Q9 sim, layer layer, Grid2D grid)
{
    propagate_layer_at<<<layer.layer_Nblock, layer.layer_block>>>(sim, layer, grid, Geometry::Ny - 1);
    collide_layer_at<<<layer.layer_Nblock, layer.layer_block>>>(sim, Geometry::Ny - 1);
}