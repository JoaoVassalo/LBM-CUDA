#include "layer_pipeline.cuh"

__global__ void seed_layer(D2Q9 sim, layer layer, Grid2D grid)
{
    propagate_layer_at(sim, layer, grid, 0);
}

__global__ void advance_layer(D2Q9 sim, layer layer, Grid2D grid, int y)
{
    propagate_layer_at(sim, layer, grid, y);
}

__global__ void end_layer(D2Q9 sim, layer layer, Grid2D grid)
{
    propagate_layer_at(sim, layer, grid, Geometry::Ny - 1);
}

__global__ void final_layers(D2Q9 sim, layer layer, Grid2D grid)
{
    propagate_layer_at(sim, layer, grid, Geometry::Ny - 1);
}
