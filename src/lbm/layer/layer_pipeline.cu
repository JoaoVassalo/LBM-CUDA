#include "layer_pipeline.cuh"

__global__ void seed_layer(D2Q9 sim, Grid2D grid)
{
    first_layer(sim);
    propagation(sim, grid);
}

__global__ void advance_layer(D2Q9 sim, Grid2D grid, int y)
{
    for (int i = 0; i < sim.layer_size; i++)
    {
        other_layers(sim, y, i);
        propagation(sim, grid);
    }
}

__global__ void end_layer(D2Q9 sim, Grid2D grid)
{
    last_layer(sim);
    propagation(sim, grid);
}