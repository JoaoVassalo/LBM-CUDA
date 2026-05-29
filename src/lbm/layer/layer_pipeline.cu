#include "layer_pipeline.cuh"

__global__ void seed_layer(D2Q9 sim, layer layer, Grid2D grid)
{
    init_layers(sim, layer);
    propagation(sim, layer, grid);
}

__global__ void advance_layer(D2Q9 sim, layer layer, Grid2D grid, int y)
{
    propagation(sim, layer, grid);
    constexpr_for<int(0), layer::LNy>([&](const auto i)
                                      { swap_layers(sim, layer); });
}

__global__ void end_layer(D2Q9 sim, layer layer, Grid2D grid)
{
    propagation(sim, layer, grid);
}