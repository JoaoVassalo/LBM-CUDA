#include "step.cuh"

__host__ void step(D2Q9 sim, layer layer, Grid2D grid)
{
    seed_layer(sim, layer, grid);

    for (int y = 1; y < Geometry::Ny - 1; y++)
    {
        advance_layer(sim, layer, grid, y);
    }

    final_layers(sim, layer, grid);
}
