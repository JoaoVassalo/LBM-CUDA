#include "step.cuh"

__host__ void step(D2Q9 sim, Grid2D grid)
{
    seed_layer<<<sim.layer_Nblock, sim.layer_block>>>(sim, grid);

    for (int y = 1; y < Geometry::Ny - 1; y++)
    {
        advance_layer<<<sim.layer_Nblock, sim.layer_block>>>(sim, grid, y);
    }

    end_layer<<<sim.layer_Nblock, sim.layer_block>>>(sim, grid);
}