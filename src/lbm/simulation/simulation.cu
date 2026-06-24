#include "simulation.cuh"

__host__ void simulation(moments sim, Grid2D grid)
{
    auto t1 = std::chrono::high_resolution_clock::now();

    for (size_t t = 0; t < time::tf; t++)
    {
        step(sim, grid);

        if (t == 1 || t % time::tInterval == 0)
            writeOutput(sim, grid, t, D2Q9::momByteSize, output::vtkPath);
    }

    auto t2 = std::chrono::high_resolution_clock::now();
}

__host__ void step(moments sim, Grid2D grid)
{
    initLayer(sim, grid);

    for (size_t y = 1; y < Geometry::NY - 1; y++)
    {
        }
}