#include "simulation.cuh"

__host__ void step(moments &sim, Grid2D &grid)
{
    initLayer(sim, grid);

    for (int y = 1; y < Geometry::NY - 1; y++)
    {
        midLayer(sim, grid, y);
    }

    lastLayer(sim, grid);
}

__host__ void simulation(moments &sim, Grid2D &grid)
{
    auto t1 = std::chrono::high_resolution_clock::now();

    for (int t = 0; t < timeConfig::tf; t++)
    {
        step(sim, grid);

        if (t == 1 || t % timeConfig::tInterval == 0)
            writeOutput(sim, grid, t, D2Q9::momByteSize, output::vtkPath);

        return;
    }

    auto t2 = std::chrono::high_resolution_clock::now();
}
