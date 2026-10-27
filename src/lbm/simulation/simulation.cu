#include "simulation.cuh"

__host__ void step(moments &sim, Grid2D &grid)
{
    if constexpr (Geometry::NX <= 1024)
    {
        constexpr int blockX = ((Geometry::NX + 31) / 32) * 32;
        streamCollideStep<<<1, blockX>>>(sim, grid);
    }
    else
    {
        initLayer(sim, grid);

        for (int y = 1; y < Geometry::NY - 1; y++)
        {
            midLayer(sim, grid, y);
        }

        lastLayer(sim, grid);
    }
}

__host__ void simulation(moments &sim, Grid2D &grid)
{
    const auto startTime = std::chrono::steady_clock::now();

    for (int t = 0; t < timeConfig::tf; t++)
    {
        step(sim, grid);

        printProgress(t + 1, timeConfig::tf, startTime);

        if (t == 1 || t % timeConfig::tInterval == 0)
            writeOutput(sim, grid, t, D2Q9::momByteSize, output::vtkPath);

    }

    std::cout << '\n';
}
