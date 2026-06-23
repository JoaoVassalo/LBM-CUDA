#include "config/stencilConfig.cuh"
#include "config/gridConfig.cuh"
#include "lbm/initialization/initDomain.cuh"
#include "lbm/initialization/initGrid.cuh"
#include "lbm/initialization/initialization.cuh"
#include "lbm/simulation/simulation.cuh"

int main()
{
    moments sim;

    Grid2D grid;

    initialization(sim, grid);

    simulation(sim, grid);
}