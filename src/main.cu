#include <iostream>
#include <string>
#include <chrono>
#include <fstream>
#include <cstdint>

#include "config/stencil.cuh"
#include "config/geometry.h"
#include "config/mom_config.cuh"
#include "config/physics.h"

#include "lbm/init/init_domain.cuh"
#include "lbm/step.cuh"
#include "lbm/build/build_mom.cuh"
#include "lbm/build/build_grid.cuh"

#include "core/calc_tke.cuh"
#include "core/indexing.cuh"

#include "io/vtk.cuh"

int main()
{
    D2Q9 sim;

    layer layer;

    cudaMalloc((void **)&sim.mom, sim.size);

    for (int i = 0; i < layer::LNy; ++i)
    {
        cudaMalloc((void **)&layer.buffer[i], layer::buffer_bytesize);
    }

    initDomain<<<sim.N_block, sim.block>>>(sim.mom);

    Grid2D grid;
    build_grid(sim, grid);

    cudaError_t err = cudaGetLastError();
    printf("Kernel launch error: %s\n", cudaGetErrorString(err));

    cudaDeviceSynchronize();

    float *mom_host = (float *)malloc(sim.size);

    std::ofstream file("animation/tke.csv");
    file << "time,tke\n";

    auto t1 = std::chrono::high_resolution_clock::now();

    for (int t = 0; t < tf; t++)
    {
        step(sim, layer, grid);

        if (t % t_interval == 0)
        {
            cudaDeviceSynchronize();
            cudaMemcpy(mom_host, sim.mom, sim.size, cudaMemcpyDeviceToHost);

            std::cout << "Iteration " << t << std::endl;

            std::cout << std::endl;

            std::string path = "./plot";

            write_vti(t, path, mom_host);

            calc_tke(file, mom_host, t);
        }
    }

    auto t2 = std::chrono::high_resolution_clock::now();

    auto tempo = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1);

    double MLUPS = (double)(Geometry::Nx * Geometry::Ny) * double(tf) / ((double)tempo.count() * 1e6);
    printf("MLUPS: %f\n", MLUPS);

    file.close();

    std::ofstream file2("animation/vel.csv");
    file2 << "x,ux,uy\n";

    float ux_plot;
    float uy_plot;
    float xplot;

    int mid_up = Geometry::Nx / 2;
    int mid_down = (Geometry::Nx / 2) - 1;
    for (int i = 0; i < Geometry::Nx; i++)
    {
        ux_plot = (mom_host[grid_plot(mid_up, i) + 1] + mom_host[grid_plot(mid_down, i) + 1]) / 2.f; // ux
        // ux_plot /= u_max;
        uy_plot = (mom_host[grid_plot(i, mid_up) + 2] + mom_host[grid_plot(i, mid_down) + 2]) / 2.f; // uy
        // uy_plot /= u_max;

        xplot = (float)(i) / (float)Geometry::Nx;

        file2 << xplot << "," << ux_plot << "," << uy_plot << "\n";
    }
    file2.close();
}