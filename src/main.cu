#include <iostream>
#include <string>
#include <chrono>
#include <fstream>

#include "presets/stencil.cuh"
#include "presets/geometry.h"
#include "presets/config.h"
#include "presets/physics.h"

#include "lbm_config/mom_config.cuh"
#include "lbm_config/var.cuh"

#include "functions/init_domain.cuh"
#include "functions/lbm_step.cuh"
#include "functions/calc_tke.cuh"
#include "functions/grid_id.cuh"
#include "functions/grid_plot.cuh"

#include "vtk.cuh"

int main()
{

    D2Q9 sim;

    cudaMalloc((void **)&sim.mom, sim.size);
    cudaMalloc((void **)&sim.layer, sim.layer_size);

    initDomain<<<sim.N_block, sim.block>>>(sim.mom, sim.layer);

    cudaError_t err = cudaGetLastError();
    printf("Kernel launch error: %s\n", cudaGetErrorString(err));

    cudaDeviceSynchronize();

    float *mom_host = (float *)malloc(sim.size);

    std::ofstream file("animation/tke.csv");
    file << "time,tke\n";

    auto t1 = std::chrono::high_resolution_clock::now();

    for (int t = 0; t < tf; t++)
    {

        if (t & 1)
        {
            lbm_step<<<sim.N_block, sim.block>>>(sim.momA, sim.momB);
        }
        else
        {
            lbm_step<<<sim.N_block, sim.block>>>(sim.momB, sim.momA);
        }

        if (t % t_interval == 0)
        {
            cudaDeviceSynchronize();
            cudaMemcpy(mom_host, sim.momA, sim.size, cudaMemcpyDeviceToHost);

            std::cout << "Iteration " << t << std::endl;

            std::cout << std::endl;

            std::string path = "./plot";

            write_vti(t, path, mom_host);

            calc_tke(file, mom_host, t);
        }
    }

    auto t2 = std::chrono::high_resolution_clock::now();

    auto tempo = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1);

    double MLUPS = (double)(Nx * Ny) * double(tf) / ((double)tempo.count() * 1e6);
    printf("MLUPS: %f\n", MLUPS);

    file.close();

    std::ofstream file2("animation/vel.csv");
    file2 << "x,ux,uy\n";

    float ux_plot;
    float uy_plot;
    float xplot;

    int mid_up = Nx / 2;
    int mid_down = (Nx / 2) - 1;
    for (int i = 0; i < Nx; i++)
    {
        ux_plot = (mom_host[grid_plot(mid_up, i) + 1] + mom_host[grid_plot(mid_down, i) + 1]) / 2.f; // ux
        ux_plot /= u_max;
        uy_plot = (mom_host[grid_plot(i, mid_up) + 2] + mom_host[grid_plot(i, mid_down) + 2]) / 2.f; // uy
        uy_plot /= u_max;

        xplot = (float)(i) / (float)Nx;

        file2 << xplot << "," << ux_plot << "," << uy_plot << "\n";
    }
    file2.close();
}