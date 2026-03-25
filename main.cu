#include <iostream>
#include <string>
#include <chrono>
#include <fstream>

#include "presets/stencil.cuh"
#include "presets/geometry.h"
#include "presets/config.h"
#include "presets/physics.h"

#include "functions/init_domain.cuh"
#include "functions/lbm_step.cuh"
#include "functions/calc_tke.cuh"
#include "functions/grid_id.cuh"
#include "vtk.cuh"

int main()
{

    dim3 block(32, 16);
    dim3 N_block((Nx + block.x - 1) / block.x, (Ny + block.y - 1) / block.y);

    float *rhoA,
        *rhoB;
    float *uxA;
    float *uxB;
    float *uyA;
    float *uyB;
    float *mxxA;
    float *mxxB;
    float *mxyA;
    float *mxyB;
    float *myyA;
    float *myyB;

    int size = Nx * Ny * sizeof(float);

    cudaMalloc((void **)&rhoA, size);
    cudaMalloc((void **)&rhoB, size);
    cudaMalloc((void **)&uxA, size);
    cudaMalloc((void **)&uxB, size);
    cudaMalloc((void **)&uyA, size);
    cudaMalloc((void **)&uyB, size);
    cudaMalloc((void **)&mxxA, size);
    cudaMalloc((void **)&mxxB, size);
    cudaMalloc((void **)&mxyA, size);
    cudaMalloc((void **)&mxyB, size);
    cudaMalloc((void **)&myyA, size);
    cudaMalloc((void **)&myyB, size);

    initDomain<<<N_block, block>>>(rhoA, uxA, uyA, mxxA, mxyA, myyA);

    cudaError_t err = cudaGetLastError();
    printf("Kernel launch error: %s\n", cudaGetErrorString(err));

    cudaDeviceSynchronize();

    cudaMemcpy(rhoB, rhoA, size, cudaMemcpyDeviceToDevice);
    cudaMemcpy(uxB, uxA, size, cudaMemcpyDeviceToDevice);
    cudaMemcpy(uyB, uyA, size, cudaMemcpyDeviceToDevice);
    cudaMemcpy(mxxB, mxxA, size, cudaMemcpyDeviceToDevice);
    cudaMemcpy(mxyB, mxyA, size, cudaMemcpyDeviceToDevice);
    cudaMemcpy(myyB, myyA, size, cudaMemcpyDeviceToDevice);

    float *rho_host = (float *)malloc(size);
    float *ux_host = (float *)malloc(size);
    float *uy_host = (float *)malloc(size);

    std::ofstream file("animation/tke.csv");
    file << "time,tke\n";

    auto t1 = std::chrono::high_resolution_clock::now();

    for (int t = 0; t < tf; t++)
    {

        if (t & 1)
        {
            lbm_step<<<N_block, block>>>(rhoA, uxA, uyA, mxxA, mxyA, myyA,
                                         rhoB, uxB, uyB, mxxB, mxyB, myyB);
        }
        else
        {
            lbm_step<<<N_block, block>>>(rhoB, uxB, uyB, mxxB, mxyB, myyB,
                                         rhoA, uxA, uyA, mxxA, mxyA, myyA);
        }

        if (t % t_interval == 0)
        {
            cudaDeviceSynchronize();
            cudaMemcpy(rho_host, rhoA, size, cudaMemcpyDeviceToHost);
            cudaMemcpy(ux_host, uxA, size, cudaMemcpyDeviceToHost);
            cudaMemcpy(uy_host, uyA, size, cudaMemcpyDeviceToHost);
            std::cout << "Iteration " << t << std::endl;

            std::cout << std::endl;

            std::string path = "./plot";

            write_vti(t, path, rho_host, ux_host, uy_host);

            calc_tke(file, ux_host, uy_host, t);
        }
    }

    auto t2 = std::chrono::high_resolution_clock::now();

    auto tempo = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1);

    double MLUPS = (double)(Nx * Ny) * double(tf) / ((double)tempo.count() * 1e6);
    printf("MLUPS: %f\n", MLUPS);

    file.close();

    int mid_up = Nx / 2;
    int mid_down = (Nx / 2) - 1;

    float ux_plot = 0;
    float uy_plot = 0;

    std::ofstream file2("animation/vel.csv");
    file2 << "x,ux,uy\n";

    for (int i = 0; i < Nx; i++)
    {

        ux_plot = (ux_host[grid_id(mid_up, i)] + ux_host[grid_id(mid_down, i)]) / 2.f;
        ux_plot /= u_max;
        uy_plot = (uy_host[grid_id(i, mid_up)] + uy_host[grid_id(i, mid_down)]) / 2.f;
        uy_plot /= u_max;

        float xplot = (float)(i) / (float)Nx;

        file2 << xplot << "," << ux_plot << "," << uy_plot << "\n";
    }
}