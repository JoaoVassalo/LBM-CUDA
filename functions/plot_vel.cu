#include "../presets/geometry.h"
#include "../presets/physics.h"

#include "grid_plot.cuh"

__global__ void plot_vel(float *ux_host, float *uy_host, float *xplot, float *ux_plot, float *uy_plot)
{

    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= Nx)
        return;

    int mid_up = Nx / 2;
    int mid_down = (Nx / 2) - 1;

    ux_plot[i] = (ux_host[grid_plot(mid_up, i)] + ux_host[grid_plot(mid_down, i)]) / 2.f;
    ux_plot[i] /= u_max;
    uy_plot[i] = (uy_host[grid_plot(i, mid_up)] + uy_host[grid_plot(i, mid_down)]) / 2.f;
    uy_plot[i] /= u_max;

    xplot[i] = (float)(i) / (float)Nx;
}