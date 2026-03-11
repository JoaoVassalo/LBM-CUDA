#pragma once

__device__ void calc_density(int x, int y, float* rho_in, float* ux_in, float* uy_in, float* mxx_in, float* mxy_in,  float* myy_in,
                             float* rho_out);