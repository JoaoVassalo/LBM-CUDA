#pragma once

__device__ void calc_momentum(int x, int y, float* rho_in, float* ux_in, float* uy_in, float* mxx_in, float* mxy_in,  float* myy_in,
                                             float* rho_out, float* ux_out, float* uy_out, float* mxx_out, float* mxy_out,  float* myy_out);