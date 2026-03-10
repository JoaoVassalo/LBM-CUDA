#include "calc_velocity.cuh"

#include "../presets/stencil.h"

#include "grid_id.cuh"

__device__ float calc_velocity(int x, int y, float* rho, float* u_x, float* u_y, float* mxx, float* mxy, float* myy){

    int index = grid_id(x, y);
;
    u_x[index] = 0.f;
    u_y[index] = 0.f;

    for (int i = 0; i<4; i++){
        float f_i = rho[index]*w[i] * ( 1+ 
                a_s2*u_x[index]*c_ix[i] + a_s2*u_y[index]*c_iy[i] + 
                a_s4*0.5f*mxx[index]*(c_ix[i]*c_ix[i] - inv_as2) + 
                a_s4*0.5f*mxy[index]*c_ix[i]*c_iy[i] + 
                a_s4*0.5f*myy[index]*(c_iy[i]*c_iy[i] - inv_as2));
        u_x[index] += f_i*c_ix[i];
        u_y[index] += f_i*c_iy[i];
    }

    u_x[index] *= 1/rho[index];
    u_y[index] *= 1/rho[index];

}