#include "calc_velocity.cuh"

#include "../presets/stencil.h"

__device__ float calc_velocity(float mxx, float mxy, float myy, float rho){
    float u_x = 0.f;
    float u_y = 0.f;

    for (int i = 0; i<4; i++){
        float f_i = rho*w[i] * ( 1+ 
                a_s2*u_x*c_ix[i] + a_s*a_s*u_y*c_iy[i] + 
                a_s4*0.5f*mxx*(c_ix[i]*c_ix[i] - inv_as2) + 
                a_s4*0.5f*mxy*c_ix[i]*c_iy[i] + 
                a_s4*0.5f*myy*(c_iy[i]*c_iy[i] - inv_as2));
        u_x += f_i*c_ix[i];
        u_y += f_i*c_iy[i];
    }

    u_x *= 1/rho;
    u_y *= 1/rho;


    return u_x, u_y;
}