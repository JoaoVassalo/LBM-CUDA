#include "calc_velocity.cuh"

#include "../presets/stencil.cuh"
#include "../presets/physics.h"

#include "grid_id.cuh"
#include "from_id.cuh"

__device__ void calc_velocity(int x, int y, float* rho_in, float* ux_in, float* uy_in, float* mxx_in, float* mxy_in,  float* myy_in,
                             float* rho_out, float* ux_out, float* uy_out){

    int index = grid_id(x, y);

    if(x==0 && y==0){ //Sudoeste
        ux_out[index] = 0.f;
        uy_out[index] = 0.f;
    }
    else if (x==Nx-1 && y==0){ //Sudeste
        ux_out[index] = 0.f;
        uy_out[index] = 0.f;
    }
    else if (x==0 && y==Ny-1){ //Noroeste
        ux_out[index] = u_max;
        uy_out[index] = 0.f;
    }
    else if (x==Nx-1 && y==Ny-1){ //Nordeste
        ux_out[index] = u_max;
        uy_out[index] = 0.f;
    }
    else if (y==0){ //Sul
        ux_out[index] = 0.f;
        uy_out[index] = 0.f;
    }
    else if (x==0 && y==Ny-1){ //Norte
        ux_out[index] = u_max;
        uy_out[index] = 0.f;
    }
    else if (x==0){ //Oeste
        ux_out[index] = 0;
        uy_out[index] = 0.f;
    }
    else if (x==Nx-1){ //Leste
        ux_out[index] = 0.f;
        uy_out[index] = 0.f;
    }
    else { //Centro
        ux_out[index] = 0.f;
        uy_out[index] = 0.f;
        for (int i = 0; i<4; i++){

            int index_from = from_id(x, y, i);

            float f_i = rho_in[index_from]*w[i] * ( 1+ 
                    a_s2*ux_in[index_from]*c_ix[i] + a_s2*uy_in[index_from]*c_iy[i] + 
                    a_s4*0.5f*mxx_in[index_from]*(c_ix[i]*c_ix[i] - inv_as2) + 
                    a_s4*0.5f*mxy_in[index_from]*c_ix[i]*c_iy[i] + 
                    a_s4*0.5f*myy_in[index_from]*(c_iy[i]*c_iy[i] - inv_as2));
            ux_out[index] += f_i*c_ix[i];
            uy_out[index] += f_i*c_iy[i];
        }

        ux_out[index] *= 1/rho_out[index];
        uy_out[index] *= 1/rho_out[index];
    }

}