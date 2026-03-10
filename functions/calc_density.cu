#include "calc_density.cuh"
#include "../presets/stencil.h"
#include "../presets/geometry.h"
#include "../presets/physics.h"

__device__ float calc_density(float x, float y, float rho, float u_x, float u_y, float mxx, float mxy, float myy){

    float density = 0.f;

    
    if (x==0 && y==0){ //Sudoeste
        constexpr int I_s[4] = {0, 3, 4, 7};

        float density_I = 0.f;
        
        #pragma unroll
        for (int k = 0; k<4; k++){
            int i = I_s[k];

            density_I += rho*w[i] * ( 1+ 
                a_s2*u_x*c_ix[i] + a_s*a_s*u_y*c_iy[i] + 
                a_s4*0.5f*mxx*(c_ix[i]*c_ix[i] - inv_as2) + 
                a_s4*0.5f*mxy*c_ix[i]*c_iy[i] + 
                a_s4*0.5f*myy*(c_iy[i]*c_iy[i] - inv_as2));

        }

        float sum_wi = w[0] + w[1] + w[2] + w[5];

        density = density_I/sum_wi;

    }


    else if (x==Nx-1 && y==0){ //Sudeste   
        constexpr int I_s[4] = {0, 1, 4, 8};

        float density_I = 0.f;
        
        #pragma unroll
        for (int k = 0; k<4; k++){
            int i = I_s[k];

            density_I += rho*w[i] * ( 1+ 
                        a_s2*u_x*c_ix[i] + a_s*a_s*u_y*c_iy[i] + 
                        a_s4*0.5f*mxx*(c_ix[i]*c_ix[i] - inv_as2) + 
                        a_s4*0.5f*mxy*c_ix[i]*c_iy[i] + 
                        a_s4*0.5f*myy*(c_iy[i]*c_iy[i] - inv_as2));

        }

        float sum_wi = w[0] + w[2] + w[3] + w[6];

        density = density_I/sum_wi;
    }
    else if (x==0 && y==Ny-1){ //Noroeste
        constexpr int I_s[4] = {0, 2, 3, 6};
        constexpr int O_s[4] = {0, 1, 4, 8};

        float density_I_density = 0.f;
        float density_I = 0.f;

        #pragma unroll
        for (int k = 0; k<4; k++){
            int i = O_s[k];

            density_I_density += w[i]*(1 + 
                                a_s2*u_x*c_ix[i] + 
                                a_s4*0.5f*u_x*u_x*(c_ix[i]*c_ix[i] - inv_as2)) + 
                                w[i]*(1-omega)*a_s4*mxy*c_ix[i]*c_iy[i];

        }

        #pragma unroll
        for (int k = 0; k<4; k++){
            int i = I_s[k];

            density_I += rho*w[i] * ( 1+ 
                        a_s2*u_x*c_ix[i] + a_s*a_s*u_y*c_iy[i] + 
                        a_s4*0.5f*mxx*(c_ix[i]*c_ix[i] - inv_as2) + 
                        a_s4*0.5f*mxy*c_ix[i]*c_iy[i] + 
                        a_s4*0.5f*myy*(c_iy[i]*c_iy[i] - inv_as2));
        }

        density = density_I/density_I_density;

    }
    else if (x==Nx-1 && y==Ny-1){ //Nordeste
        constexpr int I_s[4] = {0, 1, 2, 5};
        constexpr int O_s[4] = {0, 3, 4, 7};

        float density_I_density = 0.f;
        float density_I = 0.f;

        #pragma unroll
        for (int k = 0; k<4; k++){
            int i = O_s[k];

            density_I_density += w[i]*(1 + 
                                a_s2*u_x*c_ix[i] + 
                                a_s4*0.5f*u_x*u_x*(c_ix[i]*c_ix[i] - inv_as2)) + 
                                w[i]*(1-omega)*a_s4*mxy*c_ix[i]*c_iy[i];

        }

        #pragma unroll
        for (int k = 0; k<4; k++){
            int i = I_s[k];

            density_I += rho*w[i] * ( 1+ 
                        a_s2*u_x*c_ix[i] + a_s*a_s*u_y*c_iy[i] + 
                        a_s4*0.5f*mxx*(c_ix[i]*c_ix[i] - inv_as2) + 
                        a_s4*0.5f*mxy*c_ix[i]*c_iy[i] + 
                        a_s4*0.5f*myy*(c_iy[i]*c_iy[i] - inv_as2));
        }

        density = density_I/density_I_density;
    }
    else if (y==0){ //Sul
        constexpr int I_s[6] = {0, 1, 3, 4, 7, 8};
        constexpr int O_s[6] = {0, 1, 2, 3, 5, 6};

        float density_I_density = 0.f;
        float density_I = 0.f;

        #pragma unroll
        for (int k = 0; k<6; k++){
            int i = O_s[k];

            density_I_density += w[i] + 
                                w[i]*(1-omega)*a_s4*mxy*c_ix[i]*c_iy[i];

        }

        #pragma unroll
        for (int k = 0; k<6; k++){
            int i = I_s[k];

            density_I += rho*w[i] * ( 1+ 
                        a_s2*u_x*c_ix[i] + a_s*a_s*u_y*c_iy[i] + 
                        a_s4*0.5f*mxx*(c_ix[i]*c_ix[i] - inv_as2) + 
                        a_s4*0.5f*mxy*c_ix[i]*c_iy[i] + 
                        a_s4*0.5f*myy*(c_iy[i]*c_iy[i] - inv_as2));
        }

        density = density_I/density_I_density;
    }
    else if (y==Ny-1){ //Norte
        constexpr int I_s[6] = {0, 1, 2, 3, 5, 6};
        constexpr int O_s[6] = {0, 1, 3, 4, 7, 8};

        float density_I_density = 0.f;
        float density_I = 0.f;

        #pragma unroll
        for (int k = 0; k<6; k++){
            int i = O_s[k];

            density_I_density += w[i]*(1 + 
                                a_s2*u_x*c_ix[i] + 
                                a_s4*0.5f*u_x*u_x*(c_ix[i]*c_ix[i] - inv_as2)) + 
                                w[i]*(1-omega)*a_s4*mxy*c_ix[i]*c_iy[i];

        }

        #pragma unroll
        for (int k = 0; k<6; k++){
            int i = I_s[k];

            density_I += rho*w[i] * ( 1+ 
                        a_s2*u_x*c_ix[i] + a_s*a_s*u_y*c_iy[i] + 
                        a_s4*0.5f*mxx*(c_ix[i]*c_ix[i] - inv_as2) + 
                        a_s4*0.5f*mxy*c_ix[i]*c_iy[i] + 
                        a_s4*0.5f*myy*(c_iy[i]*c_iy[i] - inv_as2));
        }

        density = density_I/density_I_density;
    }
    else if (x==0){ //Oeste
        constexpr int I_s[6] = {0, 2, 3, 4, 6, 7};
        constexpr int O_s[6] = {0, 1, 2, 4, 5, 8};

        float density_I_density = 0.f;
        float density_I = 0.f;

        #pragma unroll
        for (int k = 0; k<6; k++){
            int i = O_s[k];

            density_I_density += w[i] + 
                                w[i]*(1-omega)*a_s4*mxy*c_ix[i]*c_iy[i];

        }

        #pragma unroll
        for (int k = 0; k<6; k++){
            int i = I_s[k];

            density_I += rho*w[i] * ( 1+ 
                        a_s2*u_x*c_ix[i] + a_s*a_s*u_y*c_iy[i] + 
                        a_s4*0.5f*mxx*(c_ix[i]*c_ix[i] - inv_as2) + 
                        a_s4*0.5f*mxy*c_ix[i]*c_iy[i] + 
                        a_s4*0.5f*myy*(c_iy[i]*c_iy[i] - inv_as2));
        }

        density = density_I/density_I_density;
    }
    else if (x==Nx-1){ //Leste
        constexpr int I_s[6] = {0, 1, 2, 4, 5, 8};
        constexpr int O_s[6] = {0, 2, 3, 4, 6, 7};

        float density_I_density = 0.f;
        float density_I = 0.f;

        #pragma unroll
        for (int k = 0; k<6; k++){
            int i = O_s[k];

            density_I_density += w[i] + 
                                w[i]*(1-omega)*a_s4*mxy*c_ix[i]*c_iy[i];

        }

        #pragma unroll
        for (int k = 0; k<6; k++){
            int i = I_s[k];

            density_I += rho*w[i] * ( 1+ 
                        a_s2*u_x*c_ix[i] + a_s*a_s*u_y*c_iy[i] + 
                        a_s4*0.5f*mxx*(c_ix[i]*c_ix[i] - inv_as2) + 
                        a_s4*0.5f*mxy*c_ix[i]*c_iy[i] + 
                        a_s4*0.5f*myy*(c_iy[i]*c_iy[i] - inv_as2));
        }

        density = density_I/density_I_density;
    }
    else{ //Centro
        #pragma unroll
        for (int i = 0; i < Q; i++) {

            density += rho*w[i] * ( 1+ 
                a_s2*u_x*c_ix[i] + a_s*a_s*u_y*c_iy[i] + 
                a_s4*0.5f*mxx*(c_ix[i]*c_ix[i] - inv_as2) + 
                a_s4*0.5f*mxy*c_ix[i]*c_iy[i] + 
                a_s4*0.5f*myy*(c_iy[i]*c_iy[i] - inv_as2));

        }
    }


    return density;
}