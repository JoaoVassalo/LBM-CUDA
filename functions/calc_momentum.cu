#include "calc_momentum.cuh"

#include "../presets/stencil.h"
#include "../presets/geometry.h"
#include "../presets/physics.h"

__device__ float calc_momentum(int x, int y, float rho, float u_x, float u_y, float mxx, float mxy, float myy){

    if (x==0 && y==0){ //Sudoeste

        float mxx = 0.f;
        float myy = 0.f;
        float mxy = 0.f;

    }
    else if (x==Nx-1 && y==0){ //Sudeste

        float mxx = 0.f;
        float myy = 0.f;
        float mxy = 0.f;

    }
    else if (x==0 && y==Ny-1){ //Noroeste
        constexpr int I_s[4] = {0, 2, 3, 6};
        constexpr int O_s[4] = {0, 1, 4, 8};


        float mxx = 0.f;
        float myy = u_max*u_max;
        float mxy = 0.f;

        float num = 0;
        float div = 0;

        float Is_up = 0;
        float Is_down = 0;
        float Os_up = 0;
        float Os_down = 0;

        #pragma unroll
        for (int k=0; k<4; k++){
            int i = I_s[k];
            int j = O_s[k];

            num += rho*w[i]*(1 + 
                    a_s2*u_x*c_ix[i] + a_s2*u_y*c_iy[i] + 
                    a_s4*0.5f*mxx*(c_ix[i]*c_ix[i] - inv_as2) + 
                    a_s4*0.5f*mxy*c_ix[i]*c_iy[i] + 
                    a_s4*0.5f*myy*(c_iy[i]*c_iy[i] - inv_as2))*c_ix[i]*c_iy[i];

            div += rho*w[i]*(1 + 
                    a_s2*u_x*c_ix[i] + a_s2*u_y*c_iy[i] + 
                    a_s4*0.5f*mxx*(c_ix[i]*c_ix[i] - inv_as2) + 
                    a_s4*0.5f*mxy*c_ix[i]*c_iy[i] + 
                    a_s4*0.5f*myy*(c_iy[i]*c_iy[i] - inv_as2));


            Is_up += w[i]*c_ix[i]*c_iy[i]*(1 + 
                         a_s2*u_max*c_ix[i] + 
                         a_s4*0.5f*u_max*(c_ix[i]*c_ix[i] - inv_as2));

            Is_down += w[i]*a_s4*c_ix[i]*c_ix[i]*c_iy[i]*c_iy[i];

            Os_up += w[i]*(1 + a_s2*u_max*c_ix[i] + a_s4*0.5f*u_max*u_max*(c_ix[i]*c_ix[i] - inv_as2));

            Os_down += w[i]*a_s4*c_ix[i]*c_iy[i];

        }

        float mxy_I = num/div;

        mxy = (Is_up - mxy_I*Os_up)/(mxy_I*(1-omega)*Os_down - Is_down);


    }
    else if (x==Nx-1 && y==0){ //Nordeste
        constexpr int I_s[4] = {0, 1, 2, 5};
        constexpr int O_s[4] = {0, 3, 4, 7};


        float mxx = 0.f;
        float myy = u_max*u_max;
        float mxy = 0.f;

        float num = 0;
        float div = 0;

        float Is_up = 0;
        float Is_down = 0;
        float Os_up = 0;
        float Os_down = 0;

        #pragma unroll
        for (int k=0; k<4; k++){
            int i = I_s[k];
            int j = O_s[k];

            num += rho*w[i]*(1 + 
                    a_s2*u_x*c_ix[i] + a_s2*u_y*c_iy[i] + 
                    a_s4*0.5f*mxx*(c_ix[i]*c_ix[i] - inv_as2) + 
                    a_s4*0.5f*mxy*c_ix[i]*c_iy[i] + 
                    a_s4*0.5f*myy*(c_iy[i]*c_iy[i] - inv_as2))*c_ix[i]*c_iy[i];

            div += rho*w[i]*(1 + 
                    a_s2*u_x*c_ix[i] + a_s2*u_y*c_iy[i] + 
                    a_s4*0.5f*mxx*(c_ix[i]*c_ix[i] - inv_as2) + 
                    a_s4*0.5f*mxy*c_ix[i]*c_iy[i] + 
                    a_s4*0.5f*myy*(c_iy[i]*c_iy[i] - inv_as2));


            Is_up += w[i]*c_ix[i]*c_iy[i]*(1 + 
                         a_s2*u_max*c_ix[i] + 
                         a_s4*0.5f*u_max*(c_ix[i]*c_ix[i] - inv_as2));

            Is_down += w[i]*a_s4*c_ix[i]*c_ix[i]*c_iy[i]*c_iy[i];

            Os_up += w[i]*(1 + a_s2*u_max*c_ix[i] + a_s4*0.5f*u_max*u_max*(c_ix[i]*c_ix[i] - inv_as2));

            Os_down += w[i]*a_s4*c_ix[i]*c_iy[i];

        }

        float mxy_I = num/div;

        mxy = (Is_up - mxy_I*Os_up)/(mxy_I*(1-omega)*Os_down - Is_down);


    }
    else if (x==Nx-1 && y==0){ //Sul
        constexpr int I_s[6] = {0, 1, 3, 4, 7, 8};
        constexpr int O_s[6] = {0, 1, 2, 3, 5, 6};


        float mxx = 0.f;
        float myy = 0.f;
        float mxy = 0.f;

        float num = 0;
        float div = 0;

        float Is_up = 0;
        float Is_down = 0;
        float Os_up = 0;
        float Os_down = 0;

        #pragma unroll
        for (int k=0; k<6; k++){
            int i = I_s[k];
            int j = O_s[k];

            num += rho*w[i]*(1 + 
                    a_s2*u_x*c_ix[i] + a_s2*u_y*c_iy[i] + 
                    a_s4*0.5f*mxx*(c_ix[i]*c_ix[i] - inv_as2) + 
                    a_s4*0.5f*mxy*c_ix[i]*c_iy[i] + 
                    a_s4*0.5f*myy*(c_iy[i]*c_iy[i] - inv_as2))*c_ix[i]*c_iy[i];

            div += rho*w[i]*(1 + 
                    a_s2*u_x*c_ix[i] + a_s2*u_y*c_iy[i] + 
                    a_s4*0.5f*mxx*(c_ix[i]*c_ix[i] - inv_as2) + 
                    a_s4*0.5f*mxy*c_ix[i]*c_iy[i] + 
                    a_s4*0.5f*myy*(c_iy[i]*c_iy[i] - inv_as2));


            Is_up += w[i]*c_ix[i]*c_iy[i];

            Is_down += w[i]*a_s4*c_ix[i]*c_ix[i]*c_iy[i]*c_iy[i];

            Os_up += w[i];

            Os_down += w[i]*a_s4*c_ix[i]*c_iy[i];

        }

        float mxy_I = num/div;

        mxy = (Is_up - mxy_I*Os_up)/(mxy_I*(1-omega)*Os_down - Is_down);


    }
    else if (x==Nx-1 && y==0){ //Norte
        constexpr int I_s[6] = {0, 1, 2, 3, 5, 6};
        constexpr int O_s[6] = {0, 1, 3, 4, 7, 8};


        float mxx = 0.f;
        float myy = u_max*u_max;
        float mxy = 0.f;

        float num = 0;
        float div = 0;

        float Is_up = 0;
        float Is_down = 0;
        float Os_up = 0;
        float Os_down = 0;

        #pragma unroll
        for (int k=0; k<6; k++){
            int i = I_s[k];
            int j = O_s[k];

            num += rho*w[i]*(1 + 
                    a_s2*u_x*c_ix[i] + a_s2*u_y*c_iy[i] + 
                    a_s4*0.5f*mxx*(c_ix[i]*c_ix[i] - inv_as2) + 
                    a_s4*0.5f*mxy*c_ix[i]*c_iy[i] + 
                    a_s4*0.5f*myy*(c_iy[i]*c_iy[i] - inv_as2))*c_ix[i]*c_iy[i];

            div += rho*w[i]*(1 + 
                    a_s2*u_x*c_ix[i] + a_s2*u_y*c_iy[i] + 
                    a_s4*0.5f*mxx*(c_ix[i]*c_ix[i] - inv_as2) + 
                    a_s4*0.5f*mxy*c_ix[i]*c_iy[i] + 
                    a_s4*0.5f*myy*(c_iy[i]*c_iy[i] - inv_as2));


            Is_up += w[i]*c_ix[i]*c_iy[i]*(1 + 
                         a_s2*u_max*c_ix[i] + 
                         a_s4*0.5f*u_max*(c_ix[i]*c_ix[i] - inv_as2));

            Is_down += w[i]*a_s4*c_ix[i]*c_ix[i]*c_iy[i]*c_iy[i];

            Os_up += w[i]*(1 + a_s2*u_max*c_ix[i] + a_s4*0.5f*u_max*u_max*(c_ix[i]*c_ix[i] - inv_as2));

            Os_down += w[i]*a_s4*c_ix[i]*c_iy[i];

        }

        float mxy_I = num/div;

        mxy = (Is_up - mxy_I*Os_up)/(mxy_I*(1-omega)*Os_down - Is_down);

    }
    else if (x==Nx-1 && y==0){ //Oeste
        constexpr int I_s[6] = {0, 2, 3, 4, 6, 7};
        constexpr int O_s[6] = {0, 1, 2, 4, 5, 8};


        float mxx = 0.f;
        float myy = 0.f;
        float mxy = 0.f;

        float num = 0;
        float div = 0;

        float Is_up = 0;
        float Is_down = 0;
        float Os_up = 0;
        float Os_down = 0;

        #pragma unroll
        for (int k=0; k<6; k++){
            int i = I_s[k];
            int j = O_s[k];

            num += rho*w[i]*(1 + 
                    a_s2*u_x*c_ix[i] + a_s2*u_y*c_iy[i] + 
                    a_s4*0.5f*mxx*(c_ix[i]*c_ix[i] - inv_as2) + 
                    a_s4*0.5f*mxy*c_ix[i]*c_iy[i] + 
                    a_s4*0.5f*myy*(c_iy[i]*c_iy[i] - inv_as2))*c_ix[i]*c_iy[i];

            div += rho*w[i]*(1 + 
                    a_s2*u_x*c_ix[i] + a_s2*u_y*c_iy[i] + 
                    a_s4*0.5f*mxx*(c_ix[i]*c_ix[i] - inv_as2) + 
                    a_s4*0.5f*mxy*c_ix[i]*c_iy[i] + 
                    a_s4*0.5f*myy*(c_iy[i]*c_iy[i] - inv_as2));


            Is_up += w[i]*c_ix[i]*c_iy[i];

            Is_down += w[i]*a_s4*c_ix[i]*c_ix[i]*c_iy[i]*c_iy[i];

            Os_up += w[i];

            Os_down += w[i]*a_s4*c_ix[i]*c_iy[i];

        }

        float mxy_I = num/div;

        mxy = (Is_up - mxy_I*Os_up)/(mxy_I*(1-omega)*Os_down - Is_down);

    }
    else if (x==Nx-1 && y==0){ //Leste
        constexpr int I_s[6] = {0, 1, 2, 4, 5, 8};
        constexpr int O_s[6] = {0, 2, 3, 4, 6, 7};


        float mxx = 0.f;
        float myy = 0.f;
        float mxy = 0.f;

        float num = 0;
        float div = 0;

        float Is_up = 0;
        float Is_down = 0;
        float Os_up = 0;
        float Os_down = 0;

        #pragma unroll
        for (int k=0; k<6; k++){
            int i = I_s[k];
            int j = O_s[k];

            num += rho*w[i]*(1 + 
                    a_s2*u_x*c_ix[i] + a_s2*u_y*c_iy[i] + 
                    a_s4*0.5f*mxx*(c_ix[i]*c_ix[i] - inv_as2) + 
                    a_s4*0.5f*mxy*c_ix[i]*c_iy[i] + 
                    a_s4*0.5f*myy*(c_iy[i]*c_iy[i] - inv_as2))*c_ix[i]*c_iy[i];

            div += rho*w[i]*(1 + 
                    a_s2*u_x*c_ix[i] + a_s2*u_y*c_iy[i] + 
                    a_s4*0.5f*mxx*(c_ix[i]*c_ix[i] - inv_as2) + 
                    a_s4*0.5f*mxy*c_ix[i]*c_iy[i] + 
                    a_s4*0.5f*myy*(c_iy[i]*c_iy[i] - inv_as2));


            Is_up += w[i]*c_ix[i]*c_iy[i];

            Is_down += w[i]*a_s4*c_ix[i]*c_ix[i]*c_iy[i]*c_iy[i];

            Os_up += w[i];

            Os_down += w[i]*a_s4*c_ix[i]*c_iy[i];

        }

        float mxy_I = num/div;

        mxy = (Is_up - mxy_I*Os_up)/(mxy_I*(1-omega)*Os_down - Is_down);


    }
    else if (x==Nx-1 && y==0){ //Centro

        float mxx_new = 0.f;
        float myy_new = 0.f;
        float mxy_new = 0.f;

        for (int i=0; i<Q; i++){
            float f_i = rho*w[i]*(1 + a_s2*u_x*c_ix[i] + a_s2*u_y*c_iy[i] + a_s4*0.5f*mxx*(c_ix[i]*c_ix[i]-inv_as2) + a_s4*0.5f*mxy*c_ix[i]*c_iy[i] + a_s4*0.5f*myy*(c_ix[i]*c_ix[i]-inv_as2)); 
            mxx_new += f_i*(c_ix[i]*c_ix[i] - inv_as2);
            myy_new += f_i*(c_iy[i]*c_iy[i] - inv_as2);
            mxy_new += f_i*(c_ix[i]*c_iy[i]);
        }

        mxx_new *= 1/rho;
        mxy_new *= 1/rho;
        myy_new *= 1/rho;

        mxx = mxx_new;
        mxy = mxy_new;
        myy = myy_new;

    }

    return mxx, mxy, myy;

}