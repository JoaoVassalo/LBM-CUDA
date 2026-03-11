#include "calc_momentum.cuh"

#include "../presets/stencil.cuh"
#include "../presets/geometry.h"
#include "../presets/physics.h"

#include "grid_id.cuh"
#include "from_id.cuh"

__device__ void calc_momentum(int x, int y, float *rho_in, float *ux_in, float *uy_in, float *mxx_in, float *mxy_in, float *myy_in,
                              float *rho_out, float *ux_out, float *uy_out, float *mxx_out, float *mxy_out, float *myy_out)
{

    int index = grid_id(x, y);

    if (x == 0 && y == 0)
    { // Sudoeste

        mxx_out[index] = 0.f;
        myy_out[index] = 0.f;
        mxy_out[index] = 0.f;
    }
    else if (x == Nx - 1 && y == 0)
    { // Sudeste

        mxx_out[index] = 0.f;
        myy_out[index] = 0.f;
        mxy_out[index] = 0.f;
    }
    else if (x == 0 && y == Ny - 1)
    { // Noroeste
        constexpr int I_s[4] = {0, 2, 3, 6};
        constexpr int O_s[4] = {0, 1, 4, 8};

        mxx_out[index] = 0.f;
        myy_out[index] = u_max * u_max;
        mxy_out[index] = 0.f;

        float num = 0;
        float div = 0;

        float Is_up = 0;
        float Is_down = 0;
        float Os_up = 0;
        float Os_down = 0;

#pragma unroll
        for (int k = 0; k < 4; k++)
        {
            int i = I_s[k];
            int j = O_s[k];

            int index_from = from_id(x, y, i);

            num += rho_in[index_from] * w[i] *
                   (1 + a_s2 * ux_in[index_from] * c_ix[i] + a_s2 * uy_in[index_from] * c_iy[i] +
                    a_s4 * 0.5f * mxx_in[index_from] * (c_ix[i] * c_ix[i] - inv_as2) + a_s4 * mxy_in[index_from] * c_ix[i] * c_iy[i] +
                    a_s4 * 0.5f * myy_in[index_from] * (c_iy[i] * c_iy[i] - inv_as2)) *
                   c_ix[i] * c_iy[i];

            div += rho_in[index_from] * w[i] *
                   (1 + a_s2 * ux_in[index_from] * c_ix[i] + a_s2 * uy_in[index_from] * c_iy[i] +
                    a_s4 * 0.5f * mxx_in[index_from] * (c_ix[i] * c_ix[i] - inv_as2) + a_s4 * mxy_in[index_from] * c_ix[i] * c_iy[i] +
                    a_s4 * 0.5f * myy_in[index_from] * (c_iy[i] * c_iy[i] - inv_as2));

            Is_up += w[i] * c_ix[i] * c_iy[i] * (1 + a_s2 * u_max * c_ix[i] + a_s4 * 0.5f * u_max * (c_ix[i] * c_ix[i] - inv_as2));

            Is_down += w[i] * a_s4 * c_ix[i] * c_ix[i] * c_iy[i] * c_iy[i];

            Os_up += w[j] * (1 + a_s2 * u_max * c_ix[j] + a_s4 * 0.5f * u_max * u_max * (c_ix[j] * c_ix[j] - inv_as2));

            Os_down += w[j] * a_s4 * c_ix[j] * c_iy[j];
        }

        float mxy_I = num / div;

        mxy_out[index] = (Is_up - mxy_I * Os_up) / (mxy_I * (1 - omega) * Os_down - Is_down);
    }
    else if (x == Nx - 1 && y == Ny - 1)
    { // Nordeste
        constexpr int I_s[4] = {0, 1, 2, 5};
        constexpr int O_s[4] = {0, 3, 4, 7};

        mxx_out[index] = 0.f;
        myy_out[index] = u_max * u_max;
        mxy_out[index] = 0.f;

        float num = 0;
        float div = 0;

        float Is_up = 0;
        float Is_down = 0;
        float Os_up = 0;
        float Os_down = 0;

#pragma unroll
        for (int k = 0; k < 4; k++)
        {
            int i = I_s[k];
            int j = O_s[k];

            int index_from = from_id(x, y, i);

            num += rho_in[index_from] * w[i] *
                   (1 + a_s2 * ux_in[index_from] * c_ix[i] + a_s2 * uy_in[index_from] * c_iy[i] +
                    a_s4 * 0.5f * mxx_in[index_from] * (c_ix[i] * c_ix[i] - inv_as2) + a_s4 * mxy_in[index_from] * c_ix[i] * c_iy[i] +
                    a_s4 * 0.5f * myy_in[index_from] * (c_iy[i] * c_iy[i] - inv_as2)) *
                   c_ix[i] * c_iy[i];

            div += rho_in[index_from] * w[i] *
                   (1 + a_s2 * ux_in[index_from] * c_ix[i] + a_s2 * uy_in[index_from] * c_iy[i] +
                    a_s4 * 0.5f * mxx_in[index_from] * (c_ix[i] * c_ix[i] - inv_as2) + a_s4 * mxy_in[index_from] * c_ix[i] * c_iy[i] +
                    a_s4 * 0.5f * myy_in[index_from] * (c_iy[i] * c_iy[i] - inv_as2));

            Is_up += w[i] * c_ix[i] * c_iy[i] * (1 + a_s2 * u_max * c_ix[i] + a_s4 * 0.5f * u_max * (c_ix[i] * c_ix[i] - inv_as2));

            Is_down += w[i] * a_s4 * c_ix[i] * c_ix[i] * c_iy[i] * c_iy[i];

            Os_up += w[j] * (1 + a_s2 * u_max * c_ix[j] + a_s4 * 0.5f * u_max * u_max * (c_ix[j] * c_ix[j] - inv_as2));

            Os_down += w[j] * a_s4 * c_ix[j] * c_iy[j];
        }

        float mxy_I = num / div;

        mxy_out[index] = (Is_up - mxy_I * Os_up) / (mxy_I * (1 - omega) * Os_down - Is_down);
    }
    else if (y == 0)
    { // Sul
        constexpr int I_s[6] = {0, 1, 3, 4, 7, 8};
        constexpr int O_s[6] = {0, 1, 2, 3, 5, 6};

        mxx_out[index] = 0.f;
        myy_out[index] = 0.f;
        mxy_out[index] = 0.f;

        float num = 0;
        float div = 0;

        float Is_up = 0;
        float Is_down = 0;
        float Os_up = 0;
        float Os_down = 0;

#pragma unroll
        for (int k = 0; k < 6; k++)
        {
            int i = I_s[k];
            int j = O_s[k];

            int index_from = from_id(x, y, i);

            num += rho_in[index_from] * w[i] *
                   (1 + a_s2 * ux_in[index_from] * c_ix[i] + a_s2 * uy_in[index_from] * c_iy[i] +
                    a_s4 * 0.5f * mxx_in[index_from] * (c_ix[i] * c_ix[i] - inv_as2) + a_s4 * mxy_in[index_from] * c_ix[i] * c_iy[i] +
                    a_s4 * 0.5f * myy_in[index_from] * (c_iy[i] * c_iy[i] - inv_as2)) *
                   c_ix[i] * c_iy[i];

            div += rho_in[index_from] * w[i] *
                   (1 + a_s2 * ux_in[index_from] * c_ix[i] + a_s2 * uy_in[index_from] * c_iy[i] +
                    a_s4 * 0.5f * mxx_in[index_from] * (c_ix[i] * c_ix[i] - inv_as2) + a_s4 * mxy_in[index_from] * c_ix[i] * c_iy[i] +
                    a_s4 * 0.5f * myy_in[index_from] * (c_iy[i] * c_iy[i] - inv_as2));

            Is_up += w[i] * c_ix[i] * c_iy[i];

            Is_down += w[i] * a_s4 * c_ix[i] * c_ix[i] * c_iy[i] * c_iy[i];

            Os_up += w[j];

            Os_down += w[j] * a_s4 * c_ix[j] * c_iy[j];
        }

        float mxy_I = num / div;

        mxy_out[index] = (Is_up - mxy_I * Os_up) / (mxy_I * (1 - omega) * Os_down - Is_down);
    }
    else if (y == Ny - 1)
    { // Norte
        constexpr int I_s[6] = {0, 1, 2, 3, 5, 6};
        constexpr int O_s[6] = {0, 1, 3, 4, 7, 8};

        mxx_out[index] = 0.f;
        myy_out[index] = u_max * u_max;
        mxy_out[index] = 0.f;

        float num = 0;
        float div = 0;

        float Is_up = 0;
        float Is_down = 0;
        float Os_up = 0;
        float Os_down = 0;

#pragma unroll
        for (int k = 0; k < 6; k++)
        {
            int i = I_s[k];
            int j = O_s[k];

            int index_from = from_id(x, y, i);

            num += rho_in[index_from] * w[i] *
                   (1 + a_s2 * ux_in[index_from] * c_ix[i] + a_s2 * uy_in[index_from] * c_iy[i] +
                    a_s4 * 0.5f * mxx_in[index_from] * (c_ix[i] * c_ix[i] - inv_as2) + a_s4 * mxy_in[index_from] * c_ix[i] * c_iy[i] +
                    a_s4 * 0.5f * myy_in[index_from] * (c_iy[i] * c_iy[i] - inv_as2)) *
                   c_ix[i] * c_iy[i];

            div += rho_in[index_from] * w[i] *
                   (1 + a_s2 * ux_in[index_from] * c_ix[i] + a_s2 * uy_in[index_from] * c_iy[i] +
                    a_s4 * 0.5f * mxx_in[index_from] * (c_ix[i] * c_ix[i] - inv_as2) + a_s4 * mxy_in[index_from] * c_ix[i] * c_iy[i] +
                    a_s4 * 0.5f * myy_in[index_from] * (c_iy[i] * c_iy[i] - inv_as2));

            Is_up += w[i] * c_ix[i] * c_iy[i] * (1 + a_s2 * u_max * c_ix[i] + a_s4 * 0.5f * u_max * (c_ix[i] * c_ix[i] - inv_as2));

            Is_down += w[i] * a_s4 * c_ix[i] * c_ix[i] * c_iy[i] * c_iy[i];

            Os_up += w[j] * (1 + a_s2 * u_max * c_ix[j] + a_s4 * 0.5f * u_max * u_max * (c_ix[j] * c_ix[j] - inv_as2));

            Os_down += w[j] * a_s4 * c_ix[j] * c_iy[j];
        }

        float mxy_I = num / div;

        mxy_out[index] = (Is_up - mxy_I * Os_up) / (mxy_I * (1 - omega) * Os_down - Is_down);
    }
    else if (x == 0)
    { // Oeste
        constexpr int I_s[6] = {0, 2, 3, 4, 6, 7};
        constexpr int O_s[6] = {0, 1, 2, 4, 5, 8};

        mxx_out[index] = 0.f;
        myy_out[index] = 0.f;
        mxy_out[index] = 0.f;

        float num = 0;
        float div = 0;

        float Is_up = 0;
        float Is_down = 0;
        float Os_up = 0;
        float Os_down = 0;

#pragma unroll
        for (int k = 0; k < 6; k++)
        {
            int i = I_s[k];
            int j = O_s[k];

            int index_from = from_id(x, y, i);

            num += rho_in[index_from] * w[i] *
                   (1 + a_s2 * ux_in[index_from] * c_ix[i] + a_s2 * uy_in[index_from] * c_iy[i] +
                    a_s4 * 0.5f * mxx_in[index_from] * (c_ix[i] * c_ix[i] - inv_as2) + a_s4 * mxy_in[index_from] * c_ix[i] * c_iy[i] +
                    a_s4 * 0.5f * myy_in[index_from] * (c_iy[i] * c_iy[i] - inv_as2)) *
                   c_ix[i] * c_iy[i];

            div += rho_in[index_from] * w[i] *
                   (1 + a_s2 * ux_in[index_from] * c_ix[i] + a_s2 * uy_in[index_from] * c_iy[i] +
                    a_s4 * 0.5f * mxx_in[index_from] * (c_ix[i] * c_ix[i] - inv_as2) + a_s4 * mxy_in[index_from] * c_ix[i] * c_iy[i] +
                    a_s4 * 0.5f * myy_in[index_from] * (c_iy[i] * c_iy[i] - inv_as2));

            Is_up += w[i] * c_ix[i] * c_iy[i];

            Is_down += w[i] * a_s4 * c_ix[i] * c_ix[i] * c_iy[i] * c_iy[i];

            Os_up += w[j];

            Os_down += w[j] * a_s4 * c_ix[j] * c_iy[j];
        }

        float mxy_I = num / div;

        mxy_out[index] = (Is_up - mxy_I * Os_up) / (mxy_I * (1 - omega) * Os_down - Is_down);
    }
    else if (x == Nx - 1)
    { // Leste
        constexpr int I_s[6] = {0, 1, 2, 4, 5, 8};
        constexpr int O_s[6] = {0, 2, 3, 4, 6, 7};

        mxx_out[index] = 0.f;
        myy_out[index] = 0.f;
        mxy_out[index] = 0.f;

        float num = 0;
        float div = 0;

        float Is_up = 0;
        float Is_down = 0;
        float Os_up = 0;
        float Os_down = 0;

#pragma unroll
        for (int k = 0; k < 6; k++)
        {
            int i = I_s[k];
            int j = O_s[k];

            int index_from = from_id(x, y, i);

            num += rho_in[index_from] * w[i] *
                   (1 + a_s2 * ux_in[index_from] * c_ix[i] + a_s2 * uy_in[index_from] * c_iy[i] +
                    a_s4 * 0.5f * mxx_in[index_from] * (c_ix[i] * c_ix[i] - inv_as2) + a_s4 * mxy_in[index_from] * c_ix[i] * c_iy[i] +
                    a_s4 * 0.5f * myy_in[index_from] * (c_iy[i] * c_iy[i] - inv_as2)) *
                   c_ix[i] * c_iy[i];

            div += rho_in[index_from] * w[i] *
                   (1 + a_s2 * ux_in[index_from] * c_ix[i] + a_s2 * uy_in[index_from] * c_iy[i] +
                    a_s4 * 0.5f * mxx_in[index_from] * (c_ix[i] * c_ix[i] - inv_as2) + a_s4 * mxy_in[index_from] * c_ix[i] * c_iy[i] +
                    a_s4 * 0.5f * myy_in[index_from] * (c_iy[i] * c_iy[i] - inv_as2));

            Is_up += w[i] * c_ix[i] * c_iy[i];

            Is_down += w[i] * a_s4 * c_ix[i] * c_ix[i] * c_iy[i] * c_iy[i];

            Os_up += w[j];

            Os_down += w[j] * a_s4 * c_ix[j] * c_iy[j];
        }

        float mxy_I = num / div;

        mxy_out[index] = (Is_up - mxy_I * Os_up) / (mxy_I * (1 - omega) * Os_down - Is_down);
    }
    else
    { // Centro

        float mxx_new = 0.f;
        float myy_new = 0.f;
        float mxy_new = 0.f;

        for (int i = 0; i < Q; i++)
        {

            int index_from = from_id(x, y, i);

            float f_i = rho_in[index_from] * w[i] *
                        (1 + a_s2 * ux_in[index_from] * c_ix[i] + a_s2 * uy_in[index_from] * c_iy[i] +
                         a_s4 * 0.5f * mxx_in[index_from] * (c_ix[i] * c_ix[i] - inv_as2) + a_s4 * mxy_in[index_from] * c_ix[i] * c_iy[i] +
                         a_s4 * 0.5f * myy_in[index_from] * (c_ix[i] * c_ix[i] - inv_as2));

            mxx_new += f_i * (c_ix[i] * c_ix[i] - inv_as2);
            myy_new += f_i * (c_iy[i] * c_iy[i] - inv_as2);
            mxy_new += f_i * (c_ix[i] * c_iy[i]);
        }

        mxx_new /= rho_out[index];
        mxy_new /= rho_out[index];
        myy_new /= rho_out[index];

        mxx_out[index] = mxx_new;
        mxy_out[index] = mxy_new;
        myy_out[index] = myy_new;
    }
}