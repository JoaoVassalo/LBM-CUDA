#include "calc_density.cuh"

#include "../presets/stencil.cuh"
#include "../presets/geometry.h"
#include "../presets/physics.h"

#include "grid_id.cuh"
#include "from_id.cuh"

__device__ void calc_density(int x, int y, float *rho_in, float *ux_in, float *uy_in, float *mxx_in, float *mxy_in, float *myy_in,
                             float *rho_out)
{

    int index = grid_id(x, y);

    if (x == 0 && y == 0)
    { // Sudoeste
        constexpr int I_s[4] = {0, 3, 4, 7};

        float density_I = 0.f;

#pragma unroll
        for (int k = 0; k < 4; k++)
        {
            int i = I_s[k];

            int index_from = from_id(x, y, i);

            density_I += rho_in[index_from] * w[i] * (1 + a_s2 * ux_in[index_from] * c_ix[i] + a_s2 * uy_in[index_from] * c_iy[i] + a_s4 * 0.5f * mxx_in[index_from] * (c_ix[i] * c_ix[i] - inv_as2) + a_s4 * mxy_in[index_from] * c_ix[i] * c_iy[i] + a_s4 * 0.5f * myy_in[index_from] * (c_iy[i] * c_iy[i] - inv_as2));
        }

        float sum_wi = w[0] + w[1] + w[2] + w[5];

        rho_out[index] = density_I / sum_wi;
    }
    else if (x == Nx - 1 && y == 0)
    { // Sudeste
        constexpr int I_s[4] = {0, 1, 4, 8};

        float density_I = 0.f;

#pragma unroll
        for (int k = 0; k < 4; k++)
        {
            int i = I_s[k];

            int index_from = from_id(x, y, i);

            density_I += rho_in[index_from] * w[i] * (1 + a_s2 * ux_in[index_from] * c_ix[i] + a_s2 * uy_in[index_from] * c_iy[i] + a_s4 * 0.5f * mxx_in[index_from] * (c_ix[i] * c_ix[i] - inv_as2) + a_s4 * mxy_in[index_from] * c_ix[i] * c_iy[i] + a_s4 * 0.5f * myy_in[index_from] * (c_iy[i] * c_iy[i] - inv_as2));
        }

        float sum_wi = w[0] + w[2] + w[3] + w[6];

        rho_out[index] = density_I / sum_wi;
    }
    else if (x == 0 && y == Ny - 1)
    { // Noroeste
        constexpr int I_s[4] = {0, 2, 3, 6};
        constexpr int O_s[4] = {0, 1, 4, 8};

        float density_I_density = 0.f;
        float density_I = 0.f;

#pragma unroll
        for (int k = 0; k < 4; k++)
        {
            int i = O_s[k];

            int index_from = from_id(x, y, i);

            density_I_density += w[i] * (1 +
                                         a_s2 * ux_in[index] * c_ix[i] +
                                         a_s4 * 0.5f * ux_in[index] * ux_in[index] * (c_ix[i] * c_ix[i] - inv_as2)) +
                                 w[i] * (1 - omega) * a_s4 * mxy_in[index] * c_ix[i] * c_iy[i];

            i = I_s[k];

            density_I += rho_in[index] * w[i] * (1 + a_s2 * ux_in[index] * c_ix[i] + a_s2 * uy_in[index] * c_iy[i] + a_s4 * 0.5f * mxx_in[index] * (c_ix[i] * c_ix[i] - inv_as2) + a_s4 * mxy_in[index] * c_ix[i] * c_iy[i] + a_s4 * 0.5f * myy_in[index] * (c_iy[i] * c_iy[i] - inv_as2));
        }

        rho_out[index] = density_I / density_I_density;
    }
    else if (x == Nx - 1 && y == Ny - 1)
    { // Nordeste
        constexpr int I_s[4] = {0, 1, 2, 5};
        constexpr int O_s[4] = {0, 3, 4, 7};

        float density_I_density = 0.f;
        float density_I = 0.f;

#pragma unroll
        for (int k = 0; k < 4; k++)
        {
            int i = O_s[k];

            int index_from = from_id(x, y, i);

            density_I_density += w[i] * (1 +
                                         a_s2 * ux_in[index_from] * c_ix[i] +
                                         a_s4 * 0.5f * ux_in[index_from] * ux_in[index_from] * (c_ix[i] * c_ix[i] - inv_as2)) +
                                 w[i] * (1 - omega) * a_s4 * mxy_in[index_from] * c_ix[i] * c_iy[i];

            i = I_s[k];
            index_from = from_id(x, y, i);

            density_I += rho_in[index_from] * w[i] * (1 + a_s2 * ux_in[index_from] * c_ix[i] + a_s2 * uy_in[index_from] * c_iy[i] + a_s4 * 0.5f * mxx_in[index_from] * (c_ix[i] * c_ix[i] - inv_as2) + a_s4 * mxy_in[index_from] * c_ix[i] * c_iy[i] + a_s4 * 0.5f * myy_in[index_from] * (c_iy[i] * c_iy[i] - inv_as2));
        }

        rho_out[index] = density_I / density_I_density;
    }
    else if (y == 0)
    { // Sul
        constexpr int I_s[6] = {0, 1, 3, 4, 7, 8};
        constexpr int O_s[6] = {0, 1, 2, 3, 5, 6};

        float density_I_density = 0.f;
        float density_I = 0.f;

#pragma unroll
        for (int k = 0; k < 6; k++)
        {
            int i = O_s[k];
            int index_from = from_id(x, y, i);

            density_I_density += w[i] +
                                 w[i] * (1 - omega) * a_s4 * mxy_in[index_from] * c_ix[i] * c_iy[i];

            i = I_s[k];
            index_from = from_id(x, y, i);

            density_I += rho_in[index_from] * w[i] * (1 + a_s2 * ux_in[index_from] * c_ix[i] + a_s2 * uy_in[index_from] * c_iy[i] + a_s4 * 0.5f * mxx_in[index_from] * (c_ix[i] * c_ix[i] - inv_as2) + a_s4 * mxy_in[index_from] * c_ix[i] * c_iy[i] + a_s4 * 0.5f * myy_in[index_from] * (c_iy[i] * c_iy[i] - inv_as2));
        }

        rho_out[index] = density_I / density_I_density;
    }
    else if (y == Ny - 1)
    { // Norte
        constexpr int I_s[6] = {0, 1, 2, 3, 5, 6};
        constexpr int O_s[6] = {0, 1, 3, 4, 7, 8};

        float density_I_density = 0.f;
        float density_I = 0.f;

#pragma unroll
        for (int k = 0; k < 6; k++)
        {
            int i = O_s[k];

            int index_from = from_id(x, y, i);

            density_I_density += w[i] * (1 +
                                         a_s2 * ux_in[index_from] * c_ix[i] +
                                         a_s4 * 0.5f * ux_in[index_from] * ux_in[index_from] * (c_ix[i] * c_ix[i] - inv_as2)) +
                                 w[i] * (1 - omega) * a_s4 * mxy_in[index_from] * c_ix[i] * c_iy[i];
            i = I_s[k];

            index_from = from_id(x, y, i);

            density_I += rho_in[index_from] * w[i] * (1 + a_s2 * ux_in[index_from] * c_ix[i] + a_s2 * uy_in[index_from] * c_iy[i] + a_s4 * 0.5f * mxx_in[index_from] * (c_ix[i] * c_ix[i] - inv_as2) + a_s4 * mxy_in[index_from] * c_ix[i] * c_iy[i] + a_s4 * 0.5f * myy_in[index_from] * (c_iy[i] * c_iy[i] - inv_as2));
        }

        rho_out[index] = density_I / density_I_density;
    }
    else if (x == 0)
    { // Oeste
        constexpr int I_s[6] = {0, 2, 3, 4, 6, 7};
        constexpr int O_s[6] = {0, 1, 2, 4, 5, 8};

        float density_I_density = 0.f;
        float density_I = 0.f;

#pragma unroll
        for (int k = 0; k < 6; k++)
        {
            int i = O_s[k];

            int index_from = from_id(x, y, i);

            density_I_density += w[i] +
                                 w[i] * (1 - omega) * a_s4 * mxy_in[index_from] * c_ix[i] * c_iy[i];

            i = I_s[k];

            index_from = from_id(x, y, i);

            density_I += rho_in[index_from] * w[i] * (1 + a_s2 * ux_in[index_from] * c_ix[i] + a_s2 * uy_in[index_from] * c_iy[i] + a_s4 * 0.5f * mxx_in[index_from] * (c_ix[i] * c_ix[i] - inv_as2) + a_s4 * mxy_in[index_from] * c_ix[i] * c_iy[i] + a_s4 * 0.5f * myy_in[index_from] * (c_iy[i] * c_iy[i] - inv_as2));
        }

        rho_out[index] = density_I / density_I_density;
    }
    else if (x == Nx - 1)
    { // Leste
        constexpr int I_s[6] = {0, 1, 2, 4, 5, 8};
        constexpr int O_s[6] = {0, 2, 3, 4, 6, 7};

        float density_I_density = 0.f;
        float density_I = 0.f;

#pragma unroll
        for (int k = 0; k < 6; k++)
        {
            int i = O_s[k];
            int index_from = from_id(x, y, i);

            density_I_density += w[i] +
                                 w[i] * (1 - omega) * a_s4 * mxy_in[index_from] * c_ix[i] * c_iy[i];

            i = I_s[k];
            index_from = from_id(x, y, i);

            density_I += rho_in[index_from] * w[i] * (1 + a_s2 * ux_in[index_from] * c_ix[i] + a_s2 * uy_in[index_from] * c_iy[i] + a_s4 * 0.5f * mxx_in[index_from] * (c_ix[i] * c_ix[i] - inv_as2) + a_s4 * mxy_in[index_from] * c_ix[i] * c_iy[i] + a_s4 * 0.5f * myy_in[index_from] * (c_iy[i] * c_iy[i] - inv_as2));
        }

        rho_out[index] = density_I / density_I_density;
    }
    else
    { // Centro

        float density = 0.f;
#pragma unroll
        for (int i = 0; i < Q; i++)
        {
            int index_from = from_id(x, y, i);

            density += rho_in[index_from] * w[i] *
                       (1 +
                        a_s2 * ux_in[index_from] * c_ix[i] + a_s2 * uy_in[index_from] * c_iy[i] +
                        a_s4 * 0.5f * mxx_in[index_from] * (c_ix[i] * c_ix[i] - inv_as2) + a_s4 * mxy_in[index_from] * c_ix[i] * c_iy[i] +
                        a_s4 * 0.5f * myy_in[index_from] * (c_iy[i] * c_iy[i] - inv_as2));
        }

        rho_out[index] = density;
    }
}