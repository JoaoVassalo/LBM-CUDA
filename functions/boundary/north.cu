#include "../equations/f_i.cuh"
#include "../../presets/config.h"
#include "../../presets/physics.h"
#include "../from_id.cuh"
#include "../grid_id.cuh"

__device__ void north(int size, CInt *I_s, CInt *O_s, int x, int y, float *rho_in, float *ux_in, float *uy_in, float *mxx_in, float *mxy_in, float *myy_in,
                      float *rho_out, float *mxy_out)
{
    float sum_fi = 0;
    float mxy_I = 0;
    float rho_I_rho = 0;

    float Is_up = 0;
    float Is_down = 0;
    float Os_up = 0;
    float Os_down = 0;

    int index = grid_id(x, y);

#pragma unroll
    for (int k = 0; k < size; k++)
    {
        int i = I_s[k];
        int index_from = from_id(x, y, i);

        float fi = f_i(index_from, i, rho_in, ux_in, uy_in, mxx_in, mxy_in, myy_in);
        sum_fi += fi;

        mxy_I += fi * c_ix[i] * c_iy[i];

        Is_up += w[i] *
                 (1 +
                  a_s2 * u_max * c_ix[i] +
                  a_s4 * 0.5f * u_max * (c_ix[i] * c_ix[i] - inv_as2)) *
                 c_ix[i] * c_iy[i];

        Is_down += w[i] * a_s4 * c_ix[i] * c_ix[i] * c_iy[i] * c_iy[i];

        /*---------------------------------------------------------------------*/
        i = O_s[k];
        int index_from = from_id(x, y, i);

        rho_I_rho += w[i] *
                         (1 +
                          a_s2 * u_max * c_ix[i] +
                          a_s4 * 0.5f * u_max * u_max * (c_ix[i] * c_ix[i] - inv_as2)) +
                     w[i] * (1 - omega) * a_s4 * mxy_in[index_from] * (c_ix[i] * c_iy[i]);

        Os_up += w[i] *
                 (1 +
                  a_s2 * u_max * c_ix[i] +
                  a_s4 * 0.5f * u_max * (c_ix[i] * c_ix[i] - inv_as2));

        Os_down += w[i] * a_s4 * c_ix[i] * c_ix[i] * c_iy[i] * c_iy[i];
    }
    mxy_I /= sum_fi;

    rho_out[index] = sum_fi / rho_I_rho;
    mxy_out[index] = (Is_up - mxy_I * Os_up) / (mxy_I * (1 - omega) * Os_down - Is_down);
}