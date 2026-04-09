#include "../../presets/config.h"
#include "../../presets/physics.h"

#include "../grid_id.cuh"
#include "../from_id.cuh"
#include "../equations/f_i.cuh"

__device__ void wall(CInt *I_s, CInt *O_s, int x, int y,
                     float *mom_in,
                     float *mom_out)
{
    float sum_fi = 0.f;
    float mxy_I = 0.f;
    float rho_I_rho = 0.f;

    float Is_up = 0.f;
    float Is_down = 0.f;
    float Os_up = 0.f;
    float Os_down = 0.f;

    int index = grid_id();

#pragma unroll
    for (int k = 0; k < 6; k++)
    {
        int i = I_s[k];
        int index_from = from_id(x, y, i);

        float fi = f_i(index_from, i, mom_in);
        sum_fi += fi;

        mxy_I += fi * c_ix[i] * c_iy[i];

        Is_up += w[i] * c_ix[i] * c_iy[i];

        Is_down += w[i] * a_s4 * c_ix[i] * c_ix[i] * c_iy[i] * c_iy[i];

        /*---------------------------------------------------------------------*/
        i = O_s[k];
        index_from = from_id(x, y, i);

        rho_I_rho += w[i] +
                     w[i] * (1.f - omega) * a_s4 * mom_in[index_from + 5] * (c_ix[i] * c_iy[i]); // mxy

        Os_up += w[i];

        Os_down += w[i] * a_s4 * c_ix[i] * c_iy[i];
    }
    mxy_I /= sum_fi;

    mom_out[index] = sum_fi / rho_I_rho;                                                      // rho
    mom_out[index + 5] = (Is_up - mxy_I * Os_up) / (mxy_I * (1 - omega) * Os_down - Is_down); // mxy
}