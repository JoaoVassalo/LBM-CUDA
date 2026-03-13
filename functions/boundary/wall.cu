#include "../../presets/config.h"
#include "../../presets/physics.h"

#include "../grid_id.cuh"
#include "../from_id.cuh"
#include "../equations/f_i.cuh"

__device__ void wall(CInt *I_s, CInt *O_s, int x, int y, float *rho_in, float *ux_in, float *uy_in, float *mxx_in, float *mxy_in, float *myy_in,
                     float *rho_out, float *mxy_out)
{
    float sum_fi = 0.f;
    float mxy_I = 0.f;
    float rho_I_rho = 0.f;

    float Is_up = 0.f;
    float Is_down = 0.f;
    float Os_up = 0.f;
    float Os_down = 0.f;

    int index = grid_id(x, y);

#pragma unroll
    for (int k = 0; k < 6; k++)
    {
        int i = I_s[k];
        int index_from = from_id(x, y, i);

        float fi = f_i(index_from, i, rho_in, ux_in, uy_in, mxx_in, mxy_in, myy_in);
        sum_fi += fi;

        mxy_I += fi * c_ix[i] * c_iy[i];

        Is_up += w[i] * c_ix[i] * c_iy[i];

        Is_down += w[i] * a_s4 * c_ix[i] * c_ix[i] * c_iy[i] * c_iy[i];

        /*---------------------------------------------------------------------*/
        i = O_s[k];
        index_from = from_id(x, y, i);

        rho_I_rho += w[i] +
                     w[i] * (1.f - omega) * a_s4 * mxy_in[index_from] * (c_ix[i] * c_iy[i]);

        Os_up += w[i];

        Os_down += w[i] * a_s4 * c_ix[i] * c_iy[i];
    }
    mxy_I /= sum_fi;

    rho_out[index] = sum_fi / rho_I_rho;
    mxy_out[index] = (Is_up - mxy_I * Os_up) / (mxy_I * (1 - omega) * Os_down - Is_down);
}