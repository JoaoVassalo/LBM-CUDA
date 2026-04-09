#include "../equations/f_i.cuh"
#include "../../presets/config.h"
#include "../../presets/physics.h"
#include "../from_id.cuh"
#include "../grid_id.cuh"

__device__ void north(int size, CInt *I_s, CInt *O_s, int x, int y,
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
    for (int k = 0; k < size; k++)
    {
        int i = I_s[k];
        int index_from = from_id(x, y, i);

        float fi = f_i(index_from, i, mom_in);
        sum_fi += fi;

        mxy_I += fi * c_ix[i] * c_iy[i];

        Is_up += w[i] *
                 (1.f +
                  a_s2 * u_max * c_ix[i] +
                  a_s4 * 0.5f * u_max * u_max * (c_ix[i] * c_ix[i] - inv_as2)) *
                 c_ix[i] * c_iy[i];

        Is_down += w[i] * a_s4 * c_ix[i] * c_ix[i] * c_iy[i] * c_iy[i];

        /*---------------------------------------------------------------------*/
        i = O_s[k];
        index_from = from_id(x, y, i);

        rho_I_rho += w[i] *
                         (1.f +
                          a_s2 * u_max * c_ix[i] +
                          a_s4 * 0.5f * u_max * u_max * (c_ix[i] * c_ix[i] - inv_as2)) +
                     w[i] * (1.f - omega) * a_s4 * mom_in[momIdx<MomentId::mxy>(index)] * (c_ix[i] * c_iy[i]); // mxy

        Os_up += w[i] *
                 (1.f +
                  a_s2 * u_max * c_ix[i] +
                  a_s4 * 0.5f * u_max * u_max * (c_ix[i] * c_ix[i] - inv_as2));

        Os_down += w[i] * a_s4 * c_ix[i] * c_iy[i];
    }
    mxy_I /= sum_fi;

    mom_out[momIdx<MomentId::rho>(index)] = sum_fi / rho_I_rho;                                                    // rho
    mom_out[momIdx<MomentId::mxy>(index)] = (Is_up - mxy_I * Os_up) / (mxy_I * (1.f - omega) * Os_down - Is_down); // mxy
}