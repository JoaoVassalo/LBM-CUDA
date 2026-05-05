#include "../equations/f_i.cuh"
#include "../../presets/config.h"
#include "../from_id.cuh"
#include "../grid_id.cuh"

#include <iostream>

__device__ void corner(CInt *I_s, CInt *O_s, int x, int y, float *mom_in,
                       float *mom_out)
{
    int index = grid_id();

    float rho_I = 0.f;
    float sum_wi = 0.f;

#pragma unroll
    for (int k = 0; k < 4; k++)
    {
        int i = I_s[k];

        int index_from = from_id(x, y, i);

        rho_I += f_i(index_from, i, mom_in);

        /*---------------------------------------------------------------------*/
        i = O_s[k];

        sum_wi += w[i];
    }

    mom_out[momIdx<MomentId::rho>(index)] = rho_I / sum_wi;
}