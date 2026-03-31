#include "../equations/f_i.cuh"
#include "../../presets/config.h"
#include "../from_id.cuh"
#include "../grid_id.cuh"

#include <iostream>

__device__ void corner(CInt *I_s, CInt *O_s, int x, int y, float *rho_in, float *ux_in, float *uy_in, float *mxx_in, float *mxy_in, float *myy_in,
                       float *rho_out)
{
    int index = grid_id();

    float rho_I = 0.f;
    float sum_wi = 0.f;

#pragma unroll
    for (int k = 0; k < 4; k++)
    {
        int i = I_s[k];

        int index_from = from_id(x, y, i);

        rho_I += f_i(index_from, i, rho_in, ux_in, uy_in, mxx_in, mxy_in, myy_in);

        /*---------------------------------------------------------------------*/
        i = O_s[k];

        sum_wi += w[i];
    }

    rho_out[index] = rho_I / sum_wi;
}