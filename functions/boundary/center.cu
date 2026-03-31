#include "../../presets/config.h"
#include "../../presets/stencil.cuh"

#include "../equations/f_i.cuh"
#include "../grid_id.cuh"
#include "../from_id.cuh"

__device__ void center(int x, int y,
                       float *rho_in, float *ux_in, float *uy_in, float *mxx_in, float *mxy_in, float *myy_in,
                       float *rho_out, float *ux_out, float *uy_out, float *mxx_out, float *mxy_out, float *myy_out)
{
    int index = grid_id();
    float rho = 0.f;
    float ux = 0.f;
    float uy = 0.f;
    float mxx = 0.f;
    float mxy = 0.f;
    float myy = 0.f;

#pragma unroll
    for (int i = 0; i < Q; i++)
    {
        int index_from = from_id(x, y, i);

        float fi = f_i(index_from, i, rho_in, ux_in, uy_in, mxx_in, mxy_in, myy_in);

        rho += fi;

        ux += fi * (float)c_ix[i];
        uy += fi * (float)c_iy[i];

        mxx += fi * (c_ix[i] * c_ix[i] - inv_as2);
        mxy += fi * (c_ix[i] * c_iy[i]);
        myy += fi * (c_iy[i] * c_iy[i] - inv_as2);
    }
    rho_out[index] = rho;

    ux_out[index] = ux / rho;
    uy_out[index] = uy / rho;

    mxx_out[index] = mxx / rho;
    mxy_out[index] = mxy / rho;
    myy_out[index] = myy / rho;
}