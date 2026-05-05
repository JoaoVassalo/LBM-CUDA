#include "../config/config.h"
#include "../config/stencil.cuh"

#include "../lbm/f_i.cuh"
#include "../core/grid_id.cuh"
#include "../core/from_id.cuh"

#include "../config/var.cuh"

__device__ void center(int x, int y,
                       float *mom_in,
                       float *mom_out)
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

        float fi = f_i(index_from, i, mom_in);

        rho += fi;

        ux += fi * (float)c_ix[i];
        uy += fi * (float)c_iy[i];

        mxx += fi * (c_ix[i] * c_ix[i] - inv_as2);
        mxy += fi * (c_ix[i] * c_iy[i]);
        myy += fi * (c_iy[i] * c_iy[i] - inv_as2);
    }
    mom_out[momIdx<MomentId::rho>(index)] = rho;

    mom_out[momIdx<MomentId::ux>(index)] = ux / rho;
    mom_out[momIdx<MomentId::uy>(index)] = uy / rho;

    mom_out[momIdx<MomentId::mxx>(index)] = mxx / rho;
    mom_out[momIdx<MomentId::myy>(index)] = myy / rho;
    mom_out[momIdx<MomentId::mxy>(index)] = mxy / rho;
}