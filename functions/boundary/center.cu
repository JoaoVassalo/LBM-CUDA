#include "../../presets/config.h"
#include "../../presets/stencil.cuh"

#include "../equations/f_i.cuh"
#include "../grid_id.cuh"
#include "../from_id.cuh"

#include "../../cuda_config/var.cuh"

__device__ void center(int x, int y,
                       varUnit *mom_in,
                       varUnit *mom_out)
{
    int index = grid_id();
    varUnit rho = 0.f;
    varUnit ux = 0.f;
    varUnit uy = 0.f;
    varUnit mxx = 0.f;
    varUnit mxy = 0.f;
    varUnit myy = 0.f;

#pragma unroll
    for (int i = 0; i < Q; i++)
    {
        int index_from = from_id(x, y, i);

        varUnit fi = f_i(index_from, i, mom_in);

        rho += fi;

        ux += fi * (varUnit)c_ix[i];
        uy += fi * (varUnit)c_iy[i];

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