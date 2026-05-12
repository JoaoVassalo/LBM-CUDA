#include "boundary.cuh"

#include "../../config/mom_config.cuh"
#include "../build/build_grid.cuh"

#include "../core/grid_id.cuh"
#include "../core/to_u8.cuh"

__device__ void applyBoundary(uint8_t *node, uint8_t *mask, float *mom_in, float *mom_out)
{
    const int x = blockIdx.x * blockDim.x + threadIdx.x;
    const int y = blockIdx.y * blockDim.y + threadIdx.y;

    const int index = grid_id();

    if (node[index] & to_u8(domainTags::Boundary))
    {
        boundary(mask[index], x, y, mom_in, mom_out);
    }
    else
    {
        center(x, y, mom_in, mom_out);
    }
}

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

__device__ void boundary(uint8_t mask_in, int x, int y, float *mom_in,
                         float *mom_out)
{
}
