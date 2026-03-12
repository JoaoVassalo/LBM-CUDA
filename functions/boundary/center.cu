#include "../../presets/config.h"
#include "../../presets/stencil.cuh"

#include "../equations/f_i.cuh"
#include "../grid_id.cuh"
#include "../from_id.cuh"

__device__ void center(CInt *I_s, CInt *O_s, int x, int y,
                       float *rho_in, float *ux_in, float *uy_in, float *mxx_in, float *mxy_in, float *myy_in,
                       float *rho_out, float *ux_out, float *uy_out, float *mxx_out, float *mxy_out, float *myy_out)
{
    int index = grid_id(x, y);
    float rho = 0;

#pragma unroll
    for (int i = 0; i < Q; i++)
    {
        int index_from =
            rho += f_i(index)
    }
}