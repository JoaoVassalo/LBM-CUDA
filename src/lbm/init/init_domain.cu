#include "../../config/stencil.cuh"
#include "../../config/geometry.h"
#include "../../config/mom_config.cuh"

#include "../build/build_mom.cuh"

#include "../../core/indexing.cuh"

__global__ void initDomain(float *mom)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    if (x >= Geometry::Nx || y >= Geometry::Ny)
        return;

    mom[momIdx<MomentId::rho>(x, y)] = 1.f; // rho
    mom[momIdx<MomentId::ux>(x, y)] = 0.f;  // ux
    mom[momIdx<MomentId::uy>(x, y)] = 0.f;  // uy
    mom[momIdx<MomentId::mxx>(x, y)] = 0.f; // mxx
    mom[momIdx<MomentId::mxy>(x, y)] = 0.f; // mxy
    mom[momIdx<MomentId::myy>(x, y)] = 0.f; // myy
}
