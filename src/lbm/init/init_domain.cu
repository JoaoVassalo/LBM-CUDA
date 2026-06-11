#include "init_domain.cuh"

__global__ void initDomain(float *mom)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    if (x >= Geometry::Nx || y >= Geometry::Ny)
        return;

    int index = grid_id(x, y);

    mom[momIdx<MomentId::rho>(index)] = 1.f; // rho
    mom[momIdx<MomentId::ux>(index)] = 0.f;  // ux
    mom[momIdx<MomentId::uy>(index)] = 0.f;  // uy
    mom[momIdx<MomentId::mxx>(index)] = 0.f; // mxx
    mom[momIdx<MomentId::mxy>(index)] = 0.f; // mxy
    mom[momIdx<MomentId::myy>(index)] = 0.f; // myy

    if (x == 1)
    {
        mom[momIdx<MomentId::ux>(index)] = u_max; // ux
    }
}