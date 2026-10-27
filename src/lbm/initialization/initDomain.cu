#include "initDomain.cuh"

__global__ void initDomain(varType *mom)
{
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;

    if (x >= Geometry::NX || y >= Geometry::NY)
        return;

    mom[momIdx<momId::rho>(x, y)] = 1.;
    mom[momIdx<momId::ux>(x, y)] = 0.;
    mom[momIdx<momId::uy>(x, y)] = 0.;
    mom[momIdx<momId::mxx>(x, y)] = 0.;
    mom[momIdx<momId::mxy>(x, y)] = 0.;
    mom[momIdx<momId::myy>(x, y)] = 0.;
}