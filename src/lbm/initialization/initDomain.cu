#include "initDomain.cuh"

__global__ void initDomain(float *mom)
{
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;

    if (x >= Geometry::NX || y >= Geometry::NY)
        return;

    mom[momIdx<momId::rho>(x, y)] = 1.f;
    mom[momIdx<momId::ux>(x, y)] = 0.f;
    mom[momIdx<momId::uy>(x, y)] = 0.f;
    mom[momIdx<momId::mxx>(x, y)] = 0.f;
    mom[momIdx<momId::mxy>(x, y)] = 0.f;
    mom[momIdx<momId::myy>(x, y)] = 0.f;
}