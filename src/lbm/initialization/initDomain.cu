#include "initDomain.cuh"

__global__ void initDomain(float *mom)
{
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;

    if (x >= Geometry::NX || y >= Geometry::NY)
        return;

    mom[momIdx<momIdD2Q9::rho>(x, y)] = 1.f;
    mom[momIdx<momIdD2Q9::ux>(x, y)] = 0.f;
    mom[momIdx<momIdD2Q9::uy>(x, y)] = 0.f;
    mom[momIdx<momIdD2Q9::mxx>(x, y)] = 0.f;
    mom[momIdx<momIdD2Q9::mxy>(x, y)] = 0.f;
    mom[momIdx<momIdD2Q9::myy>(x, y)] = 0.f;
}