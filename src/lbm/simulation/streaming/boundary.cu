#include "boundary.cuh"

__device__ void center(moments sim)
{
    int x = threadIdx.x + blockDim.x * blockIdx.x;
    int y = threadIdx.y + blockDim.y * blockIdx.y;

    int index = gridId(x, y);

    float rho = 0.f;
    float ux = 0.f;
    float uy = 0.f;
    float mxx = 0.f;
    float mxy = 0.f;
    float myy = 0.f;

    float f = 0.f;

    for (int i = 1; i < D2Q9::Q; i++)
    {
        f += fi(x, (int)i, sim.layer);

        rho += f;
        ux += f * D2Q9::cx(i);
        uy += f * D2Q9::cy(i);
        mxx += f * (D2Q9::cx(i) * D2Q9::cx(i) - D2Q9::inv_as2);
        mxy += f * (D2Q9::cx(i) * D2Q9::cy(i));
        myy += f * (D2Q9::cy(i) * D2Q9::cy(i) - D2Q9::inv_as2);
    }
}