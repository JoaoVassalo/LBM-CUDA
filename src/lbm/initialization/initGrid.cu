#include "initGrid.cuh"

#include <cstdint>

__global__ void initGrid(uint8_t *mask, uint8_t *node)
{
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;

    if (x >= Geometry::NX || y >= Geometry::NY)
        return;

    const int index = gridId(x, y);

    node[index] = to_u8(defBC(x, y));

    uint8_t m = 0u;

    for (int i = 1; i < D2Q9::Q; i++)
    {
        int xn = x + c_ix[i];
        int yn = y + c_iy[i];

        if (xn < 0 || xn >= Geometry::NX ||
            yn < 0 || yn >= Geometry::NY)
        {
            continue;
        }
    }
}