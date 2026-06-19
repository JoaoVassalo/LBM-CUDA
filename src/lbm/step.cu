#include "step.cuh"

__global__ void step(float **mom, float *layer, uint8_t *mask, uint8_t *node, int step_i)
{

    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    if (x == 0 && y == 0)
        printf("x = %d, y = %d, ", x, y);

    if (x >= Nx || y >= Ny)
        return;

    if (x == 1 && y == 45)
        printf("step?, ");

    propagation(x, y, mom, layer, mask, node, step_i);

    collision(mom[step_i]);
}