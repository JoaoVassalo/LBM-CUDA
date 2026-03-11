#include "../presets/geometry.h"
#include "../presets/physics.h"

#include "collision.cuh"

__device__ void collision(int index, float *ux, float *uy, float *mxx, float *mxy, float *myy)
{
    mxx[index] = ((1 - omega) * mxx[index] + omega * ux[index] * ux[index]);
    mxy[index] = ((1 - omega) * mxy[index] + omega * ux[index] * uy[index]);
    myy[index] = ((1 - omega) * myy[index] + omega * uy[index] * uy[index]);
}