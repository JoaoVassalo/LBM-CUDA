#include "../presets/geometry.h"
#include "../presets/physics.h"

#include "collision.cuh"
#include "cuda_config/var.cuh"

__device__ void collision(int index, float *mom)
{
    const float ux = mom[momIdx<MomentId::ux>(index)];
    const float uy = mom[index + MomentId::uy];

    const float mxx = mom[index + MomentId::mxx];
    const float myy = mom[index + MomentId::myy];
    const float mxy = mom[index + MomentId::mxy];

    mom[index + 3] = ((1.0f - omega) * mxx + omega * ux * ux); // mxx, mxx, ux, ux
    mom[index + 4] = ((1.0f - omega) * myy + omega * uy * uy); // myy, myy, uy, uy
    mom[index + 5] = ((1.0f - omega) * mxy + omega * ux * uy); // mxy, mxy, ux, uy
}