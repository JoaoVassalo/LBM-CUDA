#include "collision.cuh"

__device__ void collision(float *mom)
{
    int index = grid_id();

    const float ux = mom[momIdx<MomentId::ux>(index)];
    const float uy = mom[momIdx<MomentId::uy>(index)];

    const float mxx = mom[momIdx<MomentId::mxx>(index)];
    const float myy = mom[momIdx<MomentId::myy>(index)];
    const float mxy = mom[momIdx<MomentId::mxy>(index)];

    mom[momIdx<MomentId::mxx>(index)] = ((1.0f - omega) * mxx + omega * ux * ux); // mxx, mxx, ux, ux
    mom[momIdx<MomentId::myy>(index)] = ((1.0f - omega) * myy + omega * uy * uy); // myy, myy, uy, uy
    mom[momIdx<MomentId::mxy>(index)] = ((1.0f - omega) * mxy + omega * ux * uy); // mxy, mxy, ux, uy
}