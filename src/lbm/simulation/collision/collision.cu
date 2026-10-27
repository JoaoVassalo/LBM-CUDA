#include "collision.cuh"

__device__ void collide(varType *mom, int index)
{
    const varType ux = mom[momIdx<momId::ux>(index)];
    const varType uy = mom[momIdx<momId::uy>(index)];
    const varType mxx = mom[momIdx<momId::mxx>(index)];
    const varType mxy = mom[momIdx<momId::mxy>(index)];
    const varType myy = mom[momIdx<momId::myy>(index)];

    mom[momIdx<momId::mxx>(index)] = ((1. - physics::omega) * mxx + physics::omega * ux * ux);
    mom[momIdx<momId::mxy>(index)] = ((1. - physics::omega) * mxy + physics::omega * ux * uy);
    mom[momIdx<momId::myy>(index)] = ((1. - physics::omega) * myy + physics::omega * uy * uy);
}
