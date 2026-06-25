#include "collision.cuh"

#include <iostream>

__device__ void collide(moments sim, int y)
{
    const int x = threadIdx.x + blockDim.x * blockIdx.x;

    if (x >= Geometry::NX || x < 0)
        printf("Out of bounds"); // DEBUG
    return;

    const float ux = sim.layer[y][layerIdx<momId::ux>(x)];
    const float uy = sim.layer[y][layerIdx<momId::uy>(x)];
    const float mxx = sim.layer[y][layerIdx<momId::mxx>(x)];
    const float mxy = sim.layer[y][layerIdx<momId::mxy>(x)];
    const float myy = sim.layer[y][layerIdx<momId::myy>(x)];

    sim.layer[y][layerIdx<momId::mxx>(x)] = ((1.f - physics::omega) * mxx + ux * ux);
    sim.layer[y][layerIdx<momId::mxy>(x)] = ((1.f - physics::omega) * mxy + ux * uy);
    sim.layer[y][layerIdx<momId::myy>(x)] = ((1.f - physics::omega) * myy + uy * uy);
}