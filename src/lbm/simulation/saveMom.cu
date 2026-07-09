#include "saveMom.cuh"

#include <cstdio>

__device__ void saveMom(moments sim, int x, int y)
{
    int index = gridId(x, y);

    sim.mom[momIdx<momId::rho>(index)] = sim.layer[y][layerIdx<momId::rho>(x)];
    sim.mom[momIdx<momId::ux>(index)] = sim.layer[y][layerIdx<momId::ux>(x)];
    sim.mom[momIdx<momId::uy>(index)] = sim.layer[y][layerIdx<momId::uy>(x)];
    sim.mom[momIdx<momId::mxx>(index)] = sim.layer[y][layerIdx<momId::mxx>(x)];
    sim.mom[momIdx<momId::mxy>(index)] = sim.layer[y][layerIdx<momId::mxy>(x)];
    sim.mom[momIdx<momId::myy>(index)] = sim.layer[y][layerIdx<momId::myy>(x)];
}