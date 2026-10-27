#include "saveMom.cuh"

#include <cstdio>

__device__ void saveMom(moments sim, int x, int globalY, int layerSlot)
{
    int index = gridId(x, globalY);

    sim.mom[momIdx<momId::rho>(index)] = sim.layer[layerSlot][layerIdx<momId::rho>(x)];
    sim.mom[momIdx<momId::ux>(index)] = sim.layer[layerSlot][layerIdx<momId::ux>(x)];
    sim.mom[momIdx<momId::uy>(index)] = sim.layer[layerSlot][layerIdx<momId::uy>(x)];
    sim.mom[momIdx<momId::mxx>(index)] = sim.layer[layerSlot][layerIdx<momId::mxx>(x)];
    sim.mom[momIdx<momId::mxy>(index)] = sim.layer[layerSlot][layerIdx<momId::mxy>(x)];
    sim.mom[momIdx<momId::myy>(index)] = sim.layer[layerSlot][layerIdx<momId::myy>(x)];
}
