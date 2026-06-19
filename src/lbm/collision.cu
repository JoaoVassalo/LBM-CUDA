#include "../config/geometry.h"
#include "../config/physics.h"
#include "../core/indexing.cuh"

#include "collision.cuh"
#include "build/build_mom.cuh"

__global__ void collide_layer_at(D2Q9 sim, int y)
{
    const int x = blockIdx.x * blockDim.x + threadIdx.x;

    if (x >= Geometry::Nx)
        return;

    const int index = grid_id(x, y);

    const float ux = sim.mom[momIdx<MomentId::ux>(index)];

    const float uy = sim.mom[momIdx<MomentId::uy>(index)];

    const float mxx = sim.mom[momIdx<MomentId::mxx>(index)];
    const float myy = sim.mom[momIdx<MomentId::myy>(index)];
    const float mxy = sim.mom[momIdx<MomentId::mxy>(index)];

    sim.mom[momIdx<MomentId::mxx>(index)] = ((1.0f - omega) * mxx + omega * ux * ux); // mxx, mxx, ux, ux
    sim.mom[momIdx<MomentId::myy>(index)] = ((1.0f - omega) * myy + omega * uy * uy); // myy, myy, uy, uy
    sim.mom[momIdx<MomentId::mxy>(index)] = ((1.0f - omega) * mxy + omega * ux * uy); // mxy, mxy, ux, uy
}