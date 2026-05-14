#include "build_layer.cuh"

#include "../grid_id.cuh"
#include "../rest.cuh"
#include "../../lbm/build/build_mom.cuh"
#include "../../config/geometry.h"

__global__ void build_layer(int Nlayer, D2Q9 sim)
{

    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    if (x >= layer::LNx || y >= layer::LNy)
        return;

    int index = grid_id();

    switch (Nlayer)
    {
    case 0:
        first_layer();
        break;
    default:
        other_layers();
        break;
    }
}

__device__ __forceinline__ void first_layer(D2Q9 sim, int index)
{
    sim.layer[layerIndex<rho>(index)] = sim.mom[momIdx<rho>(index)];
    sim.layer[layerIndex<ux>(index)] = sim.mom[momIdx<ux>(index)];
    sim.layer[layerIndex<uy>(index)] = sim.mom[momIdx<uy>(index)];
    sim.layer[layerIndex<mxx>(index)] = sim.mom[momIdx<mxx>(index)];
    sim.layer[layerIndex<mxy>(index)] = sim.mom[momIdx<mxy>(index)];
    sim.layer[layerIndex<myy>(index)] = sim.mom[momIdx<myy>(index)];
}

__device__ __forceinline__ void other_layers(D2Q9 sim, int index)
{
    sim.layer[layerIndex<rho>(index)] = sim.mom[momIdx<rho>(index)];
    sim.layer[layerIndex<ux>(index)] = sim.mom[momIdx<ux>(index)];
    sim.layer[layerIndex<uy>(index)] = sim.mom[momIdx<uy>(index)];
    sim.layer[layerIndex<mxx>(index)] = sim.mom[momIdx<mxx>(index)];
    sim.layer[layerIndex<mxy>(index)] = sim.mom[momIdx<mxy>(index)];
    sim.layer[layerIndex<myy>(index)] = sim.mom[momIdx<myy>(index)];
}

template <int I>
__device__ __forceinline__ int layerIndex(int index)
{
    constexpr int layerSize = layer::LNx * layer::LNy;

    int localId = rest(index, layerSize);

    return I * layerSize + localId;
}