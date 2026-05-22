#include "build_layer.cuh"

template <int I>
__device__ __forceinline__ int grid_id_from_layer(int y)
{
    int x =
        blockIdx.x * blockDim.x +
        threadIdx.x;

    return x + (y + I - 1) * Geometry::Nx;
}

__device__ void first_layer(D2Q9 sim)
{
    int index = grid_id();

    sim.layer[0][layerIndex<rho>(index)] = sim.mom[momIdx<rho>(grid_id_from_layer<bot>(1))];
    sim.layer[0][layerIndex<ux>(index)] = sim.mom[momIdx<ux>(grid_id_from_layer<bot>(1))];
    sim.layer[0][layerIndex<uy>(index)] = sim.mom[momIdx<uy>(grid_id_from_layer<bot>(1))];
    sim.layer[0][layerIndex<mxx>(index)] = sim.mom[momIdx<mxx>(grid_id_from_layer<bot>(1))];
    sim.layer[0][layerIndex<myy>(index)] = sim.mom[momIdx<myy>(grid_id_from_layer<bot>(1))];
    sim.layer[0][layerIndex<mxy>(index)] = sim.mom[momIdx<mxy>(grid_id_from_layer<bot>(1))];

    sim.layer[1][layerIndex<rho>(index)] = sim.mom[momIdx<rho>(grid_id_from_layer<mid>(1))];
    sim.layer[1][layerIndex<ux>(index)] = sim.mom[momIdx<ux>(grid_id_from_layer<mid>(1))];
    sim.layer[1][layerIndex<uy>(index)] = sim.mom[momIdx<uy>(grid_id_from_layer<mid>(1))];
    sim.layer[1][layerIndex<mxx>(index)] = sim.mom[momIdx<mxx>(grid_id_from_layer<mid>(1))];
    sim.layer[1][layerIndex<myy>(index)] = sim.mom[momIdx<myy>(grid_id_from_layer<mid>(1))];
    sim.layer[1][layerIndex<mxy>(index)] = sim.mom[momIdx<mxy>(grid_id_from_layer<mid>(1))];

    sim.layer[2][layerIndex<rho>(index)] = sim.mom[momIdx<rho>(grid_id_from_layer<top>(1))];
    sim.layer[2][layerIndex<ux>(index)] = sim.mom[momIdx<ux>(grid_id_from_layer<top>(1))];
    sim.layer[2][layerIndex<uy>(index)] = sim.mom[momIdx<uy>(grid_id_from_layer<top>(1))];
    sim.layer[2][layerIndex<mxx>(index)] = sim.mom[momIdx<mxx>(grid_id_from_layer<top>(1))];
    sim.layer[2][layerIndex<myy>(index)] = sim.mom[momIdx<myy>(grid_id_from_layer<top>(1))];
    sim.layer[2][layerIndex<mxy>(index)] = sim.mom[momIdx<mxy>(grid_id_from_layer<top>(1))];
}

__device__ void other_layers(D2Q9 sim, int y, int i)
{
    int index = grid_id();

    sim.layer[i][layerIndex<rho>(index)] = sim.layer[i + 1][layerIndex<rho>(index)];
    sim.layer[i][layerIndex<ux>(index)] = sim.layer[i + 1][layerIndex<ux>(index)];
    sim.layer[i][layerIndex<uy>(index)] = sim.layer[i + 1][layerIndex<uy>(index)];
    sim.layer[i][layerIndex<mxx>(index)] = sim.layer[i + 1][layerIndex<mxx>(index)];
    sim.layer[i][layerIndex<myy>(index)] = sim.layer[i + 1][layerIndex<myy>(index)];
    sim.layer[i][layerIndex<mxy>(index)] = sim.layer[i + 1][layerIndex<mxy>(index)];

    if (i == sim.layer_size - 1)
    {
        sim.layer[sim.layer_size][layerIndex<rho>(index)] = sim.mom[momIdx<rho>(grid_id_from_layer<top>(y))];
        sim.layer[sim.layer_size][layerIndex<ux>(index)] = sim.mom[momIdx<ux>(grid_id_from_layer<top>(y))];
        sim.layer[sim.layer_size][layerIndex<uy>(index)] = sim.mom[momIdx<uy>(grid_id_from_layer<top>(y))];
        sim.layer[sim.layer_size][layerIndex<mxx>(index)] = sim.mom[momIdx<mxx>(grid_id_from_layer<top>(y))];
        sim.layer[sim.layer_size][layerIndex<myy>(index)] = sim.mom[momIdx<myy>(grid_id_from_layer<top>(y))];
        sim.layer[sim.layer_size][layerIndex<mxy>(index)] = sim.mom[momIdx<mxy>(grid_id_from_layer<top>(y))];
    }
}

__device__ void last_layer(D2Q9 sim, int y)
{
}

template <int I>
__device__ __forceinline__ int layerIndex(int index)
{
    constexpr int layerSize = layer::LNx * layer::LNy;

    int localId = rest(index, layerSize);

    return I * layerSize + localId;
}
