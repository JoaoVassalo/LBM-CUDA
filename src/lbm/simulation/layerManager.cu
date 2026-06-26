#include "layerManager.cuh"

__global__ void defLayerMem(moments sim)
{
    const int x = threadIdx.x + blockDim.x * blockIdx.x;

    if (x >= Geometry::NX)
        return;

    constexprEnumFor<momId, static_cast<int>(momId::rho), static_cast<int>(momId::count)>([&](auto mom)
                                                                                          {
        constexpr int moment = static_cast<int>(decltype(mom)::value);

        constexprFor<0, Geometry::LNY>([&](auto ly)
                                   {
            constexpr int y = decltype(ly)::value;

            sim.layer[y][layerIdx<moment>(x)]=sim.mom[momIdx<moment>(x, y)]; 
        }); });
}

__global__ void swapLayerMem(moments sim, int y)
{
    const int x = threadIdx.x + blockDim.x * blockIdx.x;

    if (x >= Geometry::NX)
        return;

    constexprEnumFor<momId, static_cast<int>(momId::rho), static_cast<int>(momId::count)>([&](auto mom)
                                                                                          {
            constexpr int moment = static_cast<int>(decltype(mom)::value);

            const int layerIndex = layerIdx<moment>(x); 
        
            const float mid = sim.layer[1][layerIndex];
            const float top = sim.layer[2][layerIndex];

            sim.layer[0][layerIndex] = mid;
            sim.layer[1][layerIndex] = top;
            sim.layer[2][layerIndex] = sim.mom[momIdx<moment>(x, y)]; });
}