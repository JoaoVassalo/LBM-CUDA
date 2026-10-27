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
            constexpr int slot = decltype(ly)::value;
            constexpr int sourceY = slot == 0 ? 0 : slot - 1;

            sim.layer[slot][layerIdx<moment>(x)] = sim.mom[momIdx<moment>(x, sourceY)];
        }); });
}

__global__ void swapLayerMem(moments sim, int sourceY)
{
    const int x = threadIdx.x + blockDim.x * blockIdx.x;

    if (x >= Geometry::NX)
        return;

    constexprEnumFor<momId, static_cast<int>(momId::rho), static_cast<int>(momId::count)>([&](auto mom)
                                                                                          {
            constexpr int moment = static_cast<int>(decltype(mom)::value);

            const int layerIndex = layerIdx<moment>(x); 
        
            const varType mid = sim.layer[1][layerIndex];
            const varType top = sim.layer[2][layerIndex];

            sim.layer[0][layerIndex] = mid;
            sim.layer[1][layerIndex] = top;
            sim.layer[2][layerIndex] = sim.mom[momIdx<moment>(x, sourceY)]; });
}
