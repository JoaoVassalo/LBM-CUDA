#include "layer_manager.cuh"

__global__ void init_layers(D2Q9 sim, layer current_layer)
{
    const int x = blockIdx.x * blockDim.x + threadIdx.x;

    if (x >= layer::LNx)
        return;

    constexpr_enum_for<MomentId, static_cast<int>(MomentId::rho), static_cast<int>(MomentId::Count)>([&](auto mom)
    {
        constexpr int moment = static_cast<int>(decltype(mom)::value);

        constexpr_for<0, layer::LNy>([&](auto ly)
        {
            constexpr int slot = decltype(ly)::value;
            const int y = wrap_y(current_layer.yref + slot - 1);
            const int domain_index = grid_index_from_xy(x, y);

            layer_moment_write<moment>(
                current_layer,
                slot,
                x,
                sim.mom[momIdx<moment>(domain_index)]);
        });
    });
}

__global__ void swap_layers(D2Q9 sim, layer current_layer)
{
    const int x = blockIdx.x * blockDim.x + threadIdx.x;

    if (x >= layer::LNx)
        return;

    constexpr_enum_for<MomentId, static_cast<int>(MomentId::rho), static_cast<int>(MomentId::Count)>([&](auto mom)
    {
        constexpr int moment = static_cast<int>(decltype(mom)::value);

        const int layer_index = layer_moment_index<moment>(x);
        const int next_y = wrap_y(current_layer.yref + 2);
        const int next_domain_index = grid_index_from_xy(x, next_y);

        const float middle = current_layer.buffer[1][layer_index];
        const float front = current_layer.buffer[2][layer_index];

        current_layer.buffer[0][layer_index] = middle;
        current_layer.buffer[1][layer_index] = front;
        current_layer.buffer[2][layer_index] = sim.mom[momIdx<moment>(next_domain_index)];
    });
}
