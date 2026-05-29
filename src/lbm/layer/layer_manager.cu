#include "layer_manager.cuh"

__global__ void init_layers(D2Q9 sim, layer layer)
{

    constexpr_enum_for<MomentId, static_cast<int>(MomentId::rho), static_cast<int>(MomentId::Count)>([&](auto M)
                                                                                                     {
        constexpr MomentId mom = M;

        constexpr_for<(int) 0, layer::LNy>([&] (auto ly)
        {

            constexpr int i = ly;

            layer.layer[ly][layerIndex<]

        }); });

    // constexpr_enum_for<MomentId, static_cast<int>(MomentId::rho), static_cast<int>(MomentId::Count)>([&](const auto mom)
    //                                                                                                  {
    //     constexpr MomentId momtype = decltype(mom)::value;
    //             constexpr_for< 0, layer::LNy>([&](const auto ly)
    //         {
    //             constexpr int i = decltype(ly - layer::LNy)::value;

    //             layer.layer[i][layerIndex<static_cast<int>(decltype(mom)::value)>(grid_id_from_layer(i))] =
    //             sim.mom[
    //                 momIdx<static_cast<int>(decltype(mom)::value)>(
    //                     grid_id_from_layer(i)
    //                 )
    //             ];
    //         }); });
}

__global__ void swap_layers(D2Q9 sim, layer layer)
{
    constexpr_enum_for<MomentId, static_cast<int>(MomentId::rho), static_cast<int>(MomentId::Count)>([&](auto mom)
                                                                                                     {
        constexpr MomentId momtype = decltype(mom)::value;
                constexpr_for< 0, layer::LNy>([&](auto ly)
            {
                constexpr int i = decltype(ly)::value;

                layer.layer[i][layerIndex<static_cast<int>(decltype(mom)::value)>(grid_id_from_layer(i))] = 
                layer.layer[i+1][layerIndex<static_cast<int>(decltype(mom)::value)>(grid_id_from_layer(i))];

                if constexpr (i==layer::LNy)
                {
                    layer.layer[i][layerIndex<static_cast<int>(decltype(mom)::value)>(grid_id_from_layer(i))] = 
                    sim.mom[
                    momIdx<static_cast<int>(decltype(mom)::value)>(
                        grid_id_from_layer(i)
                    )];
                }
            }); });
}

__global__ void final_layers(D2Q9 sim, layer layer)
{
}
