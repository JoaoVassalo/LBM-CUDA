#pragma once

#include <cstdint>
#include "../lbm/build/build_mom.cuh"

__device__ __forceinline__ uint8_t to_u8(auto num)
{
    return static_cast<uint8_t>(num);
}