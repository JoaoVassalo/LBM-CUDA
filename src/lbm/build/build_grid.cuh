#pragma once

#include <cstdint>

struct Grid2D;

enum class domainTags : uint8_t
{
    Fluid = 0,
    Boundary = 1,
};

__host__ __device_builtin__ __forceinline__ int income(int i);

__host__ void build_grid(D2Q9 sim);

__device__ void init_grid(uint8_t *node, uint8_t *mask);